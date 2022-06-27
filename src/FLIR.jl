using Images, FileIO
using Dates, TimeZones
using DataFrames, CSV
using ProgressMeter



# see duo-pro-r_manual.pdf page 39
"""
    TIFFtoTemp(img)

Convert a FLIR thermal image to temperature units in °C
"""
function TIFFtoTemp(img::AbstractArray)
    #NOTE: Images.jl interprets Uint16 values to N0f16 for us automatically
    # this means values are scaled to 0-1
    return rawview(channelview(img)) .* 0.04 .- 273.15
end





"""
    function exiftool(imgpath)

Run the command line tool `exiftool` on the image at `imgpath`. Returns a dictionary of metadata from the image file. 
"""
function exiftool(imgpath)
    cmd = `exiftool $(imgpath)`
    txt = read(cmd, String)


    ret_dict = Dict{String, String}()
    for line ∈ split(txt, "\n")
        splitline = split(line, ":")
        if size(splitline, 1) >= 2
            key = splitline[1]
            val = join(splitline[2:end], ":")

            # remove whitespace
            key = strip(key)
            val = strip(val)

            # add to dict
            ret_dict[key] = val
        end
    end

    return ret_dict
end





"""
    function getExifDate(dict)

Given a dictionary `dict` returned from `exiftool`, return the correct time field as a `DateTime` object.
"""
function getExifDate(dict::Dict)
    date = dict["Date/Time Original"]
    res = DateTime(date, "Y:m:d H:M:S.s")
end



"""
    function getExifDate(tiffpath::String)

Dispatch that returns the time of creation of a FLIR tiff file using `exiftool`
"""
function getExifDate(tiffpath::String)
    getExifDate(exiftool(tiffpath))
end






"""
    function getFlirTimes(dirPath, tz)

Given directory `dirPath`, generate a list of all FLIR files and their associated time of capture. Assume the time used is in timezone `tz`. See the [TimeZones.jl docs](https://juliatime.github.io/TimeZones.jl/dev/) for more information.
"""
function getFlirTimes(dirPath, tz)
    flist = filter(x->endswith(x, ".TIFF"), readdir(dirPath))
    basenames = [join(split(f, "_")[1:end-1], "_") for f ∈ flist]
    unique!(basenames)

    tiff_end = "_IR.TIFF"
    jpg_end = "_8b.JPG"


    out_df = DataFrame()

    out_df.basename = basenames
    out_df.tiffpath = [joinpath(dirPath, b * tiff_end) for b ∈ basenames]
    out_df.jpgpath = [joinpath(dirPath, b * jpg_end) for b ∈ basenames]

    # preallocate times
    out_df.times_central = Array{DateTime}(undef, size(out_df, 1))
    out_df.basename_times_central = Array{DateTime}(undef, size(out_df, 1))
    # generate UTC versions from the Central Standard Time. Note: the conversion depends on whether or not we're in daylight savings time
    out_df.times_utc = Array{DateTime}(undef, size(out_df, 1))
    out_df.times_basename_utc = Array{DateTime}(undef, size(out_df, 1))




    # process time in parallel
    p = Progress(size(out_df, 1); dt=1, desc="Getting image capture DateTimes...")
    Threads.@threads for i ∈ 1:size(out_df,1)
        out_df.times_central[i] = getExifDate(out_df.tiffpath[i])
        out_df.basename_times_central[i] = DateTime(out_df.basename[i], "YYYYmmdd_HHMMSS_sss")

        time_central = out_df.times_central[i]
        time_basename_central = out_df.basename_times_central[i]

        # make into zoned time and then convert to utc and get back a regular datetime
        out_df.times_utc[i] = DateTime(ZonedDateTime(time_central, tz), UTC)
        out_df.times_basename_utc[i] = DateTime(ZonedDateTime(time_basename_central, tz), UTC)
        next!(p)
    end


    sort!(out_df, :times_utc)

    return out_df
end





"""
    function getFlirTimes(dirPath)

Calls `getFlirTimes` on `dirPath` with the timezone assumed to be US/Central
"""
function getFlirTimes(dirPath)
    tz_central = TimeZone("US/Central", TimeZones.Class(:LEGACY))
    getFlirTimes(dirPath, tz_central)
end




"""
    function generateFlirSummaries(basepath::String)

Loop through directories in `basepath` and generate dataframes with capture times for all FLIR images in each folder.
"""
function generateFlirSummaries(basepath::String)
    dirlist = filter(x -> !endswith(x, ".csv"), readdir(basepath))
    for dir ∈ dirlist
        if isdir(joinpath(basepath, dir))
            println("working on $(dir)...")
            try
                df = getFlirTimes(joinpath(basepath, dir))
                CSV.write(joinpath(basepath, dir*".csv"), df)
            catch e
                println("Didn't work...\n")
                println(e)
            end
        end
    end
end



"""
    function DuringFlight(flir_row, df_lcf)

Check if a row from a Flir Summary csv, `flir_row` occurred during flight detailed in `df_lcf`
"""
function DuringFlight(flir_row::DataFrameRow, df_lcf::DataFrame)
    flight_start = df_lcf.tstart[1]
    flight_end = df_lcf.tend[end]

    tmin = minimum([flir_row.times_utc, flir_row.times_basename_utc])
    tmax = maximum([flir_row.times_utc, flir_row.times_basename_utc])

    (flight_start <= tmin && tmax <= flight_end) ? true : false
end



"""
    function DuringFlight(df_flir, df_lcf)

Given a master lcf list `df_lcf`, return a all rows of `df_flir` that occur during the flight
"""
function DuringFlight(df_flir::DataFrame, df_lcf::DataFrame)
    idx = [DuringFlight(df_flir[i,:], df_lcf) for i ∈ 1:size(df_flir, 1)]
    return df_flir[idx, :]
end





"""
    function matchLCfile(df_flir, df_lcf)

Given a table of flir files and the master lcf file, add a column specifying which lcf_path goes with each row in the flir_df. Rows that do not have a match are given type `:missing`.
"""
function matchLCFfile!(df_flir, df_lcf)
    # add new column for paths to df_flir
    # NOTE: when the Drone goes out-of-bounds, we don't collect imagery so there are some images do not have an LCF path
    df_flir.lcf_path = Array{Union{String, Missing}}(missing, size(df_flir, 1))
    df_flir.lcf_path_basename = Array{Union{String, Missing}}(missing, size(df_flir, 1))

    # df_flir.lcf_path = ["N/A" for _ ∈ 1:size(df_flir, 1)]
    # df_flir.lcf_path_basename = ["N/A" for _ ∈ 1:size(df_flir, 1)]


    # brute force search for the matching lcf file
    for row_flir ∈ eachrow(df_flir)
        for row_lcf ∈ eachrow(df_lcf)
            if (row_lcf.tstart <= row_flir.times_utc && row_flir.times_utc <= row_lcf.tend)
                row_flir.lcf_path = row_lcf.files
            end
        end
    end

    for row_flir ∈ eachrow(df_flir)
        for row_lcf ∈ eachrow(df_lcf)
            if (row_lcf.tstart <= row_flir.times_basename_utc && row_flir.times_basename_utc <= row_lcf.tend)
                row_flir.lcf_path_basename = row_lcf.files
            end
        end
    end

    basepath = join(split(df_flir.tiffpath[1], "/")[1:end-2], "/")
    basename = split(df_flir.tiffpath[1], "/")[end-1]
    basename = basename * "_matched.csv"
    CSV.write(joinpath(basepath, basename), df_flir)

end

