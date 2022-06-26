using Images, FileIO
using Dates, TimeZones
using Plots
using ImageMagick: magickinfo
using DataFrames, CSV
using ProgressMeter
using BenchmarkTools


# 20201123_113254_869_IR.TIFF

# some test paths
basepath = "/media/john/HSDATA/FLIR/20201123_112658"
ispath(basepath)

jpgPath = "/media/john/HSDATA/FLIR/20201123_112658/20201123_113418_852_8b.JPG"
tiffPath = "/media/john/HSDATA/FLIR/20201123_112658/20201123_113418_852_IR.TIFF"

jpgPath = "/media/john/HSDATA/FLIR/20201123_112658/20201123_113254_869_8b.JPG"
tiffPath = "/media/john/HSDATA/FLIR/20201123_112658/20201123_113254_869_IR.TIFF"

isfile(jpgPath)
isfile(tiffPath)

jpg = load(jpgPath)
tiff = load(tiffPath)

p1 = plot(jpg, title="visible", aspect_ratio=:equal)
p2 = heatmap(channelview(tiff), title="thermal", aspect_ratio=:equal, yflip=true)
plot(p1, p2)

# get all metadata from jpg file
p = magickinfo(jpgPath, magickinfo(jpgPath))
q = magickinfo(tiffPath, magickinfo(tiffPath))



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

D = exiftool(tiffPath)
date = getExifDate(D)
date2 = getExifDate(tiffPath)


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


df = getFlirTimes(basepath)

names(df)

df.times_central[1]
df.times_utc[1]

# @btime getFlirTimes(basepath)


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


flirpath = "/media/john/HSDATA/FLIR"
generateFlirSummaries(flirpath)

test_csv = joinpath(flirpath, "20201123_112658.csv")
df = CSV.File(test_csv) |> DataFrame

df.times[end]



using georectification
df_dye1 = masterLCF("/media/john/HSDATA/raw/12-09", "Dye_1")
df_dye2 = masterLCF("/media/john/HSDATA/raw/12-09", "Dye_2")
df_nodye1 = masterLCF("/media/john/HSDATA/raw/12-09", "NoDye_1")
df_nodye2 = masterLCF("/media/john/HSDATA/raw/12-09", "NoDye_2")

size(df_dye1)
size(df_dye2)
size(df_nodye1)
size(df_nodye2)

df_nodye1.tstart[1]
df_nodye1.tend[end]


df.times[1]
df.times[end]

"""
    DuringFlight(thermalPath::String, df::DataFrame)

Check each row of the dataframe (coming from masterLCF() in IMU.jl) and determine if the FLIR image was taken during the HSI acquisition.
"""
function DuringFlight(thermalPath::String, df::DataFrame)
    t, date_time = FLIRtime(thermalPath)
    for row ∈ eachrow(df)
        if row.tstart <= t && t <= row.tend
            for f ∈ readdir(row.paths)
                if endswith(f, ".lcf")
                    return true, joinpath(row.paths, f), row.paths
                end
            end
        end
    end
    return false, "", ""
end




function getMasterFlightPath(basepath)
    # loop through all HSI folders to generate dataframe of combined LCF values
end



function wasDuringFlight(time, flightPath_df)
    # check if a given time occurred during the flight (we need to filter for ground photos)
end


function duringFlight(flir_df, flightPath_df)
    # given the data_frame with image times, return the subset which occur during the flight
end

function getPositionOrientation(flir_df, flightPath_df)
    # update the flir_df with interpolated position/orientation data from flightPath_df
end


getImageTime(tiffPath)
FLIRtime(tiffPath)
tiff_temp = TIFFtoTemp(tiff)

p1 = plot(jpg, title="visible", aspect_ratio=:equal)
p2 = heatmap(tiff_temp, title="thermal", aspect_ratio=:equal, yflip=true)
plot(p1, p2)






