using Images, FileIO
using Dates
using Plots
using ImageMagick: magickinfo


20201123_113254_869_IR.TIFF


# some test paths
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




"""
    FLIRtime(path::String)

Examine metadata for a thermal FLIR image and return the time when the photo was taken.
"""
function FLIRtime(thermal_path::String)
    cmd_therm = `exiftool $(thermal_path)`
    txt = read(cmd_therm, String)
    lines = split(txt, "\n")

    i = 1
    for line ∈ lines
        if occursin("Date/Time Original", line)
            if i == 2
                dt = split(line, ":")[2:end]
                year = parse(Int, dt[1])
                month = parse(Int, dt[2])
                day, hour = parse.(Int, split(dt[3], " "))
                hour += 6 # NOTE texas is central daylight time. this is +5 hr ahead of UTC time
                min = parse(Int, dt[4])
                sec, ms = parse.(Int, split(dt[5], "."))

                # sec += 18 #GPS time is offset from UTC by 'leapseconds'

                res = DateTime(year, month, day, hour, min, sec, ms*10)
                res = res + Dates.Second(18)  #GPS time is offset from UTC by 'leapseconds'
                return datetime2unix(res), res
            end
            i += 1
        end
    end
end



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


# see duo-pro-r_manual.pdf page 39
"""
    TIFFtoTemp(img)

Convert a FLIR thermal image to temperature units in °C
"""
function TIFFtoTemp(img)
    #NOTE: Images.jl interprets Uint16 values to N0f16 for us automatically
    # this means values are scaled to 0-1
    return rawview(channelview(img)) .* 0.04 .- 273.15
end







function getImageList(basepath)
    # generate a list of all (JPG, TIFF) pairs for a given directory
    # return a dataframe with columns JPG_path, TIFF_path
end


function getImageTime(tiffPath)
    p = magickinfo(tiffPath, magickinfo(tiffPath))
    T = p["exif:DateTimeOriginal"]
end

function addTimesToDF(flir_df)
    # generate new column of image times using the TIFF_path column
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
