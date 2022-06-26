using Dates

"""
    function leap_count(year::Int)

Determine the number of leap seconds introduced before the given date. See the links below:
- [Stack Overflow Post](https://stackoverflow.com/questions/33415475/how-to-get-current-date-and-time-from-gps-unsegment-time-in-python)
- [Official Leap Seconds Table](https://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en)
"""
function leap_count(year::Int)
    leap_years =[1980,
                 1981,
                 1982,
                 1983,
                 1985,
                 1988,
                 1990,
                 1991,
                 1992,
                 1993,
                 1994,
                 1996,
                 1997,
                 1999,
                 2006,
                 2009,
                 2012,
                 2015,
                 2017,
                 ];

    leap_secs = [19,
                 20,
                 21,
                 22,
                 23,
                 24,
                 25,
                 26,
                 27,
                 28,
                 29,
                 30,
                 31,
                 32,
                 33,
                 34,
                 35,
                 36,
                 37,];


    idx = findfirst(year .< leap_years)
    if idx == nothing
        return leap_secs[end]
    else
        return leap_secs[idx - 1]
    end

end


"""
    function gpsToUTC(gps_sec, year)

Given the gps_time in seconds and the current year, return the UTC time accounting for leap seconds. See this [Stack Overflow Post](https://stackoverflow.com/questions/33415475/how-to-get-current-date-and-time-from-gps-unsegment-time-in-python) for more details.
"""
function gpsToUTC(gps_sec, year)
    sec = round(gps_sec)
    ms = (gps_sec - sec)
    # convert to DateTime object
    sec = Second(Int(sec))
    ms = Millisecond(Int(round(1000*ms, digits=0)))

    gps = sec + ms

    utc = DateTime(1980, 1, 6) + (gps - Second(leap_count(year) - leap_count(1980)) )
end


# try it and compare with https://www.andrews.edu/~tzs/timeconv/timeconvert.php
gpsToUTC(1290187878.260, 2019)




"""
    getIMUdata(pathToLCF::String)

Read a .lcf file and return a DataFrame `df` with corresponding IMU data and the value `start_time` indicating the start time of the clock in linux epoch time.

# Data Returned
- **time:** system time in seconds starting at zero
- **roll:** angle in radians.
- **pitch:** angle in radians.
- **heading:** aka yaw. angle in radians.
- **longitude:** GPS longitude in degrees, negative for west longitudes
- **latitude:** GPS latitude in degrees
- **altitude**: height in meters above the WGS-54 ellipsoid
- **x**: x coordinate in meters in UTMz coordinate system
- **y**: y coordinate in meters in UTMz coordinate system
- **z**: z coordinate in meters in UTMz coordinate system
- **zone**: the UTM zone for the local coordinate system
- **isnorth**: boolean for whether or not position is above equator

# NOTES
The conversion between GPS time and Unix Epoch time was a little confusing. It is good to double check that this is correct (how should we deal with leap seconds?). This is relevant for parsing data from the FLIR as well.
"""
function getIMUdata(pathToLCF::String)

    df = DataFrame(CSV.File(pathToLCF, header=[:time, :roll, :pitch, :heading, :longitude, :latitude, :altitude, :placeHolder1, :placeHolder2, :placeHolder3, :placeHolder4]))
    select!(df, 1:7)

    # now we want to compute the positions
    x = Float64[]
    y = Float64[]
    z = Float64[]
    isnorth = Bool[]
    zone = UInt8[]
    for row in eachrow(df[!, [:latitude, :longitude, :altitude]])
        ϕ, θ, alt = row
        x_lla = LLA(ϕ, θ, alt)
        utmz = UTMZ(x_lla, wgs84)
        push!(x, utmz.x)
        push!(y, utmz.y)
        push!(z, utmz.z)
        push!(isnorth, utmz.isnorth)
        push!(zone, utmz.zone)
    end

    df.x = x
    df.y = y
    df.z = z
    df.isnorth = isnorth
    df.zone = zone


    # use 2019 for now. May need to update once the leap_years table changes

    # start_time = df.time[1]  # note this is gps time (i.e. the number of seconds since 0h 6-Jan-1980, UTC)
    # start_time += 315964800.00  # convert to linux epoch time to make working with DateTime object easier
    start_time = gpsToUTC(df.time[1], 2019)


    # now we want to set the times to start at 0.0 seconds to make it easier to understand
    df.time .= (df.time .- df.time[1])


    # correct heading for meridian deviation as in Muller paper
    utm_zones = range(-180, stop=180, step=6)  # utm zones are every 6 degrees


    Δ = atan.(tand.(df.longitude .- (utm_zones[df.zone] .+ 3.0)).*sind.(df.latitude))
    df.heading_correct = df.heading  .- Δ  # <-- this is what's used in the paper by Muller
    # df.heading_correct = df.heading


    # final assignments
    df.roll .= -df.roll  # <-- due to opposite convention used by GPS
    df.pitch .= df.pitch
    df.heading_correct .= df.heading_correct

    return df, start_time
end




# define rotations as in paper by Yeh "Direct georeferencing of airborne pushbroom images"
# Rz(κ) = [cos(κ) -sin(κ) 0; sin(κ) cos(κ) 0; 0 0 1]
# Ry(ϕ) = [cos(ϕ) 0 sin(ϕ); 0 1 0; -sin(ϕ) 0 cos(ϕ)]
# Rx(ω) = [1 0 0; 0 cos(ω) -sin(ω); 0 sin(ω) cos(ω)]
# Rotate(heading, pitch, roll) = Rz(heading)*Ry(pitch)*Rx(roll)

# use static arrays implementation from Rotations.jl

"""
    Rotate(heading, pitch, roll)

Return a rotation matrix corresponding to the given orientation angles.
"""
Rotate(heading, pitch, roll) = RotZ(heading)*RotY(pitch)*RotX(roll)
T_n_E = @SMatrix [0 1 0; 1 0 0; 0 0 -1 ]  # matrix to convert navigation system (IMU) to the Earth System (UTM) from Baumker paper page 7





"""
    masterLCF(folder::String)

Loop through all .lcf files in `folder` and return a DataFrame containing the start and end times for each datacube aquisition.

# Fields
- **paths**: Paths to each HSI directory
- **files**: file name for .lcf files
- **tstart**: starting time for each datacube acquisition
- **tend**: ending time for each datacube acquisition
"""
function masterLCF(folder::String, basename::String)
    lcf_paths = []
    lcf_files = []
    for (root, dirs, files) ∈ walkdir(folder, topdown=true)
        for file ∈ files
            if endswith(file, ".lcf") && occursin(basename, split(file, "-")[1])
                push!(lcf_files, joinpath(root, file))
                push!(lcf_paths, root)
            end
        end
    end


    t_start = []
    t_end = []
    for f ∈ lcf_files
        println("Processing $(f)")
        df = DataFrame(CSV.File(f, header=[:time, :roll, :pitch, :heading, :longitude, :latitude, :altitude, :placeHolder1, :placeHolder2, :placeHolder3, :placeHolder4]))


        # start_time = df.time[1]  # note this is gps time (i.e. the number of seconds since 0h 6-Jan-1980, UTC)
        # start_time += 315964800.00  # convert to linux epoch time to make working with DateTime object easier
        start_time = gpsToUTC(df.time[1], 2019)


        df.time .= (df.time .- df.time[1])

        t_s = start_time

        sec = round(df.time[end])
        ms = (df.time[end] - sec)
        # convert to DateTime object
        sec = Second(Int(sec))
        ms = Millisecond(Int(round(1000*ms, digits=0)))
        Δt = sec + ms
        t_e = start_time + Δt

        push!(t_start, t_s)
        push!(t_end, t_e)
        println("  $(t_e-t_s)")

    end

    master_lcf_df = DataFrame()
    master_lcf_df."paths" = lcf_paths
    master_lcf_df."files" = lcf_files
    master_lcf_df."tstart" = t_start
    master_lcf_df."tend" = t_end


    sort!(master_lcf_df, :tstart)


    CSV.write(joinpath(folder, "master_lcf.csv"), master_lcf_df)

    return master_lcf_df
end
