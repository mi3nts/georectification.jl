using georectification
using Images, FileIO
using Dates, TimeZones
using Plots
using DataFrames, CSV
using ProgressMeter
using BenchmarkTools



include("../../mintsRobotTeam/config.jl")

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


tiff_temp = TIFFtoTemp(tiff)

p1 = plot(jpg, title="visible", aspect_ratio=:equal)
p2 = heatmap(tiff_temp, title="thermal", aspect_ratio=:equal, yflip=true)
plot(p1, p2)



D = exiftool(tiffPath)
date = getExifDate(D)
date2 = getExifDate(tiffPath)



df = getFlirTimes(basepath)

names(df)

df.times_central[1]
df.times_utc[1]

# @btime getFlirTimes(basepath)


flirpath = "/media/john/HSDATA/FLIR"
# generateFlirSummaries(flirpath)



# test_csv = joinpath(flirpath, "20201123_112658.csv")
test_csv = joinpath(flirpath, "20201209_110416.csv")
df = CSV.File(test_csv) |> DataFrame

names(df)
df.times_utc[1]
df.times_utc[end]



df_dye1 = masterLCF("/media/john/HSDATA/raw/12-09", "Dye_1")
df_dye2 = masterLCF("/media/john/HSDATA/raw/12-09", "Dye_2")
df_nodye1 = masterLCF("/media/john/HSDATA/raw/12-09", "NoDye_1")
df_nodye2 = masterLCF("/media/john/HSDATA/raw/12-09", "NoDye_2")

size(df_dye1)
size(df_dye2)
size(df_nodye1)
size(df_nodye2)

println("FLIR times: $(df.times_utc[1]) to $(df.times_utc[end])")
println("NoDye_1 times: $(df_nodye1.tstart[1]) to $(df_nodye1.tend[end])")
println("NoDye_2 times: $(df_nodye2.tstart[1]) to $(df_nodye2.tend[end])")
println("Dye_1 times: $(df_dye1.tstart[1]) to $(df_dye1.tend[end])")
println("Dye_2 times: $(df_dye2.tstart[1]) to $(df_dye2.tend[end])")

size(df)
size(df_nodye2)



df_sub = DuringFlight(df, df_nodye2)

names(df_sub)
names(df_nodye2)

matchLCFfile!(df_sub, df_nodye2)

df_sub.tiffpath[1]

df_sub.lcf_path[1:20]

df_nodye2[1:10, [:tstart, :tend]]
df_nodye2.tend[1]
df_nodye2.tstart[2]


df_sub.tiffpath[1]

res = exiftool(df_sub.tiffpath[1])
for key ∈ keys(res)
    println(key)
end

res["Yaw"]
res["Roll"]
res["Pitch"]

res["File Modification Date/Time"]
res["File Access Date/Time"]
res["File Inode Change Date/Time"]

df_sub.times_utc[1]
df_sub.times_basename_utc[1]
df_nodye2.tstart[1]

df_sub.lcf_path[1:20]

img_vis, img_therm, visCoords, visTimes, thermCoords, thermTimes = georectifyFLIR(df_sub.tiffpath[2],
                                                                                  df_sub.jpgpath[2],
                                                                                  df_sub.lcf_path_basename[2],
                                                                                  df_sub.times_basename_utc[2],
                                                                                  location_data["scotty"]["z"],
                                                                                  10.0,
                                                                                  )


rgb_tiff, longitudes, latitudes = getBackgroundTile(location_data["scotty"]["w"],
                                                    location_data["scotty"]["n"],
                                                    location_data["scotty"]["e"], location_data["scotty"]["s"],
                                                    "scotty_background",
                                                    )


p1 = plot_background(rgb_tiff,
                     longitudes,
                     latitudes;
                     size=(3*600, 3*400),
                     xlim=(location_data["scotty"]["w"],location_data["scotty"]["e"]),
                     ylim=(location_data["scotty"]["s"],location_data["scotty"]["n"]),
                     left_margin = 10*Plots.mm,
                     guidefontsize=18,
                     tickfontsize=13,
                     tick_direction=:out,
                     )


size(thermCoords)  # lat, lon, alt
size(img_vis)

plot!(p1, visCoords[2,1:10:end,1:10:end], visCoords[1,1:10:end,1:10:end], seriestype=:scatter, color = img_vis[1:10:end,1:10:end], alpha=0.7, ms = 1, markerstrokewidth=0, label="")
plot!(p1, thermCoords[2,:,:], thermCoords[1, :, :], seriestype=:scatter,  zcolor=img_therm, clims=(0,30), ms=1, markerstrokewidth=0, label="")

# savefig("12__utc_time_basename__3π_4_heading_pitchAdjust_10.png")
# savefig("utc_time_basename__3π_4_heading.png")
# savefig("utc_time__3π_4_heading.png")
# savefig("utc_time__π_4_heading.png")
# savefig("utc_time__π_2_heading.png")
# savefig("utc_time.png")

display(p1)


names(df_sub)
imu_df, start_time = getIMUdata(df_sub.lcf_path_basename[2])
plot(imu_df.x, imu_df.y, imu_df.z, line_z=imu_df.time, xlabel="x", ylabel="y", zlabel="z", label="")

t_cap = Dates.value(df_sub.times_basename_utc[3] - start_time)/1000

pa = plot(imu_df.time, imu_df.heading_correct)
vline!([t_cap])

pb = plot(imu_df.time, imu_df.roll)
vline!([t_cap])

pc = plot(imu_df.time, imu_df.pitch)
vline!([t_cap])

plot(pa, pb, pc, layout=(3,1))

t_cap
