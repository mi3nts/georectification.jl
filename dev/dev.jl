using georectification
using DataFrames
using Plots
using Interpolations
using Images

# set up paths to test files
testPath = "/media/john/HSDATA/raw/11-23/Scotty_1-1"
testBilPath = joinpath(testPath, "Scotty_1_Pika_XC2_1-radiance.bil")
testBilHdrPath = joinpath(testPath, "Scotty_1_Pika_XC2_1-radiance.bil.hdr")
testTimesPath = joinpath(testPath, "Scotty_1_Pika_XC2_1.bil.times")
testSpecPath = joinpath(testPath, "Scotty_1_downwelling_1_pre.spec")
testSpecHdrPath = joinpath(testPath, "Scotty_1_downwelling_1_pre.spec.hdr")
lcf_path = joinpath(testPath, "Scotty_1_Pika_XC2_1.lcf")
# grab the default values from config
include("../../mintsRobotTeam/config.jl")

# overwrite the calibrationPath value
calibrationPath = "../../mintsRobotTeam/calibration/"


df = readToDataFrame(testBilPath, testBilHdrPath, testTimesPath)

generateReflectance!(df, testSpecPath, testSpecHdrPath, calibrationPath, wavelengths)


println("attempting georectification")

df = georectify(testBilPath,
                testBilHdrPath,
                testTimesPath,
                testSpecPath,
                testSpecHdrPath,
                calibrationPath,
                lcf_path,
                wavelengths,
                location_data["scotty"]["z"],
                θ_view,
                true,
                6
                )



cols = ["λ_$(i)" for i ∈ 1:length(wavelengths)]
refHSI = Array(df[:, cols])


sRGB = HSI_2_RGB(wavelengths, refHSI, 65)

RGB_final = RGB.(sRGB[1, :], sRGB[2,:], sRGB[3,:])

rgb_tiff, longitudes, latitudes = getBackgroundTile(location_data["scotty"]["w"],
                                                    location_data["scotty"]["n"],
                                                    location_data["scotty"]["e"], location_data["scotty"]["s"],
                                                    "scotty_background"
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

plot!(df.longitude, df.latitude, seriestype=:scatter, color = RGB_final, alpha=0.7, ms = 1, markerstrokewidth=0, label="")

savefig("./no_correction.png")


