using georectification
using DataFrames


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
                false,
                6
                )

println(typeof(df))

println(describe(df))


