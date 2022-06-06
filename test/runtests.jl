using georectification
using Test

@testset "georectification.jl" begin


    # set up paths to test files
    testPath = "/media/john/HSDATA/raw/11-23/Scotty_1-1"
    @test ispath(testPath) == true

    testBilPath = joinpath(testPath, "Scotty_1_Pika_XC2_1-radiance.bil")
    @test isfile(testBilPath) == true

    testBilHdrPath = joinpath(testPath, "Scotty_1_Pika_XC2_1-radiance.bil.hdr")
    @test isfile(testBilHdrPath) == true

    testTimesPath = joinpath(testPath, "Scotty_1_Pika_XC2_1.bil.times")
    @test isfile(testTimesPath) == true

    testSpecPath = joinpath(testPath, "Scotty_1_downwelling_1_pre.spec")
    @test isfile(testSpecPath)

    testSpecHdrPath = joinpath(testPath, "Scotty_1_downwelling_1_pre.spec.hdr")
    @test isfile(testSpecHdrPath)

    lcf_path = joinpath(testPath, "Scotty_1_Pika_XC2_1.lcf")
    @test isfile(lcf_path)

    # grab the default values from config
    include("../../mintsRobotTeam/config.jl")

    # overwrite the calibrationPath value
    calibrationPath = "../../mintsRobotTeam/calibration/"
    @test ispath(calibrationPath)

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

    println(describe(df))


end
