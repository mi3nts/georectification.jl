module georectification

using CSV
using DataFrames
using Geodesy
using Rotations
using StaticArrays
using DataInterpolations
# using Interpolations
using Statistics
using Images
using Plots
using Dates
using ProgressMeter

include("HSI.jl")
include("IMU.jl")
include("satplot.jl")
include("HSI2RGB.jl")
include("FLIR.jl")

export getHdrFile
export getTimes
export readToDataFrame
export getHdrFile

export leap_count
export gpsToUTC
export getIMUdata
export Rotate
export masterLCF

export generateReflectance!
export generateDerivedMetrics!
export getHSIFlightData
export imagecoords
export imagecoordsFlipped
export Coords!
export ilat_ilon!
export georectify

export getBackgroundTile
export plot_background

export Spec_2_RGB
export HSI_2_RGB

# FLIR stuff
export TIFFtoTemp
export exiftool
export getExifDate
export getFlirTimes
export generateFlirSummaries
export DuringFlight
export matchLCFfile!
export georectifyFLIR




"""
    generateReflectance!(df::DataFrame, specPath::String, specHdrPath::String, calibrationPath::String, wavelengths::Array{Float64, 1})

Given a dataframe `df` containing the radiance data, generate reflectance data using collected irradiance spectrum and Resonon-supplied calibration data.

# NOTE
The conversion assumes a Lambertian surface. This is clearly *not* the case for water that isn't perfectly still.

# NOTE Also add the corrected radiance spectrum values as additional features.
"""
function generateReflectance!(df::DataFrame, specPath::String, specHdrPath::String, calibrationPath::String, wavelengths::Array{Float64, 1})

    # Grab the calibration gain and offset HSIs for the
    # downwelling sensor
    #----------------------------------------------------------------------
    calibrationGainPath = joinpath(calibrationPath, "gain.spec")
    calibrationGainHdrPath = joinpath(calibrationPath, "gain.spec.hdr")

    calibrationOffsetPath = joinpath(calibrationPath, "offset.spec")
    calibrationOffsetHdrPath = joinpath(calibrationPath, "offset.spec.hdr")

    gain_df = readToDataFrame(calibrationGainPath, calibrationGainHdrPath)
    gain_hdr = getHdrFile(calibrationGainHdrPath)


    offset_df = readToDataFrame(calibrationOffsetPath, calibrationOffsetHdrPath)
    offset_hdr = getHdrFile(calibrationOffsetHdrPath)

    spec_df = readToDataFrame(specPath, specHdrPath)
    spec_hdr = getHdrFile(specHdrPath)


    # calculate shutter differences
    #----------------------------------------------------------------------
    cal_shutter = gain_hdr["shutter"]
    spec_shutter = spec_hdr["shutter"]
    gain_factor = cal_shutter/spec_shutter # should automatically convert to type float

    # for easier looping
    nλ_gain = size(gain_hdr["wavelength"], 1)
    nλ_offset = size(offset_hdr["wavelength"], 1)
    nλ_spec = size(spec_hdr["wavelength"], 1)


    # produce the correction frames
    #----------------------------------------------------------------------
    adjusted_gain = gain_factor .* vec(Matrix(gain_df[!, 1:nλ_gain]))
    frame = vec(Matrix(spec_df[!, 1:nλ_spec])) .- vec(Matrix(offset_df[!, 1:nλ_offset]))

    # calculate correction
    #----------------------------------------------------------------------
    correction = π .* adjusted_gain .* frame  # need to figure out what this π is for... I suspect this is from the cosine corrector.
    clamp!(correction, 0, typemax(eltype(frame)))


    # interpolate the correction to match the datacube's wavelengths
    #----------------------------------------------------------------------
    interp = CubicSpline(correction[:], spec_hdr["wavelength"])
    adjustedSpec = interp.(wavelengths)


    # setup outgoing cube
    # Note: R = π(sensor radiance / downwelling irradiance) as per spectronon docs
    # http://docs.resonon.com/spectronon/SpectrononUserManual/html/plugin_docs.html#reflectance-from-radiance-data-and-downwelling-irradiance-spectrum
    #----------------------------------------------------------------------
    npixels = size(df)[1]
    for i ∈ 1:462
        @inbounds df[!, "λ_$(i)"] = clamp.(π .* df[!, "λ_$(i)_rad"] ./ adjustedSpec[i], 0.0, typemax(Float64))
    end

    # # Also copy over the "correction" values i.e. the calibrated spectrum from the downwelling irradiance spectrometer with cosine corrector
    # for i ∈ 1:nλ_spec
    #     df[!, "λ_$(i)_downwelling"] = correction[i] .* ones(nrow(df))
    # end

    return  correction
end




"""
    generateDerivedMetrics!(df::DataFrame, λs::Array{Float64})

Given a dataframe `df` with HSI data, compute a variety of derived wavelength indices such as the NDVI.
"""
function generateDerivedMetrics!(df::DataFrame, λs::Array{Float64})
    # compute total radiance at each point
    rads = ["λ_$(i)_rad" for i ∈ 1:462]
    refs = ["λ_$(i)" for i ∈ 1:462]


    # mNDWI
    df.mNDWI = -1 .* ones(size(df)[1])
    df.rad_mNDWI  = -1 .* ones(size(df)[1])

    jgreen = findfirst(λs .> 550)
    jswir = findfirst(λs .> 770)
    green = df[!, "λ_$(jgreen)"]
    swir = df[!, "λ_$(jswir)"]
    numer = (green .- swir)
    denom = (green .+ swir)
    df.mNDWI[denom .> 0] = (numer ./ denom)[denom .> 0]

    green = df[!, "λ_$(jgreen)_rad"]
    swir = df[!, "λ_$(jswir)_rad"]
    numer = (green .- swir)
    denom = (green .+ swir)
    df.rad_mNDWI[denom .> 0] = (numer ./ denom)[denom .> 0]

    # # NDVI "normalized difference vegetative index ∈ [-1, 1]"
    df.NDVI = -2 .* ones(size(df)[1])
    df.rad_NDVI  = -2 .* ones(size(df)[1])

    jir = findfirst(λs .> 800)
    jred = findfirst(λs .> 680)

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]

    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)

    df.NDVI[denom .> 0] .= (numer ./ denom)[denom .> 0]


    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)
    df.rad_NDVI[denom .> 0] .= (numer ./ denom)[denom .> 0]


    # # SR "simple ratio ∈ [0, 30]"
    df.SR = -1 .* ones(size(df)[1])
    df.rad_SR = -1 .* ones(size(df)[1])

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    numer = ir_band
    denom = red_band
    df.SR[denom .> 0] .= (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    numer = ir_band
    denom = red_band
    df.rad_SR[denom .> 0] .= (numer ./ denom)[denom .> 0]


    # # EVI "enhanced vegetative index ∈ [-1, 1]"
    jblue = findfirst(λs .> 450)
    df.EVI = -2 .* ones(size(df)[1])
    df.rad_EVI = -2 .* ones(size(df)[1])

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    blue_band = df[!, "λ_$(jblue)"]

    numer = 2.5 .* (ir_band .- red_band)
    denom = ir_band .+ 6 .* red_band .- 7.5 .* blue_band
    df.EVI[denom .> 0] .= (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    blue_band = df[!, "λ_$(jblue)_rad"]

    numer = 2.5 .* (ir_band .- red_band)
    denom = ir_band .+ 6 .* red_band .- 7.5 .* blue_band
    df.rad_EVI[denom .> 0] .= (numer ./ denom)[denom .> 0]



    # # AVRI "Atmospherical Resistant Vegitative Indes"
    df.AVRI = -2 .* ones(size(df)[1])
    df.rad_AVRI = -2 .* ones(size(df)[1])

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    blue_band = df[!, "λ_$(jblue)"]
    numer = (ir_band .- 2 .* red_band .+ blue_band)
    denom = (ir_band .+ 2 .* red_band .- blue_band)
    df.AVRI[denom .> 0] .=  (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    blue_band = df[!, "λ_$(jblue)_rad"]
    numer = (ir_band .- 2 .* red_band .+ blue_band)
    denom = (ir_band .+ 2 .* red_band .- blue_band)
    df.rad_AVRI[denom .> 0] .=  (numer ./ denom)[denom .> 0]



    # # NDVI_705 "Red Edge Normalized Difference Vegetation Index"S
    df.NDVI_705 = -2 .* ones(size(df)[1])
    df.rad_NDVI_705  = -2 .* ones(size(df)[1])

    jir = findfirst(λs .> 750)
    jred = findfirst(λs .> 705)
    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)
    df.NDVI_705[denom .> 0] .= (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)
    df.rad_NDVI_705[denom .> 0] .= (numer ./ denom)[denom .> 0]




    # # MSR_705 "Modified Red Edge Simple Ratio Index"
    df.MSR_705 = -1 .* ones(size(df)[1])
    df.rad_MSR_705 = -1 .* ones(size(df)[1])

    jblue = findfirst(λs .> 445)
    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    blue_band = df[!, "λ_$(jblue)"]
    numer = (ir_band .- blue_band)
    denom = (red_band .- blue_band)
    df.MSR_705[denom .> 0] .= (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    blue_band = df[!, "λ_$(jblue)_rad"]
    numer = (ir_band .- blue_band)
    denom = (red_band .- blue_band)
    df.rad_MSR_705[denom .> 0] .= (numer ./ denom)[denom .> 0]



    # # MNDVI "modified red edge normalized vegetation index"
    df.MNDVI = -2 .* ones(size(df)[1])
    df.rad_MNDVI = -2 .* ones(size(df)[1])

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    blue_band = df[!, "λ_$(jblue)"]
    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band .- 2 .* blue_band)
    df.MNDVI[denom .> 0] .=  (numer ./ denom)[denom .> 0]


    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    blue_band = df[!, "λ_$(jblue)_rad"]
    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band .- 2 .* blue_band)
    df.rad_MNDVI[denom .> 0] .=  (numer ./ denom)[denom .> 0]


    # # VOG1 "vogelmann red edge index"
    df.VOG1 = -1 .* ones(size(df)[1])
    df.rad_VOG1 = -1 .* ones(size(df)[1])
    jir = findfirst(λs .> 740)
    jred = findfirst(λs .> 720)

    ir_band = df[!, "λ_$(jir)"]
    red_band = df[!, "λ_$(jred)"]
    numer = (ir_band)
    denom = (red_band)
    df.VOG1[denom .> 0] .= (numer ./ denom)[denom .> 0]

    ir_band = df[!, "λ_$(jir)_rad"]
    red_band = df[!, "λ_$(jred)_rad"]
    numer = (ir_band)
    denom = (red_band)
    df.rad_VOG1[denom .> 0] .= (numer ./ denom)[denom .> 0]


    # # VOG2 "vogelmann red edge index 2"
    df.VOG2 = -1 .* ones(size(df)[1])
    df.rad_VOG2 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 734)
    j2 = findfirst(λs .> 747)
    j3 = findfirst(λs .> 715)
    j4 = findfirst(λs .> 726)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    band4 = df[!, "λ_$(j4)"]
    numer = (band1 .- band2)
    denom = (band3 .+ band4)
    df.VOG2[denom .> 0] .=  (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    band4 = df[!, "λ_$(j4)_rad"]
    numer = (band1 .- band2)
    denom = (band3 .+ band4)
    df.rad_VOG2[denom .> 0] .=  (numer ./ denom)[denom .> 0]



    # # VOG3 "vogelmann red edge index 3 ∈ [0, 20]"
    df.VOG3 = -1 .* ones(size(df)[1])
    df.rad_VOG3 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 734)
    j2 = findfirst(λs .> 747)
    j3 = findfirst(λs .> 715)
    j4 = findfirst(λs .> 720)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    band4 = df[!, "λ_$(j4)"]
    numer = (band1 .- band2)
    denom = (band3 .+ band4)
    df.VOG3[denom .> 0] .=  (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    band4 = df[!, "λ_$(j4)_rad"]
    numer = (band1 .- band2)
    denom = (band3 .+ band4)
    df.rad_VOG3[denom .> 0] .=  (numer ./ denom)[denom .> 0]



    # # PRI "photochemical reflectance index" ∈[-1, 1]
    df.PRI = -2 .* ones(size(df)[1])
    df.rad_PRI = -2 .* ones(size(df)[1])
    j1 = findfirst(λs .> 531)
    j2 = findfirst(λs .> 570)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    numer = (band1 .- band2)
    denom = (band1 .+ band2)
    df.PRI[denom .> 0] .= (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    numer = (band1 .- band2)
    denom = (band1 .+ band2)
    df.rad_PRI[denom .> 0] .= (numer ./ denom)[denom .> 0]


    # # SIPI "structure intensive pigment index" ∈[0, 2]
    df.SIPI = -1 .* ones(size(df)[1])
    df.rad_SIPI = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 800)
    j2 = findfirst(λs .> 445)
    j3 = findfirst(λs .> 680)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    numer = (band1 .- band2)
    denom = (band1 .+ band3)
    df.SIPI[denom .> 0] .= (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    numer = (band1 .- band2)
    denom = (band1 .+ band3)
    df.rad_SIPI[denom .> 0] .= (numer ./ denom)[denom .> 0]



    # # PSRI "Plant Senescence Reflectance Index"
    df.PSRI = -2 .* ones(size(df)[1])
    df.rad_PSRI = -2 .* ones(size(df)[1])
    j1 = findfirst(λs .> 680)
    j2 = findfirst(λs .> 500)
    j3 = findfirst(λs .> 750)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    numer = (band1 .- band2)
    denom = band3
    df.PSRI[denom .> 0] .= (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    numer = (band1 .- band2)
    denom = band3
    df.rad_PSRI[denom .> 0] .= (numer ./ denom)[denom .> 0]


    # # CRI1 "carotenoid reflectance index"
    df.CRI1 = -1 .* ones(size(df)[1])
    df.rad_CRI1 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 510)
    j2 = findfirst(λs .> 550)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    df.CRI1[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    df.rad_CRI1[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]


    # # CRI2 "carotenoid reflectance index 2"
    df.CRI2 = -1 .* ones(size(df)[1])
    df.rad_CRI2 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 510)
    j2 = findfirst(λs .> 700)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[:, "λ_$(j2)"]
    df.CRI2[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[:, "λ_$(j2)_rad"]
    df.rad_CRI2[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]


    # # ARI1 "anthocyanin reflectance index"
    df.ARI1 = -1 .* ones(size(df)[1])
    df.rad_ARI1 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 550)
    j2 = findfirst(λs .> 700)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    df.ARI1[(band1 .> 0) .& (band2 .> 0)] = ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    df.rad_ARI1[(band1 .> 0) .& (band2 .> 0)] = ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]




    # # ARI2 "anthocyanin reflectance index 2"
    df.ARI2 = -1 .* ones(size(df)[1])
    df.rad_ARI2 = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 550)
    j2 = findfirst(λs .> 700)
    j3 = findfirst(λs .> 800)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[:, "λ_$(j2)"]
    band3 = df[:, "λ_$(j3)"]
    df.ARI2[(band1 .> 0) .& (band2 .> 0)] .= (band3.*((1 ./ band1) .- (1 ./ band2)))[(band1 .> 0) .& (band2 .> 0)]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[:, "λ_$(j2)_rad"]
    band3 = df[:, "λ_$(j3)_rad"]
    df.rad_ARI2[(band1 .> 0) .& (band2 .> 0)] .= (band3.*((1 ./ band1) .- (1 ./ band2)))[(band1 .> 0) .& (band2 .> 0)]


    # # WBI "water band index"
    df.WBI = -1 .* ones(size(df)[1])
    df.rad_WBI = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 900)
    j2 = findfirst(λs .> 970)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    numer = band1
    denom = band2
    df.WBI[denom .> 0] .= (numer ./ denom)[denom .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    numer = band1
    denom = band2
    df.rad_WBI[denom .> 0] .= (numer ./ denom)[denom .> 0]



    # # MCRI "Modified Chlorophyll Absorption Reflectance Index"
    df.MCRI = -1 .* ones(size(df)[1])
    df.rad_MCRI  = -1 .* ones(size(df)[1])
    j1 = findfirst(λs .> 550)
    j2 = findfirst(λs .> 670)
    j3 = findfirst(λs .> 701)
    j4 = findfirst(λs .> 780)

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    band4 = df[!, "λ_$(j4)"]
    df.MCRI[band2 .> 0] .= (((band3 .-band2) .- 0.2 .* (band3 .- band1)) .* (band3 ./ band2))[band2 .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    band4 = df[!, "λ_$(j4)_rad"]
    df.rad_MCRI[band2 .> 0] .= (((band3 .-band2) .- 0.2 .* (band3 .- band1)) .* (band3 ./ band2))[band2 .> 0]



    # # TCARI "transformed chlorophyll absorption reflectance index"
    df.TCARI = -1 .* ones(size(df)[1])
    df.rad_TCARI  = -1 .* ones(size(df)[1])

    band1 = df[!, "λ_$(j1)"]
    band2 = df[!, "λ_$(j2)"]
    band3 = df[!, "λ_$(j3)"]
    df.TCARI[band2 .> 0] .= (3 .* ((band3 .- band2) .- 0.2 .* (band3 .- band1) .* (band3 ./band2)))[band2 .> 0]

    band1 = df[!, "λ_$(j1)_rad"]
    band2 = df[!, "λ_$(j2)_rad"]
    band3 = df[!, "λ_$(j3)_rad"]
    df.rad_TCARI[band2 .> 0] .= (3 .* ((band3 .- band2) .- 0.2 .* (band3 .- band1) .* (band3 ./band2)))[band2 .> 0]

end






"""
    getHSIFlightData(df::DataFrame, lcf_path::String)

Given IMU data found at `lcf_path`, compute the position and orientation of the Aerial Vehicle at time of capture for each pixel. Return's a dataframe and the corresponding start time of the acquisition reported in UTC.
"""
function getHSIFlightData(df::DataFrame, lcf_path::String)
    imu_df, start_time = getIMUdata(lcf_path)  # note start time is in unix epoch time

    # for each of the variables, set up a cubic spline interpolation
    α_interp = CubicSpline(imu_df.heading_correct, imu_df.time)
    β_interp = CubicSpline(imu_df.pitch, imu_df.time)
    γ_interp = CubicSpline(imu_df.roll, imu_df.time)
    x_interp = CubicSpline(imu_df.x, imu_df.time)
    y_interp = CubicSpline(imu_df.y, imu_df.time)
    z_interp = CubicSpline(imu_df.z, imu_df.time)


    # generate new dataframe with interpolated results
    times = unique(df.times)

    α = α_interp.(times)
    β = β_interp.(times)
    γ = γ_interp.(times)
    x = x_interp.(times)
    y = y_interp.(times)
    z = z_interp.(times)

    # estimate the  forward pixel resolution
    ΔLs = sqrt.((x[2:end].-x[1:end-1]).^2 .+ (y[2:end] .- y[1:end-1]).^2) # copy first value to get correct length
    ΔLs = vcat(ΔLs[1], ΔLs...)

    isnorth = [imu_df[1, :isnorth] for i∈1:length(times)]
    zone = [imu_df[1, :zone] for i∈1:length(times)]

    hsi_df = DataFrame()
    hsi_df.times = times
    hsi_df.heading = α
    hsi_df.pitch = β
    hsi_df.roll = γ
    hsi_df.x = x
    hsi_df.y = y
    hsi_df.z = z
    hsi_df.isnorth = isnorth
    hsi_df.zone = zone
    hsi_df.forward_res = ΔLs

    return hsi_df, start_time
end





"""
    imagecoords(i, j, N, f)

Given pixel indices `i` and `j`, compute the coordinates of a HSI pixel in the image coordinate system. The height is the focal length `f`, and the pixels are assumed to lie along the y-axis.
"""
function imagecoords(i, j, N, f)
    if i == 1
        return 0.0
    elseif i == 2
        return (N-1)/2 - (j-1)
    else
        return f
    end
end

"""
    imagecoordsFlipped(i, j, N, f)

Same as `imagecoord()` but with the y axis flipped. This is useful in case the camera settings for direction of flight are backwards.
"""
function imagecoordsFlipped(i, j, N, f)
    if i == 1
        return 0.0
    elseif i == 2
        return -(N-1)/2 + (j-1)
    else
        return f
    end
end





"""
    Coords!(df::DataFrame,
            lcf_path::String,
            z_ground::Float64,
            θ::Float64,
            flipped::Bool)

Compute the ground coordinates for each pixel of an HSI `df`. Update the original `df` to include
- latitude
- longitude
- altitude
- pixeltimes
- roll
- pitch
- heading
and return `start_time` in UTC corresponding to `pixeltime = 0.0`.
"""
function Coords!(df::DataFrame,
                lcf_path::String,
                z_ground::Float64,
                θ::Float64,
                flipped::Bool)
    lines = size(unique(df.row_index))[1]
    samples = Int(size(df)[1]/lines)

    # collect interpolated flight data
    hsi_df, start_time = getHSIFlightData(df, lcf_path)
    # set up return arrays

    ## we may be able to get rid of T_n_E below by swapping samples and lines.
    ## the row index is effectively the y-value in image coordinates so
    ## we may have these swapped... 
    pixelCoords = Array{Float64}(undef, 3, lines, samples)  # val, row, col
    # pixelRes = Array{Float64}(undef, lines, samples) # estimates of pixel resolution
    pixelTimes = Array{Float64}(undef, lines, samples) # time each frame was captured
    # pixelTimes = Array{DateTime}(undef, lines, samples) # time each frame was captured

    viewingGeom = Array{Float64}(undef, 3, lines, samples) # store roll, pitch, yaw

    # compute pixel locations in
    N = samples
    f = ((N-1)/2)/tand(θ/2) # compute the focal length in units of "pixels". tand is tangent in degrees

    if flipped
        # rs_pixel_sensor = SMatrix{3,N}([imagecoordsFlipped(i,j,N,f) for i∈1:3, j∈1:N])
        rs_pixel_sensor = [imagecoordsFlipped(i,j,N,f) for i∈1:3, j∈1:N]
    else
        # rs_pixel_sensor = SMatrix{3,N}([imagecoords(i,j,N,f) for i∈1:3, j∈1:N])
        rs_pixel_sensor = [imagecoords(i,j,N,f) for i∈1:3, j∈1:N]
    end

    i = 1
    for frame ∈ eachrow(hsi_df)
        s = (frame.z-z_ground)/f  # scale factor

        rs_object_utm = [frame.x; frame.y; frame.z] .+ s .* T_n_E*Rotate(frame.heading, frame.pitch, frame.roll)*rs_pixel_sensor

        # convert back to latitude, longitude, altitude
        rs_object_lla = [LLAfromUTMZ(wgs84)(UTMZ(rs_object_utm[:,j]..., frame.zone, frame.isnorth)) for j∈1:N]
        @inbounds pixelCoords[1, i, :] .= [lla.lat for lla ∈ rs_object_lla]
        @inbounds pixelCoords[2, i, :] .= [lla.lon for lla ∈ rs_object_lla]
        @inbounds pixelCoords[3, i, :] .= [lla.alt for lla ∈ rs_object_lla]


        @inbounds pixelTimes[i, :] .= [frame.times for i ∈ 1:N]
        @inbounds viewingGeom[1, i, :] .= [frame.roll for i∈1:N]
        @inbounds viewingGeom[2, i, :] .= [frame.pitch for i∈1:N]
        @inbounds viewingGeom[3, i, :] .= [frame.heading for i∈1:N]

        # Turn off resolution for speed up

        # compute distance between pixels in a sample
        # Δxs_sample = [rs_object_utm[1,j]-rs_object_utm[1, j-1] for j∈2:N]
        # Δys_sample = [rs_object_utm[2,j]-rs_object_utm[2, j-1] for j∈2:N]
        # ΔL_sample = sqrt.(Δxs_sample.^2 .+ Δys_sample.^2)

        # pixelRes[i, 1] = maximum([ΔL_sample[1], frame.forward_res])
        # pixelRes[i, 2:end] .= [maximum([dl, frame.forward_res]) for dl ∈ ΔL_sample] # for each pixel, save maximum possible res between sample and line


        i += 1
    end

    # boundary = cat(pixelCoords[1:2,:,1], pixelCoords[1:2,end:-1:1,end], pixelCoords[1:2,1,1], dims=2)


    # update df with results
    # lat, lon, alt
    pxcoords = reshape(PermutedDimsArray(pixelCoords, (1, 3, 2)), (3, lines*samples))
#    pxres = reshape(pixelRes', (lines*samples))
    pxtimes = reshape(pixelTimes', (lines*samples))
    # roll, pitch, heading
    viewgeom = reshape(PermutedDimsArray(viewingGeom, (1, 3, 2)), (3, lines*samples))

    df[!, :latitude] = vec(pxcoords[1,:])
    df[!, :longitude] = vec(pxcoords[2,:])
    df[!, :altitude] = vec(pxcoords[3,:])
    # df[!, :resolution] = pxres
    df[!, :pixeltimes] = pxtimes
    df[!, :roll] = vec(viewgeom[1,:])
    df[!, :pitch] = vec(viewgeom[2,:])
    df[!, :heading] = vec(viewgeom[3,:])

    return start_time
end





"""
 ilat_ilon!(df::DataFrame, ndigits)

Update a dataframe to include the rounded coordinates `ilat` and `ilon` with `ndigits` of precision.
"""
function ilat_ilon!(df::DataFrame, ndigits)
    df[!, :ilat] = round.(df[!,:latitude], digits=ndigits)
    df[!, :ilon] = round.(df[!,:longitude], digits=ndigits)
end


"""
    ilat_ilon!(df::DataFrame)

Dispatch of `ilat_ilon()` to use default value of `6` for `ndigits`.
"""
function ilat_ilon!(df::DataFrame)
    ilat_ilon!(df, 6)
end




"""
    georectify(bilpath::String,
               bilhdrpath::String,
               timespath::String,
               specpath::String,
               spechdrpath::String,
               calibrationpath::String,
               lcfpath::String,
               λs::Array{Float64},
               z_ground::Float64,
               θ_view::Float64,
               isFlipped::Bool,
               ndigits::Int,
               )

Generate a georectified dataframe including radiance, reflectance, derived metrics, `utc_time`, and position data.
"""
function georectify(bilpath::String,
                    bilhdrpath::String,
                    timespath::String,
                    specpath::String,
                    spechdrpath::String,
                    calibrationpath::String,
                    lcfpath::String,
                    λs::Array{Float64},
                    z_ground::Float64,
                    θ_view::Float64,
                    isFlipped::Bool,
                    ndigits::Int,
                    )
    println("\treading in the BIL file")
    df = readToDataFrame(bilpath, bilhdrpath, timespath)
    println("\tgenerating reflectance data")
    irrad_df = generateReflectance!(df,
                                    specpath,
                                    spechdrpath,
                                    calibrationpath,
                                    λs
                                    )

    println("\tgenerating derived metrics")
    generateDerivedMetrics!(df, λs)

    start_time = Coords!(df,
                         lcfpath,
                         z_ground,
                         θ_view,
                         isFlipped
                         )


    println("\tgenerating ilat and ilon")
    ilat_ilon!(df, ndigits)

    println("\tgrouping by ilat and ilon")
    gdf = groupby(df, [:ilat, :ilon])

    println("\trecombining via average")
    res = combine(gdf, Not([:ilat, :ilon]) .=> mean; renamecols=false)


    res.utc_times = Array{DateTime}(undef, size(res, 1)) # time each frame was captured

    for i ∈ 1:size(res,1)
        sec = round(res.pixeltimes[i])
        ms = (res.pixeltimes[i] - sec)
        sec = Second(Int(sec))
        ms = Millisecond(Int(round(1000*ms, digits=0)))
        Δt = sec + ms

        res.utc_times[i] = start_time + Δt
    end

    # add new column with utc_time



    # hocus pocus for freeing locked memory
    # df = nothing
    # gdf = nothing

    # GC.gc()
    # ccall(:malloc_trim, Cvoid, (Cint,), 0)
    # GC.gc()

    return res, irrad_df
end







"""
    function georectifyFLIR(thermal_path::String, visible_path::String, lcf_path::String, capture_time::DateTime, z_ground::Float64)
"""
function georectifyFLIR(thermal_path::String, visible_path::String,  lcf_path::String, capture_time::DateTime, z_ground::Float64, pitch_adjustment::Float64)
    # df is dataframe with start and end times for each lcf file

    # read in the image
    img_vis = load(visible_path)
    img_therm = load(thermal_path)
    img_therm = TIFFtoTemp(img_therm) # convert to temperature units

    Nx_vis, Ny_vis = size(img_vis)
    Nx_therm, Ny_therm = size(img_therm)

    θx_therm = 20
    θy_therm = 25
    θx_vis = 45
    θy_vis = 56

    imu_df, start_time = getIMUdata(lcf_path)

    t_cap = Dates.value(capture_time - start_time)/1000

    # interpolate the lcf data to match the FLIR time
    α_interp = CubicSpline(imu_df.heading_correct, imu_df.time)
    β_interp = CubicSpline(imu_df.pitch, imu_df.time)
    γ_interp = CubicSpline(imu_df.roll, imu_df.time)
    x_interp = CubicSpline(imu_df.x, imu_df.time)
    y_interp = CubicSpline(imu_df.y, imu_df.time)
    z_interp = CubicSpline(imu_df.z, imu_df.time)

    α = α_interp(t_cap)
    β = β_interp(t_cap)
    γ = γ_interp(t_cap)
    x = x_interp(t_cap)
    y = y_interp(t_cap)
    z = z_interp(t_cap)

    println("Yaw: $(α*180/π)")
    println("Pitch: $(β*180/π)")
    println("Roll: $(γ*180/π)")


    visCoords = Array{Float64}(undef, 3, Nx_vis, Ny_vis)  # val, row, col
    visTimes = Array{Float64}(undef, Nx_vis, Ny_vis)

    thermCoords = Array{Float64}(undef, 3, Nx_therm, Ny_therm)
    thermTimes = Array{Float64}(undef, Nx_therm, Ny_therm)


    f_vis = ((Ny_vis-1)/2)/tand(θy_vis/2)
    f_therm = ((Ny_therm-1)/2)/tand(θy_therm/2)

    # i,j ordering: 43.412 s with 13.36 GiB allocations
    # j,i ordering: 42.117 s with 13.35 GiB allocations




    s = (z-z_ground)/f_vis  # scale factor
    # loop through pixels and apply the transformation
    # @showprogress 1 "Georectifying the visible image..." for i∈1:Nx_vis, j∈1:Ny_vis
    @showprogress 1 "Georectifying the visible image..." for j∈1:Ny_vis, i∈1:Nx_vis

        # NOTE: we need an extra rotation since FLIR is on the side of the HSI
        # rs_utm = [x; y; z] .+ s .*Rotate(-π,0,0) * T_n_E*Rotate(α, β, γ)*[(Nx_vis-1)/2 - (i-1); -(Ny_vis-1)/2 + (j-1); f_vis]
        rs_utm = [x; y; z] .+ s .*Rotate(3π/4,pitch_adjustment*π/180,0) * T_n_E*Rotate(α, β, γ)*[(Nx_vis-1)/2 - (i-1); -(Ny_vis-1)/2 + (j-1); f_vis]

        rs_object_lla = LLAfromUTMZ(wgs84)(UTMZ(rs_utm..., imu_df.zone[1], imu_df.isnorth[1]))

        @inbounds visCoords[1, i, j] = rs_object_lla.lat
        @inbounds visCoords[2, i, j] = rs_object_lla.lon
        @inbounds visCoords[3, i, j] = rs_object_lla.alt
        @inbounds visTimes[i,j] = t_cap
    end


    s = (z-z_ground)/f_therm  # scale factor
    # @showprogress 1 "Georectifying the thermal image..." for i∈1:Nx_therm, j∈1:Ny_therm
    @showprogress 1 "Georectifying the thermal image..." for j∈1:Ny_therm, i∈1:Nx_therm
        # NOTE: we need an extra rotation since FLIR is on the side of the HSI
        # rs_utm = [x; y; z] .+ s .* Rotate(-π,0,0) * T_n_E*Rotate(α, β, γ)*[(Nx_therm-1)/2 - (i-1); -(Ny_therm-1)/2 + (j-1); f_therm]
        rs_utm = [x; y; z] .+ s .*Rotate(3π/4,pitch_adjustment*π/180,0) * T_n_E*Rotate(α, β, γ)*[(Nx_therm-1)/2 - (i-1); -(Ny_therm-1)/2 + (j-1); f_therm]

        rs_object_lla = LLAfromUTMZ(wgs84)(UTMZ(rs_utm..., imu_df.zone[1], imu_df.isnorth[1]))

        @inbounds thermCoords[1, i, j] = rs_object_lla.lat
        @inbounds thermCoords[2, i, j] = rs_object_lla.lon
        @inbounds thermCoords[3, i, j] = rs_object_lla.alt
        @inbounds thermTimes[i, j] = t_cap
    end

    return img_vis, img_therm, visCoords, visTimes, thermCoords, thermTimes
end



end
