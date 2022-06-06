
"""
    getHdrFile(pathToHdrFile::String)

Read an ENVI-style header file and return a dictionary with image metadata.
"""
function getHdrFile(path::String)
    lines = readlines(path)
    splitlines = [split(line, " = ") for line in lines]
    hdr_dict = Dict()

    for splitline ∈ splitlines
        if length(splitline) > 1
            K = join([word for word ∈ split(splitline[1], " ") if word != ""], "_") 
            if K == "data_type"
                V = parse(Int, splitline[2])
            elseif K == "lines"
                V = parse(Int, splitline[2])
            elseif K == "samples"
                V = parse(Int, splitline[2])
            elseif K == "bands"
                V = parse(Int, splitline[2])
            elseif K == "ceiling"
                V = parse(Int, splitline[2])
            elseif K == "sample_binning"
                V = parse(Int, splitline[2])
            elseif K == "spectral_binning"
                V = parse(Int, splitline[2])
            elseif K == "shutter"
                V = parse(Float64, splitline[2])
            elseif K == "gain"
                V = parse(Int, splitline[2])
            elseif K == "framerate"
                V = parse(Int, splitline[2])
            elseif K == "byte_order"
                V = parse(Int, splitline[2])
            elseif K == "header_offset"
                V = parse(Int, splitline[2])
            elseif K == "flip_radiometric_calibration"
                V = parse(Bool, lowercase(splitline[2]))

            elseif K == "wavelength"
                V = [parse(Float64, val) for val ∈ split(splitline[2][2:end-1], ",")]
            else
                V = splitline[2]
            end
            hdr_dict[K] = V
            #println("$(K)...$(V)")
        end
    end
    return hdr_dict
end



"""
    getTimes(pathToTimes::String)

Read a .times and return values starting at `t=0`. Times are assumed to be aligned with the corresponding .lcf file.
"""
function getTimes(pathToTimes::String)
    lines = readlines(pathToTimes)
    times = [parse(Float64, line) for line in lines]
    times = times .- times[1] # assume .lcf and .times start at the same time
    return times
end





"""
    readToDataFrame(pathToBinaryFile::String, pathToHdrFile::String, pathToTimesFile::String)

Read ENVI file with corresponding header (.hdr) and times (.times) file and return a DataFrame containing the image data flattened to (nbands, ncols*nrows).

"""
function readToDataFrame(pathToBinaryFile::String, pathToHdrFile::String, pathToTimes::String="")
    typeDict = Dict()
    # fill in the dict with the type translation
    typeDict[2] = Int16
    typeDict[3] = Int32
    typeDict[4] = Float32
    typeDict[5] = Float64 # Double
    typeDict[9] = ComplexF64 # 2xDouble complex number
    typeDict[12] = UInt16

    info = getHdrFile(pathToHdrFile)
    ncols = info["samples"]
    nbands = info["bands"]
    nrows = info["lines"]

    # set up arrays to read data into
    if info["interleave"] == "bil"
        data = Array{typeDict[info["data_type"]]}(undef, info["samples"], info["bands"], info["lines"])
        read!(pathToBinaryFile, data)
        data = reshape(PermutedDimsArray(data, (2,1,3)), (nbands, ncols*nrows))
        info["data_shape"] = "band, col×row"

    elseif info["interleave"] == "bip"
        data = Array{typeDict[info["data_type"]]}(undef, info["bands"], info["samples"], info["lines"])
        read!(pathToBinaryFile, data)
        data = reshape(data), (nbands, ncols*nrows)
        info["data_shape"] = "band, col×row"
    elseif info["interleave"] == "bsq"
        data = Array{typeDict[info["data_type"]]}(undef, info["samples"], info["lines"], info["bands"])
        read!(pathToBinaryFile, data)
        data = reshape(PermutedDimsArray(data, (3,1,2)), (nbands, ncols*nrows))
        info["data_shape"] = "band, col×row"

    else
        println("Interleave unrecognized... Leaving data empty")
        data = []
    end

    # fill array with the data

    if pathToTimes == ""
        times = []
    else
        times = getTimes(pathToTimes)
    end

    times_array = vcat([times' for _ ∈ 1:ncols]...)
    times_final = reshape(times_array, ncols*nrows)

    row_indices = vcat([collect(1:nrows)' for _ ∈1:ncols]...)
    rows = reshape(row_indices, ncols*nrows)


    names = ["λ_$(i)_rad" for i ∈ 1:462]
    df = DataFrame(data', names)  # include float64 cast

    # generate Σrad now for efficiency
    df[!, "Σrad"] = vec(sum(data, dims=1))


    df[!, "times"] = times_final

    df[!, "row_index"] = rows
    return df
end




# for now this is necessary for dealing with the calibration files

struct HSI{S<:AbstractArray, T<:Dict, U<:AbstractArray}
    data::S
    info::T
    times::U
end

