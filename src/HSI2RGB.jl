# using DataInterpolations
using LinearAlgebra
using Interpolations
using Trapz


include("./D_illuminants.jl")

"""
    Spec_2_RGB(λs_spec::AbstractVector, Spectrum::AbstractVector, d::Int)

Convert a pixel spectrum from an HSI to an RGB representation.
# Arguments
- `λs_spec`::AbstractVector`: Wavelength values in nm of the spectrum
- `Spectrum::AbstractVector`: Vector of reflectance or radiance values
- `d::Int`: 50, 55, 65, 75. Determines the illuminant used. If in doubt, use d=65.

Source:
*M. Magnusson, J. Sigurdsson, S. E. Armansson, M. O. Ulfarsson, H. Deborah and J. R. Sveinsson, "Creating RGB Images from Hyperspectral Images Using a Color Matching Function," IGARSS 2020 - 2020 IEEE International Geoscience and Remote Sensing Symposium, 2020, pp. 2045-2048, doi: 10.1109/IGARSS39084.2020.9323397.*

"""
function Spec_2_RGB(λs_spec::AbstractVector, Spectrum::AbstractMatrix, d::Int)
    wxyz = D_illuminants["wxyz"]
    D = D_illuminants["D"]
    λs_cmf = wxyz[1]
    λs_cmf_range = range(λs_cmf[1], λs_cmf[end], length=size(λs_cmf, 1))
    x = wxyz[2]
    y = wxyz[3]
    z = wxyz[4]


    i = Dict(50 => 3,
             55 => 4,
             65 => 2,
             75 => 5
             )


    λs_illuminant = D[1]  # wavelengths of the illuminant reference spectrum
    λs_illuminant_range = range(λs_illuminant[1], λs_illuminant[end], length=size(λs_illuminant, 1))
    I = D[i[d]]


    I_itp = CubicSplineInterpolation(λs_illuminant_range, I, extrapolation_bc=Linear())
    Î = I_itp.(λs_spec)

    x_itp = CubicSplineInterpolation(λs_cmf_range, x, extrapolation_bc=Linear())
    x̄ = x_itp.(λs_spec)

    y_itp = CubicSplineInterpolation(λs_cmf_range, y, extrapolation_bc=Linear())
    ȳ = y_itp.(λs_spec)

    z_itp = CubicSplineInterpolation(λs_cmf_range, z, extrapolation_bc=Linear())
    z̄ = z_itp.(λs_spec)



    # clip HSI wavelengths to fall within the visible range i.e. up to 780 nm
    idx = findfirst(λs_cmf[end] .< λs_spec) - 1
    λs_spec = λs_spec[1:idx]
    Spectrum = Spectrum[1:idx]
    Î = Î[1:idx]
    x̄ = x̄[1:idx]
    ȳ = ȳ[1:idx]
    z̄ = z̄[1:idx]



    # compute k := 1/N, N = ∫I(λ)ȳ(λ)dλ
    k = 1.0/trapz(λs_spec, ȳ .* Î)

    # compute X,Y, and Z
    X = k * trapz(λs_spec, Spectrum .* Î .* x̄)  # X = ∫S(λ)I(λ)x̄(λ)dλ
    Y = k * trapz(λs_spec, Spectrum .* Î .* ȳ)  # Y = ∫S(λ)I(λ)ȳ(λ)dλ
    Z = k * trapz(λs_spec, Spectrum .* Î .* z̄)  # Z = ∫S(λ)I(λ)z̄(λ)dλ

    # collect into 3 x n matrix
    XYZ = vcat(X', Y', Z')

    # XYZ → RGB Matrix
    M = [
          3.2404542 -1.5371385 -0.4985314;
          -0.9692660 1.8760108 0.0415560;
          0.0556434 -0.2040259 1.0572252
        ]

    sRGB = M * XYZ   # get an RGB x n matrix

    # gamma correction
    idx_high = sRGB .> 0.0031308
    idx_low = sRGB .<= 0.0031308

    sRGB[idx_high] .= 1.055 .* sRGB[idx_high].^(1.0/2.4) .- 0.055
    sRGB[idx_low] .= 12.92 .* sRGB[idx_low]

    # clip RGB values
    sRGB[sRGB .> 1] .= 1
    sRGB[sRGB .< 0] .= 0


    return sRGB
end



"""
    HSI_2_RGB(λs_spec::AbstractVector, HSI::AbstractMatrix, d::Int)

Convert reflectance data from an HSI to an RGB representation.
# Arguments
- `λs_spec`::AbstractVector`: Wavelength values in nm of the spectrum
- `HSI::AbstractMatrix`:  Matrix of reflectance or radiance values with shape npixels x nbands
- `d::Int`: 50, 55, 65, 75. Determines the illuminant used. If in doubt, use d=65.

Source:
*M. Magnusson, J. Sigurdsson, S. E. Armansson, M. O. Ulfarsson, H. Deborah and J. R. Sveinsson, "Creating RGB Images from Hyperspectral Images Using a Color Matching Function," IGARSS 2020 - 2020 IEEE International Geoscience and Remote Sensing Symposium, 2020, pp. 2045-2048, doi: 10.1109/IGARSS39084.2020.9323397.*

"""
function HSI_2_RGB(λs_spec::AbstractVector, HSI::AbstractMatrix, d::Int)
    wxyz = D_illuminants["wxyz"]
    D = D_illuminants["D"]
    λs_cmf = wxyz[1]
    λs_cmf_range = range(λs_cmf[1], λs_cmf[end], length=size(λs_cmf, 1))
    x = wxyz[2]
    y = wxyz[3]
    z = wxyz[4]


    i = Dict(50 => 3,
             55 => 4,
             65 => 2,
             75 => 5
             )


    λs_illuminant = D[1]  # wavelengths of the illuminant reference spectrum
    λs_illuminant_range = range(λs_illuminant[1], λs_illuminant[end], length=size(λs_illuminant, 1))
    I = D[i[d]]


    I_itp = CubicSplineInterpolation(λs_illuminant_range, I, extrapolation_bc=Linear())
    Î = I_itp.(λs_spec)

    x_itp = CubicSplineInterpolation(λs_cmf_range, x, extrapolation_bc=Linear())
    x̄ = x_itp.(λs_spec)

    y_itp = CubicSplineInterpolation(λs_cmf_range, y, extrapolation_bc=Linear())
    ȳ = y_itp.(λs_spec)

    z_itp = CubicSplineInterpolation(λs_cmf_range, z, extrapolation_bc=Linear())
    z̄ = z_itp.(λs_spec)



    # clip HSI wavelengths to fall within the visible range i.e. up to 780 nm
    idx = findfirst(λs_cmf[end] .< λs_spec) - 1
    λs_spec = λs_spec[1:idx]
    HSI = HSI[:, 1:idx]
    Î = Î[1:idx]
    x̄ = x̄[1:idx]
    ȳ = ȳ[1:idx]
    z̄ = z̄[1:idx]



    # compute k := 1/N, N = ∫I(λ)ȳ(λ)dλ
    k = 1.0/trapz(λs_spec, ȳ .* Î)

    # compute X,Y, and Z
    X = k * trapz((:,λs_spec), HSI * diagm(Î .* x̄))  # X = ∫S(λ)I(λ)x̄(λ)dλ
    Y = k * trapz((:,λs_spec), HSI * diagm(Î .* ȳ))  # Y = ∫S(λ)I(λ)ȳ(λ)dλ
    Z = k * trapz((:,λs_spec), HSI * diagm(Î .* z̄))  # Z = ∫S(λ)I(λ)z̄(λ)dλ

    # collect into 3 x n matrix
    XYZ = vcat(X', Y', Z')

    # XYZ → RGB Matrix
    M = [
          3.2404542 -1.5371385 -0.4985314;
          -0.9692660 1.8760108 0.0415560;
          0.0556434 -0.2040259 1.0572252
        ]

    sRGB = M * XYZ   # get an RGB x n matrix

    # gamma correction
    idx_high = sRGB .> 0.0031308
    idx_low = sRGB .<= 0.0031308

    sRGB[idx_high] .= 1.055 .* sRGB[idx_high].^(1.0/2.4) .- 0.055
    sRGB[idx_low] .= 12.92 .* sRGB[idx_low]

    # clip RGB values
    sRGB[sRGB .> 1] .= 1
    sRGB[sRGB .< 0] .= 0


    return sRGB
end


