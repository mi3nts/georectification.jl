# do this if you don't have contextily installed
# using Conda
# Conda.pip_interop(true)
# Conda.pip("install", ["toml", "packaging","contextily"])


using Plots
using PyCall
using Images


"""
    getBackgroundTile(w::Float64, n::Float64, e::Float64, s::Float64, outName::String)


Grab Esri World Imagery tiles for region with corners (w,n), (e,s) in longitude and latitude.
Saves resulting image to `outpath`

**Note:** result images are saved in Web-Mercator projection by default. See `WebMercatorfromLLA` and `LLAfromWebMercator` from `Geodesy.jl` for conversion details.
"""
function getBackgroundTile(w::Float64, n::Float64, e::Float64, s::Float64, outName::String)
    ctx = pyimport("contextily")

    # look at contextily docs for other provider options
    tiff, ext = ctx.bounds2raster(w, s, e, n, outName*".tiff", source=ctx.providers["Esri"]["WorldImagery"], ll=true)
    warped_tiff, warped_ext = ctx.warp_tiles(tiff, ext, "EPSG:4326")

    transformed_tiff = reinterpret(N0f8, warped_tiff)
    rgb_tiff = RGBA.(transformed_tiff[:,:,1],
                     transformed_tiff[:,:,2],
                     transformed_tiff[:,:,3],
                     transformed_tiff[:,:,4])
    longitudes = range(warped_ext[1], stop=warped_ext[2], length=size(rgb_tiff)[2])
    latitudes = range(warped_ext[4], stop=warped_ext[3], length=size(rgb_tiff)[1])

    # let's do we we mentioned before
    rgb_tiff = collect(rawview(channelview(reverse(rgb_tiff, dims=1)))) ./ 255
    longitudes = collect(longitudes)
    latitudes = collect(latitudes)

   return RGB.(rgb_tiff[1,:,:], rgb_tiff[2,:,:], rgb_tiff[3,:,:]), longitudes, latitudes
end



"""
    plot_background(sat_tile, longitudes, latitudes)

Given a satellite background tile `sat_tile` and associated longitude and latitude ranges, construct a plot 
"""
function plot_background(sat_tile, longitudes, latitudes; args...)
    Plots.set_theme(:dao)  # I like this theme currently
    # p = plot([longitudes[1], longitudes[end]],
    #          [latitudes[1], latitudes[end]],
    #          sat_tile,
    #          yflip=false,
    #          xlabel="Longitude",
    #          ylabel="Latitude",
    #          )

    p = plot(longitudes,
             latitudes,
             sat_tile,
             yflip=false,
             xlabel="Longitude",
             ylabel="Latitude";
             args...
             )



    return p
end

