using georectification
using Plots

scotty = Dict("w" => -97.717472,
              "n" => 33.703572,
              "s" => 33.700797,
              "e" => -97.712413,
              )




# Scotty's ranch
w = scotty["w"]
n = scotty["n"]
s = scotty["s"]
e = scotty["e"]
out = "scotty"

rgb_tiff, longitudes, latitudes = getBackgroundTile(w, n, e, s, out)


typeof(rgb_tiff)
typeof(longitudes)
typeof(latitudes)


plot_background(rgb_tiff, longitudes, latitudes)


rgb_tiff


f(x, a) = begin
    0.0
end
xs = collect(0.1:0.05:2.0)
as = collect(0.2:0.1:2.0)
x_grid = [x for x = xs for y = as]
a_grid = [y for x = xs for y = as]
plot(x_grid, a_grid, f.(x_grid, a_grid), st = :surface xlabel = "longer xlabel", ylabel = "longer ylabel", zlabel = "longer zlabel")

plotattr(:Plot)

# we want a 3d plot with the bottom being the basemap image
# revisit this later in Makie http://juliaplots.org/MakieReferenceImages/gallery/index.html
# also: http://juliaplots.org/MakieReferenceImages/gallery//surface_with_image/index.html
