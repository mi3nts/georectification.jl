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

