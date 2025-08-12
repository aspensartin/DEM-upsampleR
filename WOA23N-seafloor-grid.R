# generate 2D seafloor grid from WOA23N NetCDF

library(dplyr)
library(reshape2)
library(RNetCDF)

WOA23N_nc <- open.nc("WOA23N/woa23_decav71A0_o00_01.nc")
lat <- var.get.nc(WOA23N_nc, "lat")
lon <- var.get.nc(WOA23N_nc, "lon")
depth <- var.get.nc(WOA23N_nc, "depth")
o2 <- var.get.nc(WOA23N_nc, "o_an")
close.nc(WOA23N_nc)

WOA23N <- as.data.frame(
    cbind(
        rep(lon, times = length(lat), each = 1),
        rep(lat, times = 1, each = length(lon)),
        as.data.frame(melt(depth))$value,
        as.data.frame(melt(o2))$value
    )
)
names(WOA23N) <- c("lon", "lat", "depth", "o2")

    group_by(lon, lat) %>%
    slice_max(order(depth), n = 1) %>%
    ungroup() %>%






