#' cGENIE DEM Interpolation
#'
#' This function interpolates a cGENIE 3D output field to the points of a
#' digital elevation model (DEM). If a point is outside the range of cGENIE
#' mid-layer depths, it is assigned the value from the closest valid top/bottom
#' layer cell. If a point is within that range, a value is interpolated using
#' the two closest valid cells immediately above and below it. This function
#' returns a dataframe 'df_DEM' and saves a plot.
#' 
#' @param vari A string specifying the variable name to extract from the NetCDF
#'              file. Examples include "ocn_temp", "ocn_O2", etc.
#' @param expt A string specifying the path to the folder or file prefix
#'                   where the experiment's data is stored.
#' @param dem  A string specifying the path to the DEM NetCDF file.
#' @param dist A float specifying maximum nearest-neighbour distance when matching
#'             DEM points to the centre of cGENIE cells in km (recommended: 1000)
#'
#' @return A dataframe (`df.sum`) containing:
#'   - `lon`: Longitude of DEM points.
#'   - `lat`: Latitude of DEM points.
#'   - `depth`: Depth at DEM points.
#'   - `var`: Interpolated tracer values (e.g., [O2])
#'
#' @import RNetCDF, dplyr, reshape2, ggplot2, RNetCDF, sf
#' @export

# testing
expt <- "../L23j/muffin.CB.L23j_025_k32_bathy+clim.BASES-1.config"
dem <- "../Scotese/Map08_PALEOMAP_1deg_Late_Oligocene_25Ma.nc"
vari <- "ocn_O2"
dist <- 1000
zlim <- 2
hlim <- 9

cGENIE.DEM.interp <- function(expt, dem, vari, hlim, zlim) {

    library(dplyr)
    library(PaleoClimR)
    library(geosphere)
    library(RNetCDF)

    # extract cGENIE data
    df_cGENIE <- cGENIE.data.3D(vari, expt, "default", "biogem")
    depths_cGENIE <- as.vector(sort(unique(df_cGENIE$depth)))
    slices_cGENIE <- vector("list", length = length(depths_cGENIE))
    for (i in seq_along(depths_cGENIE)) {
        slices_cGENIE[[i]] <- filter(df_cGENIE, depth.level == i)
    }

    # extract DEM data and compile to dataframe
    nc <- open.nc(dem)
    lat <- var.get.nc(nc, "lat")
    lon <- var.get.nc(nc, "lon")
    depth <- var.get.nc(nc, "z")
    close.nc(nc)
    depth <- depth * -1
    df_DEM <- as.data.frame(
        cbind(
            rep(lon, times = length(lat), each = 1),
            rep(lat, times = 1, each = length(lon)),
            as.data.frame(melt(depth))$value
        )
    )
    names(df_DEM) <- c("lon", "lat", "depth")

    # filter out land from DEM dataframe
    df_DEM <- df_DEM %>%
        filter(depth >= 0)

    # adjust DEM longitudes
    if (mean(between(df_DEM$lon, -180, 180)) < 1) {
        df_DEM$lon[df_DEM$lon <= -180] <- df_DEM$lon[df_DEM$lon <= -180] + 360                            
    }

    # create cGENIE-DEM distances lookup matrix
    mat_DEM <- as.matrix(df_DEM[, c("lon", "lat")])
    mat_cGENIE <- as.matrix(slices_cGENIE[[1]][, c("lon.mid", "lat.mid")])
    mat_dist <- distm(mat_DEM, mat_cGENIE, fun = distCosine)

    df_DEM <- df_DEM %>% mutate(row_index = row_number()) %>%
        rowwise() %>%
        mutate(
            # find indices of hlim nearest cGENIE cells
            ord = list(order(mat_dist[row_index, ])),
            order_index = list(ord[seq_len(min(hlim, length(ord)))]),
            var_interp = case_when(
                # shallow case:
                depth <= min(depths_cGENIE) ~ {
                    df_layer <- slices_cGENIE[[1]]
                    # find nearest valid cell
                    first_valid <- which(!is.na(df_layer$var[order_index]))
                    if (length(first_valid) > 0) {
                        df_layer$var[order_index[first_valid[1]]]
                    } else {
                        NA_real_
                    }
                },
                # deep case:
                depth >= max(depths_cGENIE) ~ {
                    df_layer <- slices_cGENIE[[length(slices_cGENIE)]]
                    # find nearest valid cell
                    first_valid <- which(!is.na(df_layer$var[order_index]))
                    if (length(first_valid) > 0) {
                        df_layer$var[order_index[first_valid[1]]]
                    }
                    # if zlim > 0 pull from layers above
                    else if (zlim > 0) {
                        found_var <- NA_real_
                        for (i in seq(1, zlim)) {
                            layer_index <- length(slices_cGENIE) - i
                            if (layer_index < 1) {
                                break
                            }
                            df_layer <- slices_cGENIE[[layer_index]]
                            first_valid <- which(!is.na(df_layer$var[order_index]))
                            if (length(first_valid) > 0) {
                                found_var <- df_layer$var[order_index[first_valid[1]]]
                                break
                            }
                        }
                        found_var  
                    }
                    else {
                        NA_real_
                    }
                },
                # intermediate case:
                TRUE ~ {
                    # find layers immediately above and below depth
                    index_above <- findInterval(
                        x = depth,
                        vec = depths_cGENIE,
                        all.inside = TRUE
                    )
                    index_below <- index_above + 1
                    layer_above <- slices_cGENIE[[index_above]]
                    layer_below <- slices_cGENIE[[index_below]]
                    # find nearest valid cells for each layer
                    first_valid_above <- which(!is.na(layer_above$var[order_index]))
                    first_valid_below <- which(!is.na(layer_below$var[order_index]))
                    # interpolate if both cells are valid
                    if (length(first_valid_above) > 0 && length(first_valid_below) > 0) {
                        depth_above <- depths_cGENIE[index_above]
                        depth_below <- depths_cGENIE[index_below]
                        var_above <- layer_above$var[order_index[first_valid_above[1]]]
                        var_below <- layer_below$var[order_index[first_valid_below[1]]]
                        approx(
                            x = c(depth_above, depth_below),
                            y = c(var_above, var_below),
                            xout = depth
                        )$y
                    }
                    # if zlim > 0 pull from layers above
                    else if (zlim > 0) {
                        found_var <- NA_real_
                        for (i in seq(1, zlim)) {
                            layer_index <- index_above - i
                            if (layer_index < 1) {
                                break
                            }
                            df_layer <- slices_cGENIE[[layer_index]]
                            first_valid <- which(!is.na(df_layer$var[order_index]))
                            if (length(first_valid) > 0) {
                                found_var <- df_layer$var[order_index[first_valid[1]]]
                                break
                            }
                        }
                        found_var
                    }
                    else {
                        NA_real_
                    }
                }
            )
        ) %>%
        ungroup() %>%
        select(-row_index, -ord, -order_index)

    return(df_DEM)
}

