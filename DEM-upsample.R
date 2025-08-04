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
expt <- "cGENIE/muffin.CB.L23z_000_z32.BASES-1.config/biogem/fields_biogem_3d.nc"
dem <- "PALEOMAP/Map01_PALEOMAP_1deg_Holocene_0Ma.nc"
vari <- "ocn_O2"
hlim <- 12
zlim <- 1

DEM.upsample <- function(expt, dem, vari, hlim, zlim) {

    library(dplyr)
    library(reshape2)
    library(geosphere)
    library(RNetCDF)

    # variable names within model output and DEM .nc files
    # configured for cGENIE and PALEOMAP by default, modify as needed
    latn_model <- "lat"
    lonn_model <- "lon"
    depthn_model <- "zt"
    timen_model <- "time"
    latn_DEM <- "lat"
    lonn_DEM <- "lon"
    depthn_DEM <- "z"

    # extract cGENIE data
    nc_model <- open.nc(expt)
    lat_model <- var.get.nc(nc_model, latn_model)
    lon_model <- var.get.nc(nc_model, lonn_model)
    depth_model <- var.get.nc(nc_model, depthn_model)
    var_model <- var.get.nc(nc_model, vari)
    time_model <- var.get.nc(nc_model, timen_model)
    close.nc(nc_model)
    # find final timestep
    timestep_model <- length(time_model)
    # adjust longitudes
    if (mean(between(lon_model, -180, 180)) < 1) {
        lon_model[lon_model <= -180] <- lon_model[lon_model <= -180] + 360
    }
    # create dataframes for each depth level
    slices_model <- vector("list", length = length(unique(depth_model)))
    for (i in seq_along(slices_model)) {
        slices_model[[i]] <- as.data.frame(
            cbind(
                rep(lon_model, times = length(lat_model), each = 1),
                rep(lat_model, times = 1, each = length(lon_model)),
                as.data.frame(melt(var_model[, , i, timestep_model]))$value
            )
        )
        names(slices_model[[i]]) <- c("lon", "lat", "var")
    }

    # extract DEM data and compile to dataframe
    nc_DEM <- open.nc(dem)
    lat_DEM <- var.get.nc(nc_DEM, latn_DEM)
    lon_DEM <- var.get.nc(nc_DEM, lonn_DEM)
    depth_DEM <- var.get.nc(nc_DEM, depthn_DEM)
    close.nc(nc_DEM)
    depth_DEM <- depth_DEM * -1
    df_DEM <- as.data.frame(
        cbind(
            rep(lon_DEM, times = length(lat_DEM), each = 1),
            rep(lat_DEM, times = 1, each = length(lon_DEM)),
            as.data.frame(melt(depth_DEM))$value
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
    mat_model <- as.matrix(slices_model[[1]][, c("lon", "lat")])
    mat_dist <- distm(mat_DEM, mat_model, fun = distCosine)
    depths_model <- as.vector(sort(unique(depth_model)))

    df_DEM <- df_DEM %>% mutate(row_index = row_number()) %>%
        rowwise() %>%
        mutate(
            # find indices of hlim nearest cGENIE cells
            ord = list(order(mat_dist[row_index, ])),
            order_index = list(ord[seq_len(min(hlim, length(ord)))]),
            var_interp = case_when(
                # shallow case:
                depth <= min(depths_model) ~ {
                    df_layer <- slices_model[[1]]
                    # find nearest valid cell
                    first_valid <- which(!is.na(df_layer$var[order_index]))
                    if (length(first_valid) > 0) {
                        df_layer$var[order_index[first_valid[1]]]
                    } else {
                        NA_real_
                    }
                },
                # deep case:
                depth >= max(depths_model) ~ {
                    df_layer <- slices_model[[length(slices_model)]]
                    # find nearest valid cell
                    first_valid <- which(!is.na(df_layer$var[order_index]))
                    if (length(first_valid) > 0) {
                        df_layer$var[order_index[first_valid[1]]]
                    }
                    # if zlim > 0 pull from layers above
                    else if (zlim > 0) {
                        found_var <- NA_real_
                        for (i in seq(1, zlim)) {
                            layer_index <- length(slices_model) - i
                            if (layer_index < 1) {
                                break
                            }
                            df_layer <- slices_model[[layer_index]]
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
                        vec = depths_model,
                        all.inside = TRUE
                    )
                    index_below <- index_above + 1
                    layer_above <- slices_model[[index_above]]
                    layer_below <- slices_model[[index_below]]
                    # find nearest valid cells for each layer
                    first_valid_above <- which(!is.na(layer_above$var[order_index]))
                    first_valid_below <- which(!is.na(layer_below$var[order_index]))
                    # interpolate if both cells are valid
                    if (length(first_valid_above) > 0 && length(first_valid_below) > 0) {
                        depth_above <- depths_model[index_above]
                        depth_below <- depths_model[index_below]
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
                            df_layer <- slices_model[[layer_index]]
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

