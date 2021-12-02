########################################################################
### Title: Run GRTS
### Author: Will Vieira
### December 9, 2020
########################################################################


library(tidyverse)
library(spatstat)
library(sf)


hexas <- readRDS('data/hexa_complete.RDS')



# calculate the inclusion probability p for each hexagon

    # First filter for hexagons that have at least 50% of samplable pixel habitats within the hexagon (no NAs)
    hexas <- subset(hexas, propNA <= 0.8)
    
    # calculate p
    hexas$p <- (hexas$hab_prob * hexas$cost_prob) / sum(hexas$hab_prob * hexas$cost_prob)
    
    # remove zeros
    hexas <- subset(hexas, p != 0)
    
    # sequencial row names
    row.names(hexas) <- seq(1, nrow(hexas))

#




# Get sample size for each ecoregion

    sampleSize <- readRDS('data/nbHexa_ecoregion.RDS')

    # Sample 2% of hexagons available hexagons
    N_tmp <- sampleSize$hexagonProp80 * 0.02
    
    # Sample size in function of size proportion between ecoregions
    eco_ignore <- which(sampleSize$ecoregion %in% gsub('N', '', grep('N', sampleSize$ecoregion, value = TRUE)))
    totalN <- sum(N_tmp[-eco_ignore])
    sampleSize$N <- unclass(round(totalN * sampleSize$sizeProp, 0))
    
    # Remove ecoregions with N == 0
    sampleSize <- subset(sampleSize, N != 0)
    
    # Total sample size
    N <- sum(sampleSize$N[-eco_ignore])

#
 


## Inclusion probability depening on simulation

    ## SIMULATION 2:
    # Weight p according to the number of legacy site
    we = function(legacySite, lower, upper, mid, beta)
    {
        if(legacySite == 0) {
            return( 1 )
        }else{
            return( lower + ((upper - lower) / (1 + exp(-beta * (legacySite - mid)))) )
        }
    } 

    # weight = 3
    hexas$p_sim2 <- hexas$p * sapply(hexas$legacySite, we, lower = 1, upper = 3, mid = 4, beta = 0.9)
    hexas$p_sim2 <- hexas$p_sim2/sum(hexas$p_sim2)

#




# GRTS simulations

    # number of replications for each simulation
    nb_rep <- 1

    # Ecoregions
    eco_sim <- c('7', '28', '30', '31', '46', '47', '48', '49', '73', '77', '78', '86',
                 '72', '74', '75', '76',
                 '96N', #'96', '96S',
                 '99',
                 '100N', #'100', '100S',
                 '101N', #'101', '101S',
                 '102',
                 '103N', #'103', '103S',
                 '117',
                 '216',
                 '217N')#, '217', '217S')
    

    # Design list
    overSampleSize <- 0.2
 
    Stratdsgn <- list()
    for(eco in eco_sim)
        Stratdsgn[[paste0('eco_', eco)]] <- list(panel = c(PanelOne = subset(sampleSize, ecoregion == eco)$N), over = round(subset(sampleSize, ecoregion == eco)$N * overSampleSize, 0), seltype = 'Continuous')

    # Sample frame for specific ecoregions
    sample_frame <- subset(hexas, ecoregion %in% eco_sim)

    # Get x and y from centroid of hexagon
    sample_frame$geometry <- sample_frame %>% sf::st_centroid() %>% sf::st_geometry()
    sample_frame[c('X', 'Y')] <- sf::st_coordinates(sample_frame)
    
    # rename ecoregion name to match design name
    sample_frame$eco_name <- paste0('eco_', sample_frame$ecoregion)

    # Get only attributes table (remove spatial information)
    attframe <- sf::st_drop_geometry(sample_frame)



    # RUN GRTS
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% eco_sim)$N) * attframe$p/sum(attframe$p)

    set.seed(2)
    out_sim <- list()
    for (i in 1:nb_rep)
    {
        out_sample <- spsurvey::grts(design = Stratdsgn,
                                            DesignID = "ET_ID",
                                            type.frame = "finite",
                                            src.frame = "att.frame",
                                            att.frame = attframe,
                                            xcoord = 'X',
                                            ycoord = 'Y',
                                            stratum = "eco_name",
                                            mdcaty = "mdcaty",
                                            shapefile = FALSE,
                                            out.shape = "Cost_design")

        out_sim[[i]] <- setNames(out_sample$ET_Index, out_sample$panel)
    }

    # save output
    hexas_sample <- subset(hexas, ecoregion %in% eco_sim)

    for(i in 1:nb_rep)
    {
        hexas_sample[[paste0('rep', i, '_over')]] <- hexas_sample[[paste0('rep', i, '_main')]] <- rep(0, nrow(hexas_sample))
        
        # assign 1 to the selected hexagons (main sample)
        hexas_sample[[paste0('rep', i, '_main')]][which(hexas_sample$ET_Index %in% out_sim[[i]][which(names(out_sim[[i]]) == 'PanelOne')])] <- 1

        # assign 1 to the selected hexagons (over sample)
        hexas_sample[[paste0('rep', i, '_over')]][which(hexas_sample$ET_Index %in% out_sim[[i]][which(names(out_sim[[i]]) == 'OverSamp')])] <- 1

    }

    # TODO FIX REPRODUCIBILITY
    write_sf(hexas_sample, '../../ownCloud/BMS_Bruno/sample/sample_hypothese2.shp')


    # SAVE OUTPUT V2
    hexas %>%
        filter(ET_Index %in% out_sim[[1]]) %>%
        mutate(
            panel = if_else(ET_Index %in% out_sim[[1]][names(out_sim[[1]]) == 'PanelOne'], 'main', 'over')
        ) %>%
        write_sf('../../ownCloud/BMS_Bruno/sample/sample_quebec.shp')
        
#




#   # SELET THE BEST REPETITION BASED ON THE LOWEST SAMPLING COST
    ################################################

    cheapest_rep <- hexas_sample %>%
        st_drop_geometry() %>%
        select(ecoregion, ET_Index, costSum, starts_with('rep')) %>%
        pivot_longer(
            cols = starts_with('rep'),
            names_to = 'Rep',
            values_to = 'select'
        ) %>%
        mutate(rep_total = parse_number(Rep)) %>%
        group_by(rep_total) %>%
        filter(select == 1) %>%
        summarise(costTotal = sum(costSum)) %>%
        filter(costTotal == min(costTotal)) %>%
        select(rep_total) %>%
        as.numeric()




    # SELECT SECONDATY SAMPLE UNITS WITHIN HEXAGON
    ################################################

    hexas_sample <- read_sf('../../ownCloud/BMS_Bruno/sample/sample_hypothese2.shp')
    

    # Get cheapest repetition only 
    selected_hexas <- hexas_sample[which(
                        hexas_sample[[paste0('rep', cheapest_rep, '_main')]] == 1 |
                        hexas_sample[[paste0('rep', cheapest_rep, '_over')]] == 1), ]


    
    # FIRST APPROACH:
    # Complete random points (independant of habitat and availability of the unit)
    ########################
    
    out_ls <- list(); Count = 1
    for(eco in unique(hexas_sample$ecoregion))
    {
        hexas_eco <- subset(selected_hexas, ecoregion == eco)

        for(sim in c('main', 'over'))
        {
            
            hexas_eco_sim <- hexas_eco[which(hexas_eco[[paste0('rep', cheapest_rep, '_', sim)]] == 1), ]

            for(hex in 1:nrow(hexas_eco_sim))
            {
                hex_eco_sim <- hexas_eco_sim[hex, ]
                point_grid <- hex_eco_sim %>%
                    st_make_grid(cellsize = 300, what = 'centers') %>%
                    st_intersection(hex_eco_sim)

                point_ID <- 1:length(point_grid)
                for(Rep in 1:10)
                {
                    # Sample first point
                    sample_1 <- sample(point_ID, size = 1)

                    # Remove all points around the already sampled point
                    # using a buffer of 1200 x 1200 meters
                    toKeep <- !st_intersects(point_grid, st_buffer(point_grid[sample_1], dist = 855), sparse = FALSE)
                    sample_2 <- sample(point_ID[toKeep], size = 1)
                    
                    sample_grid <- point_grid[c(sample_1, sample_2)]
                    
                    # add attribute table
                    att_df <- data.frame(
                        ecoregion = eco,
                        ET_Index = hexas_eco_sim$ET_Index[hex],
                        simulation = sim,
                        repetition = Rep,
                        SSU_Index = c(sample_1, sample_2),
                        point = 1:2
                    )

                    out_ls[[Count]] <- st_sf(att_df, geometry = sample_grid)
                    Count = Count + 1
                }
            }
        }
    }
    
    out_point <- do.call(rbind, out_ls)

    out_point %>%
        st_drop_geometry() %>%
        group_by(ET_Index, point) %>%
        summarise(
            tb = sum(table(SSU_Index) > 1)
        ) %>%
        ggplot(aes(tb)) + geom_histogram()


    write_sf(out_point, '../../ownCloud/BMS_Bruno/sample/SSU.shp')


    # SECOND APPROACH:
    # Calculate probability following habitat proportion from ecoregion AND
    # Filter for squares with at least o
    ########################

    library(raster)
    library(exactextractr)

    land_ca <- raster("data/landcover_ca_30m.tif")
    prev_all <- readRDS('data/prev_all.RDS')    

    # Function to return squares polygons for a specific hexagon
    # Return sf object with square polygons for a specific hexagon
    get_squares <- function(hexa, cellSize, roads_bf)
    {
        # Get grid centroid of cellSize
        hex_cent <- hexa %>%
            st_make_grid(cellsize = cellSize, what = 'centers')
    
        # Get squares around grid centroid and filter the hexagons
        # in which the centroid are inside the hexagon
        hex_squares <- hexa %>%
            st_make_grid(cellsize = cellSize, square = TRUE)
        hex_squares <- hex_squares[st_intersects(hex_cent, hexa, sparse = FALSE), ]

        # Check if the squares have roads inside
        if(hexa$haveRoads == 1)
        {
            hvRoads <- st_intersects(hex_squares, st_crop(roads_bf, hexa), sparse = FALSE)
            hvRoads <- apply(hvRoads, 1, function(x) ifelse(sum(x) > 0, 1, 0))
        }else{
            hvRoads <- rep(0, length(hex_squares))
        }
        
        # squares with attribute table
        hex_squares <- st_sf(
            data.frame(
                ecoregion = hexa$ecoregion,
                ET_Index = hexa$ET_Index,
                SSU_Index = 1:length(hex_squares),
                haveRoads = hvRoads
            ),
            geometry = hex_squares)
    }


    # Function to calculate inclusion probability for each square within hexagon
    # Return a vector of inclusion probability for each square
    calc_habProb <- function(squares, landUse, incProb)
    {
        # extract pixels for each square
        hab_pixels <- exactextractr::exact_extract(landUse, squares, progress = FALSE)
        
        # get frequence of each class of habitat
        count_hab <- Map(
                    function(x, y) {
                        freq <- table(x$value)
                        if(length(freq) > 0) {
                            data.frame(freq, ecoregion = y)
                        }else{
                            NA
                        }
                    },
                    x = hab_pixels,
                    y = squares$ecoregion
                )
        

        # merge with inclusion probability
        # and calculate inclusion probbaility for each NON empty square
        squares$incl_prob <- unlist(
                lapply(
                    count_hab,
                    function(x) {
                        if(is.data.frame(x)) {
                            mg_df <- merge(x, subset(incProb, ID_ecoregion == x$ecoregion[1]), by.x = "Var1", by.y = "code" , all.x = TRUE)
                            sum(mg_df$Freq * mg_df$incl_prob)
                        }else{
                            NA
                        }
                    }
                )
            )

        squares <- squares %>%
            st_drop_geometry() %>%
            group_by(ET_Index) %>%
            mutate(
                inclProb = incl_prob/sum(incl_prob, na.rm = TRUE)
            )
        
        return(squares$inclProb)
    }

    # Function to sample SSU within selected hexagons
    sample_SSU <- function(hex_cent, repetitions)
    {
        # hex_cent index to later samples
        point_ID <- hex_cent$SSU_Index
        
        # create a second index vector to store the points in which
        # have a road AND have inclusion probability > 0
        roads_ID <- ifelse(hex_cent$haveRoads == 1, TRUE, FALSE)
        if(all(!roads_ID))
            roads_ID <- rep(TRUE, nrow(hex_cent))
        
        # If roads + incProb are not enough to sample (>1), forget about it
        if(sum(roads_ID & !is.na(hex_cent$incProb)) > 1)
        {
            roadsProb_ID <- roads_ID & !is.na(hex_cent$incProb)
        }else{
            roadsProb_ID <- !is.na(hex_cent$incProb)
        }
        
        
        # If all IDs are FALSE, we stop the function here as none of the points have habitat
        if(all(!roadsProb_ID))
            stop(paste('Inclusion probability returns a vector of NAs for Hexagon:', unique(hex_cent$ET_Index)))

        out_ls <- list(); Count = 1
        for(Rep in 1:repetitions)
        {
            # Sample first point
            sample_1 <- sample(point_ID[roadsProb_ID], size = 1, prob = hex_cent$incProb[roadsProb_ID])
            
            # Remove all points around the already sampled point
            # using a buffer of 1200 x 1200 meters
            toKeep <- !st_intersects(
                    hex_cent,
                    st_buffer(hex_cent[sample_1, ], dist = 855),
                    sparse = FALSE
                )

            if(sum(toKeep & roadsProb_ID) > 1)
            {
                sample_2 <- sample(point_ID[toKeep & roadsProb_ID], size = 1, prob = hex_cent$incProb[toKeep & roadsProb_ID])
            }else{
                sample_2 <- sample(point_ID[toKeep & !is.na(hex_cent$incProb)], size = 1, prob = hex_cent$incProb[toKeep & !is.na(hex_cent$incProb)])
            }
            
            sample_grid <- hex_cent[c(sample_1, sample_2), ]
                    
            # add attribute table
            sample_grid$Rep = Rep
            sample_grid$point = 1:2
                            
            out_ls[[Count]] <- sample_grid
            Count = Count + 1
        }
        hexa_points <- do.call(rbind, out_ls)
    }


    # Get which hexagons cross the buffer of 1km around the roads
    # For the hexagons with roads, we will try to sample inside the buffer only
    # If it's not possible for both points to be inside the buffer, then the second will be outside
    roads_bf <- roads %>%
        st_transform(st_crs(selected_hexas)) %>%
        st_buffer(dist = 1000)

    haveRoads <- st_intersects(selected_hexas, roads_bf, sparse = FALSE)     
    selected_hexas$haveRoads <- apply(haveRoads, 1, function(x) ifelse(sum(x) > 0, 1, 0))




    # Get squares for all selected hexagons
    selHx_squares <- map_dfr(
        seq_len(nrow(selected_hexas)),
        function(x) get_squares(selected_hexas[x, ], cellSize = 300, roads_bf)
    )
    
    saveRDS(selHx_squares, file = 'selHx_squares.RDS') 
    
    # Get inclusion probability for each square
    selHx_squares$incProb = calc_habProb(
        squares = selHx_squares,
        landUse = land_ca,
        incProb = prev_all
    )


    # Now the habitat is extracted, we can keep the sampling point only (centroid)
    selHx_cent <- selHx_squares %>%
        st_centroid()


    # Sample points within hexagon
    out_sample <- map_dfr(
        unique(selHx_cent$ET_Index),
        function(x)
            sample_SSU(
                 hex_cent = subset(selHx_cent, ET_Index == x),
                 repetitions = 10
            )
        )

## Get the points around the selected SSU for Bruno, and save the X and Y coordinates

# list shp files
filesBruno <- list.files('../../Downloads/SSU_sélectionnés_Bruno/', pattern = "\\.shp$")
dir.create('SSU_selectionnes_Bruno')

for(File in filesBruno)
{
    # load hexagons
    couche_file <- read_sf(paste0('../../Downloads/SSU_sélectionnés_Bruno/', File))
    
    hexas_ID <- unique(couche_file$ET_Index)
    hexas_neighbors <- list(); Count = 1
    for(hx_ID in hexas_ID)
    {
        all_points_hx <- hexas %>%
            filter(ET_Index == hx_ID) %>%
            st_make_grid(cellsize = 300, what = 'centers')
        all_points_hx <- st_sf(
            data.frame(
                SSU_Index = NA
            ),
            geometry = all_points_hx
            
        )

        # Which are inside the hexagon?
        toKeep <- st_intersects(all_points_hx, subset(hexas, ET_Index == hx_ID), sparse = FALSE)
        all_points_hx$SSU_Index[toKeep] <- 1:sum(toKeep)
        all_points_hx$SSU_Index[!toKeep] <- seq(sum(toKeep) + 1, sum(toKeep) + sum(!toKeep))

        ref_points_hx <- subset(couche_file, ET_Index == hx_ID)

        all_points_hx$ecoregion = unique(ref_points_hx$ecoregion)
        all_points_hx$ET_Index = hx_ID

        neighbors_ls <- list()
        for(i in 1:nrow(ref_points_hx))
        {
            neighbors_i = all_points_hx[which(st_intersects(st_buffer(ref_points_hx[i, ], dist = 555), all_points_hx, sparse = FALSE)), ]    

            # adjust attribute table
            neighbors_i$repetition = ref_points_hx$repetition[i]
            neighbors_i$group <- ref_points_hx$point[i]
            neighbors_i$point <- 0
            neighbors_i$point[!ref_points_hx$SSU_Index[i] == neighbors_i$SSU_Index] = 1:sum(!ref_points_hx$SSU_Index[i] == neighbors_i$SSU_Index)

            neighbors_ls[[i]] <- neighbors_i
        }
        hexas_neighbors[[Count]] <- do.call(rbind, neighbors_ls)
        Count = Count + 1
    }

    couche_neighbors <- do.call(rbind, hexas_neighbors)

    # Organize
    couche_neighbors <- couche_neighbors[, 
                match(
                    c('ecoregion', 'ET_Index', 'repetition', 'SSU_Index', 'group', 'point', 'geometry'),
                    names(couche_neighbors)
                )
    ]

    # Get X and Y
    coords <- couche_neighbors %>%
        sf::st_transform(4326) %>%
        sf::st_coordinates() %>%
        as.data.frame()

    couche_neighbors$X <- coords$X
    couche_neighbors$Y <- coords$Y

    # save
    write_sf(couche_neighbors, paste0('SSU_selectionnes_Bruno/',File))

}
