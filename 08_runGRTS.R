########################################################################
### Title: Run GRTS
### Author: Will Vieira
### December 9, 2020
########################################################################

library(raster)
library(exactextractr)
library(tidyverse)
library(spatstat)
library(sf)
library(spsurvey)
library(MBHdesign)

set.seed(0.0)



# parameters

    # Maximum proportion of NA pixels (non-habitat) in a hexagon
    # if 0.8, it means an hexagon has to have at least 20% of habitat pixels
    prop_na = 0.8

    # Percentage of hexagons to cover a region
    sample_effort = 0.02

    # Total sample size for SSU (Main + Over)
    # It must be a even number
    ssu_N = 4

    # Number of replications when running the GRTS
    nb_rep = 15

    # Ecoregions
    eco_sim <- c(
        '7', '28', '30', '31',
        '46', '47', '48', '49',
        '73', '77', '78', '86',
        '72', '74', '75', '76',
        '96', #'96N', '96S',
        '99',
        '100', #'100N', '100S',
        '101', #'101N', '101S',
        '102',
        '103', #'103N', '103S',
        '117',
        '216',
        '217'#, '217N', '217S'
    )
    
    # Buffer size (in Km) to adjust sample size given nb of legacy sites
    bufferSize_N = 18

    # Buffer size (in Km) to adjust inclusion probability around legacy sites
    bufferSize_p = 10

    # Distance between SSU centroid (in meters)
    ssu_dist = 294

#



# Prepare hexagons

    hexas <- readRDS('data/hexa_complete.RDS') %>%
        filter(propNA <= prop_na) %>%
        filter(ecoregion %in% eco_sim) %>%
        mutate(
            p = (hab_prob * cost_prob) / sum(hab_prob * cost_prob)
        ) %>%
        filter(p != 0)

#



# legacy sites

    # load data
    legacySites <- read_csv('data/SitesLegacy_GRTS20220314.csv') %>%
        group_by(Hexagone_Num) %>%
        summarise(legacyNew = n()) %>%
        rename(ET_Index = Hexagone_Num)

    # merge to hexagons
    hexas <- hexas %>%
        left_join(legacySites) %>%
        mutate(
            legacyNew = replace_na(legacyNew, 0),
            legacy = legacySite + legacyNew
        )

#



# Calculate sample size given number of hexagons and legacy sites

    get_sampleSize <- function(eco, hx, bf_N, sample_e)
    {
        hexa_eco <- hx %>%
            filter(ecoregion == eco) %>%
            st_centroid()

        hexa_legacy_bf <- hexa_eco %>%
            filter(legacy > 0) %>%
            st_buffer(bf_N * 1000) %>%
            st_union()

        nbHexas_legacy <- hexa_eco %>%
            st_intersects(hexa_legacy_bf) %>%
            unlist() %>%
            sum()
        
        round((nrow(hexa_eco) - nbHexas_legacy) * sample_e, 0)
    }


    sampleSize <- map_dbl(
        setNames(eco_sim, paste0('eco_', eco_sim)),
        get_sampleSize,
        hx = hexas,
        bf_N = bufferSize_N,
        sample_e = sample_effort
    )

#
 


# GRTS simulations

    # Sample size by stratum (ecoregion)
    Stratdsgn  <- sampleSize[sampleSize > 0]

    # Prepare sample frame
    sampleFrame <- hexas %>%
        mutate(    
            eco_name = paste0('eco_', ecoregion), # to match design name
            mdcaty  = sum(Stratdsgn) * p/sum(p),
            geometry = sf::st_geometry(sf::st_centroid(geometry))
        )
        
    # Coordinates of all hexagons in matrix format for MBHdesign
    coord_mt <- sampleFrame %>%
        st_coordinates()

    # Coordinates of legacy hexagons
    legacySites <- sampleFrame %>%
        filter(legacy > 0) %>%
        st_coordinates()

    sampleFrame$adj_p <- MBHdesign::alterInclProbs(
        legacy.sites = legacySites,
        potential.sites = coord_mt,
        inclusion.probs = sampleFrame$mdcaty,
        sigma = bufferSize_p * 1000
    )

    # run GRTS    
    grts_out <- list()
    for(Rep in 1:nb_rep)
    {
        out_sample <- spsurvey::grts(
            sframe = sampleFrame,
            n_base = Stratdsgn,
            stratum_var = 'eco_name',
            aux_var = 'adj_p'
        )

        grts_out[[paste0('Rep_', Rep)]] <- out_sample$sites_base$ET_Index
    }

#



# Select cheapest replication

    cheapest_rep <- map_df(
        grts_out,
        ~ hexas %>%
            st_drop_geometry() %>%
            filter(ET_Index %in% .x) %>%
            summarise(totalCost = sum(costSum))
    ) %>%
    pull(totalCost) %>%
    which.min()


    selected_hexas <- hexas %>%
        filter(ET_Index %in% grts_out[[cheapest_rep]])

#



# Select the closest hexagon with the higher inclusion probability

    # Extract neighbours hexagons
    neigh_mt <- hexas %>%
        st_centroid() %>%
        st_intersects(
            y = st_buffer(selected_hexas, dist = 3500),
            sparse = FALSE
        )

    # Remove the focus hexagon (the selected one)
    rr = match(selected_hexas$ET_Index, hexas$ET_Index)
    cc = seq_along(rr)
    neigh_mt[rr + nrow(neigh_mt) * (cc - 1)] <- FALSE


    # Select extra hexagons based on the higher p
    best_p <- apply(
        neigh_mt,
        2,
        function(x)
            hexas$ET_Index[x][which.max(hexas$p[x])]
    )
    
    selected_extra <- hexas[match(best_p, hexas$ET_Index), ]

#



# Select SSU

    # Calculate probability following habitat proportion from ecoregion
    ########################

    land_ca <- raster("data/landcover_ca_30m.tif")
    prev_all <- readRDS('data/prev_all.RDS')    

    # function to generate SSU points (from: https://github.com/dhope/BASSr)
    genSSU <- function(h, spacing)
    {
        ch <- as_tibble(st_coordinates(h))
        top_point <- ch[which.max(ch$Y),]
        bottom_point <- ch[which.min(ch$Y),]
        gridsize <- 2*floor(abs(top_point$Y-bottom_point$Y)/spacing)+3
        rowAngle <- tanh((top_point$X-bottom_point$X)/(top_point$Y-bottom_point$Y))


        cent <- st_centroid(h) %>%
            bind_cols(as_tibble(st_coordinates(.))) %>%
            st_drop_geometry %>%
            dplyr::select(ET_Index, X, Y)


        genRow <- function(cX, cY, sp,...){
            tibble(rowid = seq(-gridsize,gridsize)) %>%
            mutate(X = sin(60*pi/180+rowAngle) *sp*rowid + {{cX}},
                    Y = cos(60*pi/180+rowAngle) *sp*rowid  + {{cY}})
        }

        centroids <- tibble(crowid=seq(-gridsize,gridsize)) %>%
            mutate(cY = cos(rowAngle) *spacing*crowid + cent$Y,
                #spacing/2*crowid + cent$Y,
                cX =  sin(rowAngle) *spacing*crowid + cent$X) %>%
            #cent$X + crowid* sqrt(spacing**2-(spacing/2)**2)) %>%
            rowwise() %>%
            mutate(row = list(genRow(cX = cX,cY = cY,sp = spacing))) %>%
            unnest(row) %>%
            dplyr::select(X,Y) %>%
            st_as_sf(coords = c("X", "Y"), crs = st_crs(h)) %>%
            st_filter(h) %>%
            mutate(
                ET_Index = h$ET_Index,
                ecoregion = h$ecoregion,
                ssuID = row_number()
            )
        return(centroids)
    }


    # function to sample SSU
    sample_SSU <- function(ssuid, prob, geom, filtered, ssuDist, N)
    {
        # check if N is even
        if(N %% 2 != 0)
            stop('`ssu_N` must be a even number.')

        filtered_1 <- filtered

        # loop to sample 4 
        for(i in 1:N)
        {
            # sample point
            assign(
                paste0('sample_', i),
                sample(
                    ssuid[get(paste0('filtered_', i))],
                    size = 1,
                    prob = prob[get(paste0('filtered_', i))]
                )
            )

            # remove points around the first sample for second point
            toKeep <- !st_intersects(
                geom,
                st_buffer(
                    geom[which(get(paste0('sample_', i)) == ssuid)],
                    dist = ssuDist * 2 + ssuDist * 0.1),
                    sparse = FALSE
            )[, 1]

            # update available points
            assign(
                paste0('filtered_', i + 1),
                get(paste0('filtered_', i)) & toKeep
            )
        }

        # create vector of 0, 1 (main), and 2 (over) codes
        # number of main OR over given total N
        n_sub <- N/2
        
        out <- rep(0, length(filtered))
        out[
            ssuid %in% map_dbl(1:n_sub, ~ get(paste0('sample_', .x)))
        ] <- 1      
        out[
            ssuid %in% map_dbl(n_sub + (1:n_sub), ~ get(paste0('sample_', .x)))
        ] <- 2

        return( out )

    }



    # Generate SSU points
    SSUs <- map_dfr(
        seq_len(nrow(selected_hexas)),
        ~ genSSU(selected_hexas[.x, ], spacing = ssu_dist)
    )

    # Buffer of half `ssu_dist` to compute habitat inclusion prob
    SSU_bf <- st_buffer(SSUs, dist = ssu_dist/2)

    # extract pixels for each SSU polygon
    hab_pixels <- exactextractr::exact_extract(
        land_ca,
        SSU_bf,
        progress = FALSE
    )
    rm(SSU_bf)

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
                y = SSUs$ecoregion
            )

    # merge with inclusion probability
    # and calculate inclusion probbaility for each NON empty polygon
    SSUs$incl_prob <- unlist(
            lapply(
                count_hab,
                function(x) {
                    if(is.data.frame(x)) {
                        mg_df <- merge(x, subset(prev_all, ID_ecoregion == x$ecoregion[1]), by.x = "Var1", by.y = "code" , all.x = TRUE)
                        sum(mg_df$Freq * mg_df$incl_prob)
                    }else{
                        NA
                    }
                }
            )
        )
    rm(count_hab)

    SSUs <- SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            incl_prob = incl_prob/sum(incl_prob, na.rm = TRUE)
        ) %>%
        ungroup()


    # Calculate proportion of NA
    SSUs$propNA <- map_dbl(
        hab_pixels,
        ~ sum(is.na(.x$value))/nrow(.x)
    )
    rm(hab_pixels)

    # neighbours matrix
    neighbour_ls <- list()
    for(hx in unique(SSUs$ET_Index))
    {
        ssuhx <- subset(SSUs, ET_Index == hx)

        neighbour_ls[[hx]] <- st_intersects(
            ssuhx,
            st_buffer(ssuhx, dist = ssu_dist + ssu_dist * 0.1),
            sparse = FALSE
        )
    }


    # These are the following rules to a SSU be available to be sampled:
    # - Must have 6 neighbours (less than that means it's a border SSU)
    # - Must have at least 1 - `prop_na` of non empty pixels
    # - 4 out 6 neighbours must also respect the above rule
    SSU_dispo = SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            nbNeighb = map_dbl(
                ssuID,
                ~ sum(neighbour_ls[[unique(ET_Index)]][, .x]) - 1
            ),
            nbPropNA = map_int(
                ssuID,
                .f = function(x, pNA, ETI)
                    sum(
                        pNA[setdiff(which(neighbour_ls[[ETI]][, x]), x)] <= prop_na
                    ),
                pNA = propNA,
                ETI = unique(ET_Index)
            ),
            sampled = sample_SSU(
                ssuid = ssuID,
                prob = incl_prob,
                geom = geometry,
                filtered = propNA <= prop_na & nbNeighb == 6  & nbPropNA >= 4,
                ssuDist = ssu_dist,
                N = ssu_N
            )
        ) %>%
        ungroup()
    
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
