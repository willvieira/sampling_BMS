########################################################################
### Simulations to test how implement information from legacy sites
### Author: Will Vieira
### March 8, 2021
########################################################################


library(tidyverse)
library(sf)
load('data/spatialVectors.rda')
set.seed(0.0)



# Merge all ecoregions into one single sf object

    # load hexagons
    hexa_ls <- list()
    for (id in unique(districts$ECOREGION))
    {
    hexa_ls[[as.character(id)]] <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == id)$REGION_NAM))), '_', id, '/hexa_', id, '.shp'), quiet = TRUE)
    hexa_ls[[as.character(id)]]$ecoregion <- id
    }

    # rbind all shapefiles into one sf oject
    hexas <- do.call(rbind, hexa_ls)
    rm(hexa_ls)

#



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
    sampleSize$N2 <- sampleSize$hexagonProp80 * 0.02
    
    # Sample size in function of size proportion between ecoregions
    sampleSize$N <- unclass(round(sum(sampleSize$N2) * sampleSize$sizeProp, 0))
    
    # Total sample size
    N <- sum(sampleSize$N)

#
 


# Weight p according to the number of legacy site

    we = function(legacySite, threshold, weight)
    {
        if(legacySite >= threshold) {
            return( 2 )
        }else{
            return( (legacySite * (weight - 1)/threshold) + 1 )
        }
    } 

    # weight = 2
    hexas$p_sim2 <- hexas$p * sapply(hexas$legacySite, we, threshold = 15, weight = 2)
    hexas$p_sim2 <- hexas$p_sim2/sum(hexas$p_sim2)

    # weight = 3
    hexas$p_sim3 <- hexas$p * sapply(hexas$legacySite, we, threshold = 15, weight = 3)
    hexas$p_sim3 <- hexas$p_sim3/sum(hexas$p_sim3)

#




# GRTS simulations

    # number of replications for each simulation
    nb_rep <- 50

    # design
    Stratdsgn <- list(
        eco_101 = list(panel = c(PanelOne = subset(sampleSize, ecoregion == 101)$N), over = 0, seltype = "Continuous")
    )

    sample_frame <- subset(hexas, ecoregion %in% 101)

    # Get x and y from centroid of hexagon
    sample_frame$geometry <- sample_frame %>% sf::st_centroid() %>% sf::st_geometry()
    sample_frame[c('X', 'Y')] <- sf::st_coordinates(sample_frame)
    
    # rename ecoregion name to match design name
    sample_frame$eco_name <- paste0('eco_', sample_frame$ecoregion)

    # Get only attributes table (remove spatial information)
    attframe <- sf::st_drop_geometry(sample_frame)


    ################################################
    # Simulation 1: habitat + cost probability only
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% 101)$N) * attframe$p/sum(attframe$p)

    sim1_samp <- list()
    for (i in 1:nb_rep)
    {
        out_sample_sim1 <- spsurvey::grts(design = Stratdsgn,
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

        sim1_samp[[i]] <- out_sample_sim1$ET_Index
    }

    saveRDS(sim1_samp, file = 'sim1_samp.RDS')


    ################################################
    # Simulation 2: habitat + cost probability Weighted by nb of legacy sites
    # Weight = 2
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% 101)$N) * attframe$p_sim2/sum(attframe$p_sim2)

    sim2_samp <- list()
    for (i in 1:nb_rep)
    {
        out_sample_sim2 <- spsurvey::grts(design = Stratdsgn,
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

        sim2_samp[[i]] <- out_sample_sim2$ET_Index
    }

    saveRDS(sim2_samp, file = 'sim2_samp.RDS')

    ################################################
    # Simulation 3: habitat + cost probability Weighted by nb of legacy sites
    # Weight = 3
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% 101)$N) * attframe$p_sim3/sum(attframe$p_sim3)

    sim3_samp <- list()
    for (i in 1:nb_rep)
    {
        out_sample_sim3 <- spsurvey::grts(design = Stratdsgn,
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

        sim3_samp[[i]] <- out_sample_sim3$ET_Index                                       
    }
    
    saveRDS(sim3_samp, file = 'sim3_samp.RDS')


    ################################################
    # Simulation 4: split ecoregion in two according to the north limit of legacy site
    ################################################

    sample_frame <- subset(hexas, ecoregion %in% 101)

    # Split ecoregion in two parts
    coords <- sample_frame %>%
                sf::st_transform(4326) %>%
                sf::st_centroid() %>%
                sf::st_coordinates() %>%
                as.data.frame()

    sample_frameA <- sample_frame[coords$Y >= 50.3, ]
    sample_frameB <- sample_frame[!(coords$Y >= 50.3), ]

    # calculate sample size
    N_a <- as.numeric(round(sum(sampleSize$N2) * st_area(st_union(sample_frameA))/sum(sampleSize$area), 0))
    N_b <- as.numeric(round(sum(sampleSize$N2) * st_area(st_union(sample_frameB))/sum(sampleSize$area), 0))
    
    # design
    Stratdsgn <- list(
        eco_101a = list(panel = c(PanelOne = N_a), over = 0, seltype = "Continuous"),
        eco_101b = list(panel = c(PanelOne = N_b), over = 0, seltype = "Continuous")
    )

    # Get x and y from centroid of hexagon
    sample_frameA$geometry <- sample_frameA %>% sf::st_centroid() %>% sf::st_geometry()
    sample_frameA[c('X', 'Y')] <- sf::st_coordinates(sample_frameA)
    
    sample_frameB$geometry <- sample_frameB %>% sf::st_centroid() %>% sf::st_geometry()
    sample_frameB[c('X', 'Y')] <- sf::st_coordinates(sample_frameB)

    # rename ecoregion name to match design name
    sample_frameA$eco_name <- paste0('eco_', sample_frameA$ecoregion, 'a')
    sample_frameB$eco_name <- paste0('eco_', sample_frameB$ecoregion, 'b')

    # Get only attributes table (remove spatial information)
    sample_frame <- rbind(sample_frameA, sample_frameB)
    attframe <- sf::st_drop_geometry(sample_frame)

    attframe$mdcaty <- sum(N_a, N_b) * attframe$p/sum(attframe$p)

    sim4_samp <- list()
    for (i in 1:nb_rep)
    {
        out_sample_sim4 <- spsurvey::grts(design = Stratdsgn,
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

        sim4_samp[[i]] <- setNames(out_sample_sim4$ET_Index, out_sample_sim4$stratum)
    }

    saveRDS(sim4_samp, file = 'sim4_samp.RDS')

#





###############################################################
# Compare the different simulations
###############################################################

# Code to define zoom in mapview
m <- mapview(hx, zcol = 'sel')
mapshot(m, file = "italy.png", vwidth = 700, vheight = 744)


# Load simulations
sim1 <- readRDS('sim1_samp.RDS')
sim2 <- readRDS('sim2_samp.RDS')
sim3 <- readRDS('sim3_samp.RDS')
sim4 <- readRDS('sim4_samp.RDS')



# Extract hexagon info from selected samples
sumary_df <- data.frame()
prob_sample_df <- data.frame()
for(sim in paste0('sim', 1:4))
{
    # get hexagons ID from simulation
    sim_i <- get(sim)
    
    for(Rep in 1:nb_rep)
    {
        # get selected hexagons only
        if(sim == 'sim4')
        {
            hexagons_ID <- sim_i[[Rep]][which(names(sim_i[[Rep]]) %in% 'eco_101a')]
        }else{
            hexagons_ID <- sim_i[[Rep]]
        }

        hexas_rep <- subset(hexas, ET_Index %in% hexagons_ID)

        # save probability of inclusion of selected hexagons
        prob_sample_df <- rbind(prob_sample_df,
                                data.frame(sim = rep(sim, nrow(hexas_rep)),
                                           prob = hexas_rep[[ifelse(sim %in% c('sim1', 'sim4'), 'p', ifelse(sim == 'sim2', 'p_sim2', 'p_sim3'))]]))

        # total legacy site poitns within selected hexagons
        totalLegacy <- sum(hexas_rep$legacySite)

        # Total of hexagons with 15 or more legacy sites
        filledHexas <- sum(hexas_rep$legacySite >= 15)
        
        # remove hexagons with 15 >= legacySites
        hexas_rep <- subset(hexas_rep, legacySite < 15)
        
        # total cost of ecoregion
        totalCost <- sum(hexas_rep$costSum)

        sumary_df <- rbind(sumary_df,
                        data.frame(sim = sim,
                                   Rep = Rep,
                                   totalLegacy = totalLegacy,
                                   filledHexas = filledHexas,
                                   totalCost = totalCost))
    }
}


library(ggplot2)
library(viridis)
library(ggpubr)

# Total cost of ecoregion
p1 <- ggplot(sumary_df, aes(x = totalCost, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Coût total de l'écorégion") +
                scale_fill_viridis(discrete=TRUE) +
                theme_classic()


# Total number of legacy sites
p2 <- ggplot(sumary_df, aes(x = totalLegacy, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Nombre total de station d'ecoute") +
                scale_fill_viridis(discrete=TRUE) +
                theme_classic()

# Total number of hexagons above threshold
p3 <- ggplot(sumary_df, aes(x = filledHexas, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Nombre d'hexagon avec 15 ou plus station d'ecoute") +
                ylim(0, 1) +
                scale_fill_viridis(discrete=TRUE) +
                theme_classic()


# Frequency
ggarrange(p1, p2, p3, common.legend = TRUE)


# Inclusion probability
hexas_eco <- subset(hexas, ecoregion == 101)

prob_sim1 <- rbind(subset(prob_sample_df, sim == 'sim1'),
                          data.frame(sim = rep('All eco', nrow(hexas_eco)),
                                     prob = hexas_eco$p))
prob_sim2 <- rbind(subset(prob_sample_df, sim == 'sim2'),
                          data.frame(sim = rep('All eco', nrow(hexas_eco)),
                                     prob = hexas_eco$p_sim2))
prob_sim3 <- rbind(subset(prob_sample_df, sim == 'sim3'),
                          data.frame(sim = rep('All eco', nrow(hexas_eco)),
                                     prob = hexas_eco$p_sim3))
prob_sim4 <- rbind(subset(prob_sample_df, sim == 'sim4'),
                          data.frame(sim = rep('All eco', nrow(hexas_eco)),
                                     prob = hexas_eco$p))


p1 <- ggplot(subset(prob_sim1, prob < 0.002), aes(x = prob, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Probabilité d'inclusion (hab + coût)") +
                ggtitle('sim1') +
                scale_fill_viridis(discrete=TRUE) +
                scale_color_hue(labels = c("Prob ecoregion", "Prob selection")) +
                theme_classic()

p2 <- ggplot(subset(prob_sim2, prob < 0.002), aes(x = prob, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Probabilité d'inclusion (hab + coût)") +
                ggtitle('sim2') +
                scale_fill_viridis(discrete=TRUE) +
                scale_color_hue(labels = c("Prob ecoregion", "Prob selection")) +
                theme_classic()

p3 <- ggplot(subset(prob_sim3, prob < 0.002), aes(x = prob, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Probabilité d'inclusion (hab + coût)") +
                ggtitle('sim3') +
                scale_fill_viridis(discrete=TRUE) +
                scale_color_hue(labels = c("Prob ecoregion", "Prob selection")) +
                theme_classic()

p4 <- ggplot(subset(prob_sim4, prob < 0.002), aes(x = prob, fill = sim)) +
                geom_density(alpha = 0.6, color = NA) +
                xlab("Probabilité d'inclusion (hab + coût)") +
                ggtitle('sim4') +
                scale_fill_viridis(discrete=TRUE) +
                scale_color_hue(labels = c("Prob ecoregion", "Prob selection")) +
                theme_classic()

ggarrange(p1, p2, p3, p4, common.legend = TRUE)
