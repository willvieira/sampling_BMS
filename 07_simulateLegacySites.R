########################################################################
### Simulations to test how implement information from legacy sites
### Author: Will Vieira
### March 8, 2021
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

    # SIMULATION 3:
    # Remove the already filled hexagons with 4 or more legacy sites
    hexas$p_sim3 <- hexas$p
    hexas$p_sim3[which(hexas$legacySite >= 4)] <- 0
    hexas$p_sim3 <- hexas$p_sim3/sum(hexas$p_sim3)

    # SIMULATION 4:
    # Remove the already filled hexagons with at least 1 legacy site
    hexas$p_sim4 <- hexas$p
    hexas$p_sim4[which(hexas$legacySite > 0)] <- 0
    hexas$p_sim4 <- hexas$p_sim4/sum(hexas$p_sim4)
 
#




# GRTS simulations

    # number of replications for each simulation
    nb_rep <- 25

    # Ecoregions
    eco_sim <- c('72', '74', '75', '76',
                 '96', '96N', '96S',
                 '99',
                 '100', '100N', '100S',
                 '101', '101N', '101S',
                 '102',
                 '103', '103N', '103S',
                 '117',
                 '217', '217N', '217S')

    # Design list    
    Stratdsgn <- list()
    for(eco in eco_sim)
        Stratdsgn[[paste0('eco_', eco)]] <- list(panel = c(PanelOne = subset(sampleSize, ecoregion == eco)$N), over = 0, seltype = 'Continuous')

    # Sample frame for specific ecoregions
    sample_frame <- subset(hexas, ecoregion %in% eco_sim)

    # Get x and y from centroid of hexagon
    sample_frame$geometry <- sample_frame %>% sf::st_centroid() %>% sf::st_geometry()
    sample_frame[c('X', 'Y')] <- sf::st_coordinates(sample_frame)
    
    # rename ecoregion name to match design name
    sample_frame$eco_name <- paste0('eco_', sample_frame$ecoregion)

    # Get only attributes table (remove spatial information)
    attframe <- sf::st_drop_geometry(sample_frame)




    ################################################
    # Simulation 0 AND 1: habitat + cost probability only
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% eco_sim)$N) * attframe$p/sum(attframe$p)

    set.seed(1)
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

        sim1_samp[[i]] <- setNames(out_sample_sim1$ET_Index, out_sample_sim1$panel)
    }

    saveRDS(sim1_samp, file = 'sim1_samp.RDS')


    ################################################
    # Simulation 2: habitat + cost probability Weighted by nb of legacy sites
    # Weight = 3
    ################################################

    attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% eco_sim)$N) * attframe$p_sim2/sum(attframe$p_sim2)

    set.seed(2)
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

        sim2_samp[[i]] <- setNames(out_sample_sim2$ET_Index, out_sample_sim2$panel)
    }

    saveRDS(sim2_samp, file = 'sim2_samp.RDS')



    ################################################
    # Simulation 3 and 4: Sample size reduced by the number of legacy sites
    # - 3: Reduce sample size by the number of hexagons with at least 4 legacy sites
    # - 4: Reduce sample size by the total number of legacy site for the ecoregion: nbHexagon = nbLegacySite/10
    ################################################

    # Recalculate sample size as a function of number of legacy site
    # And filter to ecoregions that have N > 0
    sampleSize_sim3 <- hexas %>%
                            sf::st_drop_geometry() %>%
                            filter(ecoregion %in% eco_sim) %>%
                            group_by(ecoregion) %>%
                            summarise(nbHexagon = sum(legacySite >= 4), # number of hexagons that already have the minimum sample size
                                      nb_hexlegacySite = round(sum(legacySite)/10, 0)) %>% # number of hexagons if we consider the total number of legacy sites (independent of their aggregation)
                            inner_join(y = sampleSize[, c('ecoregion', 'N')], by = 'ecoregion') %>%
                            mutate(N_sim3 = if_else(nbHexagon > N, 0, N - nbHexagon),
                                   N_sim4 = if_else(nb_hexlegacySite > N, 0, N - nb_hexlegacySite)) %>%
                            as.data.frame()

    # Variation A
    #=======================
    
    # Design list    
    Stratdsgn <- list()
    for(eco in subset(sampleSize_sim3, N_sim3 > 0)$ecoregion)
        Stratdsgn[[paste0('eco_', eco)]] <- list(panel = c(PanelOne = subset(sampleSize_sim3, ecoregion == eco)$N_sim3), over = 0, seltype = 'Continuous')

    attframe$mdcaty <- sum(subset(sampleSize_sim3, ecoregion %in% eco_sim)$N_sim3) * attframe$p_sim3/sum(attframe$p_sim3)

    set.seed(3)
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

        sim3_samp[[i]] <- setNames(out_sample_sim3$ET_Index, out_sample_sim3$panel)
    }

    saveRDS(sim3_samp, file = 'sim3_samp.RDS')


    # Variation B (sim4)
    #=======================

    # Design list    
    Stratdsgn <- list()
    for(eco in subset(sampleSize_sim3, N_sim4 > 0)$ecoregion)
        Stratdsgn[[paste0('eco_', eco)]] <- list(panel = c(PanelOne = subset(sampleSize_sim3, ecoregion == eco)$N_sim4), over = 0, seltype = 'Continuous')

    attframe$mdcaty <- sum(subset(sampleSize_sim3, ecoregion %in% eco_sim)$N_sim4) * attframe$p_sim4/sum(attframe$p_sim4)

    set.seed(4)
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

        sim4_samp[[i]] <- setNames(out_sample_sim4$ET_Index, out_sample_sim4$panel)
    }

    saveRDS(sim4_samp, file = 'sim4_samp.RDS')







###############################################################
# Compare the different simulations
###############################################################

# Load simulations
sim1 <- readRDS('sim1_samp.RDS')
sim2 <- readRDS('sim2_samp.RDS')
sim3 <- readRDS('sim3_samp.RDS')
sim4 <- readRDS('sim4_samp.RDS')



# Extract hexagon info from selected samples
summary_df <- data.frame()
prob_sample_df <- data.frame()

for(sim in 1:4)
{
    # get hexagons ID from simulation
    sim_i <- get(paste0('sim', sim))
    
    for(Rep in 1:nb_rep)
    {
        # get selected hexagons only
        hexagons_ID <- sim_i[[Rep]]

        hexas_rep <- subset(hexas, ET_Index %in% hexagons_ID)

        # save probability of inclusion of selected hexagons
        prob_sample_df <- rbind(prob_sample_df,
                                data.frame(Simulation = paste0('H', sim),
                                           Rep = Rep,
                                           ecoregion = hexas_rep$ecoregion,
                                           prob = hexas_rep$p))


        # total legacy site points within selected hexagons
        # Total of hexagons with 4 or more legacy sites
        legacy_summ <- hexas_rep %>%
            sf::st_drop_geometry() %>%
            group_by(ecoregion) %>%
            summarise(
                totalLegacy = sum(legacySite),
                filledHexas = sum(legacySite >= 4)
            )

        # total cost of ecoregion
        # cost by hexagon
        sim_summ <- hexas_rep %>%
            sf::st_drop_geometry() %>%
            filter(legacySite < 4) %>%
            group_by(ecoregion) %>%
            summarise(
                totalCost = sum(costSum),
                costHex = totalCost/n(),
            ) %>%
            inner_join(legacy_summ, by = 'ecoregion') %>%
            as.data.frame()

        sim_summ$Rep = Rep
        sim_summ$Simulation = paste0('H', sim)

        # append to main data.frame
        summary_df <- rbind(summary_df, sim_summ)
    }        
}

# Calculate cost for H0 (no legacy site considered: using H1 but not filtering to remove already filled hexagons)
for(Rep in 1:nb_rep)
{
    # get selected hexagons only
    hexagons_ID <- sim1[[Rep]]

    hexas_rep <- subset(hexas, ET_Index %in% hexagons_ID)

    # total legacy site points within selected hexagons
    # Total of hexagons with 4 or more legacy sites
    legacy_summ <- hexas_rep %>%
        sf::st_drop_geometry() %>%
        group_by(ecoregion) %>%
        summarise(
            totalLegacy = sum(legacySite),
            filledHexas = sum(legacySite >= 4)
        )

    # total cost of ecoregion
    # cost by hexagon
    sim_summ <- hexas_rep %>%
        sf::st_drop_geometry() %>%
        group_by(ecoregion) %>%
        summarise(
            totalCost = sum(costSum),
            costHex = totalCost/n(),
        ) %>%
        inner_join(legacy_summ, by = 'ecoregion') %>%
        as.data.frame()

    sim_summ$Rep = Rep
    sim_summ$Simulation = 'H0'

    # append to main data.frame
    summary_df <- rbind(summary_df, sim_summ)
}


# Add inclusion probability of already filled hexagons for H3 and H4

# H3
hexas_H3 <- hexas %>%
    st_drop_geometry() %>%
    filter(ecoregion %in% eco_sim) %>%
    filter(legacySite >= 4)

prob_sample_df <- rbind(prob_sample_df,
                        data.frame(Simulation = 'H3',
                                    Rep = 0,
                                    ecoregion = hexas_H3$ecoregion,
                                    prob = hexas_H3$p))

# H4
hexas_H4 <- hexas %>%
    st_drop_geometry() %>%
    filter(ecoregion %in% eco_sim) %>%
    filter(legacySite > 0)

prob_sample_df <- rbind(prob_sample_df,
                        data.frame(Simulation = 'H4',
                                    Rep = 0,
                                    ecoregion = hexas_H4$ecoregion,
                                    prob = hexas_H4$p))





# Extract the frequency a hexagon was selected among replications along with their inclusion probability
freq_df <- data.frame()
for(sim in 1:4)
{
    # get hexagons ID from simulation
    sim_i <- get(paste0('sim', sim))

    # frequence
    sim_freq <- table(unlist(sim_i))

    # Get row index of selected hexagons
    hexas_i <- match(names(sim_freq), hexas$ET_Index)

    # Extract and append hexagons info into df
    freq_df <- rbind(freq_df,
                        data.frame(
                            Simulation = paste0('H', sim),
                            ecoregion = hexas$ecoregion[hexas_i],
                            Freq = as.numeric(sim_freq),
                            p = hexas$p[hexas_i],
                            p_sim = hexas[[ifelse(sim == 1, 'p', ifelse(sim == 3, 'p_sim2', 'p_sim3'))]][hexas_i],
                            legacySite = hexas$legacySite[hexas_i])
                    )
}



# get habitat info from the whole ecoregion and save it as a repetition nb_rep + 2
hab_df <- hexas %>%
    st_drop_geometry() %>%
    group_by(ecoregion) %>%
    summarise(
        across(starts_with('land_ca_'), ~sum(.x, na.rm = TRUE))
    ) %>%
    pivot_longer(!ecoregion, names_to = 'land', values_to = 'prop') %>%
    group_by(ecoregion) %>%
    mutate(
        prop = prop/sum(prop),
        Rep = nb_rep + 2,
        Simulation = 'H1'
    ) %>%
    ungroup() %>%
    as.data.frame()


# Generate a copy of hab_df for all simulations
hab_df_i <- hab_df
for(i in 2:4)
{
    hab_df_i$Simulation = paste0('H', i)
    hab_df <- rbind(hab_df, hab_df_i)
}


# Extract habitat info from selected hexagons
for(sim in 1:4)
{
    # get hexagons ID from simulation
    sim_i <- get(paste0('sim', sim))
    
    # if H3 or H4, add habitat of the already filled hexagons
    if(sim == 3) {
        hex_toAdd <- hexas %>%
            st_drop_geometry() %>%
            filter(ecoregion %in% eco_sim) %>%
            filter(legacySite >= 4) %>%
            select(ET_Index)
        sim_i <- lapply(sim_i, append, hex_toAdd$ET_Index)
    }else if(sim == 4) {
         hex_toAdd <- hexas %>%
            st_drop_geometry() %>%
            filter(ecoregion %in% eco_sim) %>%
            filter(legacySite > 0) %>%
            select(ET_Index)
        sim_i <- lapply(sim_i, append, hex_toAdd$ET_Index)       
    }

    for(Rep in 1:nb_rep)
    {
 
        # save habitat proportion
        land_dt <- hexas %>%
            st_drop_geometry() %>%
            filter(ET_Index %in% sim_i[[Rep]]) %>%
            group_by(ecoregion) %>%
            summarise(
                across(starts_with('land_ca_'), ~sum(.x, na.rm = TRUE))
            ) %>%
            pivot_longer(!ecoregion, names_to = 'land', values_to = 'prop') %>%
            group_by(ecoregion) %>%
            mutate(
                prop = prop/sum(prop),
                Rep = Rep,
                Simulation = paste0('H', sim)
            ) %>%
            ungroup() %>%
            as.data.frame()

        # Save
        hab_df <- rbind(hab_df, land_dt)
    }
}



# Ripley's K
ripleysK <- data.frame()
Count <- 1
for(eco in eco_sim)
{
    # Get Window (ecoregion's boundary)
    owinWindow_eco <- as.owin(sf::st_union(subset(hexas, ecoregion == eco)))

    for(sim in 1:4)
    {
        sim_i <- get(paste0('sim', sim))

        # if H3 or H4, add habitat of the already filled hexagons
        if(sim == 3) {
            hex_toAdd <- hexas %>%
                st_drop_geometry() %>%
                filter(ecoregion %in% eco) %>%
                filter(legacySite >= 4) %>%
                select(ET_Index)
            sim_i <- lapply(sim_i, append, hex_toAdd$ET_Index)
        }else if(sim == 4) {
            hex_toAdd <- hexas %>%
                st_drop_geometry() %>%
                filter(ecoregion %in% eco) %>%
                filter(legacySite > 0) %>%
                select(ET_Index)
            sim_i <- lapply(sim_i, append, hex_toAdd$ET_Index)       
        }

        for(Rep in 1:nb_rep)
        {
            # Get coordinates of selected points
            coords <- hexas %>%
                        filter(ET_Index %in% sim_i[[Rep]]) %>%
                        filter(ecoregion == eco) %>%
                        sf::st_centroid() %>%
                        sf::st_coordinates() %>%
                        as.data.frame()

            # Transform in a point pattern obj
            sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
            
            # Get K
            K_sim_rep <- as.data.frame(Kest(sample_ppp, correction = 'iso'))
            
            K_sim_rep$ecoregion <- eco
            K_sim_rep$Simulation <- paste0('H', sim)
            K_sim_rep$Rep <- Rep

            ripleysK <- rbind(ripleysK, K_sim_rep)

            # progress
            cat('   Simulation ', round(Count/(length(eco_sim) + 4 + nb_rep) * 100, 0), '%\r')
            Count <- Count + 1
        }
    }
}





# PLOT
####################

library(ggplot2)
library(viridis)
library(ggpubr)

habitat_colors <- c(land_ca_1 = rgb(0, 61, 0, maxColorValue = 255),
                    land_ca_2 = rgb(148, 156, 112, maxColorValue = 255),
                    land_ca_5 = rgb(20, 140, 61, maxColorValue = 255),
                    land_ca_6 = rgb(92, 117, 43, maxColorValue = 255),
                    land_ca_8 = rgb(179, 138, 51, maxColorValue = 255),
                    land_ca_10 = rgb(225, 207, 138, maxColorValue = 255),
                    land_ca_11 = rgb(156, 117, 84, maxColorValue = 255),
                    land_ca_12 = rgb(186, 212, 143, maxColorValue = 255),
                    land_ca_13 = rgb(64, 138, 112, maxColorValue = 255),
                    land_ca_14 = rgb(107, 163, 138, maxColorValue = 255),
                    land_ca_16 = rgb(168, 171, 174, maxColorValue = 255),
                    land_ca_21 = '#03045e',
                    land_ca_22 = '#023e8a',
                    land_ca_23 = '#0077b6',
                    land_ca_24 = '#0096c7',
                    land_ca_25 = '#00b4d8',
                    land_ca_26 = '#48cae4',
                    land_ca_27 = '#90e0ef',
                    land_ca_28 = '#ade8f4',
                    land_ca_29 = '#caf0f8',
                    land_ca_30 = '#E6FFFF')


pdf('sampling_by_ecoregion.pdf', width = 10, height = 12)
for(eco in eco_sim)
{
    summary_eco <- subset(summary_df, ecoregion == eco)
    prob_sample_eco <- subset(prob_sample_df, ecoregion == eco)
    freq_eco <- subset(freq_df, ecoregion == eco)

    # Cost per hexagon
    p1 <- ggplot(summary_eco, aes(x = costHex, fill = Simulation)) +
                    geom_density(alpha = 0.6, color = NA) +
                    xlab("Coût par hexagone") +
                    scale_fill_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic()

    # Total cost of ecoregion
    p2 <- ggplot(summary_eco, aes(x = totalCost, fill = Simulation)) +
                    geom_density(alpha = 0.6, color = NA) +
                    xlab("Coût total de l'écorégion") +
                    scale_fill_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic()

    # Total number of legacy sites
    p3 <- ggplot(subset(summary_eco, Simulation %in% c('H1', 'H2')), aes(x = totalLegacy, fill = Simulation)) +
                    geom_density(alpha = 0.6, color = NA) +
                    xlab("Nombre total de station d'ecoute") +
                    scale_fill_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic() +
                    annotate(
                        'text',
                        x = Inf, y = Inf,
                        label = paste("Echantillon station d'ecoute:", subset(sampleSize_sim3, ecoregion == eco)$N * 10),
                        vjust = 1, hjust = 1)

    # Total number of hexagons above threshold
    p4 <- ggplot(subset(summary_eco, Simulation %in% c('H1', 'H2')), aes(x = filledHexas, fill = Simulation)) +
                    geom_density(alpha = 0.6, color = NA) +
                    xlab("Nombre d'hexagon avec 4 ou plus station d'ecoute") +
                    scale_fill_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic() +
                    annotate(
                        'text',
                        x = Inf, y = Inf,
                        label = paste("Echantillon hexagone:", subset(sampleSize_sim3, ecoregion == eco)$N),
                        vjust = 1, hjust = 1)

    # Inclusion probability without weight of legacy sites
    p5 <- ggplot(subset(prob_sample_eco, prob < 0.0018 & Simulation %in% paste0('H', 1:4)), aes(x = prob, fill = Simulation)) +
                    geom_density(alpha = 0.6, color = NA) +
                    xlab("Probabilité d'inclusion p (sans poids de station d'ecoute)") +
                    scale_fill_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic()

    # Inclusion probability in function of number of times a hexagon was choosen
    p6 <- ggplot(subset(freq_eco, p < 0.0018 & Simulation %in% paste0('H', 1:4)), aes(as.factor(Freq), p, fill = Simulation, color = Simulation)) +
                    geom_boxplot(outlier.size = 0.1, lwd = 0.2) +
                    xlab("Fréquence un hexagone a été choisi parmi 50 répétitions") +
                    ylab("Probabilité d'inclusion p") +
                    scale_fill_viridis(discrete = TRUE, , limits = c('H1', 'H2', 'H3', 'H4'), alpha = 0.6) +
                    scale_color_viridis(discrete = TRUE, limits = c('H0', 'H1', 'H2', 'H3', 'H4')) +
                    theme_classic()

    print(
        annotate_figure(ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3, common.legend = TRUE, legend = 'right'),
                        top = text_grob(paste('Ecoregion', eco), size = 14))
    )
}
dev.off()


pdf('sampling_habitat.pdf', width = 10, height = 12)
for(eco in eco_sim)
{
    hab_eco <- subset(hab_df, ecoregion == eco)

    # Habitat frequence of choosen hexagons
    p1 <- ggplot(hab_eco, aes(y = Rep, x = prop)) +
                    geom_col(aes(fill = land), orientation = 'y', width = 1) +
                    scale_fill_manual(values = habitat_colors) +
                    xlab("Proportion") +
                    ylab("Repetition") +
                    facet_grid(~Simulation) +
                    theme_classic()
    
    print(
        annotate_figure(ggarrange(p1, ncol = 1, nrow = 1, common.legend = TRUE, legend = 'bottom'),
                        top = text_grob(paste('Ecoregion', eco), size = 14))
    )
}
dev.off()


# Get max and min of each simulation
ripley_range <- ripleysK %>%
    group_by(ecoregion, Simulation, r) %>%
    summarise(
        mean_iso = mean(iso),
        min_iso = min(iso),
        max_iso = max(iso)
    ) %>%
    as.data.frame()

# Simulation colors
sim_col <- viridis::viridis_pal(alpha = 1)(4)
sim_col_alpha <- viridis::viridis_pal(alpha = 0.5)(4)


pdf('sampling_RipleysK.pdf', width = 8, height = 6)
for(eco in eco_sim)
{
    theo_ripleys <- subset(ripleysK, ecoregion == eco & Simulation == 'H1' & Rep == 1)
    range_eco <- subset(ripley_range, ecoregion == eco)
    xLim <- c(0, 60000)
    yLim <- range_eco %>%
        filter(r < xLim[2]) %>%
        summarise(range(c(min_iso, max_iso)))

    par(mar = c(2.5, 2.5, 1, 0.5), oma = c(0, 0, 1, 0), mgp = c(1.4, 0.2, 0), tck = -.008)
    plot(0, pch = '', xlab = 'r', ylab = 'K(r)', xlim = xLim, ylim = yLim[, 1])
    for(sim in 1:4)
    {       
        range_sim <- subset(range_eco, Simulation == paste0('H', sim))
        polygon(c(range_sim$r, rev(range_sim$r)), c(range_sim$min_iso, rev(range_sim$max_iso)),
                col = sim_col_alpha[sim],
                border = NA)
        points(range_sim$r, range_sim$mean_iso, type = 'l', lwd = 1.8, col = sim_col[sim])
    }
    points(theo_ripleys$r, theo_ripleys$theo, type = 'l', lty = 2, col = rgb(0, 0, 0, 0.8), lwd = 2)
    mtext(paste('Ecoregion', eco), 3, line = 0)
    legend('topleft', legend = c(paste0('H', 1:4), 'Theorique'), lty = c(rep(1, 4), 2), col = c(sim_col, 'black'), bty = 'n')
}
dev.off()







################################################################################################
# Simulate the effect of increasing weight given by legacy Site on cost and habitat representatives
################################################################################################

# Folder to save sim output
dir.create('sim_weight')


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


# Define simulation parameters, ecoregion and their respective seed

    p_weight <- seq(2, 40, 2)
    eco_sim <- c('99', '102', '117', '100S', '101S')
    c <- 15

    sim_df <- expand.grid(
        ecoregion = eco_sim,
        weight = p_weight,
        Rep = 1:nb_rep
    )

    set.seed(0)
    sim_df$seed <- sample(1:5000, nrow(sim_df))

#


# Function to weight p according to the number of legacy site

    we = function(legacySite, lower, upper, mid, beta)
    {
        if(legacySite == 0) {
            return( 1 )
        }else{
            return( lower + ((upper - lower) / (1 + exp(-beta * (legacySite - mid)))) )
        }
    } 

#


# Function to get selected hexagons given the row id (to match with `sim_df`)
    run_simulation <- function(id)
    {
        print(paste('Running for id', id))

        hexas_id <- subset(hexas, ecoregion == sim_df[id, 'ecoregion'])
        
        # Calculate new probability of inclusion for simulation ID
        hexas_id$p_sim <- hexas_id$p *
                        sapply(hexas_id$legacySite, we, lower = 1, upper = sim_df[id, 'weight'], mid = 4, beta = 0.9)

        # Design list    
        Stratdsgn <- list()
        for(eco in sim_df[id, 'ecoregion'])
            Stratdsgn[[paste0('eco_', eco)]] <- list(panel = c(PanelOne = subset(sampleSize, ecoregion == eco)$N), over = 0, seltype = 'Continuous')
        

        sample_frame <- hexas_id

        # Get x and y from centroid of hexagon
        sample_frame$geometry <- sample_frame %>%
                    sf::st_centroid() %>%
                    sf::st_geometry()
        sample_frame[c('X', 'Y')] <- sf::st_coordinates(sample_frame)
        
        # rename ecoregion name to match design name
        sample_frame$eco_name <- paste0('eco_', sample_frame$ecoregion)

        # Get only attributes table (remove spatial information)
        attframe <- sf::st_drop_geometry(sample_frame)

        # standardize p by N
        attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% sim_df[id, 'ecoregion'])$N) * attframe$p_sim/sum(attframe$p_sim)

        # Run GRTS with specified seed
        set.seed(sim_df[id, 'seed'])
        out <- spsurvey::grts(design = Stratdsgn,
                                DesignID = "ET_ID",
                                type.frame = "finite",
                                src.frame = "att.frame",
                                att.frame = attframe,
                                xcoord = 'X',
                                ycoord = 'Y',
                                stratum = "eco_name",
                                mdcaty = "mdcaty",
                                shapefile = FALSE)
        # Save hexagons ID 
        saveRDS(
            out$ET_Index,
            file = paste0(
                'sim_weight/sim_',
                sim_df[id, 'ecoregion'], '_',
                sim_df[id, 'weight'], '_',
                sim_df[id, 'Rep'], '.RDS'
            )
        )
    }

#


# Run
sapply(1:1500, run_simulation)





# Plot
#######################

# Load simulations and get info from selected hexagons
out_df <- data.frame()
for(i in 1:nrow(sim_df))
{
    # load selected hexagons
    ETIndex_id <- readRDS(
        paste0(
            'sim_weight/sim_',
            sim_df[i, 'ecoregion'], '_',
            sim_df[i, 'weight'], '_',
            sim_df[i, 'Rep'], '.RDS'
        )
    )
    hexas_id <- subset(hexas, ET_Index %in% ETIndex_id)

    # bind hexas info to main data.frame
    out_df <- rbind(
        out_df,
        data.frame(
            ecoregion = sim_df[i, 'ecoregion'],
            weight = sim_df[i, 'weight'],
            Rep = sim_df[i, 'Rep'],
            cost = hexas_id$costSum,
            p = hexas_id$p,
            hab_prob = hexas_id$hab_prob,
            cost_prob = hexas_id$cost_prob,
            legacySite = hexas_id$legacySite
        )
    )
}


# load simulations and calculate Ripley's K 

# Ripley's K
ripleysK <- data.frame()
Count <- 1
for(eco in eco_sim)
{
    # Get Window (ecoregion's boundary)
    owinWindow_eco <- as.owin(sf::st_union(subset(hexas, ecoregion == eco)))

    for(wd in p_weight)
    {
        for(Rep in 1:nb_rep)
        {
            sim_i <- readRDS(paste0('sim_weight/sim_', eco, '_', wd, '_', Rep, '.RDS'))

            # Get X and Y coordinates for selected hexagons
            coords <- hexas %>%
                    filter(ET_Index %in% sim_i) %>%
                    sf::st_centroid() %>%
                    sf::st_coordinates() %>%
                    as.data.frame()

            # Transform in a point pattern obj
            sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
            
            # Get K
            K_sim_rep <- as.data.frame(Kest(sample_ppp, correction = 'iso'))
            
            K_sim_rep$ecoregion <- eco
            K_sim_rep$weight <- wd
            K_sim_rep$Rep <- Rep
            K_sim_rep$test <- 1

            ripleysK <- rbind(ripleysK, K_sim_rep, K_sim_rep2)

            # progress
            cat('   Simulation ', round(Count/nrow(sim_df) * 100, 0), '%\r')
            Count <- Count + 1

        }
    }
}





library(ggplot2)
library(ggpubr)

# FIGURE 1
###############################

p1 <- out_df %>%
    filter(legacySite < 4) %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(totalCost = sum(cost)) %>%
    ggplot(aes(x = as.factor(weight), y = totalCost, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab('Coût totale')

text_df <- sampleSize %>%
    filter(ecoregion %in% eco_sim) %>%
    select(ecoregion, N) %>%
    mutate(
        label_hex = paste('Sample size:', N, 'hexas'),
        label_station = paste('Sample size:', N * 10, 'stations')
    )

p2 <- out_df %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(totalHexagons = sum(legacySite >= 4)) %>%
    ggplot(aes(x = as.factor(weight), y = totalHexagons, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab("Nombre d'hexagons avec 4 ou plus stations d'ecoute") +
        geom_text(
            data = text_df,
            mapping = aes(x = Inf, y = -Inf, label = label_hex),
            hjust = 1,
            vjust = -0.5
        )

p3 <- out_df %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(totalLegacy = sum(legacySite)) %>%
    ggplot(aes(x = as.factor(weight), y = totalLegacy, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab("Station d'ecoute totale") +
        geom_text(
            data = text_df,
            mapping = aes(x = Inf, y = -Inf, label = label_station),
            hjust = 1,
            vjust = -0.5
        )

pdf(file = 'simulation_legacy.pdf', width = 20, height = 12)
    print(ggarrange(p1, p2, p3, ncol = 1, nrow = 3))
dev.off()


# FIGURE 2
###############################

p1 <- out_df %>%
    filter(legacySite < 4) %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(meanP = mean(cost_prob)) %>%
    ggplot(aes(x = as.factor(weight), y = meanP, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab("Probabilité d'inclusion coût (cost_prob)")


p2 <- out_df %>%
    filter(legacySite < 4) %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(meanP = mean(hab_prob)) %>%
    ggplot(aes(x = as.factor(weight), y = meanP, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab("Probabilité d'inclusion habitat (hab_prob)")

p3 <- out_df %>%
    filter(legacySite < 4) %>%
    group_by(ecoregion, weight, Rep) %>%
    summarise(meanP = mean(p)) %>%
    ggplot(aes(x = as.factor(weight), y = meanP, ecoregion)) +
        geom_boxplot(fill = 'grey') +
        geom_jitter(alpha = 0.5, size = 0.8) +
        facet_wrap(~ ecoregion, scales = 'free_y', nrow = 1) +
        theme_classic() +
        xlab('Poids (effet max)') +
        ylab("Probabilité d'inclusion p (cost_prob * hab_prob)")

pdf(file = 'simulation_legacy2.pdf', width = 20, height = 12)
    print(ggarrange(p1, p2, p3, ncol = 1, nrow = 3))
dev.off()



# FIGURE 3
###############################

# Get max and min of each simulation
ripley_range <- ripleysK %>%
    group_by(ecoregion, weight, r) %>%
    summarise(
        mean_iso = mean(iso),
        min_iso = min(iso),
        max_iso = max(iso)
    ) %>%
    as.data.frame()


# Simulation colors
sim_col <- viridis::viridis_pal(alpha = 1)(length(p_weight))
sim_col_alpha <- viridis::viridis_pal(alpha = 0.2)(length(p_weight))


pdf('simulation_legacy_RipleysK.pdf', width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(2.5, 2.5, 1, 0.5), oma = c(0, 0, 1, 0), mgp = c(1.4, 0.2, 0), tck = -.008)
for(eco in eco_sim)
{

    theo_ripleys <- subset(ripleysK, ecoregion == eco & weight == 2 & Rep == 1)
    range_eco <- subset(ripley_range, ecoregion == eco)
    xLim <- c(0, 50000)
    yLim <- range_eco %>%
        filter(r < xLim[2]) %>%
        summarise(range(c(min_iso, max_iso), na.rm = TRUE))

    plot(0, pch = '', xlab = 'r', ylab = 'K(r)', xlim = xLim, ylim = yLim[, 1])
    for(wd in p_weight)
    {       
        range_sim <- subset(range_eco, weight == wd)
        # polygon(c(range_sim$r, rev(range_sim$r)), c(range_sim$min_iso, rev(range_sim$max_iso)),
        #         col = sim_col_alpha[which(wd == p_weight)],
        #         border = NA)
        points(range_sim$r, range_sim$mean_iso, type = 'l', lwd = 1.8, col = sim_col[which(wd == p_weight)])
    }
    points(theo_ripleys$r, theo_ripleys$theo, type = 'l', lty = 2, col = rgb(0, 0, 0, 0.8), lwd = 1.2)
    mtext(paste('Ecoregion', eco), 3, line = 0)
}

# legend
plot(0, pch = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend('center', legend = p_weight, lty = 1, lwd = 1.5, col = sim_col, bty = 'n')
dev.off()
