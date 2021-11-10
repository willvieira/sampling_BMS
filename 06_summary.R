########################################################################
### Summary of all information by ecoregion
### Author: Will Vieira
### March, 1 2021
########################################################################


#####################################################################
# Steps
# - Load all ecoregion's hexa and merge them together
# - For each ecoregion (set of hexagons)
#   - Print the frequency of habitat classes
#   - Print habitat probability distribution
#   - Print cost distribution
#   - Print cost probability distribution
#   - Print proportion of NA within each hexagon 
#   - Print the frequency of legacy sites
#####################################################################



library(sf)
library(spatstat)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggrepel)


ecoregions <- readRDS('data/ecoregions.RDS')
hexas <- readRDS('data/hexa_complete.RDS')



# For each ecoregion, print all info
#####################################################################
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



pdf('images/summary_by_ecoregion.pdf', width = 8, height = 15)
for(id in ecoregions)
{
    # Select hexagons from ecoregion
    hexa_ecoregion <- subset(hexas, ecoregion == id)

    # start plot
    par(mfrow = c(5, 2), mar = c(2.5, 2.5, 2, 3.5), oma = c(0, 0, 1, 0), mgp = c(1.4, 0.2, 0), tck = -.008)

    # Proportion of habitat classes
    land_dt <- st_drop_geometry(hexa_ecoregion[, paste0('land_ca_', c(1:2, 5:6, 8, 10:14, 16, 21:30))])
    landTotal <- apply(land_dt, 2, sum, na.rm = TRUE)
    landProp <- as.matrix(landTotal/sum(landTotal))
    
    # Proportion of habitat classes for hexagons with inclusion probability > 75 / 95% quantile
    prob <- quantile(hexa_ecoregion$hab_prob, na.rm = TRUE, probs = c(.75, 0.95))
    land_dt <- st_drop_geometry(subset(hexa_ecoregion, hab_prob > prob[1])[, paste0('land_ca_', c(1:2, 5:6, 8, 10:14, 16, 21:30))])
    landTotal <- apply(land_dt, 2, sum, na.rm = TRUE)
    landProp75 <- as.matrix(landTotal/sum(landTotal))
    land_dt <- st_drop_geometry(subset(hexa_ecoregion, hab_prob > prob[2])[, paste0('land_ca_', c(1:2, 5:6, 8, 10:14, 16, 21:30))])
    landTotal <- apply(land_dt, 2, sum, na.rm = TRUE)
    landProp95 <- as.matrix(landTotal/sum(landTotal))

    # barplot    
    landProp <- cbind(landProp, landProp75, landProp95)
    colnames(landProp) <- c('Tous les\nhexagones', 'Hexag avec\nprob > 75%', 'Hexag avec\nprob > 95%')
    plot(0, pch = '', xaxt = 'n', yaxt = 'n', cex = .8, xlim = c(0, 6), ylim = c(0, 1), bty = 'n', xlab = '', ylab = "Proportion des classes d'habitat")
    barplot(landProp, add = TRUE, col = habitat_colors, border = FALSE, cex.names = 0.8)
    legend(4, 1, legend = paste0('land_', c(1:2, 5:6, 8, 10:14, 16, 21:30)), lty = 2, lwd = 5, col = habitat_colors, bty = 'n', cex = 0.8)

    # print ecoregion id
    mtext(paste('Écorégion', id, '-', gsub('_', ' ', names(ecoregions[ecoregions == id]))), 3, line = -0.5, outer = TRUE, cex = 0.75)
    
    # histogram for habitat probability
    xLim <- range(hexa_ecoregion$hab_prob, na.rm = TRUE)
    h <- hist(hexa_ecoregion$hab_prob, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Probabilité d'inclusion habitat de l'hexagone", breaks = seq(xLim[1], xLim[2], length.out = 30), col = 'grey')
    

    # Histogram for costSum
    xLim <- range(hexa_ecoregion$costSum, na.rm = TRUE)
    hist(hexa_ecoregion$costSum, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Coû totale de l'hexagone ($)", breaks = seq(xLim[1], xLim[2], length.out = 30), col = 'grey')

    # Histogram for cost probability
    xLim <- range(hexa_ecoregion$cost_prob, na.rm = TRUE)
    hist(hexa_ecoregion$cost_prob, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Probabilité d'inclusion coût de l'hexagone", breaks = seq(xLim[1], xLim[2], length.out = 30), col = 'grey')

    # Proportion of NA
    xLim <- range(hexa_ecoregion$propNA * 100, na.rm = TRUE)
    hist(hexa_ecoregion$propNA * 100, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Proportion de pixels NA dans l'hexagone (%)", breaks = seq(xLim[1], xLim[2], length.out = 30), col = 'grey')

    # Legacy sites
    if(any(hexa_ecoregion$legacySite > 0))
    {
      # remove hexagons with zero legacy sites
      legacySites <- hexa_ecoregion$legacySite[hexa_ecoregion$legacySite != 0]
      # Define xlim and histogram breaks
      rangeLegacySites <- range(legacySites)
      if(rangeLegacySites[2] > 1) {
        breaks <- seq(1, rangeLegacySites[2], by = 1)
      }else{
        breaks <- seq(1, 2, 1)
      }

      # Color depending on the number of legacy Sites
      colLegacySite <- setNames(viridis_pal()(rangeLegacySites[2]), seq(rangeLegacySites[1], rangeLegacySites[2]))
      
      hist(legacySites, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Nombre de station d'ecoute dans l'hexagone", col = colLegacySite[-1], breaks = breaks)
    }else{
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
    }

    # Spatial aggregation of legacy sites
    if(any(hexa_ecoregion$legacySite > 0))
    {
      # Get Window (ecoregion's boundary)
      owinWindow_eco <- as.owin(sf::st_union(hexa_ecoregion))
      
      # All hexagons with at least 1 legacy site
      coords <- hexa_ecoregion %>%
                      filter(legacySite > 0) %>%
                      sf::st_centroid() %>%
                      sf::st_coordinates() %>%
                      as.data.frame()

      rangeLegacySites <- range(subset(hexa_ecoregion, legacySite > 0)$legacySite)
      
      # Transform in a point pattern obj
      sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
      
      # plot rippley's K
      plot(Kest(sample_ppp, correction = "iso"), main = '', , ylab = "Ripley's K stations d'ecoute", bty = 'l')    
      
      # plot polygon and legacy sites
      plot(owinWindow_eco, main = '')
      points(coords$X, coords$Y, cex = 0.3, pch = 19, col = colLegacySite[legacySites])
      mtext("Hexagones avec au moins une station d'ecoute", 3, line = -69, outer = TRUE, cex = 0.7)

      # All hexagons with >= 4 hexagons
      coords <- hexa_ecoregion %>%
                      filter(legacySite >= 4) %>%
                      sf::st_centroid() %>%
                      sf::st_coordinates() %>%
                      as.data.frame()

      rangeLegacySites <- range(subset(hexa_ecoregion, legacySite >= 4)$legacySite)
      
      # Transform in a point pattern obj
      sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
      
      legacySites <- hexa_ecoregion$legacySite[hexa_ecoregion$legacySite >= 4]

      # plot rippley's K
      plot(Kest(sample_ppp, correction = "iso"), main = '', , ylab = "Ripley's K stations d'ecoute", bty = 'l')    
      
      # plot polygon and legacy sites
      plot(owinWindow_eco, main = '')
      points(coords$X, coords$Y, cex = 0.3, pch = 19, col = colLegacySite[legacySites])
      mtext("Hexagones avec au moins quatre station d'ecoute", 3, line = -92, outer = TRUE, cex = 0.7)

    }else{
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
    }
}
dev.off()



# Print Ripley's K for all ecoregions
######################################################

# Get largest ecoregion
    sampleSize <- readRDS('data/nbHexa_ecoregion.RDS')

    # Sample 2% of hexagons available hexagons
    sampleSize$N2 <- sampleSize$hexagonProp80 * 0.02
    
    # Sample size in function of size proportion between ecoregions
    sampleSize$N <- unclass(round(sum(sampleSize$N2) * sampleSize$sizeProp, 0))
    
    # Total sample size
    N <- sum(sampleSize$N)

    large_eco <- sampleSize$ecoregion[which.max(sampleSize$N)]
#
 
# Prepare Ripley's K data
ripley <- data.frame()
for(id in ecoregions)
{
    # Select hexagons from ecoregion
    hexa_ecoregion <- subset(hexas, ecoregion == id)

    if(any(hexa_ecoregion$legacySite > 1))
    {
      # Get Window (ecoregion's boundary)
      owinWindow_eco <- as.owin(sf::st_union(hexa_ecoregion))
      
      # All hexagons with at least 1 legacy site
      coords <- hexa_ecoregion %>%
                      filter(legacySite > 0) %>%
                      sf::st_centroid() %>%
                      sf::st_coordinates() %>%
                      as.data.frame()
      
      # Transform in a point pattern obj
      sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
      
      # plot rippley's K
      RipleyK <- as.data.frame(Kest(sample_ppp, correction = "iso"))

      # Append to data.frame
      RipleyK$ecoregion <- id
      RipleyK$legacy <- 0
      ripley <- rbind(ripley, RipleyK)
    }

    if(any(hexa_ecoregion$legacySite >= 4))
    {
      # Get Window (ecoregion's boundary)
      owinWindow_eco <- as.owin(sf::st_union(hexa_ecoregion))
      
      # All hexagons with at least 1 legacy site
      coords <- hexa_ecoregion %>%
                      filter(legacySite >= 4) %>%
                      sf::st_centroid() %>%
                      sf::st_coordinates() %>%
                      as.data.frame()
      
      # Transform in a point pattern obj
      sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
      
      # plot rippley's K
      RipleyK <- as.data.frame(Kest(sample_ppp, correction = "iso"))

      # Append to data.frame
      RipleyK$ecoregion <- as.factor(id)
      RipleyK$legacy <- as.factor(4)
      ripley <- rbind(ripley, RipleyK)

    }

    cat(' Ecoregion', which(id == ecoregions), 'of', length(ecoregions), '\r')
}


# GRTS output from the largest ecoregion as as a reference line

  Stratdsgn <- setNames(list(
                            list(panel = c(PanelOne = subset(sampleSize, ecoregion == large_eco)$N), over = 0, seltype = "Continuous")
                            ),
                        paste0('eco_', large_eco)) 

  hexas <- subset(hexas, propNA <= 0.8)
  hexas$p <- (hexas$hab_prob * hexas$cost_prob) / sum(hexas$hab_prob * hexas$cost_prob)
  hexas <- subset(hexas, p != 0)

  sample_frame <- subset(hexas, ecoregion %in% large_eco)

  # Get x and y from centroid of hexagon
  sample_frame$geometry <- sample_frame %>% sf::st_centroid() %>% sf::st_geometry()
  sample_frame[c('X', 'Y')] <- sf::st_coordinates(sample_frame)
    
  # rename ecoregion name to match design name
  sample_frame$eco_name <- paste0('eco_', sample_frame$ecoregion)

  # Get only attributes table (remove spatial information)
  attframe <- sf::st_drop_geometry(sample_frame)

  # Standard p
  attframe$mdcaty <- sum(subset(sampleSize, ecoregion %in% c(101, 103, 78))$N) * attframe$p/sum(attframe$p)
  
  # Run GRTS
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


  # Calculate Ripley's K for GRTS output
  hexa_ecoregion <- subset(hexas, ecoregion %in% large_eco)
  owinWindow_eco <- as.owin(sf::st_union(hexa_ecoregion))
  
  # Centroid coords from selected hexagons
  coords <- hexa_ecoregion %>%
                  filter(ET_Index %in% out$ET_Index) %>%
                  sf::st_centroid() %>%
                  sf::st_coordinates() %>%
                  as.data.frame()

  # Transform in a point pattern obj
  sample_ppp <- ppp(x = coords$X, y = coords$Y, window = owinWindow_eco)
  
  # plot rippley's K
  RipleyGRTS <- as.data.frame(Kest(sample_ppp, correction = "iso"))

#

p1 <- ripley %>%
        group_by(ecoregion) %>%
        mutate(label = if_else(r == max(r), ecoregion, NA_character_)) %>%
        ungroup() %>%
        ggplot(aes(x = r, y = iso, color = ecoregion)) +
          geom_line() +
          geom_label_repel(aes(label = label),
                  nudge_x = 0,
                  na.rm = TRUE,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3,
                  segment.alpha = 0.9,
                  size = 1.8,
                  alpha = 0.9) +
          geom_line(data = subset(ripley, ecoregion == large_eco), aes(x = r, y = theo), linetype = 'dashed', color = 'black') +
          geom_line(data = RipleyGRTS, aes(x = r, y = iso), linetype = 'dashed', color = 'red') +
          facet_grid(~ legacy, labeller = labeller(legacy = c('0' = "Au moins 1 station d'ecoute", '4' = "Au moins 4 station d'ecoute"))) +
          scale_fill_viridis(discrete = TRUE) + theme_classic() + theme_classic() + theme(legend.position = "none")
p2 <- ripley %>%
        filter(r < 50000) %>%
        group_by(ecoregion) %>%
        mutate(label = if_else(r == max(r), ecoregion, NA_character_)) %>%
        ungroup() %>%
        ggplot(aes(x = r, y = iso, color = ecoregion)) +
          geom_line() +
          geom_line(data = subset(ripley, ecoregion == large_eco & r < 50000), aes(x = r, y = theo), linetype = 'dashed', color = 'black') +
          geom_line(data = subset(RipleyGRTS, r < 50000), aes(x = r, y = iso), linetype = 'dashed', color = 'red') +
          ylim(c(0, 50000000000)) +
          geom_label_repel(aes(label = label),
                  nudge_x = 0,
                  na.rm = TRUE,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3,
                  segment.alpha = 0.9,
                  size = 1.8,
                  alpha = 0.9) +
          facet_grid(~ legacy, labeller = labeller(legacy = c('0' = "Au moins 1 station d'ecoute", '4' = "Au moins 4 station d'ecoute"))) +
          scale_fill_viridis(discrete = TRUE) + theme_classic() + theme_classic() + theme(legend.position = "none")

ggsave(
  gridExtra::grid.arrange(p1, p2, nrow = 2, ncol = 1),
  filename = 'images/summary_RipleyK.pdf',
  width = 17,
  height = 15
)



# Polygons for each ecoregion with at least 2 legacy sites
for(eco in unique(ripley$ecoregion))
{
  # Select hexagons from ecoregion
  hexa_ecoregion <- subset(hexas, ecoregion == eco)

  # All hexagons with at least 1 legacy site
  coords <- hexa_ecoregion %>%
                  filter(legacySite > 0) %>%
                  sf::st_centroid() %>%
                  sf::st_coordinates() %>%
                  as.data.frame()  
  
  gg <- ggplot(sf::st_union(hexa_ecoregion)) +
               geom_sf() + 
               geom_point(data = coords, aes(X, Y), size = 0.5) +
               theme_bw() + xlab('') + ylab('') + ggtitle(paste('Ecoregion', eco))

  assign(paste0('plotEco_', eco), gg)

  rm(gg)
}

gg_main <- ripley %>%
              filter(legacy == 4) %>%
              group_by(ecoregion) %>%
              mutate(label = if_else(r == max(r), ecoregion, NA_character_)) %>%
              ungroup() %>%
              ggplot(aes(x = r, y = iso, color = ecoregion)) +
                  geom_line() +
                  geom_line(data = subset(ripley, ecoregion == large_eco), aes(x = r, y = theo), linetype = 'dashed', color = 'black') +
                  geom_line(data = RipleyGRTS, aes(x = r, y = iso), linetype = 'dashed', color = 'red') +
                  geom_label_repel(aes(label = label),
                      nudge_x = 0,
                      na.rm = TRUE,
                      max.overlaps = Inf,
                      min.segment.length = 0,
                      segment.size = 0.3,
                      segment.alpha = 0.9,
                      size = 2.2,
                      alpha = 0.9) +
                  scale_fill_viridis(discrete = TRUE) + theme_classic() + theme(legend.position = "none")


lay <- rbind(c(1, 2, 3, 4),
             c(5, 6, 6, 7),
             c(8, 6, 6, 9),
             c(10, 11, 12, 13))

ggsave(
  gridExtra::grid.arrange(plotEco_72, plotEco_74, plotEco_75, plotEco_76,
                          plotEco_96, gg_main, plotEco_99,
                          plotEco_100, plotEco_101,
                          plotEco_102, plotEco_103, plotEco_117, plotEco_217,
                          layout_matrix = lay),
  filename = 'images/summary_RipleyK_map.pdf',
  width = 17,
  height = 13)



# Table to sumarise legacy sites per ecoregion
# - Sample size (in hexagons)
# - Sample size (in stations)
# - Total of legacy sites
# - Total of hexagons with X legacy sites
# - Sample size for hypothesis 0 to 4


# Get sample size for each ecoregion

    sampleSize <- readRDS('data/nbHexa_ecoregion.RDS')

    # Sample 2% of hexagons available hexagons
    N_tmp <- sampleSize$hexagonProp80 * 0.02

    # Sample size in function of size proportion between ecoregions
    eco_ignore <- which(sampleSize$ecoregion %in% gsub('N', '', grep('N', sampleSize$ecoregion, value = TRUE)))
    totalN <- sum(N_tmp[-eco_ignore])
    sampleSize$N <- unclass(round(totalN * sampleSize$sizeProp, 0))

#
 

# Max of legacy site a hexagon have for the whole study area
maxLegacySite <- max(hexas$legacySite)

summary_legacy <- data.frame()
for(eco in unique(hexas$ecoregion))
{
    # get hexagons for specific ecoregion
    hex_eco <- subset(hexas, ecoregion == eco)

    # sample size in hexagons
    N_eco <- subset(sampleSize, ecoregion == eco)$N

    # Frequence of legacy site stations within hexagons
    legacySite_freq <- setNames(rep(0, 32), paste0('nb_', 0:31))
    legacySite_freq[paste0('nb_', names(table(hex_eco$legacySite)))] <- table(hex_eco$legacySite)

    # Sample size after legacy sites
    N_h1 <- N_h2 <- N_eco
    N_h3 <- N_eco - nrow(subset(hex_eco, legacySite >= 4))
    N_h4 <- round((N_eco * 10 - sum(hex_eco$legacySite))/10, 0)

    # assign to data.frame
    summary_legacy <- rbind(summary_legacy,
                            data.frame(
                                ecoregion = eco,
                                total_hex = nrow(hex_eco),
                                total_legacySite = sum(hex_eco$legacySite),
                                data.frame(as.list(legacySite_freq)),
                                sampleSize_hex = N_eco,
                                sampleSize_station = N_eco * 10,
                                sampleSize_h1 = N_h1,
                                sampleSize_h2 = N_h2,
                                sampleSize_h3 = ifelse(N_h3 < 0, 0, N_h3),
                                sampleSize_h4 = ifelse(N_h4 < 0, 0, N_h4)
                            ))

}

write.csv(summary_legacy, 'summary_legacySite.csv', row.names = FALSE)
