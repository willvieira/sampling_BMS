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
load('data/spatialVectors.rda')


# Load all ecoregion's hexa and merge them together
#####################################################################

hexa_ls <- list()
for (id in unique(districts$ECOREGION))
{
  hexa_ls[[as.character(id)]] <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == id)$REGION_NAM))), '_', id, '/hexa_', id, '.shp'), quiet = TRUE)
  hexa_ls[[as.character(id)]]$ecoregion <- id
}

# rbind all shapefiles into one sf oject
hexas <- do.call(rbind, hexa_ls)
rm(hexa_ls)



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
)


pdf('summary_by_ecoregion2.pdf', width = 8, height = 9)
for(id in sort(unique(districts$ECOREGION)))
{
    # Select hexagons from ecoregion
    hexa_ecoregion <- subset(hexas, ecoregion == id)

    # start plot
    par(mfrow = c(3, 2), mar = c(2.5, 2.5, 1, 3.5), oma = c(0, 0, 1, 0), mgp = c(1.4, 0.2, 0), tck = -.008)
    

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
    mtext(paste('Écorégion', id, '-', unique(subset(districts, ECOREGION == id)$REGION_NAM)), 3, line = -0.5, outer = TRUE, cex = 0.75)
    
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
      hist(legacySites, main = '', ylab = 'Frequence (# de hexagones)', xlab = "Nombre de station d'ecoute dans l'hexagone", col = 'grey')
    }else{
      plot(0, pch = '', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', bty = 'n')
      text(1, 0, "Pas de station d'ecoute dans cette écorégion")
    }
}
dev.off()
