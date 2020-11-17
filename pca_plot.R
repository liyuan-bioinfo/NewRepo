rm(list = ls())
graphics.off()

require(plyr)
require(dplyr)
require(devtools)
require(easyGgplot2)
require(factoextra) # PCA
require(FactoMineR) # PCA
require(gdata) # cbindX columns with diff no. of rows
require(ggplot2)
require(ggpubr)
require(readr)
require(reshape2)
require(RColorBrewer)
require(R.utils)
require(pls)

#----------------------#
#   1-initial params
#----------------------#
#getwd()
# Set working direct to the folder containing spectra in .txt
setwd("~/pca分析/bacteria and archaea/")
files <- list.files(path = "~/pca分析/bacteria and archaea/")

# Create a final_data_wide dataframe contaning all spectra
final_data_wide <- read.table(files[1], 
                              header = FALSE, sep = "\t")[, 1]
for (filename in files)
{
  data <- read.table(filename, header = FALSE, sep="\t")[2]
  colnames(data) <- filename
  final_data_wide <- cbind(final_data_wide, data)
}
colnames(final_data_wide)[1] <- "Wavenumber"

# Subset data to fingerprint region
#final_data_wide <- final_data_wide[c(79:444), ]

# Convert final_data_wide to final_data_long
final_data_long <- melt(final_data_wide, id.vars = "Wavenumber",
                        variable.name = "Filename",
                        value.name = "Intensity")

final_data_long["Condition"] <- final_data_long$Filename

# Fix file names
final_data_long$Condition <- gsub("_.*", "", final_data_long$Condition, perl = TRUE)

# Convert final_data_long to fianl_data_rtoc for PCA analysis
final_data_rtoc <- reshape(final_data_long, timevar = "Wavenumber", 
                           idvar = c("Filename", "Condition"), direction = "wide")
#-------------------------------------
# 2，Plot spectra
#-------------------------------------

# Stats calculation of mean, standard devation and error from single-cell data
final_data_stats1 <- melt(final_data_long, id.vars=c("Condition","Wavenumber"),
                          measure.vars="Intensity")

stat.1 <- ddply(final_data_stats1, c("Condition", "Wavenumber"), summarise,
                average = mean(value), stdev = sd(value),
                error = sd(value)/sqrt(length(value)))

#  Order condition and Modify average to modified average to scatter spectra
stat.1$Condition <- factor(stat.1$Condition, 
                           levels = c("H7","SCM1"))

stat.1$modified_average = stat.1$average - 0.005 * as.numeric(as.factor(stat.1$Condition)) 


# Plot averages and ribbons
SpectraPlot <- ggplot(stat.1, aes(x = Wavenumber, y = modified_average, group = Condition)) + 
  geom_line(aes(color = Condition)) + 
  geom_ribbon(aes(ymin = modified_average - stdev, 
                  ymax = modified_average + stdev, 
                  fill = Condition),
              alpha = 0.5) +
  labs(title="",x=expression(Wavenumber/cm^{"-1"}), y = "Intensity") +
  scale_x_continuous(breaks=seq(400,3200,400)) +
  coord_cartesian(xlim = c(400, 3200)) +
  theme_linedraw() + 
  # coord_cartesian(ylim = c(0, 5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank()
  ) + ggtitle("")
# annotate("segment", x = 2185, xend = 2195, y = 0.25, yend = 1.1, linetype = "dashed",
#         colour = "grey", alpha = 0.5) +
# annotate("text", x = 2187.5, y = 1.15, label = "C-D", colour = "black", size = 3.7)

SpectraPlot




#---------------------  3，for PCA analysis -----------------
#
#---------------------------------------------------------
final_data_rtoc <- reshape(final_data_long, timevar = "Wavenumber", idvar = c("Filename", "Condition"), direction = "wide")

final_data_rtoc$Condition <- factor(final_data_rtoc$Condition, 
                                    levels = c("H7","SCM1"))


# The variable Condition and Filename (index = 1, 2) is removed before PCA analysis
final_data_rtoc.pca <- PCA(final_data_rtoc[, -(1:2)], scale.unit = FALSE, ncp = 10, 
                           graph = FALSE)
get_eigenvalue(final_data_rtoc.pca)
eigenvalues <- final_data_rtoc.pca$eig
head(eigenvalues[, 1:2])
final_data_rtoc.pca$var$contrib[,1]

# Visualize
# Use habillage to specify groups for coloring note to use as.factor
final_data_rtoc.PCA <- fviz_pca_ind(final_data_rtoc.pca,
                                    label="none", axes = c(1,2),
                                    habillage = as.factor(final_data_rtoc$Condition),
                                    repel = TRUE,
                                    addEllipses=TRUE, ellipse.level=0.75) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_blank()) + 
  ggtitle("")


final_data_rtoc.PCA
#--------------------------------------------#
#achieve the dim1 and dim2 important features
#--------------------------------------------#

var <- get_pca_var(final_data_rtoc.pca)
var
head(var$cos2, 4)
#The larger the value of the contribution, the more the variable contributes to the component.
# Contributions of variables to PC1
#-------------------
#PC1 is 88%,so used the PC1 to get plot of which contribute more
# ------------------
fviz_contrib(final_data_rtoc.pca, choice = "var", axes = 1, top = 200)
#The total contribution to PC1 and PC2
#fviz_contrib(final_data_rtoc.pca, choice = "var", axes = 1:2, top = 100)

intensity <- as.double(gsub("Intensity.","",names(temp_plot)))
plot(var$contrib[,1] ~ intensity,type = "l",col="black",
    xlab= "Intensity",ylab="Contribute to PC1",lwd="2",
    ylim=c(0,0.1),xlim=c(600,1800))
    #text(3410,1.6,"3403")

#pca2 doesn't contirbute too much
#plot(var$contrib[,2] ~ intensity,type = "l",col="black",
#     xlab= "Intensity",ylab="Contribute to PC2",lwd="2",ylim=c(0,2))
