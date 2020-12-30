# Goal: run t-SNE on Wellderly subjects and plot ancestry admixture results
# The Input file used here is obtained using the script "ancestry_Wdy.sh"


# Load required libraries
library("Rtsne")
library("ggplot2")

# Import input file and transform into a matrix
mydata <- read.table("/home/mpinsach/Escriptori/Ancestry/Wdy/Pruned/Wdy_1000Genomes.eigenvec", sep=" ")
matrix <- as.matrix (as.data.frame(mydata[,3:8])) # We use the first 6 principal components

# Run Rtsne
set.seed(9)
tsne_model_1 <- Rtsne(matrix, dims=2, perplexity=30, theta=0.0, check_duplicates=FALSE, pca=FALSE, max_iter=1000, eta=200, exaggeration_factor=12)

# Plot Rtsne results
d_tsne_1 <- as.data.frame(tsne_model_1$Y)
Population <- mydata[,1]
plot <- ggplot(d_tsne_1, aes(x=V1, y=V2, colour=Population)) + geom_point(size=0.5, alpha=0.7) + xlab("tSNE-1") + ylab("tSNE-2") + theme_classic(base_size=11) + theme(panel.background=element_rect(fill="white", color="black")) + scale_colour_manual(values=c("#7d000f", "#659e24", "#cc3c14", "#0b5bae", "#e6a219", "#dd1c77"))
tiff("/home/mpinsach/Escriptori/Ancestry/Wdy/Wdy_ancestry_plot.tiff", width = 15, height = 10, units = 'cm', res = 300)
plot
dev.off()
