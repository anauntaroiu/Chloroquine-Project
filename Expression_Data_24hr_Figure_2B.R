
# set working directory to were the files can be found
setwd("C:/Users/Ana/Documents/Chloroquine Project/Processed Expression Data")

# Import 4hr/24hr expression data
Nonsignificant_genes_24hr <- read.csv("Nonsignificant_Genes_For_CQ_24hr.csv")
Diff_Exp_genes_24hr <- read.csv("Diff_Exp_Genes_For_CQ_24hr.csv")

# Create Scatterplot
plot(x=Nonsignificant_genes_24hr[,2],y=Nonsignificant_genes_24hr[,6], type ="p", col = "black", xlab="logFC", ylab="Adj. P-value", cex.lab=1.4)
points(x=Diff_Exp_genes_24hr[,2],y=Diff_Exp_genes_24hr[,6], type ="p", col="red")

# Add lengend and label
legend(x= 0.13, y= 1.03, legend=c("Differentially Expressed Genes","Nonsignificant Genes"), col=c("red","black"),cex=1, text.col=c("black","black"), pch=1)
text(x = -0.96, y = 0.95, labels = 'B', cex = 4)


