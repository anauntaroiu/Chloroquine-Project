
# set working directory to were the files can be found
setwd("C:/Users/Ana/Documents/Chloroquine Project/Processed Expression Data")

# Import 4hr/24hr expression data
Nonsignificant_genes_4hr <- read.csv("Nonsignificant_Genes_For_CQ_4hr.csv")
Diff_Exp_genes_4hr <- read.csv("Diff_Exp_Genes_For_CQ_4hr.csv")

# Create Scatterplot
plot(x=Nonsignificant_genes_4hr[,2],y=Nonsignificant_genes_4hr[,6], type ="p", col = "black", xlab="logFC", ylab="Adj. P-value", cex.lab=1.4)
points(x=Diff_Exp_genes_4hr[,2],y=Diff_Exp_genes_4hr[,6], type ="p", col="red")

# Add lengend and label
legend(x= -1.6, y= 0.8, legend=c("Differentially Expressed Genes","Nonsignificant Genes"), col=c("red","black"),cex=1, text.col=c("black","black"), pch=1)
text(x = -1.5, y = 0.95, labels = 'A', cex = 4)

