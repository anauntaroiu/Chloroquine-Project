
# set working directory to were the file can be found
setwd("C:/Users/Ana/Documents/Chloroquine Project/Processed Expression Data")

# Load and Save Differentially Expressed Genes
Diff_Exp_4hr <- as.matrix(read.csv("Diff_Exp_Genes_For_CQ_4hr.csv"))[,1]
Diff_Exp_24hr <- as.matrix(read.csv("Diff_Exp_Genes_For_CQ_24hr.csv"))[,1]

index <- 0 # Object to hold number of overlapping genes

for (i in 1:length(Diff_Exp_4hr)){ # Loop through all 4hr Diff Exp Genes
    
    # Check for matched between 4hr and 24hr Diff Exp Genes 
    a <- agrep(Diff_Exp_4hr[i], Diff_Exp_24hr, max.distance = 0)
    
    if (length(a) > 0){ # If match found
      index <- index + 1 # Note it
    }
}


