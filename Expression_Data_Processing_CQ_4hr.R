
## set working directory to were input files can be found
setwd("C:/Users/Ana/Documents/Chloroquine Project/Input Data")

## Load 4 hour Chloroquine-treated Expression Data
Drug_list_1 <- as.matrix(read.csv("CQ 4h, replicate 1.csv"))
Drug_list_2 <- as.matrix(read.csv("CQ 4h, replicate 2.csv"))
Drug_list_3 <- as.matrix(read.csv("CQ 4h, replicate 3.csv"))

###########################################################
# Create function to replace missing expression values with 
# lowest recorded values across replicates

## INPUTS: Three replicates from one sample

## OUTPUTS: Replicate with replaced missing values; index positions
## where no data was measured for all replicates

replace_missing_values <- function(Rep1, Rep2, Rep3){
  
  missing_data_probes <- matrix() # Create variable to hold gene IDs w/ missing data
  k <- 1 # Index for 'missing_data_probes' object
  
  for (i in 1:length(Rep1)){ # Loop through gene IDs in replicate 1
    
    if (is.na(Rep1[i])){ # If expression data is missing
      
      if (is.na(Rep2[i]) && is.na(Rep3[i])){ # If exprs data is missing for all Reps
        
        missing_data_probes[k] <- i # Save position of gene ID
        k <- k + 1 # Increase index for 'missing_data_probes' by one
        
      }else{ # Otherwise, other replicates have recorded values
        
        if (is.na(Rep2[i]) || is.na(Rep3[i])){ # If one replicate is missing
          
          if (is.na(Rep2[i])){ # Rep2 value is missing
            
            Rep1[i] <- Rep3[i] # Save recorded value
            
          }else{ # Rep3 value is missing
            
            Rep1[i] <- Rep2[i] # Save recorded value
            
          }
          
        }else{ # Otherwise, values for Rep2 and Rep3 are recorded
          
          if (Rep2[i] >= Rep3[i]){ # Check if Rep 2 has greater exprs value than rep 3
            
            Rep1[i] <- Rep3[i] # Save lower exprs value for Rep 1
          }
          
          if (Rep3[i] >= Rep2[i]){ # Otherwise, assume value for Rep 2 is greater than Rep 3
            
            Rep1[i] <- Rep2[i] # Save lower exprs value for Rep 1
          }
        }
      }
    }
  }
  return(list(Rep1, missing_data_probes))
}

###########################################################

## Call replace missing value function for each replicate

Drug_Replicate_1 <- replace_missing_values(Drug_list_1, Drug_list_2, Drug_list_3)
Drug_Replicate_2 <- replace_missing_values(Drug_list_2, Drug_list_1, Drug_list_3)  
Drug_Replicate_3 <- replace_missing_values(Drug_list_3, Drug_list_2, Drug_list_1)

###########################################################

## Load Control Expression Data

No_Drug_list_1 <- as.matrix(read.csv("Trophozoite stage, replicate 1.csv"))
No_Drug_list_2 <- as.matrix(read.csv("Trophozoite stage, replicate 2.csv"))
No_Drug_list_3 <- as.matrix(read.csv("Trophozoite stage, replicate 3.csv"))

###########################################################

## Call replace missing value function for each replicate

No_Drug_Replicate_1 <- replace_missing_values(No_Drug_list_1, No_Drug_list_2, No_Drug_list_3)
No_Drug_Replicate_2 <- replace_missing_values(No_Drug_list_2, No_Drug_list_1, No_Drug_list_3)  
No_Drug_Replicate_3 <- replace_missing_values(No_Drug_list_3, No_Drug_list_2, No_Drug_list_1)

###########################################################
## Remove missing data probes from Drug/No Drug Data

cut <- c(union(Drug_Replicate_1[[2]], No_Drug_Replicate_1[[2]]), 
         setdiff(Drug_Replicate_1[[2]], No_Drug_Replicate_1[[2]]),
         setdiff(No_Drug_Replicate_1[[2]], Drug_Replicate_1[[2]]))

No_Drug_Replicate_1 <- No_Drug_Replicate_1[[1]][-c(cut),]
No_Drug_Replicate_2 <- No_Drug_Replicate_2[[1]][-c(cut),]
No_Drug_Replicate_3 <- No_Drug_Replicate_3[[1]][-c(cut),]

Drug_Replicate_1 <- Drug_Replicate_1[[1]][-c(cut),]
Drug_Replicate_2 <- Drug_Replicate_2[[1]][-c(cut),]
Drug_Replicate_3 <- Drug_Replicate_3[[1]][-c(cut),]

###########################################################
# Load Gene IDs for each probe
Probe_names <- as.matrix(read.csv("Probe Information.csv"))
Probe_names <- Probe_names[-c(cut),] # Remove missing data probes from Drug/No Drug Data

###
# Replace multiple repeating probes with the probe with the largest variance
hits <- matrix() # Object to hold
count <- 1 # Index for 'hits' object

for (i in 1:dim(Probe_names)[1]){ # Loop through each Gene ID
  
  if ((i-1) == dim(Probe_names)[1]){ # If index exceeded matrix length
    break # Exit loop
  }
  
  hits <- matrix()
  count <- 1
  
  # This loop finds the probes that repeat in the list of probe names
  for (j in 1:dim(Probe_names)[1]){ # Loop through each Gene ID
    
    if (Probe_names[i,2] == Probe_names[j,2] && Probe_names[i,2] != ""){ # If a repeat is found
      
     hits[count] <- j # Save index of repeated probe
     count <- count + 1 # Increase index by one
    }
  }
  
  if (length(hits) > 1){ # If a repeat was found
  
  # Create matrices to hold the variance of each repeated probe
  Drug_variances <- matrix(nrow = length(hits), ncol = 1)
  No_Drug_variances <- matrix(nrow = length(hits), ncol = 1)
  
  # Find the variance of each probe that is repeated 
  for (k in 1:length(hits)){ # Loop through repeated probes
    
    # Save expression values from the replicates into one object
    x <- c(Drug_Replicate_1[k], Drug_Replicate_2[k], Drug_Replicate_3[k])
    Drug_variances[k,1] <- var(x) # Find variance
    
    # Save expression values from the replicates into one object
    x <- c(No_Drug_Replicate_1[k], No_Drug_Replicate_2[k], No_Drug_Replicate_3[k])
    No_Drug_variances[k,1] <- var(x) # Find variance 
  }
  
  # Search the vector with the variances for each probe and find the one with the largest variance
  largest_var <- Drug_variances[1,1] # Set first probe as having the largest variance be default
  index <- 1 # Index for 'largest_var' object
  
  for (k in 1:length(Drug_variances)){ # Loop through each repeated probe
    
    if (largest_var < Drug_variances[k,1]){ # If a larger variance if found
      
      largest_var <- Drug_variances[k,1] # Save it
      index <- k # Save it's index
    }
  }
  
  # Save the values of the largest variance probe
  Drug_Replicate_1[hits[1]] <- Drug_Replicate_1[hits[k]]
  Drug_Replicate_2[hits[1]] <- Drug_Replicate_2[hits[k]]
  Drug_Replicate_3[hits[1]] <- Drug_Replicate_3[hits[k]]
  
  # Delete the rest of the repeating probes
  Drug_Replicate_1 <- Drug_Replicate_1[-c(hits[2:length(hits)])]
  Drug_Replicate_2 <- Drug_Replicate_2[-c(hits[2:length(hits)])]
  Drug_Replicate_3 <- Drug_Replicate_3[-c(hits[2:length(hits)])]
  
  # Search the vector with the variances for each probe and find the one with the largest variance
  largest_var <- No_Drug_variances[1,1] # Set first probe as having the largest variance be default
  index <- 1 # Index for 'largest_var' object
  
  for (k in 1:length(No_Drug_variances)){ # Loop through each repeated probe
    
    if (largest_var < No_Drug_variances[k,1]){ # If a larger variance if found
      
      largest_var <- No_Drug_variances[k,1] # Save it 
      index <- k # Save it's index
    }
  }
  
  # Save the values of the largest variance probe
  No_Drug_Replicate_1[hits[1]] <- No_Drug_Replicate_1[hits[k]]
  No_Drug_Replicate_2[hits[1]] <- No_Drug_Replicate_2[hits[k]]
  No_Drug_Replicate_3[hits[1]] <- No_Drug_Replicate_3[hits[k]]
  
  # Delete the rest of the repeating probes
  No_Drug_Replicate_1 <- No_Drug_Replicate_1[-c(hits[2:length(hits)])]
  No_Drug_Replicate_2 <- No_Drug_Replicate_2[-c(hits[2:length(hits)])]
  No_Drug_Replicate_3 <- No_Drug_Replicate_3[-c(hits[2:length(hits)])]
  Probe_names <- Probe_names[-c(hits[2:length(hits)]),]
  }
 }

#############################################################
## Identify probes with missing gene IDs to cut

cut <- c() # Create object to hold missing gene IDs

for (i in 1:dim(Probe_names)[1]){ # Loop though each probe
  
  if (nchar(Probe_names[i,2]) == 0){ # If gene ID is missing
    cut <- c(cut, i) # Add probe to list
  }
}

#############################################################

# Unlog expression data and save each sample in a matrix
Exprs_matrix <- cbind(2^Drug_Replicate_1, 2^Drug_Replicate_2, 2^Drug_Replicate_3, 
                      2^No_Drug_Replicate_1, 2^No_Drug_Replicate_2, 2^No_Drug_Replicate_3)

rownames(Exprs_matrix) <- Probe_names[,2] # Add Gene IDs and labels to expression matrix
colnames(Exprs_matrix) <- c("Drug_rep1","Drug_rep2", "Drug_rep3", "No_Drug_rep1", "No_Drug_rep2", "No_Drug_rep3")

Exprs_matrix <- Exprs_matrix[-c(cut),]

library("Biobase") # Import package to construct ExpressionSet Object
Exprs_matrix <- new("ExpressionSet", exprs = Exprs_matrix) # Save matrix as ExpressionSet Object

# An appropriate design matrix is created were each row
# of the design matrix corresponds to an array in your 
# experiment and each column corresponds to a coefficient
# that is used to describe the RNA sources in your experiment

library("limma")

# Create design matrix
design <- model.matrix(~0+factor(c("Drug", "Drug", "Drug", "Control", "Control", "Control")))
colnames(design) <- c("No_CQ_treatment", "CQ_treated")

# Use lmFit function in Limma to fit a linear model for each gene given a series of arrays
fit <- lmFit(Exprs_matrix, design)

# Create a contrast matrix to make all pair-wise comparisions between the groups
contrast.matrix <- makeContrasts(CQ_treated-No_CQ_treatment, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# A list of top genes differential expressed in Drug versus No Drug
Diff_Exprs_Genes <- topTable(fit2, coef=1, number = length(fit2$coefficients), sort.by = "logFC")

# FDR correction for p-values
pAdj <- p.adjust(Diff_Exprs_Genes$P.Value, method="fdr", n=length(Diff_Exprs_Genes$P.Value))
Diff_Exprs_Genes$adj.P.Val <- pAdj

#############################################################
# Generate spreadsheet of data to integrate into model w/ MADE
#############################################################

# replace '.' with '_' in gene IDs so they are compatible with the model
a <- agrep("MAL", rownames(Diff_Exprs_Genes), max.distance = 0) # Find all gene IDs with '.'

for(i in 1:length(a)){ # Loop through all gene IDs with '.'
  # Replace the '.' with '_'
  rownames(Diff_Exprs_Genes)[a[i]] <- sub(".", "_", rownames(Diff_Exprs_Genes[a[i],]), fixed=TRUE)
}

# Save only gene IDs, Fold Changes, and P-values
MADE_Data <- Diff_Exprs_Genes[,c(1,2,6)]
colnames(MADE_Data) <- c("ORF_old", "logFC", "P_Value")

# Save file for MADE integration
write.csv(MADE_Data, file = "Expression_Data_For_MADE_CQ_4hr.csv")

########################################################
# Determine differentially expressed genes as a FC greater than 2 or less than 0.5

isexpr <- matrix(Diff_Exprs_Genes$adj.P.Val) < 0.05 # Define expressed genes as having a adj.P.val less than 0.05

## Significant based on P-values
nonsignificant_genes <- Diff_Exprs_Genes[!(isexpr),] # Non-significant genes
Diff_Exprs_Genes <- Diff_Exprs_Genes[isexpr,] # Significant genes

## Significant based on Fold Change
FC <- matrix(2^((Diff_Exprs_Genes$logFC))) # Unlog fold changes
isexpr <- FC > 2 | FC < 0.5 # Significance is defined as FC > 2 or FC < 0.5
nonsignificant_genes_2 <- Diff_Exprs_Genes[!(isexpr),]

Diff_Exprs_Genes <- Diff_Exprs_Genes[isexpr,] # Final list of differential expressed genes
nonsignificant_genes <- rbind(nonsignificant_genes, nonsignificant_genes_2) # Final list of non-sig. genes

# Save Non-differentially Expressed Genes
write.csv(nonsignificant_genes, file = "Nonsignificant_Genes_For_CQ_4hr.csv")

# Save Differentially Expressed Genes
write.csv(Diff_Exprs_Genes, file = "Diff_Exp_Genes_For_CQ_4hr.csv")
