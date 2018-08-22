library(diagram)
par(mar = c(1,1,1,1))
openplotmat()
 elpos <- coordinates (c(1, 4, 3, 3))
 fromto <- matrix(ncol = 2, byrow = TRUE,
                  c(1, 2, 1, 3, 1, 4, 1, 5, 3, 7, 4, 7, 7, 10))
 nr <- nrow(fromto)
 arrpos <- matrix(ncol = 2, nrow = nr)
 for (i in 1:nr)
 arrpos[i, ] <- straightarrow (to = elpos[fromto[i, 2], ],
               from = elpos[fromto[i, 1], ],
               lwd = 2, arr.pos = 0.6, arr.length = 0.5)
 #textrect(elpos[2,], 0.15, 0.05, lab = "Expression Data", box.col = "white",
               #shadow.col = "darkgreen", shadow.size = 0.005, cex = 1.5)
 textrect (elpos[1,], 0.12, 0.04,lab = "Expression Data", box.col = "white",
               shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 #textrect (elpos[3,], 0.15, 0.05,lab = c("P. falciparum","Genome-Scale Metabolic", "Reconstruction"), box.col = "White",
 #             shadow.col = "darkblue", shadow.size = 0.005, cex = 1.2)
 #textrect(elpos[3,], 0.1, 0.1, lab = c("other","term"), box.col = "orange",
                #shadow.col = "darkblue", shadow.size = 0.005, cex = 1.5)
 textrect(elpos[2,], 0.08, 0.08, lab = c("Differentially","Expressed", "Genes"), box.col = "white",
               shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 textellipse(elpos[3,], 0.08, 0.08, lab = c("4hr CQ", "Condition","Model"),box.col = "white",
               shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 textellipse(elpos[4,], 0.08, 0.08, lab = c("24hr CQ", "Condition","Model"),box.col = "white",
               shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 textellipse(elpos[5,], 0.08, 0.08, lab = c("No CQ", "Condition","Model"),box.col = "white",
             shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 textrect (elpos[7,], 0.15, 0.03,lab = "Essential Reactions", box.col = "white",
           shadow.col = "darkblue", shadow.size = 0.005, cex = 1.2)
 textellipse(elpos[10,], 0.08, 0.08, lab = c("Candidate","Drug", "Targets"),box.col = "white",
             shadow.col = "darkblue", shadow.size = 0.004, cex = 1.2)
 # dd <- c(0.0, 0.025)
 text(arrpos[1, 2] - 0.47, arrpos[2, 2] + 0.037, "limma", cex = 1.2)
 text(arrpos[3, 1] + 0.25 , arrpos[2, 2] + 0.028, "MADE", cex = 1.2)
 text(arrpos[3, 1] + 0.25 , arrpos[2, 2] + 0.105, "P. falciparum", cex = 1.2)
 text(arrpos[3, 1] + 0.25 , arrpos[2, 2] + 0.080, "Genome-Scale Metabolic", cex = 1.2)
 text(arrpos[3, 1] + 0.25 , arrpos[2, 2] + 0.054, "Reconstruction,", cex = 1.2)
 text(arrpos[6, 1] - 0.205, arrpos[6, 2] + 0.02, "FBA,", cex = 1.2)
 text(arrpos[6, 1] - 0.2, arrpos[6, 2] - 0.006, "Reaction", cex = 1.2)
 text(arrpos[6, 1] - 0.2, arrpos[6, 2] - 0.032, "Deletions", cex = 1.2)
 
 #text(arrpos[5, 1] - 0.05, arrpos[5, 2] + 0.05, "no")
