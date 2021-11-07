workingDir = "E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\optimized_tissue\\";
setwd(workingDir);
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

nSets = 5;
multiExpr = vector(mode = "list", length = nSets);

brain_names = dir()
multiExpr = vector(mode = "list", length = nSets)

for (i in 1:5){
brain_expr = read.table(brain_names[i],sep="\t",header=T)
multiExpr[[i]] = list(data = as.data.frame(t(brain_expr[2:dim(brain_expr)[2]])));
names(multiExpr[[i]]$data) = brain_expr$Name;
rownames(multiExpr[[i]]$data) = names(brain_expr)[2:dim(brain_expr)[2]];

}

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2,networkType = "signed")[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
2
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(16, 12)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
setLabels = "";
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

net = blockwiseConsensusModules(
  multiExpr, power = 12, minModuleSize = 30, deepSplit = 3,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  networkType = "signed",
  TOMType = "signed",
  maxBlockSize = 26000,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5,
  networkCalibration="single quantile",
  calibrationQuantile=0.95,
  consensusQuantile=0.95,
  nThreads=20)

moduleLabels = net$colors;
moduleColors = labels2colors(moduleLabels);
for (i in 1:23){
     index<- moduleColors==all_color[i]
     gene_name <- names(net$colors[which(index)])
     gene_name <- as.data.frame(gene_name)
     write.table(gene_name,file=paste("/home/lxy/02ECS_eqtl/result/KGGSEE_maf05/WGCNA/Brain_tmm_normalized_geneExpression_consensus_network/moduleGene_5_brain_tissue/",all_color[i],"module_gene.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
  }