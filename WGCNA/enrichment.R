Args <- commandArgs(TRUE)
file <- Args[1]#all tissues' gmt file
eqtl_gene_file <- Args[2]
eqtl_trans_file <- Args[3]
pred <- Args[4]
output <- Args[5]
dise <- Args[6]

library(gprofiler2)
#gene_set = dir(folder)
#gmt = read.table(file,sep="\t",header=T)
id = upload_GMT_file(file)

eqtl_gene = read.csv(eqtl_gene_file,sep="\t",header=T)
eqtl_trans = read.csv(eqtl_trans_file,sep="\t",header=T)
eqtl_gene = eqtl_gene[eqtl_gene$CondiECSP<3.125e-7,]
eqtl_trans = eqtl_trans[eqtl_trans$CondiECSP<3.125e-7,]
eqtl_trans1 <- c()
for (i in 1:length(eqtl_trans$Gene)){
 item <- unlist(strsplit(as.character(eqtl_trans$Gene[i]),split=":"))[1]
 eqtl_trans1[i] <- item
}
pred = read.table(pred,sep=",",header=T)
pred = pred[pred$pvalue<8.3e-7,]
gene = gconvert(query = eqtl_gene$Gene, organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)
transcript = gconvert(query = eqtl_trans1, organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)
pred = gconvert(query = pred$gene_name, organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)
write.table(data.frame(gene$target),file=paste(dise,"gene_ENSG.txt",sep=""),row.names=F,quote=F)
write.table(data.frame(pred$target),file=paste(dise,"pred_ENSG.txt",sep=""),row.names=F,quote=F)
write.table(data.frame(transcript$target),file=paste(dise,"transcript_ENSG.txt",sep=""),row.names=F,quote=F)

Query = list(eqtl_gene_file=as.character(gene$target),eqtl_trans_file=transcript$target,predixcan=pred$target)
print(Query)
custom_gostres <- gost(query=Query,multi_query = TRUE,organism =id,significant = T)

result <- custom_gostres$result
print(result)
write.table(result$term_id,file=output,sep="\t",quote=F)

