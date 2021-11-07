library(ggplot2)
tiff("MAGMA_edese_dist_common_gene_GO.tiff",width=6,height=6,units="in",res=300,compression = 'lzw')
df <- read.table("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\MCGA_Dist_vs_MAGMA\\common_105_genes_gprofiler.txt",sep="\t",head=T)
df$V1 <- factor(df$source,levels = unique(df$source),order=T)
df$V2 <- factor(df$term_name,levels = unique(df$term_name),order=T)
col <- c('#ffb6b9','#fae3d9','#bbded6')
#col <- c('#FFCCCC',"#FF69B4")
#col <- c('#FFDDAA')#cc
#col <- c('#FFCCCC')#BP
#p <- ggplot(data = df, aes(x = V2, y = log, fill = V1)) + geom_bar(stat = 'identity', position = 'dodge')+coord_flip()+labs(x = "", y = "-Log10(adjusted p-value)")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text = element_text(size=6),axis.title = element_text(size=8),axis.text.x = element_text(angle = 90),legend.text=element_text(size=5),panel.background = element_blank())+scale_fill_manual(values = col)+guides(fill=guide_legend(title=''))
p <- ggplot(data = df, aes(x = V2, y = log, fill = V1)) + geom_bar(stat = 'identity', position = 'dodge')+labs(x = "", y = "-Log10(adjusted p-value)")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text = element_text(size=9),axis.title = element_text(size=8),axis.text.x = element_text(angle = 78,hjust =1, vjust =1),legend.text=element_text(size=5),panel.background = element_blank())+scale_fill_manual(values = col)+guides(fill=guide_legend(title=''))
p+theme(legend.key.size = unit(0.6,"line"))
dev.off()