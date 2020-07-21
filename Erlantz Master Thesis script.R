###Figure 2###
########Figure2A############

data_bsj <- read.delim("../../ciri_bsj_table.txt", header = T, stringsAsFactors = F, sep="\t")
data_bsj2 <- data_bsj

metadata_sort <- metadata[order(metadata$type, metadata$gleason),]
data_bsj2 <- data_bsj2[,match(metadata_sort$sample_id, colnames(data_bsj2))]
colnames(data_bsj2) == metadata_sort$sample_id

data_bsj2[data_bsj2 <4] <- NA #filtering# / data_bsj2[data_bsj2 == 0] <- NA #no filtering

circRNAs_perCol <- colSums(!is.na(data_bsj2))
circRNAs_perCol <- data.frame(circRNAs_perCol)
circRNAs_perCol$Sample <- rownames(circRNAs_perCol)
colnames(circRNAs_perCol) [1] <- "NcircRNAs"
circRNAs_perCol <- circRNAs_perCol[sort(circRNAs_perCol$Sample),]
circRNAs_perCol$type <- metadata$type
circRNAs_perCol <- circRNAs_perCol[order(circRNAs_perCol$type, circRNAs_perCol$NcircRNAs),]
circRNAs_perCol$Sample <- as.vector(circRNAs_perCol$Sample)
circRNAs_perCol$Sample = factor(circRNAs_perCol$Sample, circRNAs_perCol$Sample)
colnames(circRNAs_perCol) [3] <- "Type" 
number_circRNa <- ggplot(circRNAs_perCol, aes(x=Sample, y=NcircRNAs, color=Type)) +
  ylim(0,7500) + geom_bar(stat = "identity", fill="white") + 
  scale_colour_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue"))+
  labs(title = "# of unique BSJs per sample (filtered ≥3 reads)", x="Samples", y="# of unique BSJs")+
  theme(strip.background = element_blank(), strip.text = element_text(size=7, face="bold"), plot.title = element_text(size=14, face="bold",hjust = 0.5), axis.text=element_text(size=14, face="bold", colour = "black"),
        axis.title=element_text(size=14,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2),legend.text = element_text(size=14), legend.title = element_text(face="bold",size=14))

df.venn <- data.frame(x = c(-0.2, 0.866),
                      y = c(0.5, 0.5),
                      Type = c('NAP', 'PCa'))
ven <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1, fill = Type)) +
  geom_circle(alpha = .2, size = 1, colour = 'grey') +
  scale_fill_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue"))+
  scale_colour_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue"))+
  coord_fixed() +
  theme_void()+
  theme(legend.position = "none")

ven_bar <- number_circRNa + annotation_custom(ggplotGrob(ven), xmin = 40, xmax = 100, 
                                              ymin = 4000, ymax = 7500)
##Number of unique circRNAs per NAP and PCa
circRNAs_perRow <- rowSums(!is.na(data_bsj2))
circRNAs_perRow <- circRNAs_perRow[circRNAs_perRow >0]
str(circRNAs_perRow)
circRNAs_perRow_NAP <- rowSums(!is.na(data_bsj2[,1:40]))
circRNAs_perRow_PCa <- rowSums(!is.na(data_bsj2[,41:91]))

row_circRNAs_NAP <- circRNAs_perRow_NAP[circRNAs_perRow_NAP == 0]
row_circRNAs_PCa <- circRNAs_perRow_PCa[circRNAs_perRow_PCa >0]
row_circRNAs_NAP <- data.frame(row_circRNAs_NAP)
row_circRNAs_PCa <- data.frame(row_circRNAs_PCa)
unique_PCa <-row_circRNAs_PCa[match(rownames(row_circRNAs_NAP), rownames(row_circRNAs_PCa)),]
unique_PCa <- unique_PCa[!is.na(unique_PCa)]
str(unique_PCa)

row_circRNAs_NAP_1 <- circRNAs_perRow_NAP[circRNAs_perRow_NAP > 0]
row_circRNAs_PCa_1 <- circRNAs_perRow_PCa[circRNAs_perRow_PCa == 0]
unique_NAP <-row_circRNAs_NAP_1[match(names(row_circRNAs_PCa_1), names(row_circRNAs_NAP_1))]
unique_NAP <- unique_NAP[!is.na(unique_NAP)]
str(unique_NAP)

row_circRNAs_mix <- circRNAs_perRow_NAP[circRNAs_perRow_NAP == 0]
row_circRNAs_mix_1 <- circRNAs_perRow_PCa[circRNAs_perRow_PCa == 0]
unique_mix <-row_circRNAs_mix[match(names(row_circRNAs_mix_1), names(row_circRNAs_mix))]
unique_mix <- unique_mix[!is.na(unique_mix)]
str(unique_mix)

#######Figure 2B######

rownames(data_bsj) <- data_bsj$circRNA.id
data_bsj$circRNA.id <- NULL
names(data_bsj) <- gsub("X","",colnames(data_bsj), ignore.case = T)
colnames(data_bsj) <- gsub(".","-",colnames(data_bsj), fixed = T)
colnames(data_bsj) <- gsub("-gtf",".gtf",colnames(data_bsj), fixed=T)
colnames(data_bsj) <- gsub(".gtf","",colnames(data_bsj), fixed=T)
rownames(data_bsj) <- paste0(rownames(data_bsj),";",data_bsj$gene_name)
data_bsj [,92:93] <- NULL
head(data_bsj)

metadata$type <- gsub(" ","_",metadata$type, fixed = T)
metadata$type <- gsub("-","_",metadata$type, fixed = T)
metadata$gleason <- gsub("-","_",metadata$gleason, fixed = T)
metadata_Gleason <- metadata %>% filter(type == "organ_confined_PCa")

data_bsj <- data_bsj[,match(metadata$sample_id, colnames(data_bsj))]
colnames(data_bsj) == metadata$sample_id

design <- as.factor(metadata_order$type)
dds_bsj <- DESeqDataSetFromMatrix(data_bsj, DataFrame(design), ~design)
dds_bsj <- DESeq(dds_bsj)
res_bsj <- results(dds_bsj)
res_bsj <- res_bsj[order(res_bsj$padj),]
sig_bsj <- res_bsj[!is.na(res_bsj$padj) & res_bsj$padj < 0.01,]

data_bsj_norm <- assay(vst(dds_bsj))
bsj_norm_sig_type<- data_bsj_norm[match(rownames(sig_bsj), rownames(data_bsj_norm)),]
PCA_type_sig_BSR <- prcomp(t(bsj_norm_sig_type))
autoplot(PCA_type_sig, data= metadata, colour= "type", main="Supervised PCa vs. NAP") + labs(colour="Type")+
  scale_colour_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=10, face="bold", colour = "black"),
                                                                                       axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))

##Figure 2C###
metadata2 <- metadata
colnames(metadata2) [2] <- "Gleason_Score"
metadata2$`Gleason score` <- gsub("Agressive","Aggressive",metadata2$`Gleason score`,fixed=T)
metadata2$`Gleason score` <- gsub("Non-agressive","Non-Aggressive",metadata2$`Gleason score`,fixed=T)
metadata2$`Gleason score` <- gsub("Healthy","Non-cancer",metadata2$`Gleason score`,fixed=T)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(metadata2$`Gleason score`))))
mycolors <- newCols(length(unique(metadata2$`Gleason score`)))
names(mycolors) <- unique(metadata2$`Gleason score`)
mycolors <- list(`Gleason score` = mycolors)

ann_colors <- list(
  `Gleason score` = c(Aggressive="#240CF6", `Non-Aggressive`="#04FF7A", `Non-cancer`="#AE0000"),Type = c(NAP="#FF5454", PCa="#97D3F0"))

heatmap_plot <- pheatmap(bsj_norm_sig_type, main = "Supervised PCa vs. NAP",
                         show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",fontsize=12,fontsize_row = 1, fontsize_col = 0.1,
                         cex=1, clustering_distance_rows="euclidean", cex=1,cellwidth = 2,cellheight = 1,
                         clustering_distance_cols="correlation", clustering_method="complete", border_color=FALSE,
                         annotation_col=metadata2, treeheight_row = 20, treeheight_col = 30, annotation_colors = ann_colors)



#####Supplementary Figure 1######

data_fsj <- read.table("/mnt/data/ccbc_environment/project/NGSProtocol/processed/ciri/ciri_out/ciri_fsj_table.txt", header = T, stringsAsFactors = F, sep="\t")
metadata_metas <- read.table("/mnt/data/ccbc_environment/project/NGSProtocol/processed/ciri/erlantz-data/R_analysis_bsj/info_justmetas.txt",header = T, stringsAsFactors = F, sep="\t")

##BSR data
##BSR NAP vs PCa unsup##

data_bsj_norm <- assay(vst(dds_bsj))
PCA_type_un <- prcomp(t(data_bsj_norm))

unsup_PCApcanap <- autoplot(PCA_type_un, data= metadata_order, colour= "type", main="Unsupervised PCa vs. NAP") + labs(colour="Type") + geom_text(aes(label=ifelse(PC2>0.15,as.character(colnames(data_bsj_norm)),'')),hjust=-0.15,vjust=0.5) +
  scale_colour_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue")) +
  theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12, face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
        axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2),legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12))

##Gleason BSR###
metadata_Gleason <- metadata %>% filter(type == "organ_confined_PCa")
data_bsj_gle <- data_bsj[,match(metadata_Gleason$sample_id, colnames(data_bsj))]
design2 <- as.factor(metadata_Gleason$gleason)
dds_bsj2 <- DESeqDataSetFromMatrix(data_bsj_gle, DataFrame(design2), ~design2)
dds_bsj2 <- DESeq(dds_bsj2)
res2 <- results(dds_bsj2)
res2 <- res2[order(res2$padj),]
sig2_bsj <- res2[!is.na(res2$padj) & res2$padj < 0.01,]

bsj_norm_sig_gleason<- bsj_norm_gleason[match(rownames(sig2_bsj), rownames(data_bsj_norm)),]
bsj_norm_gleason <- data_bsj_norm[,match(metadata_Gleason$sample_id, colnames(data_bsj_norm))]
PCA_gleason_sig <- prcomp(t(bsj_norm_sig_gleason))
sup_gle_BSR <- autoplot(PCA_gleason_sig, data= metadata_Gleason, colour= "gleason", main="Supervised Ag vs. Nag (Gleason Score)") + labs(colour="Gleason Score") +
  scale_colour_manual(name = "Gleason Score", labels = c("Aggressive (4+3/>7)", "Non-Aggressive (3+3/3+4)"), values = c("red", "lightblue")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12,face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                                                                                    axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                    panel.background = element_blank(),legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12), axis.line = element_line(colour = "black", size = 1.2))


###Metastasis BSR##
metadata_metas$Met_Pcadeath <-  gsub("no", "y", metadata_metas$Met_Pcadeath, fixed=T) 
metadata_metas$Met_Pcadeath <-  gsub("yes", "n", metadata_metas$Met_Pcadeath, fixed=T)
data_bsj_met <- data_bsj[,match(metadata_metas$sample_id, colnames(data_bsj))]
design3 <- as.factor(metadata_metas$Met_Pcadeath)
dds_fsj3 <- DESeqDataSetFromMatrix(data_bsj_met, DataFrame(design3), ~design3)
dds_fsj3 <- DESeq(dds_fsj3)
res3 <- results(dds_fsj3)
res3 <- res3[order(res3$padj),]
sig3 <- res3[!is.na(res3$padj) & res3$padj < 0.01,]

bsj_norm_metas <- data_bsj_norm[,match(metadata_metas$sample_id, colnames(data_bsj_norm))]
bsj_norm_sig_metas<- bsj_norm_metas[match(rownames(sig3), rownames(bsj_norm_metas)),]
PCA_sig_metas <- prcomp(t(bsj_norm_sig_metas))
sup_PCAagg_non_met_BSR <- autoplot(PCA_sig_metas, data= metadata_metas, colour= "Met_Pcadeath", main="Supervised Ag vs. Nag (Metastasis/PCa death)") + labs(colour="Metastasis/PCa death") +
  scale_colour_manual(name = "Metastasis/PCa death",labels = c("Non-Aggressive (No)","Aggressive (Yes)"), values = c("lightblue","red")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12, face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                                                                                axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                panel.background = element_blank(), legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12),axis.line = element_line(colour = "black", size = 1.2))

###FSR data##
rownames(data_fsj) <- data_fsj$circRNA.id
data_fsj$circRNA.id <- NULL
names(data_fsj) <- gsub("X","",colnames(data_fsj), ignore.case = T)
colnames(data_fsj) <- gsub(".","-",colnames(data_fsj), fixed = T)
colnames(data_fsj) <- gsub("-gtf",".gtf",colnames(data_fsj), fixed=T)
colnames(data_fsj) <- gsub(".gtf","",colnames(data_fsj), fixed=T)
rownames(data_fsj) <- paste0(rownames(data_fsj),";",data_fsj$gene_name)
data_fsj [,92:93] <- NULL

##NAP vs. PCa FSR Supervised##
design <- as.factor(metadata$type)
dds_fsj <- DESeqDataSetFromMatrix(data_fsj, DataFrame(design), ~design)
dds_fsj <- DESeq(dds_fsj)
res_fsj <- results(dds_fsj)
res_fsj <- res_fsj[order(res_fsj$padj),]
sig_fsj <- res_fsj[!is.na(res_fsj$padj) & res_fsj$padj < 0.01,]


data_fsj_norm <- assay(vst(dds_fsj))
fsj_norm_sig_type<- data_fsj_norm[match(rownames(sig_fsj), rownames(data_fsj_norm)),]
PCA_type_sig_fsj <- prcomp(t(fsj_norm_sig_type))
sup_nappca_FSR <- autoplot(PCA_type_sig_fsj, data= metadata, colour= "type", main="Supervised PCa vs. NAP FSR") + labs(colour="Type")+
  scale_colour_manual(labels = c("NAP", "PCa"), values = c("red", "lightblue")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12,face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                       axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank(),legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12), axis.line = element_line(colour = "black", size = 1.2))


##Gleason FSR Supervised##
data_fsj_gle <- data_fsj1[,match(metadata_Gleason$sample_id, colnames(data_fsj1))]
design2 <- as.factor(metadata_Gleason$gleason)
dds_fsj2 <- DESeqDataSetFromMatrix(data_fsj_gle, DataFrame(design2), ~design2)
dds_fsj2 <- DESeq(dds_fsj2)
res2 <- results(dds_fsj2)
res2 <- res2[order(res2$padj),]
sig2 <- res2[!is.na(res2$padj) & res2$padj < 0.01,]

fsj_norm_gleason <- data_fsj_norm[,match(metadata_Gleason$sample_id, colnames(data_fsj_norm))]
fsj_norm_sig_gleason<- fsj_norm_gleason[match(rownames(sig2), rownames(fsj_norm_gleason)),]
PCA_gleason_sig_FSR <- prcomp(t(fsj_norm_sig_gleason))
sup_PCAagg_non_gle_FSR <- autoplot(PCA_gleason_sig_FSR, data= metadata_Gleason, colour= "gleason", main="Supervised Ag vs. Nag FSR (Gleason Score)") + labs(colour="Gleason Score") +
  scale_colour_manual(name = "Gleason Score", labels = c("Aggressive (4+3/>7)", "Non-Aggressive (3+3/3+4)"), values = c("red", "lightblue")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12, face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                                                                                    axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12),panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))

##Metastasis FSR Supervised##
fsj_norm_metas <- data_fsj_norm[,match(metadata_metas$sample_id, colnames(data_fsj_norm))]
PCA_metas_fsr <- prcomp(t(fsj_norm_metas))
unsup_met_FSR <- autoplot(PCA_metas_fsr, data= metadata_metas, colour= "Met_Pcadeath", main="Unsupervised Ag vs. Nag FSR (Metastasis/PCa death)") + labs(colour="Met_Pcadeath") + geom_text(aes(label=ifelse(PC2>(0.22),as.character(colnames(data_fsj_norm)),'')),hjust=-0.3,vjust=0.5)+ 
  scale_colour_manual(name = "Metastasis/PCa death",labels = c("Non-Aggressive (No)","Aggressive (Yes)"), values = c("lightblue","red")) +theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5,size=12, face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                                                                                axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(),legend.text = element_text(size=12), legend.title = element_text(face="bold",size=12), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))


legend_type <- get_legend(unsup_PCApcanap + 
                            guides(color=guide_legend(nrow=1))+
                            theme(legend.position = "bottom"))
legend_gleason <-  get_legend(sup_gle_BSR + 
                                guides(color=guide_legend(nrow=1))+
                                theme(legend.position = "bottom"))
legend_met <- get_legend(unsup_met_FSR + 
                           guides(color=guide_legend(nrow=1))+
                           theme(legend.position = "bottom"))

sup_type<-plot_grid(unsup_PCApcanap+theme(legend.position="none"),sup_nappca_FSR+theme(legend.position="none"),ncol = 2,align="h",labels = c("A","D"))
sup_type_l <- plot_grid(sup_type,legend_type,nrow = 2,align="v",rel_heights = c(10,1))


sup_gle<-plot_grid(sup_gle_BSR+theme(legend.position="none"),sup_PCAagg_non_gle_FSR+theme(legend.position="none"),ncol = 2,align="h",labels = c("B","E"))
sup_gle_l <- plot_grid(sup_gle,legend_gleason,nrow = 2,align="h",rel_heights = c(10,1))


sup_met<-plot_grid(sup_PCAagg_non_met_BSR+theme(legend.position="none"),unsup_met_FSR+theme(legend.position="none"),ncol = 2,align="h",labels = c("C","F"))
sup_met_l <- plot_grid(sup_met,legend_met,nrow = 2,align="v",rel_heights = c(10,1))

plot_grid(sup_type_l,sup_gle_l,sup_met_l, align = "v", nrow = 3)
ggsave(width =10.5, height = 8,"suppfigure1_04062020.tiff")







######Figure 2D#####
data_name_bsj_sig_PCa <- bsj_norm_sig_type [,41:91]
list_cor <- rownames(data_name_bsj_sig_PCa)
list_cor <- cbind(list_cor, "cor")
list_cor <- data.frame(list_cor)

list_cor$gene_id <- rownames(data_name_bsj_sig_PCa)

list_cor1 <- str_split_fixed(list_cor$gene_id, ";", 2)
list_cor2 <- cbind(list_cor, list_cor1)
list_cor2$candidates <- rep("FALSE",times=160)
list_cor2$candidates [c(7,2,5,26,1,27,18,133,6,3,30,48,108)] <- "TRUE" ###CircRNA candidates
colnames(list_cor2) [5] <- "gene_names"

circ_candidates <- list_cor[which(list_cor$candidates == "TRUE"),]
data_name_bsj_sig <- bsj_norm_sig_type[match(rownames(circ_candidates), rownames(bsj_norm_sig_type)),]


data_fsj_circsig <- data_fsj_norm[match(rownames(data_name_bsj_sig), rownames(data_fsj_norm)),]
data_fsj_circsig1 <-data_fsj_circsig


for(i in 1:nrow (data_name_bsj_sig)) {
  circRNA_1 = rownames(data_name_bsj_sig)  [i]
  pos_bsj = data_name_bsj_sig[match(circRNA_1,rownames(data_name_bsj_sig)),]
  pos_fsj = data_fsj_circsig1[match(circRNA_1,rownames(data_fsj_circsig1)),]
  name = gsub(":",".",circRNA_1,fixed=T)
  name = gsub("|",".",name,fixed=T)
  cor_num <- round(cor(pos_bsj, pos_fsj, method="spearman"),3)
  list_cor$cor [i] <- cor_num
  
}

bsj_data_candidates <- bsj_data[match(rownames(sig_bsj), rownames(data_bsj)),]
log_bsj_raw <- log(rowMeans(bsj_data_candidates+1))
list_cor$log_bsjread <- log_bsj_raw

meanbsj_cor <- ggplot(list_cor, aes(log_bsjread,cor))+
  geom_point()+
  geom_smooth()+
  #geom_label_repel(data=subset(corplot_Scatter_mean_ratio, candidates=="TRUE"), aes(log_mean_bsj,cor,label=gene_names), box.padding   = 0.65, point.padding = 0.5, segment.color = 'grey50')+
  stat_cor(method="pearson", label.x = 0, label.y = -0.5) +
  xlim(min(list_cor$log_bsjread),max(list_cor$log_bsjread))+
  ylab("circ-linear correlation")+
  xlab("log(mean BSJ raw reads)")+
  theme(strip.background = element_blank(), strip.text = element_text(size=7, face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=14, face="bold", colour = "black"),
        axis.title=element_text(size=14,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text( size=14, face = "bold"),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))

#####Supplementary Figure 2; outlier plot########
t_data_name_bsj_sig <- t(data_name_bsj_sig)

for(i in 1:ncol(t_data_name_bsj_sig)) {
  circRNA_1 = colnames(t_data_name_bsj_sig) [i]
  print(circRNA_1)
  
  exprofcircrna = t_data_name_bsj_sig[,i] 
  exprofcircrna = data.frame(expr=exprofcircrna)
  exprofcircrna$sid = rownames(exprofcircrna)
  exprofcircrna <- merge(exprofcircrna, metadata, by.x="sid", by.y="sample_id")
  ordeer_table <- order(exprofcircrna$type, exprofcircrna$expr)
  exprofcircrna = exprofcircrna[ordeer_table,]
  
  df_scatter <- data.frame(x=exprofcircrna$sid, y=exprofcircrna$expr, z=rep(c("NAP","PCa"),times=c(40,51)))
  df_scatter$x <- as.character( df_scatter$x)
  df_scatter$x <- factor( df_scatter$x, levels=unique( df_scatter$x))
  
  ggplot(df_scatter, aes(x=x, y=y), fill=z) + 
    geom_point(data = df_scatter, aes(x = x, y = y, color = z)) + 
    labs(y="vst(BSR)", x="", color="Type")+
    scale_color_manual("Type", values = c("PCa" = "lightblue", "NAP" = "red"))+
    ggtitle(paste0(circRNA_1," normalized expression per sample"))+
    theme(strip.background = element_blank(), strip.text = element_text(size=7, face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=10, face="bold", colour = "black"),
          axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle=90, size=6, face = "plain"),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))
  
  ggsave(width = 5.92, height = 3.39,paste0("Final RData/Outlier_ggplot_RNAseqDE/Outlier_loop_",circRNA_1,"13052020.tiff"))
  
  dev.set(dev.next())
  
}

#####Fig.3B####


sig_volcano <- sig_bsj
res_volcano <- res_bsj
library(limma)
sig_volcano <- sig_volcano[order(sig_volcano$padj, sig_volcano$log2FoldChange, decreasing=F),]
sig_volcano$logpadj <- -log10(sig_volcano$padj)
sig_volcano <- data.frame(sig_volcano)
res_volcano <- res_volcano[order(res_volcano$padj, res_volcano$log2FoldChange, decreasing=F),]
res_volcano$logpadj <- -log10(res_volcano$padj)
res_volcano <- data.frame(res_volcano)
rownames(sig_volcano) == list_cor2$gene_id
sig_volcano <- sig_volcano[match(list_cor2$list_cor, rownames(sig_volcano)),]
sig_volcano$candidate <- list_cor2$candidates
sig_volcano$gene_names <- list_cor2$gene_names
sig_volcano$gene_names <- as.character(sig_volcano$gene_names)
sig_volcano$gene_names [1] <- "AMACR"
sig_volcano$gene_names [2] <- "TARGC2-TARP"
sig_volcano$gene_names[133] <- "ABCC4_1"
sig_volcano$gene_names[6] <- "ABCC4_2"
sig_volcano$gene_names[48] <- "ABCC4_3"
sig_volcano$gene_names[27] <- "MATR3"
volcano_plot_BSJ_NAPPCa<- ggplot(res_volcano, aes(x = log2FoldChange, y = logpadj), xlim = c(-6,28), pch=20, cex=0.5, main="Volcano plot") +
  geom_point(data=subset(res_volcano), color="black") +
  geom_point(data=subset(sig_volcano, padj<.01 & log2FoldChange > 1), color="green") + 
  geom_point(data=subset(sig_volcano, padj<.01 & log2FoldChange < -1), color="red") +
  geom_label_repel(data=subset(sig_volcano,candidate=="TRUE"), aes(log2FoldChange, logpadj, label=gene_names), box.padding   = 0.65, point.padding = 0.5, segment.color = 'grey50') +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  theme_linedraw() + 
  theme(panel.grid = element_blank()) +
  xlab("Fold change (log2)") +
  ylab("P-Value (log10)") +
  annotate("text", x = 21, y = 1,
           label = "padj<0.01",
           color = "black", hjust = 0) +
  annotate("text", x = 0.7, y = 24, label = "FC > 1", color = "black", angle = 90) +
  annotate("text", x = -1.3, y = 24, label = "FC < 1", color = "black", angle = 90)+theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                          axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2),legend.text = element_text(size=12))


#####correlation plot######



cormat_Down <- round(cor(t(data_norm_bsj_sig)[,77:160],t(data_fsj_circsig1)[,77:160], method="spearman"),3)
colnames(cormat_UP) <- gsub(".*;", "",colnames(cormat_UP), fixed = F)
rownames(cormat_UP) <- gsub(".*;", "",rownames(cormat_UP), fixed = F)

p.mat <- cor_pmat(cormat_UP)
tiff("part2fig2_cor.tiff")
UP_corplot <- corrplot(cormat_UP, type="upper", order="hclust",col=brewer.pal(n=8, name="RdYlBu"),tl.col="black", tl.cex = 0.3,mar=c(0,0,2,0), method = "circle",
                       title= "Spearman correlation \nof Upregulated circRNAs",
                       p.mat = p.mat,
                       sig.level = 0.01, insig = "blank")


dev.off()
pie$group <- c("Upregulated","Downregulated")
PCa_updown_pie <- ggplot(pie, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1) +
  labs(fill="")+
  coord_polar("y", start=0) + theme_void() +
  geom_text(aes(y = c(45,115), label = value), color = "white", size=8)+
  ggtitle("Differentially expressed\ncircRNAs in prostate cancer") +
  theme(legend.text=element_text(size=12),plot.title = element_text(size=12, hjust = 0.5), legend.position="bottom",legend.title = element_text(hjust = 0.5), title = element_text(size = 8))+guides(fill=guide_legend(ncol=1,byrow=TRUE))


pop <- plot_grid(PCa_updown_pie,volcano_plot_BSJ_NAPPCa,align = "h", ncol=2,nrow = 1,labels = c("A", "B"),rel_widths = c(0.5, 1),rel_heights =c(2/4, 2/4)) 
ggsave(width =8, height = 5,"part1figure2_corheatmap_04062020.tiff")


##Supplementary figure 3; violing plots


for(i in 1:ncol(t(data_name_bsj_sig))) {
  circRNA_1 = colnames(t(data_name_bsj_sig)) [i]
  print(circRNA_1)
  exprofcircrna = t(data_name_bsj_sig)[,i] 
  exprofcircrna = data.frame(expr=exprofcircrna)
  exprofcircrna$sid = rownames(exprofcircrna)
  exprofcircrna <- merge(exprofcircrna, metadata_order, by.x="sid", by.y="sample_id")
  ordeer_table <- order(exprofcircrna$type, exprofcircrna$expr)
  exprofcircrna = exprofcircrna[ordeer_table,]
  df_vio <- data.frame(x=exprofcircrna$type, y=exprofcircrna$expr)
  po <- ggplot(df_vio, aes(x=x, y=y)) + 
    geom_violin(aes(fill=x), alpha=0.5,trim = FALSE) +
    stat_summary(fun.data="mean_sdl", geom="pointrange") +
    geom_jitter(aes(fill=x), shape=24, position=position_jitter(0.2)) +
    labs(y="vst(reads)", x="",fill="Type")+
    scale_color_discrete(name="Samples",labels = c("NAP", "PCa"))+
    scale_fill_discrete(labels = c("NAP", "PCa"))+
    scale_x_discrete(labels = c("NAP", "PCa"))+
    ggtitle(paste0("circ",circRNA_1))+
    theme(strip.background = element_blank(), strip.text = element_text(size=7, face="bold"), plot.title = element_text(size=12,hjust = 0.5,face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
          axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(size=12, face = "bold"),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2),legend.title = element_text(face="bold",size=12),legend.text = element_text(size=12))+
    stat_compare_means(comparisons = ,label = "p.signif", method = "wilcox.test",
                       ref.group = "normal_adjacent_prostate",label.x = c(0.8,7.5))
  name = gsub(":",".",circRNA_1,fixed=T)
  name = gsub("|",".",name,fixed=T)
  pl[[i]] <- po

}


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(po)
pl[[14]] <- legend

vio_plots <- cowplot::plot_grid(plotlist = pl,
                                ncol=4,labels=c("A","B","C","D","E","F","G","H","I","J","K","L","M"))
title_sup3 <- ggdraw() + draw_label("CircRNA normalized expression per type of selected circRNAs", fontface='bold')

cowplot::plot_grid(title_sup3,vio_plots,ncol=1, rel_heights = c(0.05,1))
ggsave(width = 10,height = 10, "supp fig3 vio plots 04062020.tiff")

########Supplementary Figure4; Correlation analysis 13 candidates ######


pc<-vector("list", length = ncol(data_name_bsj_sig))
pcc <- gList()
data_fsj_circsig <- data_fsj_norm[match(rownames(data_name_bsj_sig), rownames(data_fsj_norm)),]
rownames(data_fsj_circsig) <- paste0(rownames(data_fsj_circsig), "fsj")

candi_circsig_bsj_fsj <- rbind(data_name_bsj_sig, data_fsj_circsig)
metadata_order$sample_id == colnames(candi_circsig_bsj_fsj)
candi_circsig_bsj_fsj_out <- candi_circsig_bsj_fsj[,match(metadata$sample_id, colnames(candi_circsig_bsj_fsj))]
candi_circsig_bsj_fsj_out<-data.frame(candi_circsig_bsj_fsj_out)
candi_circsig_bsj_fsj_out <- subset(candi_circsig_bsj_fsj_out, select = -`7046-004-047`)

ciriids1 = colnames(t(candi_circsig_bsj_fsj))[gsub("^.+(...)$","\\1",colnames(t(candi_circsig_bsj_fsj)))  != "fsj"]

for(ciriidd in ciriids1) {
  ciriidff = paste0(ciriidd, "fsj")
  
  dff1 = data.frame(y = t(candi_circsig_bsj_fsj_out)[, colnames(t(candi_circsig_bsj_fsj_out)) == ciriidd],
                    x = t(candi_circsig_bsj_fsj_out)[, colnames(t(candi_circsig_bsj_fsj_out)) == ciriidff] )
  # dff1 = dff1[dff1$x >1.11,]
  p <- ggplot(dff1, aes(x=x, y=y)) + 
    geom_point(size=2, shape=23)+
    geom_smooth(method = "gam") +
    stat_cor(method="spearman", label.x = min(dff1$x) , label.y = max(dff1$y)+1) +
    labs(title = paste0("circ",gsub(".*;", "",ciriidd, fixed = F)), x="vst(fsj reads)", y= "vst(bsj reads)")+
    theme(strip.background = element_blank(), strip.text = element_text(size=7, face="bold"), plot.title = element_text(hjust = 0.5,size=12,face="bold"), axis.text=element_text(size=12, face="bold", colour = "black"),
          axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(size=12, face = "bold"),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))
  pc [[ciriidd]] <-p
  print(ciriidff)
  #ggsave(paste0("corplots_ggplot_circRNADE/candidates/ggplot correlation plots_del 24052020",ciriidd,".tiff"))
  #dev.set(dev.next())
}
pc



title_sup4 <- ggdraw() + draw_label("Circular-Linear correlation of selected circRNAs", fontface='bold')

cor_plots <- cowplot::plot_grid(plplotlist = pc$`chr7:38249243|38262218;TARGC2-TARP`,pc$`chr5:33998641|34005899;AMACR`,pc$`chr8:39833944|39838395;ADAM2`,pc$`chr2:180724596|180725291;SCHLAP1`,
                                pc$`chr11:24905991|24983293;LUZP2`,pc$`chr7:3619080|3642105;SDK1`,pc$`chr5:139278327|139279129;MATR3`,pc$`chr13:95115922|95170628;ABCC4_3`,pc$`chr6:31863555|31869637;SLC44A4`,pc$`chr1:213077696|213117410;RPS6KC1`,pc$`chr7:24620052|24650712;MPP6`,pc$`chr13:95115922|95188542;ABCC4_2`,pc$`chr13:95161189|95188542;ABCC4_1`, 
                                ncol=4,labels=c("A","B","C","D","E","F","G","H","I","J","K","L","M"))
cowplot::plot_grid(title_sup4,cor_plots,ncol=1, rel_heights = c(0.05,1))
ggsave(width = 10,height = 12.5, "supp fig4 correlation plots_gam 09072020.tiff")


#####GAPDH whole gene counts####

linear_reads <- read.table("../../counts/expression_matrix_gencode", sep="\t",stringsAsFactors = F, header=T)

#We find GAPDH from the linear reads dataset
grep("ENSG00000111640", linear_reads$Geneid)
GAPDH_reads <- linear_reads[34984,]
write.table(GAPDH_reads, file = "../erlantz-data/Research progress presentation_plots_17042020/GAPDH_reads.txt")

names = cbind (colnames(bsj_reads), colnames(GAPDH_reads))
linear_reads <- data.frame(linear_reads)
linear_reads %>% filter(Geneid == "ENSG00000111640.14")
linear_reads[34984,]

colnames(linear_reads)
colnames(bsj_reads)
GAPDH_reads$alignments_complete.7046.004.001.merged.bam <- NULL
GAPDH_reads[, c(2,3)] <- list(NULL)
rownames(GAPDH_reads) <- "GAPDH"
GAPDH_reads$Geneid <- NULL

colnames(GAPDH_reads) <- colnames(bsj_reads)

GAPDH_reads <- cbind(X7046.004.003 = 17522, GAPDH_reads)
#We find the validated genes from bsj reads dataset
v_schlap1 <- grep("chr2:180724596|180725291;SCHLAP1", rownames(bsj_reads))
rownames(bsj_reads[v_schlap1,])
val_circ_reads [1,] <-bsj_reads[v_schlap1[3],]
colnames(val_circ_reads) == colnames(bsj_reads)
val_circ_reads <- rbind(val_circ_reads, bsj_reads[v_schlap1[3],])
val_circ_reads <- val_circ_reads[-3,]
write.table("../erlantz-data/Research progress presentation_plots_17042020/validated_circRNAs_reads.txt")
val_circ_reads_1 <- val_circ_reads + 1


for(i in 1:nrow(val_circ_reads_1)) {
  Val_bsj_GAPDH_1 [i,] <- (val_circ_reads_1[i,]/GAPDH_reads)
}

#####Correlation analysis of the genes involved in the PCa pathway#####

Linear_reads
colnames(Linear_reads) <- gsub("alignments_complete.","",colnames(Linear_reads), fixed = T)
colnames(Linear_reads) <- gsub(".merged.bam","",colnames(Linear_reads), fixed = T)
colnames(Linear_reads) <- gsub("_","-",colnames(Linear_reads), fixed = T)

Linear_reads[,1:6] <- NULL
Linear_reads$Geneid_1 <- NULL
Linear_reads2 <- Linear_reads
rownames(Linear_reads2) <- gsub("\\..*","",rownames(Linear_reads2))

library('biomaRt')
library('EnsDb.Hsapiens.v79')
library(EuPathDB)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- Linear_reads$Geneid
Linear_reads$Geneid_1<- gsub("\\..*","",Linear_reads$Geneid)

genes_PCa_target <- c("CDKN1B","CREB3L2","IGF1","PDGFD","PIK3CD","NKX3-1","AKT2","ARAF","CREBBP","E2F2","E2F3","IGF1R","KLK3","PDGFA","PDGFRA","PDGFRB","PIK3R1","RB1","BCL2","GSK3B","HSP90B1","MAPK1","MAPK3","MDM2","PDGFD","PIK3R3")
G_list_PCa <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes_PCa_target,mart= mart)
G_lis_PCa <- G_list_PCa [-8,]

design <- as.factor(metadata$type)
dds_LR <- DESeqDataSetFromMatrix(Linear_reads, DataFrame(design), ~design)
dds_LR <- DESeq(dds_LR)
res <- results(dds_LR)
res <- res[order(res$padj),]
sig <- res[!is.na(res$padj) & res$padj < 0.01,]
sig$geneid<- gsub("\\..*","",rownames(sig))

#####Fig9A#####

sig_volcano <- sig
res_volcano <- res
sig_volcano <- sig_volcano[order(sig_volcano$padj, sig_volcano$log2FoldChange, decreasing=F),]
sig_volcano$logpadj <- -log10(sig_volcano$padj)
sig_volcano <- data.frame(sig_volcano)
res_volcano <- res_volcano[order(res_volcano$padj, res_volcano$log2FoldChange, decreasing=F),]
res_volcano$logpadj <- -log10(res_volcano$padj)
res_volcano <- data.frame(res_volcano)

sig_volcano_PCa <- res2
sig_volcano_PCa <- sig_volcano_PCa[match(G_lis_PCa$ensembl_gene_id,rownames(sig_volcano_PCa)),]
sig_volcano_PCa$gene_names <- G_lis_PCa$hgnc_symbol
sig_volcano_PCa$logpadj <- -log10(sig_volcano_PCa$padj)
sig_volcano_PCa <- data.frame(sig_volcano_PCa)

volcano_pca <- ggplot(res_volcano, aes(x = log2FoldChange, y = logpadj), xlim = c(-6,28), pch=20, cex=0.5, main="Volcano plot") +
  geom_point(data=subset(res_volcano), aes(color="black")) +
  geom_point(data=subset(sig_volcano, padj<.01 & log2FoldChange > 1), aes(color="green")) + 
  geom_point(data=subset(sig_volcano, padj<.01 & log2FoldChange < -1), aes(color="red")) +
  geom_point(data=subset(sig_volcano_PCa), aes(color="blue")) +
  scale_color_identity(name = "",
                       breaks = c("black","green", "red","blue"),
                       labels = c("FDR > 0.01; \n-1 > FC < +1", "Upregulated \nin PCa","Downregulated \nin PCa","Genes \nin PCa pathway \n(KEGG)"),guide="legend")+
  geom_label_repel(data=subset(sig_volcano_PCa,padj<.01 & log2FoldChange < -1 | log2FoldChange > 1), aes(log2FoldChange, logpadj, label=gene_names), box.padding   = 4, point.padding = 0, segment.color = 'grey50') +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = 2, linetype = "dashed") +
  theme_linedraw() + 
  theme(panel.grid = element_blank(),legend.position = "bottom",legend.text = element_text(size=12)) +guides(color=guide_legend(ncol=2,byrow=TRUE))+
  xlab("Fold change (log2)") +
  ylab("P-Value (log10)") +
  annotate("text", x = 5, y = 0.8,
           label = "padj<0.01",
           color = "black", hjust = 0) +
  annotate("text", x = 1.2, y = 55, label = "FC > 1", color = "black", angle = 90) +
  annotate("text", x = -1.3, y = 55, label = "FC < 1", color = "black", angle = 90)+theme(strip.background = element_blank(), strip.text = element_text(size=9, face="bold"), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12, face="bold", colour = "black"),
                                                                                          axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))
ggsave(width = 5.92, height = 3.39,"Final RData/volcano_linear_PCatargets_circRNas 16052020.tiff")

##circRNA BSR and linear Reads
rownames(res)<- gsub("\\..*","",rownames(res))
res2 <- res[match(G_lis_PCa$ensembl_gene_id,rownames(res)),]
G_list_PCa$ensembl_gene_id %in% rownames(res)
res2<- res2[order(res2$padj),]
Linear_reads2_PCa_target <- Linear_reads2[match(rownames(res2),rownames(Linear_reads2)),]

#Normalized reads through vst
norm_LR <- vst(assay(dds_LR))
rownames(norm_LR)<- gsub("\\..*","",rownames(norm_LR))
norm_LR_PCa <- norm_LR[match(rownames(Linear_reads2_PCa_target),rownames(norm_LR)),]



norm_circ_bsj <- read.table("table_normbsj_sig.csv", stringsAsFactors = F, header = T, sep="\t")
colnames(norm_circ_bsj) <- gsub("X","",colnames(norm_circ_bsj),fixed = T)
colnames(norm_circ_bsj) <- gsub(".","-",colnames(norm_circ_bsj),fixed = T)
norm_circ_bsj_4 <- norm_circ_bsj [c("chr6:31863555|31869637;SLC44A4",
                                    "chr2:180724596|180725291;SCHLAP1",
                                    "chr13:95161189|95188542;ABCC4",
                                    "chr11:24905991|24983293;LUZP2"),]

norm_LR_PCa
rownames(sig_volcano_PCa2) == rownames(norm_LR_PCa)
rownames(norm_LR_PCa) <- sig_volcano_PCa2$gene_names

colnames(norm_circ_bsj_4) == colnames(norm_LR_PCa)
norm_circ_bsj_4 <- norm_circ_bsj_4[match(metadata$sample_id,colnames(norm_circ_bsj_4))]
df_coR_norm<-rbind(norm_circ_bsj_4,
                   norm_LR_PCa)
metadata$sample_id == rownames(t_df_coR_norm)

t_df_coR_norm <- t(df_coR_norm)
colnames(t_df_coR_norm) <- gsub(".*;","",colnames(t_df_coR_norm))
dim(t_df_coR_norm)
colnames(t_df_coR_norm) [1:4] <- c("circSLC44A4","circSChLAP1","circABCC4_1","circLUZP2")
t_df_coR_norm <- t_df_coR_norm[,-3]
cormat_norm <- round(cor(t_df_coR_norm[41:91,1:3],t_df_coR_norm[41:91,4:28], method="spearman"),3)
t_cormat_norm <- t(cormat_norm)
reorder_cormat <- function(t_cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-t_cormat)/2)
  hc <- hclust(dd)
  t_cormat <-t_cormat[hc$order, hc$order]
}
reorder_cormat(t_cormat)
spe_cor_norm <- ggcorrplot::ggcorrplot(t_cormat_norm, type = "full",method="square",hc.order = F) + geom_tile()+
  scale_fill_gradient2(low = "red", high = "green", mid = "black", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation")+
  theme(strip.background = element_blank(),  plot.title = element_text(size=10,face ="bold", hjust = 0.5), axis.text=element_text(size=10, face="bold", colour = "black"), axis.text.x = element_text(angle=90, size=8, face = "bold"),axis.text.y = element_text( size=8, face = "bold")) +
  ggtitle("Spearman correlation\ncircRNA-PCa pathway target genes")

ggsave(width = 10, height = 5,"correlation heatmap plot 4circRNAs_target_linearPCagenes_PCasamples_NORM 04062020.tiff")

cor_up <- t_cormat_norm[t_cormat_norm > 0.39]
cor_up
t_df_coR_norm <- data.frame(t_df_coR_norm)
t_df_coR_norm_Pca <- t_df_coR_norm [41:91,]
SLC_E2F2 <- ggplot(t_df_coR_norm_Pca, aes(x=E2F2, y=SLC44A4))+ 
  geom_point(size=2, shape=23)+
  geom_smooth(method="lm") +
  stat_cor(method="spearman", label.x = min(t_df_coR_norm_Pca$E2F2), label.y = max(t_df_coR_norm_Pca$SLC44A4)+1) +
  labs(title = paste0("circSLC44A4-E2F2"), x="Vst(Linear gene reads)", y= "Vst(circRNA reads)")+
  theme(strip.background = element_blank(), strip.text = element_text(size=12, face="bold"), plot.title = element_text(hjust = 0.5,size=12), axis.text=element_text(size=12, face="bold", colour = "black"),
        axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))

SLC_E2F3 <- ggplot(t_df_coR_norm_Pca, aes(x=E2F3, y=SLC44A4))+ 
  geom_point(size=2, shape=23)+
  geom_smooth(method="lm") +
  stat_cor(method="spearman", label.x = min(t_df_coR_norm_Pca$E2F3), label.y = max(t_df_coR_norm_Pca$SLC44A4)+1) +
  labs(title = paste0("circSLC44A4-E2F3"), x="Vst(Linear gene reads)", y= "Vst(circRNA reads)")+
  theme(strip.background = element_blank(), strip.text = element_text(size=12, face="bold"), plot.title = element_text(hjust = 0.5,size=12), axis.text=element_text(size=12, face="bold", colour = "black"),
        axis.title=element_text(size=12,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1.2))
