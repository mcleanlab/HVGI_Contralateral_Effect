#The Following R Code is related to manuscript:
#Localized microbially induced inflammation influences distant healthy tissues in the human oral cavity by Kerns et al., 2023 PNAS

#Activate packages 
library("ggformula")
library("ggExtra")
library("ggpubr")
library("ggpmisc")
library("dplyr")
library("phyloseq")
library("ggplot2")
library("ggridges")
library("viridis")
library("Rmisc")
library("microbiome")
library("gridExtra")
library("grid")
library("ggplotify")
library("radiant.basics")
library("rstatix")  
library("ggpol")
library("speedyseq")
library("vegan")
library("circlize")
library("ComplexHeatmap")
library("metagMisc")
library("clustvis")
library("pheatmap")
library("factoextra")
library("FactoMineR")
library("gghighlight")
library("gtable")
library("grid")
library("gridExtra")

#All necessary raw files to recreate these analysis are publicly available under the Data folder of the dedicated GitHub Repository.
#https://github.com/kkerns85/HVGI_Contralateral_Effect

########################
#### Processing DATA ####
########################
HVGI_biom = "/HVGI_Contralateral.biom"
HVGI = import_biom(HVGI_biom) 
HVGI 

HVGI_MapFile = import_qiime_sample_data("/HVGI_Contralateral_Mapfile.csv")

RootedTree = "/HVGI_Contralateral_rooted_tree.nwk"
HVGI_Tree <- read_tree(RootedTree)
HVGI_Tree

#Add Phylogeny as Column names
colnames(tax_table(HVGI)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Generate Phyloseq Object
HVGI_po = merge_phyloseq(HVGI,HVGI_MapFile,HVGI_Tree)
HVGI_po

#Filter Unassigned
HVGI_po_f <- subset_taxa(HVGI_po, Phylum != "Unassigned")
HVGI_po_f

#Filter zero's
HVGI_po_f <- prune_taxa(taxa_sums(HVGI_po_f) > 0, HVGI_po_f)
HVGI_po_f 

#Filter Contaminants
Contam <- c("g__Ralstonia","g__Delftia","g__Stenotrophomonas","g__Bradyrhizobium","g__Brevundimonas","g__Bacillus",
            "g__Pseudomonas","g__Sphingomonas","g__Enterobacter","g__Microbacterium","g__Bosea","g__Kocuria","g__Micrococcus",
            "g__Paracoccus","g__Anoxybacillus","g__Leptothrix","g__Lawsonella","g__Gardnerella","g__Arthrospira",
            "g__Agrobacterium","g__Escherichia","g__Porphyrobacter","g__Ochrobactrum","g__Neisseriaceae_[G-1]","g__Flavitalea",
            "g__Finegoldia","g__Bdellovibrio","g__Mesorhizobium","g__Defluvibacter","g__Acinetobacter")

HVGI_po_f <- subset_taxa(HVGI_po_f, !(Genus %in% Contam))
HVGI_po_f 

#filter singletons
HVGI_po_fs <- prune_taxa(taxa_sums(HVGI_po_f) > 1, HVGI_po_f)
HVGI_po_fs 

#translate by adding pseudo count
HVGI_po_fs_t = transform_sample_counts(HVGI_po_fs, function(x) (x + 0.000001 - min(x)))

#transform to relative abundance with translated pseudo count
HVGI_po_fs_t_n = transform_sample_counts(HVGI_po_fs_t, function(x) (x)/sum(x))

#agglomerate data
HVGI_po_fs_t_n_pglom <- tax_glom(HVGI_po_fs_t_n, taxrank = "Phylum", NArm = TRUE)
HVGI_po_fs_t_n_gglom <- tax_glom(HVGI_po_fs_t_n, taxrank = "Genus", NArm = TRUE)
HVGI_po_fs_t_n_sglom <- tax_glom(HVGI_po_fs_t_n, taxrank = "Species", NArm = TRUE)

#convert data to df
HVGI_po_fs_t_n_pglom_df <- psmelt(HVGI_po_fs_t_n_pglom)

HVGI_po_fs_t_n_gglom_df <- psmelt(HVGI_po_fs_t_n_gglom)

HVGI_po_fs_t_n_sglom_df <- psmelt(HVGI_po_fs_t_n_sglom)

#########################
#### Figure 1 ###########
#########################
#Plaque Index (PI)
c1 <-ggplot(HVGI_MapFile, aes(x = Day, y = Bacterial_Load))+ xlab("Day") + ylab("Bacterial Load") + 
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))  + 
  stat_summary(na.rm=T,fun.y=mean, aes(group=Responders),colour = my_colors) +
  scale_fill_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  annotate("rect", xmin = 1.6, xmax = 6.4, ymin = -Inf, ymax = Inf, alpha=0.1) +
  scale_color_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") + geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 2.5, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "none")
c1


#Gingival Index (GI)
c2 <-ggplot(HVGI_MapFile, aes(x = Day, y = GI))+ xlab("Day") + ylab("Gingival Index") + 
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))  + 
  stat_summary(na.rm=T,fun.y=mean, aes(group=Responders),colour = my_colors) +
  scale_fill_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  annotate("rect", xmin = 1.6, xmax = 6.4, ymin = -Inf, ymax = Inf, alpha=0.1) +
  scale_color_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") + geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 2.5, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "none")
c2

#%Bleeding on Probing (BOP)
c3 <-ggplot(HVGI_MapFile, aes(x = Day, y = BOP))+ xlab("Day") + ylab("Bleeding on Probing") + 
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))  + 
  stat_summary(na.rm=T,fun.y=mean, aes(group=Responders),colour = my_colors) +
  scale_fill_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  annotate("rect", xmin = 1.6, xmax = 6.4, ymin = -Inf, ymax = Inf, alpha=0.1) +
  scale_color_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") + geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 2.5, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "none")
c3

#GCF Volume
c4 <-ggplot(HVGI_MapFile, aes(x =Day, y = GCF_Volume))+ xlab("Day") + ylab("GCF Volume") + 
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))  + 
  stat_summary(na.rm=T,fun.y=mean, aes(group=Responders),colour = my_colors) +
  scale_fill_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  annotate("rect", xmin = 1.6, xmax = 6.4, ymin = -Inf, ymax = Inf, alpha=0.1) +
  scale_color_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") + geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 2.5, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "none")
c4

#Bacterial Load
c5 <-ggplot(HVGI_MapFile, aes(x = Day, y = Bacterial_Load))+ xlab("Day") + ylab("Bacterial_Load") + 
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8))  + 
  stat_summary(na.rm=T,fun.y=mean, aes(group=Responders),colour = my_colors) +
  scale_fill_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  annotate("rect", xmin = 1.6, xmax = 6.4, ymin = -Inf, ymax = Inf, alpha=0.1) +
  scale_color_manual(values = c("#B2474599","#B24745FF", "#00A1D599", "#DF8F44FF"), name="IRT Group") + geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 2.5, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "none")
c5
#########################
#### Figure 2 ###########
#########################
#Alpha Diversity - Use Singletons
HVGI_po_f

alpha_p  <- plot_richness(HVGI_po_f, x="Day", color="IRT", measures = c("Observed", "Chao1", "Shannon","Simpson"))
Alpha_data <- alpha_p$data

a1 <- ggplot(Alpha_data, aes(x = Day , y= value, fill = IRT)) + geom_boxplot() + facet_wrap(~variable, scales = "free_y") + 
      geom_smooth(aes(group = IRT, color = IRT), size = 1, method = "loess", se = FALSE, alpha = 0.1)

#Beta Dissimilarity
beta_dis <- read.csv("/HVGI_Contralateral_b_diss")
beta_dis
levels(beta_dis$Day)
beta_dis$Day <- factor(beta_dis$Day, levels = c("Day_0_4", "Day_0_7", "Day_0_14", "Day_0_21"))
levels(beta_dis$Day)

bd_0_21 <- ggplot(beta_dis, aes(x = Day, y = b_dis, fill = IRT)) + geom_boxplot() + ylab('Bray Curtis Dissimilarity') + xlab('Responder Group') + facet_grid(~IRT) +
  ggtitle(label = "Bray Curtis Dissimilarity By Responder Group Day 0-21") + theme_minimal() + geom_smooth(aes(group = IRT, fill = IRT),method = "loess", se = FALSE) +
  theme(legend.position = "bottom") 
bd_0_21

#Beta Diversity - Filtered Singletons
HVGI_po_fs_t_n

ordu_pufw <- ordinate(HVGI_po_fs_t_n, "PCoA", "unifrac", weighted=TRUE)
pw = plot_ordination(HVGI_po_fs, ordu_pufw, type="samples", color="IRT") 
PCoA_w <- pw$data

b1 <- ggplot(data=PCoA_w,aes(x = Axis.1, y = Axis.2, fill = Day, colour = Day))
b1 <- b2 +geom_point()+ theme_minimal()
b1 <- b2 + theme(legend.position = "bottom") + ggtitle(label = "PCoA Unweighted Unifrac Over the Induction Phase - IRT")
b1 <- ggMarginal(b1,  type = c("boxplot"), margins = c("both"), size = 10, xparams = list(), yparams = list(), groupColour = TRUE, groupFill = TRUE)
b1

#########################
#### Figure 3 ###########
#########################
#Phylum 
HVGI_po_fs_t_n_pglom_df

p1 <- ggplot(HVGI_po_fs_t_n_pglom_df, aes(x = Day, y = Abundance ))+ xlab("Day") + ylab("Relative Abundance") + facet_wrap( ~ Phylum, scales = "free") +
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8)) + 
  scale_fill_manual("Response Group",values = c("#F4A582","#854C4B","#92C5DE","#8BAAEA","#DFC27D","#F4B951")) +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  scale_color_manual(values = c("#F4A582","#854C4B","#92C5DE","#8BAAEA","#DFC27D","#F4B951"), name="Response Group") + geom_boxplot(aes(fill = Responders_c)) +
  geom_smooth(aes(group = Responders_c, color = Responders_c), size = 1, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "bottom") + 
  ggtitle(label = "Relative Abundance by Phylum over Induction") 
p1 

#Genus LFC Heatmaps
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "firebrick3"))
col_fun = colorRamp2(c(0, 2), c("white", "firebrick3"))

#use respective LFC file for each IRT designation
matrix <- read.delim(file = "/HVGI_Contralateral_X_LFC.txt", sep = "\t")

rownames(matrix) <- matrix[,1] 
matrix <- matrix[,-1]
matrix <- as.matrix(matrix)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("lightblue", "white", "firebrick3"))

h1 <- Heatmap(matrix, name = "IRT ~ Day 0-21",width = unit(3, "cm"), col = col_fun, rect_gp = gpar(col = "white", lwd = 0),
              column_title = "IRT X", #Test or Control
              column_order = colnames(matrix),
              row_names_gp = gpar(fontsize = 5))

h1

#LFC Summary
#bake barchart and att the 
LFC_g <- matrix <- read.csv(file = "/HVGI_Contralateral_LFC_Genera.csv")
LFC_g

#wide to long
LFC_g_l <- melt(LFC_g, id.vars=c("Responder","Group"),
                variable.name="LFC_Comp",
                value.name="Genera_w_gt1LFC")
LFC_g_l

lfc1 <- ggbarplot(LFC_g_l, "LFC_Comp", "Genera_w_gt1LFC", facet.by = "Responder",
                 fill = "Group", color = "Group",
                 label = FALSE,
                 position = position_dodge(0.9)) + 
  #geom_line(aes(x=LFC_Comp, y=Genera_w_gt1LFC, group = Group, color = Group), stat="identity", size = 2)+
  scale_fill_manual("Response Group",values = c("#854C4B","#F4A582","#92C5DE","#8BAAEA","#DFC27D","#F4B951")) +
  scale_color_manual("Response Group",values = c("#854C4B","#F4A582","#92C5DE","#8BAAEA","#DFC27D","#F4B951")) + 
  rotate_x_text(angle = 45) +
  theme(legend.position = "bottom")
lfc1

#########################
#### Figure 4 ###########
#########################
HVGI_po_fs

#filter taxa of interest
toi <- c("g__Prevotella","g__Alloprevotella","g__Selenomonas","g__Treponema","g__Saccharibacteria_(TM7)_[G-1]", "g__Aggregatibacter","g__Tannerella","g__Porphyromonas")

HVGI_po_fs_toi <- subset_taxa(HVGI_po_fs, Genus %in% toi)
HVGI_po_fs_toi

#convert to presece absence 
HVGI_po_fs_toi.pa <- phyloseq_standardize_otu_abundance(HVGI_po_fs_toi, method = "pa")

#investigate by Subject Test and Control Sites
TAX <- as.data.frame(tax_table(HVGI_po_fs_toi.pa))
head(TAX)
TAX$Genus <- gsub(pattern = "g__", replacement = "", TAX$Genus)
TAX$Species <- gsub(pattern = "s__", replacement = "", TAX$Species)
TAX$Taxonomy <- paste(TAX$Genus,TAX$Species, sep = "_")
TAX$ID <- paste(TAX$Taxonomy,TAX$ASV, sep = "_")
TAX

SAM <- as.data.frame(sample_data(HVGI_po_fs_toi.pa))
SAM

COUNT <- as.dataframe(otu_table(HVGI_po_fs_toi.pa))
head(COUNT)

write.csv(COUNT, file = "/Count.csv")
write.csv(TAX, file = "/Tax.csv")
write.csv(SAM, file = "/Sample.csv")

#partition by Subject Test and Control Side
matrix <- read.delim(file = "/Subject_X.txt", sep = "\t")
rownames(matrix) <- matrix[,1] 
matrix <- matrix[,-c(1)]
matrix

#filter 0's
matrix <- matrix[rowSums(matrix[])>0,]
matrix

matrix <- as.matrix(matrix)

col_fun = colorRamp2(c(0,1), c("grey", "firebrick3"))

Toi_hm <- Heatmap(matrix, name = "Presence Absence",width = unit(3, "cm"), col = col_fun , rect_gp = gpar(col = "white", lwd = 0),
                  column_title = "Samples",
                  row_title = "ASV ID",
                  column_order = colnames(matrix),
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize = 3),
                  column_names_gp = gpar(fontsize = 3))
Toi_hm

#########################
#### Figure 5 ###########
#########################
#Host Mediators
HVGI_MapFile

HVGI_MapFile_long <- melt(HVGI_MapFile,
                          id.vars=c("X.SampleID","Day", "Subject","Sex","Age","Bacterial_Load","IRT"),
                          variable.name="condition",
                          value.name="measurement")
HVGI_MapFile_long

#Repeat for mediator of Interest (moi)
moi <- ggplot(HVGI_MapFile_long, aes(x = Day, y = moi)) + xlab("Day") + ylab("Chemokine Value (pg/30s)") +
  theme(text=element_text(size = 8,face="bold")) + 
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 8)) + 
  scale_fill_manual("IRT Group",values = c("#A68686","#541C1A","#92C5DE","#5D7CBD","#DFC27D","#C48B35")) +
  scale_color_manual("IRT Group",values = c("#A68686","#541C1A","#92C5DE","#5D7CBD","#DFC27D","#C48B35")) +
  labs(fill="")+  theme(legend.text = element_text(colour="Black", size = 4)) + 
  theme(legend.title = element_text(colour="Black", size = 4, face = "bold")) + 
  geom_boxplot(aes(fill = IRT)) +
  geom_smooth(aes(group = IRT, color = IRT), size = 1, method = "loess", se = FALSE, alpha = 0.1) +
  theme_classic() + theme(legend.position = "bottom") +
  facet_grid(~IRT) +
  scale_y_log10()
moi

#Heatmaps

#IRT of interest
file = "/Cont_heatmap_clustvis_High.txt"
imp = importData(file)
proc = processData(imp)
?generateHeatmap
RColorBrewer::brewer.pal.info

#Find optimal number of Clusters
nbr_clust <- proc$mat

set.seed(123)
fviz_nbclust(nbr_clust, kmeans, nstart = 25, k.max = 10, method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method - High Responder Controls")

#plot heatmap
hm = generateHeatmap(proc, showImputed = TRUE, transpose = FALSE, scale = "row",
                     clustDistRows = "correlation", clustMethodRows = "average",
                     clustDistCols = NA,
                     treeOrderingRows = NA, nbrClustersRows = 6,
                     colorAnnoRow = NA, legendColorScheme = "Set1", plotWidth = 25,
                     plotRatio = 1.5, colorRangeMin = -1, colorRangeMax = 2,
                     matrixColorScheme = "OrRd", revScheme = FALSE,
                     cellBorder = "grey60", fontSizeGeneral = 14, showNumbers = FALSE,
                     fontSizeNumbers = 10, precisionNumbers = 2, showRownames = TRUE,
                     fontSizeRownames = 8, showColnames = TRUE, fontSizeColnames = 14,
                     showAnnoTitlesRow = TRUE, showAnnoTitlesCol = TRUE,
                     maxAnnoLevels = 50)
hm

#export hm data and add metadata and then bring back in 
hm_mat <- hm$cells
hm_mat
write.csv(hm_mat, file = "/mat_out.csv")

hm_mat2 <- read.csv("/HC_mat_in.csv")
hm_mat2

#wide to long
hm_mat2_l <- melt(hm_mat2, id.vars=c("Group_ID", 'IRT', 'Assignment','Day'),variable.name="Chemokine",
                  value.name="Distance")
hm_mat2_l$Chemokine

#sort Chemokine by assigned optimal clusters
#make line plots for mediator changes by cluster
hm_mat2_l_c <- filter(hm_mat2_l,Chemokine %in% cluster_X)

c1 <- ggline(hm_mat2_l_c, x = "Day", y = "Distance", color = "blue",group = "Chemokine", facet.by = "Responders_c", 
             add = c("mean_se"), title = "Cluster X") #color = "Chemokine",group = "Chemokine",
c1

c2<- ggline(hm_mat2_l_c, x = "Day", y = "Distance", color = "Chemokine",group = "Chemokine", facet.by = "Responders_c", 
            add = c("mean_se"), title = "Cluster X ID") + theme(legend.position = "right") #color = "Chemokine",group = "Chemokine",
c2

grid.arrange(c1,c2,ncol = 2)
#########################
#### Figure 6 ###########
#########################
#Load the LFC data set for each IRT group
file = "//HVGI_Contralateral_X_Mediator_LFC.txt"
imp = importData(file)
proc = processData(imp)

hm = generateHeatmap(proc, showImputed = TRUE, transpose = FALSE, scale = "row",
                     clustDistRows = "euclidean", 
                     clustMethodRows = "average",
                     clustDistCols = NA,
                     treeOrderingRows = NA, nbrClustersRows = 6,
                     colorAnnoRow = NA, legendColorScheme = "Set1", plotWidth = 25,
                     plotRatio = 1.5, colorRangeMin = -1, colorRangeMax = 1,
                     matrixColorScheme = "RdBu", revScheme = TRUE,
                     cellBorder = "grey60", fontSizeGeneral = 14, showNumbers = FALSE,
                     fontSizeNumbers = 10, precisionNumbers = 2, showRownames = TRUE,
                     fontSizeRownames = 8, showColnames = TRUE, fontSizeColnames = 8,
                     showAnnoTitlesRow = TRUE, showAnnoTitlesCol = TRUE,
                     maxAnnoLevels = 50)
hm

#Dual Y Axis Plots Mediator of Interest (moi)
HVGI_MapFile

#select Test and Control of IRT of interest
roi <- c("IRT Test","IRT Control")
HVGI_MapFile_R <- filter(HVGI_MapFile, IRT == "High Responders Test")

#convert wide to long

HVGI_MapFile_R_long_moi1 <- filter(HVGI_MapFile_R_long, mediator =="IL8")
HVGI_MapFile_R_long_moi2 <- filter(HVGI_MapFile_R_long, mediator == c("IL6","TNFa"))

dy1 <- ggline(HVGI_MapFile_R_long_moi1, x = "Day", y = "Abundance", color = "IRT",
              add = c("mean_se"), scales = "free", ylab = "IL8",repel = TRUE, show.line.label = TRUE,
              linetype = "twodash", shape = "IRT", palette = c("blue","blue"), title = "IRT ~ Induction")

dy1 <- dy1 + annotate("rect", xmin = 0.8, xmax = 2.2, ymin =0.20, ymax = 0.55, alpha = .1, fill = "red")
dy1

dy2<- ggline(HVGI_MapFile_R_long_moi2, x = "Day", y = "Abundance", color = "IRT",
             add = c("mean_se"), scales = "free",ylab = "IL6 and TNFa",
             linetype = "twodash", shape = "IRT", palette = c("red"))

ggplot_dual_axis(dy1,dy2)

#Dual Y Axis Plots Phylum
HVGI_po_fs_t_n_pglom_df

#select Test and Control of IRT of interest
roi <- c("IRT Test","IRT Control")
HVGI_po_fs_t_n_pglom_df_R <- filter(HVGI_po_fs_t_n_pglom_df, IRT == "High Responders Test")

SB_asv2_fs_t_n_pglom_df_Bact <- filter(HVGI_po_fs_t_n_pglom_df_R, Phylum =="p__Bacteroidetes")
SB_asv2_fs_t_n_pglom_df_Firm <- filter(HVGI_po_fs_t_n_pglom_df_R, Phylum =="p__Firmicutes")

dy3 <- ggline(SB_asv2_fs_t_n_pglom_df_Firm, x = "Day", y = "Abundance", color = "IRT",
             add = c("mean_se"), scales = "free", ylab = "Relative Abundance Firmicutes",repel = TRUE, show.line.label = TRUE,
             linetype = "twodash", shape = "IRT", palette = c("blue","blue"), title = "IRT ~ Induction")

dy3 <- dy3 + annotate("rect", xmin = 0.8, xmax = 2.2, ymin =0.20, ymax = 0.55, alpha = .1, fill = "red")
dy3

dy4<- ggline(SB_asv2_fs_t_n_pglom_df_Bact, x = "Day", y = "Abundance", color = "IRT",
            add = c("mean_se"), scales = "free",ylab = "Relative Abundance Bacteroidetes",
            linetype = "twodash", shape = "IRT", palette = c("red","red"))

ggplot_dual_axis(dy3,dy4)

