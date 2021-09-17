##################################################
# For stickleback M-F data, analyse 10 SNP AF-vapeR res
##################################################
rm(list=ls())
lib <- c("cowplot","ghibli","ggtree","tidyverse","ggplot2","parallel","regioneR","bedr","pegas","biomaRt","ggridges")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")

# Fetch and analyse afvaper -----------------------------------------------
# Read in data
window_snps <- 10

# Set chrs
chrs <- paste0("group",as.roman(1:21))
chrs <- chrs[chrs != "groupXIX"]

# Read in the rds results if needs be...
chr_res <- lapply(chrs,function(chr){
  #  print(chr)
  return(readRDS(paste0("outputs/stickleback_rad_test_",chr,"_AF_eigen_res_windsize_",window_snps,".rds")))
})

# Just get the chr AF_inputs
AF_input_chrs <-  merge_eigen_res(lapply(chr_res,"[[",1))

# Pull all the eigen_res together
eigen_res_list_test <- lapply(chr_res,"[[",3)
allchr_eigen_res <- merge_eigen_res(eigen_res_list_test)

# Visualise the tree ------------------------------------------------------
# Tree comes from Magalhaes and Whiting et al. 2021 Nature Ecology & Evolution
stickleback_tree <- read.tree(text="(LICA,GARD,NORT,ERRO,BULL,KENN,DOUG,KLEN,STOW,TROU,BRAN,HOGG,KIRK,BEAV,SPRO,LILY,CRAN,(HOTE,AMBR));")

# Inspect
base_tree <- ggtree(stickleback_tree) + layout_inward_circular(xlim=10)

# Make vectors again
BC_pops <- rownames(AF_input_chrs[[1]])
vectors <- lapply(BC_pops,function(x) c("LICA",x))

# Add the vectors...
for(i in 1:length(vectors)){
  print(i)
  base_tree <- base_tree + 
    geom_taxalink(vectors[[i]][1],vectors[[i]][2],
                  arrow=arrow(length=unit(0.025, "npc"),type = "closed",ends = "last"),
                  colour="orange3",alpha=0.75,
                  curvature = -.6) 
}
base_tree <- base_tree +
  geom_tiplab(hjust=1.5)

# Set up phenotype table
phenos <- c("Marine","Freshwater")
phenotypes <- matrix(ncol=1,nrow=length(stickleback_tree$tip.label))
colnames(phenotypes) <- "phenotype"
rownames(phenotypes) <- stickleback_tree$tip.label
for(i in 1:length(phenos)){
  phenotypes[rownames(phenotypes) %in% sapply(vectors,'[[',i),"phenotype"] <- phenos[i]
}

# Add to tree
final_pheno_tree <- gheatmap(base_tree, phenotypes, width=0.1,colnames=FALSE)+
  scale_fill_manual(breaks = phenos,
                    values = ghibli_palettes$LaputaMedium[c(4,7)])+
  scale_x_ggtree()+
  theme(legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size=13))+
  labs(fill="")


# Analysis of results -----------------------------------------------------
# Distribution of window sizes...
window_sizes <- sapply(names(allchr_eigen_res),function(x){
  wind <- strsplit(x,":")[[1]][2]
  as.integer(strsplit(wind,"-")[[1]][2])-as.integer(strsplit(wind,"-")[[1]][1])
})
med_wind_size <- median(window_sizes)
med_wind_size

# Make a single null cutoff
null_input_list <- lapply(chr_res,"[[",2)
all_nulls <- unlist(null_input_list, recursive=FALSE)
all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                     cutoffs = c(0.95,0.99,0.999))

# Fetch all our significant windows...
chr_signif <- afvaper::signif_eigen_windows(allchr_eigen_res,cutoffs = all_null_cutoffs[,2])

# We want to look at eigenvector 2
eig1_signif <- chr_signif[[1]]
eig2_signif <- chr_signif[[2]]

# Also just fetch all pvalues for all windows...
all_pvals <- eigen_pvals(allchr_eigen_res,all_nulls)
eig1_pvals <- data.frame(window_id = rownames(all_pvals),
                         eig1_p = -log10(all_pvals[,1]))
eig1_pvals <- eig1_pvals[order(-eig1_pvals$eig1_p),]

# And see all that are over 2...
eig1_pvals[eig1_pvals$eig1_p > 2,]
nrow(eig1_pvals[eig1_pvals$eig1_p > 2,])

# What is happening at the eda locus?
all_pvals["groupIV:12797185-12829813",]
summarise_window_parallelism("groupIV:12797185-12829813",allchr_eigen_res,loading_cutoff = 0.2,eigenvector = 4)
allchr_eigen_res[["groupIV:12797185-12829813"]]$A_matrix

# And get the eig2 windows...
eig2_pvals <- data.frame(window_id = rownames(all_pvals),
                         eig2_p = -log10(all_pvals[,2]))
eig2_pvals <- eig2_pvals[order(-eig2_pvals$eig2_p),]

# Summarise eigenvector 1 windows...
if(!(is.na(eig1_signif[1]))){
  stickleback_eigenvec1_window_summaries <- afvaper::summarise_window_parallelism(eig1_signif,allchr_eigen_res,loading_cutoff = 0.2,eigenvector = 1)
} else {
  stickleback_eigenvec1_window_summaries <- NA
}

# Summarise eigenvector 2 windows...
if(!(is.na(eig2_signif[1]))){
  stickleback_eigenvec2_window_summaries <- afvaper::summarise_window_parallelism(eig2_signif,allchr_eigen_res,loading_cutoff = 0.2,eigenvector = 2)
} else {
  stickleback_eigenvec2_window_summaries <- NA
}

# Change group names to chr names
new_chr_names <- paste0("chr",1:21)
names(new_chr_names) <- paste0("group",as.roman(1:21))
eigen_names <- names(allchr_eigen_res)
for(i in 1:length(new_chr_names)){
  eigen_names <- gsub(paste0(names(new_chr_names)[i],":"),
                      paste0(new_chr_names[i],":"),eigen_names)
}
names(allchr_eigen_res) <- eigen_names
# Make whole genome figures...
genome_figs <- afvaper:::eigenval_plot(eigen_res=allchr_eigen_res[grep("chr",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T,keep_chr_scale = T)
genome_figs[[1]]
genome_figs[[2]]
genome_figs[[3]]
genome_figs[[4]]

# And for chrIV with eda marked...
eda_region <- c(12800220,12810446)
groupIV_fig <- afvaper:::eigenval_plot(eigen_res=allchr_eigen_res[grep("chr4:",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T,keep_chr_scale = T)
groupIV_fig <- groupIV_fig[[1]] + 
  geom_vline(xintercept = eda_region,colour="red2",linetype="dashed")+
  annotate("text", x = 10e6, y = 2.5, label = "italic(eda)~gene",parse = TRUE)


# Show how the landscape varies across chromosomes....
per_chrom_eig1 <- data.frame(rbindlist(lapply(new_chr_names[new_chr_names != "chr19"],function(chr){
  
  # Subset
  tmp_res <- allchr_eigen_res[grep(paste0(chr,":"),names(allchr_eigen_res))]
  
  # Return first eigenvalue
  eig1_vals <- sapply(tmp_res,function(x) return(x$eigenvals[1]))
  data.frame(chr=chr,
             eig1=eig1_vals)
})))

# Plot
per_chrom_medians <- data.frame(per_chrom_eig1 %>% group_by(chr) %>% summarise(mean_eig1=median(eig1)))
per_chrom_medians <- per_chrom_medians[order(-per_chrom_medians$mean_eig1),]
per_chrom_eig1$chr_F <- factor(gsub("chr","",per_chrom_eig1$chr),levels=gsub("chr","",rev(per_chrom_medians$chr)))

per_chrom_ridges <- ggplot(per_chrom_eig1,aes(y=chr_F,x=eig1))+
  # geom_violin(draw_quantiles = 0.5,fill="lightblue")+
  geom_density_ridges(quantile_lines = TRUE,quantiles=2,fill="lightblue")+
  # geom_boxplot(outlier.colour = NA)+
  theme_minimal()+
  labs(y="Chromosome",x="Eigenvalue 1")+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16))

# Make summary tables... --------------------------------------------------
# Eig1
eig1_winds_res <- stickleback_eigenvec1_window_summaries[order(-stickleback_eigenvec1_window_summaries$eigenvalue),]
write.table(eig1_winds_res,
            "tables/TableSX_Stickleback_Eig1_outliers.txt",
            row.names = F,quote = F,sep = "\t")

# Eig2
eig2_winds_res <- stickleback_eigenvec2_window_summaries[order(-stickleback_eigenvec2_window_summaries$eigenvalue_sum),]
write.table(eig2_winds_res,
            "tables/TableSX_Stickleback_Eig2_outliers.txt",
            row.names = F,quote = F,sep = "\t")

# Combine all results together for a single figure ------------------------
pdf("figs/FigureX_MxF_BC_stickleback.pdf",width=14,height=6)
plot_grid(final_pheno_tree,
plot_grid(genome_figs[[1]]+
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  strip.text = element_text(size=12)),
          plot_grid(groupIV_fig+
                      theme(title = element_blank())+
                      labs(x="Chr4 Position (Mb)")+
                      theme(axis.text.x = element_text(size=14,angle=45,hjust=1)),
                    per_chrom_ridges+
                      theme(axis.text.y = element_text(size=10)),
                    ncol=2,rel_widths = c(3,2),labels = c("C","D"),label_size = 26,hjust = 0.1),
          ncol=1,nrow=2,rel_heights=c(3,2.5),labels = c("B",""),label_size = 26),
ncol=2,rel_widths = c(1.5,3),labels = c("A",""),label_size = 26)
dev.off()