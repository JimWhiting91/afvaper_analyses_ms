######################################################
# For Guppy HP-LP data, evaluate window sizes
######################################################

rm(list=ls())
lib <- c("ape","ggtree","GenotypePlot","pbmcapply","tidyverse","ggplot2","parallel","regioneR","bedr","cowplot")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in data for all windows
wind_sizes <- c(25,50,100,200)

# Get per chrom perms for a total of N permutations
chr_sizes <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")
chr_sizes <- chr_sizes[chr_sizes$V2 > 1000000,]

# Run over the following chr
chrs <- chr_sizes$V1
chrs <- chrs[chrs != "000094F_0"]

# Plot Finch tree with comparisons ----------------------------------------
guppy_tree <- read.tree("data/five_aside_STAR_individuals_tree_outgroup_rooted_swap_TAC.nwk")

# Prune tree for these species
pruned.tree <- drop.tip(guppy_tree,guppy_tree$tip.label[guppy_tree$tip.label == "Outgroup"])

# Inspect
base_tree <- ggtree(pruned.tree,ladderize = F)+ 
  theme_tree()+
  # geom_text(aes(label=node))+
  # geom_tiplab()+
  xlim(0, 7)

# Flip
base_tree <- flip(base_tree,node1 = 9,node2 = 10)
base_tree <- flip(base_tree,node1 = 12,node2 = 17)
base_tree <- flip(base_tree,node1 = 14,node2 = 13)
base_tree <- flip(base_tree,node1 = 15,node2 = 16)
base_tree <- flip(base_tree,node1 = 1,node2 = 2)
base_tree <- flip(base_tree,node1 = 3,node2 = 4)

# Mkae vectors
vectors <- list(c("TACHP","TACLP"),
                c("GHP","GLP"),
                c("APHP","APLP"),
                c("OHP","OLP"),
                c("MADHP","MADLP"))

# Add the vectors...
for(i in 1:length(vectors)){
  base_tree <- base_tree + 
    geom_taxalink(vectors[[i]][1],vectors[[i]][2],
                  arrow=arrow(length=unit(0.025, "npc"),type = "closed",ends = "last"),
                  colour="#073a6c",alpha=1,
                  curvature = .6) 
}
base_tree

# Annotate
final_tree <- base_tree +
  geom_tiplab(align=TRUE, linesize=.5,fontface="bold.italic",size=4)

# Set up phenotype table
phenos <- c("HP","LP")
phenotypes <- matrix(ncol=1,nrow=length(pruned.tree$tip.label))
colnames(phenotypes) <- "phenotype"
rownames(phenotypes) <- pruned.tree$tip.label
for(i in 1:length(phenos)){
  phenotypes[rownames(phenotypes) %in% sapply(vectors,'[[',i),"phenotype"] <- phenos[i]
}

# Add to tree
final_pheno_tree <- gheatmap(final_tree, phenotypes,offset=2, width=0.2,colnames=FALSE, legend_title="genotype")+
  scale_fill_manual(breaks = phenos,
                    values = c("#073a6c","#fc8009"))+
  scale_x_ggtree()+
  theme(legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size=13))+
  labs(fill="Predation")


# Analyse AF-vaper --------------------------------------------------------
# Read in the rds results if needs be...
window_outliers <- pbmclapply(wind_sizes,function(window_snps){
  chr_res <- lapply(chrs,function(chr){
    #  print(chr)
    return(readRDS(paste0("outputs/five_aside_test_",chr,"_AF_eigen_res_windsize_",window_snps,".rds")))
  })
  
  # Just get the chr AF_inputs
  AF_input_chrs <-  lapply(chr_res,'[[',1)
  names(AF_input_chrs) <- chrs
  
  # Pull all the eigen_res together
  eigen_res_list_test <- lapply(chr_res,'[[',3)
  allchr_eigen_res <- merge_eigen_res(eigen_res_list_test)
  
  # Distribution of window sizes...
  window_sizes <- sapply(names(allchr_eigen_res),function(x){
    wind <- strsplit(x,":")[[1]][2]
    as.integer(strsplit(wind,"-")[[1]][2])-as.integer(strsplit(wind,"-")[[1]][1])
  })
  med_wind_size <- median(window_sizes)
  
  # Collate the null cutoffs
  null_cutoff_list <- lapply(chr_res,function(x){return(x[[4]][,1])})
  names(null_cutoff_list) <- chrs
  
  # Make a single null cutoff
  null_input_list <- lapply(chr_res,'[[',2)
  all_nulls <- unlist(null_input_list, recursive=FALSE)
  all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                       cutoffs = c(0.95,0.99,0.999,0.9999))
  
  # Fetch all our significant windows...
  if(window_snps %in% c(25,50)){
    signif <- signif_eigen_windows(allchr_eigen_res,cutoffs = all_null_cutoffs[,3])
  } else if(window_snps %in% c(100,200)) {
    signif <- signif_eigen_windows(allchr_eigen_res,cutoffs = all_null_cutoffs[,2])
  }
  
  # We want to look at eigenvector 2
  eig1_signif <- signif[[1]]
  eig2_signif <- signif[[2]]
  
  # Summarise eigenvector 1 windows...
  five_aside_eigenvec1_window_summaries <- summarise_window_parallelism(eig1_signif,allchr_eigen_res,loading_cutoff = 0.3,eigenvector = 1)
  
  # Summarise eigenvector 2 windows...
  five_aside_eigenvec2_window_summaries <- summarise_window_parallelism(eig2_signif,allchr_eigen_res,loading_cutoff = 0.3,eigenvector = 2)
  
  return(list(med_wind_size,eig1_signif,eig2_signif,five_aside_eigenvec1_window_summaries,five_aside_eigenvec2_window_summaries))
},mc.cores = 4)

# Print median window sizes...
for(i in 1:length(window_outliers)){
  print(window_outliers[[i]][[1]])
}

# For each, we get overlapping regions...
eig1_outliers <- sapply(window_outliers,"[[",2)
eig2_outliers <- sapply(window_outliers,"[[",3)

# Convert regions to bed
regions2bed <- function(regions){
  split1 <- strsplit(regions,":")
  out <- data.frame(chr=sapply(split1,"[[",1))
  split2 <-  strsplit(sapply(split1,"[[",2),"-")
  out$start <- as.integer(sapply(split2,"[[",1))
  out$end <- as.integer(sapply(split2,"[[",2))
  return(data.frame(out))
}

# Return indices where regions1 is in regions2
subset_by_regions <- function(regions1,regions2,values=T){
  
  # Get beds
  region1_bed <- data.frame(rbindlist(lapply(regions1,regions2bed)))
  region2_bed <- data.frame(rbindlist(lapply(regions2,regions2bed)))
  
  # Find overlaps
  overlaps <- data.frame(findOverlaps(GRanges(region1_bed),GRanges(region2_bed)))
  
  # Return the integers in first set with overlap..
  if(!(values)){
    sort(unique(overlaps$queryHits))
  } else {
    regions1[sort(unique(overlaps$queryHits))]
  }
}

# Loop over both sets
overlapping_outliers <- lapply(list(eig1_outliers,eig2_outliers),function(outlier_regions){
  
  # Get Beds
  eig_bed <- lapply(outlier_regions,regions2bed)
  
  # Filter null
  eig_bed[sapply(eig_bed, is.null)] <- NULL
  comps_to_make <- combn(1:length(eig_bed),2)
  
  # sort lexographically
  eig_bed_sort <- lapply(eig_bed,function(bed){
    #print(bed)
    bedr.sort.region(bed,check.chr = FALSE)
  })
  names(eig_bed_sort) <- wind_sizes
  
  # Get the intersection
  outlier_intersect <- bedr.join.multiple.region(x = eig_bed_sort,check.chr = F)
  
  eig_overlaps <- outlier_intersect[order(-as.integer(outlier_intersect$n.overlaps)),]
  eig_overlaps <- eig_overlaps[as.integer(eig_overlaps$n.overlaps) > 1,]
  # Return
  return(eig_overlaps)
})

# Observe final outlier sets
overlapping_outliers[[1]]
overlapping_outliers[[2]]

# Observe summaries for eig1
window_outliers[[1]][[4]]
window_outliers[[2]][[4]]
window_outliers[[3]][[4]]
window_outliers[[4]][[4]]

# Observe summaries for eig2
window_outliers[[1]][[5]]
window_outliers[[2]][[5]]
window_outliers[[3]][[5]]
window_outliers[[4]][[5]]


# Full plot of all window sizes -------------------------------------------
# Make separate plots for eig1 and eig2 outliers, over all window sizes and stacked...

# Loop over window sizes
wind_out <- lapply(wind_sizes,function(wind){
  
  # Fetch the results
  chr_res <- pbmclapply(chrs,function(chr){
    #  print(chr)
    return(readRDS(paste0("outputs/five_aside_test_",chr,"_AF_eigen_res_windsize_",wind,".rds")))
  },mc.cores=6)
  
  # Pull all the eigen_res together
  eigen_res_list_test <- lapply(chr_res,'[[',3)
  allchr_eigen_res <- merge_eigen_res(eigen_res_list_test)
  
  # Make a single null cutoff
  null_input_list <- lapply(chr_res,'[[',2)
  all_nulls <- unlist(null_input_list, recursive=FALSE)
  
  # Plot eigenvalue 1
  tmp_plots <- afvaper::eigenval_plot(allchr_eigen_res[grep("chr",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T)
  
  # Loop over eig1 and eig2
  eig_plots <- lapply(1:2,function(eig){
    
    # Extract focal plot
    tmp_plot <- tmp_plots[[eig]]
    
    # Extract focal outliers
    tmp_outliers <- overlapping_outliers[[eig]]
    
    # Build metadata for highlighting overlapping outliers...
    # Retain only outliers relevant for this window size
    tmp_outliers <- tmp_outliers[grep(wind,tmp_outliers$names),]
    colnames(tmp_outliers)[1:3] <- c("chr","start","end")
    
    # Plot the midpoints
    tmp_outliers$mid <- rowMeans(tmp_outliers[,c("start","end")])
    
    # Add factors
    tmp_outliers$chr_F <- factor(tmp_outliers$chr,levels=c(paste0("chr",1:23)))
    
    # Add plotting levels
    tmp_outliers$plot_y <- 4.33
    tmp_outliers[tmp_outliers$n.overlaps == 3,"plot_y"] <- 4.66
    tmp_outliers[tmp_outliers$n.overlaps == 4,"plot_y"] <- 5
    
    # Add these to the focal plot and also format to make pretty
    p1 <- tmp_plot +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=12),
            legend.position = "none") +
      geom_point(data=tmp_outliers[order(tmp_outliers$n.overlaps),],aes(x=mid,y=plot_y,colour=n.overlaps,shape=n.overlaps),size=3)+
      ggtitle(paste0(wind," SNP - Eig:",eig))+
      scale_color_brewer(palette = "Dark2")+
      scale_shape_manual(breaks = c(2,3,4),
                         values = c(16,17,15))+
      scale_y_continuous(breaks=1:4,labels = 1:4)
    
    # Remove the x axis until wind size is maxed
    if(wind != max(wind_sizes)){
      p1 + theme(strip.text = element_blank(),
                 title = element_blank(),
                 axis.title.x = element_blank()) 
    } else {
      p1 + theme(title = element_blank()) 
    }
  })
  
  return(eig_plots)
})

# Combine all of these together...
all_wind_plots <- unlist(wind_out,recursive = F)

# Make a dummy plot and steal the legend
legend_dd <- data.frame(x=c(1,1,1),
                        y=c(2,2,2),
                        overlap=as.character(c(2,3,4)))
legend <- get_legend(
  ggplot(legend_dd,aes(x,y,colour=overlap,shape=overlap))+
    geom_point(size=5)+
    labs(colour="Window\nSize\nOverlap",shape="Window\nSize\nOverlap")+
    scale_color_brewer(palette = "Dark2")+
    scale_shape_manual(breaks = c(2,3,4),
                       values = c(16,17,15))+
    theme_void()+
    theme(legend.title = element_text(size=16),
          legend.text = element_text(size=14))
)

# Plot them all
pdf("figs/FigureX_HPLP_guppy_results.pdf",width=18,height=6.5)
plot_grid(final_pheno_tree,
          plot_grid(plotlist = all_wind_plots,
                    ncol=2,rel_heights = c(1,1,1,1.6),hjust = 0.8,
                    labels=c("B","C"),label_size = 30),
          legend,
          ncol=3,rel_widths = c(1,4,0.5),labels = c("A",""),label_size = 30)
dev.off()

# Supp Table summarising eigenvector results ------------------------------
# Eig1
eig1_supp <- unique(data.frame(rbindlist(lapply(rownames(overlapping_outliers[[1]]),function(window_id){
  
  # Isolate the biggest window size
  biggest_window <- tail(strsplit(overlapping_outliers[[1]][window_id,"names"],",")[[1]],1)
  window_indices <- 1:4
  names(window_indices) <- wind_sizes
  
  # Fetch the relevant windows from the window SNP res
  window_res_tmp <- window_outliers[[window_indices[biggest_window]]][[4]]
  window_res_tmp <- window_res_tmp[window_res_tmp$window_id %in% subset_by_regions(window_res_tmp$window_id,window_id),]
  
  # Add in the window size
  window_res_tmp$window_size_snps <- biggest_window
  return(window_res_tmp)
}))))
eig1_supp <- eig1_supp[order(-eig1_supp$eigenvalue),]
eig1_supp$total_variance <- round((eig1_supp$eigenvalue/5)*100,2)
eig1_supp$eigenvalue <- round(eig1_supp$eigenvalue,3)
eig1_supp <- eig1_supp[,c("window_id","eigenvector","eigenvalue","total_variance","parallel_lineages","parallel_pops","antiparallel_pops","window_size_snps")]

# Write these
write.table(eig1_supp,
            "tables/TableSX_Guppy_Eig1_outliers.txt",
            row.names = F,quote = F,sep="\t")

# Eig2
eig2_supp <- unique(data.frame(rbindlist(pbmclapply(rownames(overlapping_outliers[[2]]),function(window_id){
  
  # Isolate the biggest window size
  biggest_window <- tail(strsplit(overlapping_outliers[[2]][window_id,"names"],",")[[1]],1)
  window_indices <- 1:4
  names(window_indices) <- wind_sizes
  
  # Fetch the relevant windows from the window SNP res
  window_res_tmp <- window_outliers[[window_indices[biggest_window]]][[5]]
  window_res_tmp <- window_res_tmp[window_res_tmp$window_id %in% subset_by_regions(window_res_tmp$window_id,window_id),]
  
  # Add in the window size
  window_res_tmp$window_size_snps <- biggest_window
  return(window_res_tmp)
},mc.cores=4))))
eig2_supp <- eig2_supp[order(-eig2_supp$eigenvalue_sum),]
eig2_supp$total_variance <- round((eig2_supp$eigenvalue/5)*100,2)
eig2_supp$eigenvalue <- round(eig2_supp$eigenvalue,3)
eig2_supp <- eig2_supp[,c("window_id","eigenvector","eigenvalue","total_variance","eigenvalue_sum","parallel_lineages","parallel_pops","antiparallel_pops","window_size_snps")]

# Write these
write.table(eig2_supp,
            "tables/TableSX_Guppy_Eig2_outliers.txt",
            row.names = F,quote = F,sep="\t")

############################################################
# Use 200 SNP windows to look at some specifics
chr_res <- mclapply(chrs,function(chr){
  #  print(chr)
  return(readRDS(paste0("outputs/five_aside_test_",chr,"_AF_eigen_res_windsize_200.rds")))
},mc.cores=6)

# Just get the chr AF_inputs
AF_input_chrs <-  lapply(chr_res,'[[',1)
#names(AF_input_chrs) <- chrs
allchr_input_chrs <- unlist(AF_input_chrs,recursive = F)

# Pull all the eigen_res together
eigen_res_list_test <- lapply(chr_res,'[[',3)
allchr_eigen_res <- merge_eigen_res(eigen_res_list_test)

# Make a single null cutoff
null_input_list <- lapply(chr_res,'[[',2)
all_nulls <- unlist(null_input_list, recursive=FALSE)
all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                     cutoffs = c(0.95,0.99,0.999,0.9999))

# Get signif windows for this window size
tmp_signif <- signif_eigen_windows(allchr_eigen_res,all_null_cutoffs[,1])[1:4]

##################################################
# In depth look at chr20
# The original scaf94 co-ords for this region are (1797025 - (1248797-836423)) and (1797025 - (1296686-836423)), so 000094F_0:1336762-1384651
chr20_plots <- afvaper::eigenval_plot(allchr_eigen_res[grep("chr20",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T)
plot_grid(chr20_plots[[1]],chr20_plots[[2]],ncol=1,align = "v",axis = "tblr")

# Genotype Plot of this region
pops <- c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
pops2 <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")

popmap_in <- data.frame(rbindlist(lapply(pops2,function(x){
  read.table(paste0("~/Exeter/VCFs/",x,".popmap"))
})))
popmap_in$pop<-NA
for(i in 1:length(pops2)){
  popmap_in[grep(pops2[i],popmap_in[,1]),"pop"] <- pops[i]
}

# Plot the region AFs
#chr20:1248797-1296686
#chr15:5028361-5066375
chr20_region2 <- genotype_plot(vcf="~/Exeter/VCFs/five_aside_STAR_chr20_scaf94_merged_shapeit_beagle.vcf.gz",
                               chr="chr20",
                               start=1248797,
                               end=1296686,
                               popmap = popmap_in,
                               plot_allele_frequency = T,
                               cluster=F,
                               snp_label_size = 1000000)
merged_chr20_region2 <- combine_genotype_plot(chr20_region2,heights = c(1,4))

# Also extract and plot the A_matrix...
chr20_A_mat <- data.frame(allchr_eigen_res[["chr20:1248797-1296686"]]$A_matrix)
chr20_A_mat$pos <- as.integer(sapply(strsplit(rownames(chr20_A_mat),"_"),'[[',2))

# Make a longform
A_mat_plot <- data.frame(score=c(chr20_A_mat[,1],chr20_A_mat[,2]),
                         eigenvector=rep(paste0("Eig ",1:2),each=nrow(chr20_A_mat)),
                         pos=rep(chr20_A_mat$pos,2))

# In each case, highlight the top 5 loading SNPs
A_mat_plot$top5 <- "None"
A_mat_plot[A_mat_plot$eigenvector=="Eig 1" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"])))[5],"top5"] <- "Eigenvector 1"
A_mat_plot[A_mat_plot$eigenvector=="Eig 2" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"])))[5],"top5"] <- "Eigenvector 2"

chr20_A_mat_plot <- ggplot(A_mat_plot,aes(x=pos,y=score,colour=top5))+
  geom_point()+
  facet_wrap(~eigenvector,ncol=1,strip.position = "right")+
  scale_color_brewer(palette="Dark2")+
  labs(x="Chr20 Position",y="Eigenvector Score",colour="Top 5 Score")+
  theme_bw()+
  xlim(c(1248797,1296686))+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        strip.text = element_text(size=14),
        legend.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))


# Combine with genotypeplot
pdf("figs/guppy_results_chr20_focus.pdf",width=10,height=8)
cowplot::plot_grid(chr20_A_mat_plot,chr20_region2$positions,chr20_region2$genotypes,
                   ncol=1,axis = "tblr",align="v",
                   rel_heights = c(2,0.75,4))
dev.off()

##### In-depth look at chromosome 15 #####
# Focal window: chr15:5028361-5066375
chr15_plots <- afvaper::eigenval_plot(allchr_eigen_res[grep("chr15",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T)
cowplot::plot_grid(chr15_plots[[1]],chr15_plots[[2]],ncol=1,align = "v",axis = "tblr")

# Genotype Plot of this region
# Plot the region AFs
chr15_region_geno <- genotype_plot(vcf="~/Exeter/VCFs/five_aside_STAR_3033083_allchr_shapeit_beagle.vcf.gz",
                                   chr="chr15",
                                   start=5028361,
                                   end=5066375,
                                   popmap = popmap_in,
                                   plot_allele_frequency = F,
                                   cluster=T,
                                   snp_label_size = 1000000,
                                   missingness = 1)

# Add HP/LP to the dendrogram
dendrogram_metadata <- data.frame(ind = chr15_region_geno$dendro_labels)
dendrogram_metadata$tip <- 1:nrow(dendrogram_metadata)
rivers <- c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas")
rivers2 <- rep(rivers,each=2)
pops <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
pred <- rep(c("HP","LP"),5)

# Build metadata
dendrogram_metadata$river <- NA
dendrogram_metadata$pred <- NA
for(i in 1:length(pops)){
  dendrogram_metadata[grep(pops[i],dendrogram_metadata$ind),"river"]<-rivers2[i]
  dendrogram_metadata[grep(pops[i],dendrogram_metadata$ind),"pred"]<-pred[i]
}

# Add tip labels
chr15_region_geno$dendrogram <- chr15_region_geno$dendrogram + 
  geom_jitter(aes(y=-2.5,x=dendrogram_metadata$tip,colour=dendrogram_metadata$river,shape=dendrogram_metadata$pred),
              size = 2.2,height = 2.2,width=0,alpha=0.75)+
  scale_colour_manual(breaks=five_aside_colour_rivers$river,
                      values=five_aside_colour_rivers$colour)+
  scale_shape_manual(breaks=c("HP","LP"),
                     values=c(19,17))+
  theme(legend.position="left",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  labs(colour="River",shape="Predation")+
  guides(colour = guide_legend(override.aes = list(size=8)),
         shape = guide_legend(override.aes = list(size=8)))


merged_chr15_region <- combine_genotype_plot(chr15_region_geno,widths=c(1,1))

# Extract the PCA and visualise the centroids...
# Fetch centroids
scores <- chr15_region_geno$cluster_pca$li[,1:2]
scores$ind <- rownames(scores)
for(i in 1:nrow(scores)){
  scores$pop[i] <- popmap_in[popmap_in$V1 == scores$ind[i],"pop"]
}

# Make a plot of PC scores given dendrogram order...
scores_ordered <- scores[chr15_region_geno$dendro_labels,]
#scores_ordered$ind <- nrow(scores_ordered$ind)
scores_ordered <- reshape2::melt(scores_ordered)
scores_ordered$ind_N <- rep(1:length(chr15_region_geno$dendro_labels),2)
scores_ordered$variable <- gsub("Axis","PC",scores_ordered$variable)
ordered_scores_fig <- ggplot(scores_ordered)+
  geom_segment(aes(x=0,xend=value,y=ind_N,yend=ind_N))+
  facet_wrap(~variable,ncol=2,strip.position = "bottom")+
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=16),
        strip.text = element_text(size=18))+
  geom_vline(xintercept = 0,linetype="solid")+
  scale_x_continuous(breaks=c(-5,0,5),position = "top")+
  labs(x="PC Score")+
  scale_y_continuous(limits=c(1,nrow(scores_ordered)/2))

# Merge with the above
empty_fig <- ggplot() + theme_void()
merged_chr15_region2 <- plot_grid(chr15_region_geno$dendrogram, 
                                  chr15_region_geno$genotypes,
                                  ordered_scores_fig,
                                  rel_widths = c(2,1,1),axis = "tb", align = "vh", ncol = 3, nrow = 1)

# Can now plot population centroids of PC space
score_centroids <- data.frame(scores %>% group_by(pop) %>% summarise(PC1=mean(Axis1),
                                                                     PC2=mean(Axis2)))

# Add river and pred
pops <- c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
river <- rep(c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas"),each=2)    
pred <- rep(c("HP","LP"),5)
for(i in 1:nrow(score_centroids)){
  score_centroids$river[i] <- river[which(pops == score_centroids$pop[i])]
  score_centroids$pred[i] <- pred[which(pops == score_centroids$pop[i])]
}

# Now plot...
geno_trajectories <- ggplot(score_centroids,aes(PC1,PC2,shape=pred,colour=river))+
  geom_point(size=5)+
  geom_path(aes(group=river),arrow = arrow(),show.legend = F)+
  scale_colour_manual(breaks=five_aside_colour_rivers$river_F,
                      values=five_aside_colour_rivers$colour)+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text = element_text(size=14))+
  geom_vline(xintercept = 0,linetype="dashed")+
  geom_hline(yintercept = 0,linetype="dashed")+
  labs(x="Genotype PC1",y="Genotype PC2")


# Get loading mat for regions
test_mat <- data.frame(allchr_eigen_res[["chr15:5028361-5066375"]]$eigenvecs)[,1:2]
test_mat$pop <- rownames(test_mat)
test_mat$pop <- gsub("Tac","Tacarigua",test_mat$pop)

loading_plot <- ggplot(test_mat,aes(Eigenvector_1,Eigenvector_2,colour=pop))+
  geom_point(size=5,shape=15)+
  theme_bw()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        legend.position = "none")+
  geom_vline(xintercept = 0,linetype="dashed")+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_colour_manual(values=five_aside_colour_rivers$colour,
                      breaks=five_aside_colour_rivers$river_F)+
  labs(x="Eigenvector 1 Loading",y='Eigenvector 2 Loading')


# chr15 plot together -----------------------------------------------------
pdf("figs/FigureSX_HPLP_guppies_chr15_example.pdf",width=12,height=8)
plot_grid(merged_chr15_region2,
          plot_grid(geno_trajectories,
                    loading_plot,
                    ncol=1,nrow=2,labels = c("B","C"),label_size = 30,align = "v",
                    hjust=0.05),
          axis = "tblr",align='h',ncol=2,labels = c("A",""),label_size = 30,rel_widths = c(1.7,1))
dev.off()


# Also extract and plot the A_matrix...
chr15_A_mat <- data.frame(allchr_eigen_res[["chr15:5028361-5066375"]]$A_matrix)
chr15_A_mat$pos <- as.integer(sapply(strsplit(rownames(chr15_A_mat),"_"),'[[',2))

# Make a longform
A_mat_plot <- data.frame(score=c(chr15_A_mat[,1],chr15_A_mat[,2]),
                         eigenvector=rep(paste0("Eig ",1:2),each=nrow(chr15_A_mat)),
                         pos=rep(chr15_A_mat$pos,2))

# In each case, highlight the top 5 loading SNPs
A_mat_plot$top5 <- "None"
A_mat_plot[A_mat_plot$eigenvector=="Eig 1" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"])))[5],"top5"] <- "Eigenvector 1"
A_mat_plot[A_mat_plot$eigenvector=="Eig 2" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"])))[5],"top5"] <- "Eigenvector 2"

chr15_A_mat_plot <- ggplot(A_mat_plot[A_mat_plot$top5 != "None",],aes(x=pos,y=abs(score),colour=top5))+
  geom_point()+
  facet_wrap(~eigenvector,ncol=1,strip.position = "right",scales = "free_y")+
  scale_color_brewer(palette="Dark2")+
  labs(x="chr15 Position",y="Score",colour="Top 5 Score")+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        strip.text = element_text(size=14),
        legend.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  scale_x_continuous(limits = c(5028361,5066375))

# View top SNPs
top_snps <- A_mat_plot[A_mat_plot$top5 != "None",]
top_snps[order(top_snps$pos),]

# In-depth look at chr1 region --------------------------------------------
# chr1:9791276-9798829
# 50 SNP windows
chr_res <- mclapply(chrs,function(chr){
  #  print(chr)
  return(readRDS(paste0("outputs/five_aside_test_",chr,"_AF_eigen_res_windsize_50.rds")))
},mc.cores=6)

# Just get the chr AF_inputs
AF_input_chrs <-  lapply(chr_res,'[[',1)
#names(AF_input_chrs) <- chrs
allchr_input_chrs <- unlist(AF_input_chrs,recursive = F)

# Pull all the eigen_res together
eigen_res_list_test <- lapply(chr_res,'[[',3)
allchr_eigen_res <- merge_eigen_res(eigen_res_list_test)

# Make a single null cutoff
null_input_list <- lapply(chr_res,'[[',2)
all_nulls <- unlist(null_input_list, recursive=FALSE)
all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                     cutoffs = c(0.95,0.99,0.999,0.9999))

# Plot chr1
chr1_plots <- afvaper::eigenval_plot(allchr_eigen_res[grep("chr1:",names(allchr_eigen_res))],null_vectors = all_nulls,plot.pvalues = T)
plot_grid(chr1_plots[[1]],chr1_plots[[2]],ncol=1,align = "v",axis = "tblr")

# Genotype Plot of this region
pops <- c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
pops2 <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")

popmap_in <- data.frame(rbindlist(lapply(pops2,function(x){
  read.table(paste0("~/Exeter/VCFs/",x,".popmap"))
})))
popmap_in$pop<-NA
for(i in 1:length(pops2)){
  popmap_in[grep(pops2[i],popmap_in[,1]),"pop"] <- pops[i]
}

# Plot the region AFs
chr1_region2 <- genotype_plot(vcf="~/Exeter/VCFs/five_aside_STAR_3033083_allchr_shapeit_beagle.vcf.gz",
                              chr="chr1",
                              start=9791276,
                              end=9798829,
                              popmap = popmap_in,
                              plot_allele_frequency = T,
                              cluster=F,
                              snp_label_size = 1000000)
merged_chr1_region2 <- combine_genotype_plot(chr1_region2,heights = c(1,4))
merged_chr1_region2

# Also extract and plot the A_matrix...
chr1_A_mat <- data.frame(allchr_eigen_res[["chr1:9791276-9798829"]]$A_matrix)
chr1_A_mat$pos <- as.integer(sapply(strsplit(rownames(chr1_A_mat),"_"),'[[',2))

# Make a longform
A_mat_plot <- data.frame(score=c(chr1_A_mat[,1],chr1_A_mat[,2]),
                         eigenvector=rep(paste0("Eig ",1:2),each=nrow(chr1_A_mat)),
                         pos=rep(chr1_A_mat$pos,2))

# In each case, highlight the top 5 loading SNPs
A_mat_plot$top5 <- "None"
A_mat_plot[A_mat_plot$eigenvector=="Eig 1" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 1","score"])))[5],"top5"] <- "Eigenvector 1"
A_mat_plot[A_mat_plot$eigenvector=="Eig 2" &
             abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"]) >= rev(sort(abs(A_mat_plot[A_mat_plot$eigenvector=="Eig 2","score"])))[5],"top5"] <- "Eigenvector 2"

chr1_A_mat_plot <- ggplot(A_mat_plot,aes(x=pos,y=score,colour=top5))+
  geom_point()+
  facet_wrap(~eigenvector,ncol=1,strip.position = "right")+
  scale_color_brewer(palette="Dark2")+
  labs(x="chr1 Position",y="Eigenvector Score",colour="Top 5 Score")+
  theme_bw()+
  xlim(c(1248797,1296686))+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        strip.text = element_text(size=14),
        legend.position = "top",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))

cowplot::plot_grid(chr1_A_mat_plot,chr1_region2$positions,chr1_region2$genotypes,
                   ncol=1,axis = "tblr",align="v",
                   rel_heights = c(2,0.75,4))
