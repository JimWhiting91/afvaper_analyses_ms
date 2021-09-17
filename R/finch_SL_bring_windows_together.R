################################################ 
# For finch S-L data, evaluate window sizes and bring data together
################################################ 
lib <- c("GenotypePlot","ghibli","pbmcapply","tidyverse","ggplot2","parallel","regioneR","bedr","pegas","ggtree","cowplot")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")
source("~/Exeter/five_aside/five_aside_functions.R")

# Read in data for all windows
window_snps <- c(50,100,200,500)

# Set comparison
comparison <- "hmga2_finch_final"

# Set chrs
finch_fai <- read.table("data/finch_genome.fai")
large_chrs <- finch_fai[finch_fai$V2 > 5000000,1]

# Read in the rds results if needs be...
window_outliers <- pbmclapply(window_snps,function(wind){
  print(wind)
  chr_res <- lapply(large_chrs,function(chr){
    #  print(chr)
    return(readRDS(paste0("outputs/",comparison,"_chr",chr,"_windowsize_",wind,".rds")))
  })
  
  # Just get the chr AF_inputs
  AF_input_chrs <-  lapply(chr_res,"[[",1)
  names(AF_input_chrs) <- large_chrs
  
  # Pull all the eigen_res together
  eigen_res_list_test <- lapply(chr_res,"[[",3)
  eigen_res_list_test <- merge_eigen_res(eigen_res_list_test)
  
  # Distribution of window sizes...
  window_sizes <- sapply(names(eigen_res_list_test),function(x){
    wind <- strsplit(x,":")[[1]][2]
    as.integer(strsplit(wind,"-")[[1]][2])-as.integer(strsplit(wind,"-")[[1]][1])
  })
  med_wind_size <- median(window_sizes)
  
  # Make a single null cutoff
  null_input_list <- lapply(chr_res,"[[",2)
  all_nulls <- unlist(null_input_list, recursive=FALSE)
  all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                       cutoffs = c(0.95,0.99,0.999,0.9999))
  
  # We want to look at eigenvector 2
  signif_windows <- afvaper::signif_eigen_windows(eigen_res_list_test,cutoffs = all_null_cutoffs[,3])
  
  # Summarise eigenvector 1 windows...
  if(is.na(signif_windows[[1]][1])){
    eigenvec1_window_summaries <- NA
  } else {
    eigenvec1_window_summaries <- afvaper::summarise_window_parallelism(signif_windows[[1]],eigen_res_list_test,loading_cutoff = 0.3,eigenvector = 1)
  }
  
  # Summarise eigenvector 2 windows...
  if(is.na(signif_windows[[2]][1])){
    eigenvec2_window_summaries <- NA
  } else {
    eigenvec2_window_summaries <- afvaper::summarise_window_parallelism(signif_windows[[2]],eigen_res_list_test,loading_cutoff = 0.3,eigenvector = 2)
  }
  
  return(list(med_wind_size,signif_windows[[1]],signif_windows[[2]],eigenvec1_window_summaries,eigenvec2_window_summaries))
},mc.cores = 4)


# Plot Finch tree with comparisons ----------------------------------------
# Base tree simple newick from Han et al 2017
finch_tree <- read.tree("data/galapagos_finch_tree.nwk")

# Fetch the vectors from afvaper res
vectors <- list(# Allopatry (Cocos)
  c("P.inornata_C","G.propinqua_G"),
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.pauper_F","C.psittacula_P"),
  c("C.parvulus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.acutirostris_G","G.magnirostris_G"),
  c("C.parvulus_Z","P.crassirostris_Z"))

# Prune tree for these species
pruned.tree <- drop.tip(finch_tree,finch_tree$tip.label[!(finch_tree$tip.label %in% unique(unlist(vectors)))])

# Inspect
base_tree <- ggtree(pruned.tree,ladderize = F)+ 
  theme_tree()+
  xlim(0, 7)

# Add the vectors...
for(i in 1:length(vectors)){
  base_tree <- base_tree + 
    geom_taxalink(vectors[[i]][1],vectors[[i]][2],
                  arrow=arrow(length=unit(0.025, "npc"),type = "closed",ends = "first",angle=90),
                  colour="red2",alpha=0.5,
                  curvature = -.6) +
    geom_taxalink(vectors[[i]][1],vectors[[i]][2],
                  arrow=arrow(length=unit(0.025, "npc"),type = "closed",ends = "last"),
                  colour="red2",alpha=0.75,
                  curvature = -.6) 
}
base_tree

# Annotate
final_tree <- base_tree +
  geom_tiplab(align=TRUE, linesize=.5,fontface="bold.italic",size=4)

# Set up phenotype table
phenos <- c("Small","Large")
phenotypes <- matrix(ncol=1,nrow=length(pruned.tree$tip.label))
colnames(phenotypes) <- "phenotype"
rownames(phenotypes) <- pruned.tree$tip.label
for(i in 1:length(phenos)){
  phenotypes[rownames(phenotypes) %in% sapply(vectors,'[[',i),"phenotype"] <- phenos[i]
}

# Add to tree
final_pheno_tree <- gheatmap(final_tree, phenotypes,offset=2, width=0.2,colnames=FALSE, legend_title="genotype")+
  scale_fill_manual(breaks = phenos,
                    values = ghibli_palettes$YesterdayMedium[5:6])+
  scale_x_ggtree()+
  theme(legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size=13))+
  labs(fill="Beak Phenotype")

# Analyse afvaper res -----------------------------------------------------

# Print median window sizes...
for(i in 1:length(window_snps)){
  print(window_outliers[[i]][[1]])
}

# For each, we get overlapping regions...
eig1_outliers <- sapply(window_outliers,"[[",2)
eig2_outliers <- sapply(window_outliers,"[[",3)

# Get N outliers
for(i in 1:length(window_outliers)){
  print(length(eig1_outliers[[i]]))
}

# Convert regions to bed
regions2bed <- function(regions){
  if(is.na(regions[1])){
    return(NULL)
  } else {
    split1 <- strsplit(regions,":")
    out <- data.frame(chr=sapply(split1,"[[",1))
    split2 <-  strsplit(sapply(split1,"[[",2),"-")
    out$start <- as.integer(sapply(split2,"[[",1))
    out$end <- as.integer(sapply(split2,"[[",2))
    return(data.frame(out))
  } 
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

# Loop over all sets
overlapping_outliers <- lapply(list(eig1_outliers,eig2_outliers),function(outlier_regions){
  
  # Get Beds
  eig_bed <- lapply(outlier_regions,regions2bed)
  
  # Filter null
  eig_bed[sapply(eig_bed, is.null)] <- NULL
  if(length(eig_bed) > 1){
    comps_to_make <- combn(1:length(eig_bed),2)
    
    # sort lexographically
    eig_bed_sort <- lapply(eig_bed,function(bed){
      bedr.sort.region(bed,check.chr = FALSE)
    })
    
    # If we actually have outliers
    if(length(eig_bed_sort) != 0){
      names(eig_bed_sort) <- window_snps
      
      # Get the intersection
      outlier_intersect <- bedr.join.multiple.region(x = eig_bed_sort,check.chr = F)
      eig_overlaps <- outlier_intersect[order(-as.integer(outlier_intersect$n.overlaps)),]
      eig_overlaps <- eig_overlaps[as.integer(eig_overlaps$n.overlaps) > 1,]
      
      # Return
      return(eig_overlaps)
    } else {
      return(NULL)
    }
  }
})

# Observe final outlier sets
overlapping_outliers[[1]]
overlapping_outliers[[2]]

# Observe summaries for eig1
eig1_summaries <- lapply(window_outliers,"[[",4)
eig1_summaries[[1]]
eig1_summaries[[2]]
eig1_summaries[[3]]
eig1_summaries[[4]]

# Observe summaries for eig2
eig2_summaries <- lapply(window_outliers,"[[",5)
eig2_summaries[[1]]
eig2_summaries[[2]]
eig2_summaries[[3]]
eig2_summaries[[4]]


# Supp Tables of eigenvector results --------------------------------------
# Eig1
eig1_supp <- unique(data.frame(rbindlist(lapply(rownames(overlapping_outliers[[1]][overlapping_outliers[[1]]$n.overlaps == 4,]),function(window_id){
  
  # Isolate the biggest window size
  biggest_window <- tail(strsplit(overlapping_outliers[[1]][window_id,"names"],",")[[1]],1)
  window_indices <- 1:4
  names(window_indices) <- window_snps
  
  # Fetch the relevant windows from the window SNP res
  window_res_tmp <- window_outliers[[window_indices[biggest_window]]][[4]]
  window_res_tmp <- window_res_tmp[window_res_tmp$window_id %in% subset_by_regions(window_res_tmp$window_id,window_id),]
  
  # Add in the window size
  window_res_tmp$window_size_snps <- biggest_window
  return(window_res_tmp)
}))))
eig1_supp <- eig1_supp[order(-eig1_supp$eigenvalue),]
eig1_supp$total_variance <- round((eig1_supp$eigenvalue/6)*100,2)
eig1_supp$eigenvalue <- round(eig1_supp$eigenvalue,3)
eig1_supp <- eig1_supp[,c("window_id","eigenvector","eigenvalue","total_variance","parallel_lineages","parallel_pops","antiparallel_pops","window_size_snps")]

# Write these
write.table(eig1_supp,
            "tables/TableSX_SL_finch_Eig1_outliers.txt",
            row.names = F,quote = F,sep="\t")

# Final finch plots -------------------------------------------------------
# Do genome-wide results for 200 SNP windows...
wind=200

# Fetch all res
chr_res <- pbmclapply(large_chrs,function(chr){
  #  print(chr)
  return(readRDS(paste0("outputs/",comparison,"_chr",chr,"_windowsize_",wind,".rds")))
},mc.cores=4)

# Pull all the eigen_res together
eigen_res_list_test <- merge_eigen_res(lapply(chr_res,"[[",3))

# Make a single null cutoff
null_input_list <- lapply(chr_res,"[[",2)
all_nulls <- unlist(null_input_list, recursive=FALSE)
all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                     cutoffs = c(0.95,0.99,0.999,0.9999))

# Only plot the 20 largest contigs...

# Plot genome-wide
SL_genome_fig <- eigenval_plot(eigen_res_list_test,null_vectors = all_nulls,plot.pvalues = T,keep_chr_scale = F)
SL_genome_fig[[1]]

# And plot over the target region
target_scaf_fig <- eigenval_plot(eigen_res_list_test[grep("JH739900",names(eigen_res_list_test))],
                                 null_vectors = all_nulls,plot.pvalues = T,keep_chr_scale = F)

# Add on annotation for HMGA2 gene
hmga2_region <- c(7003336,7122126)
hmga2_fig <- target_scaf_fig[[1]] +
  annotate("rect", xmin = hmga2_region[1], xmax = hmga2_region[2], ymin = -Inf, ymax = Inf,alpha = .5,fill="red2") +
  annotate("text", x = 6e6, y = 3.5, label = "italic(HMGA2)~gene",parse = TRUE)+
  xlab("Contig JH739900 Position (Mb)")+
  theme(title = element_blank())

# Combine them all together
pdf("figs/FigureX_SL_finch_results.pdf",width=16,height=6)
plot_grid(final_pheno_tree,
          plot_grid(SL_genome_fig[[1]]+
                      theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.minor.y = element_blank())+
                      xlab("Contig"),
                    hmga2_fig,
                    ncol=1,nrow=2,align = "v",axis = "tblr",labels=c("B","C"),label_size=30,rel_heights = c(1.8,1.3)),
          ncol=2,nrow=1,rel_widths = c(3,5),labels = c("A","",""),label_size = 30)
dev.off()

################################################################################  
# What is happening at the hmga2 gene? - genotypeplot
# hmga2 haplotypes, small > large
hmga2_vectors <-  list(# Allopatry (Cocos)
  c("P.inornata_C","G.propinqua_G"),
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.pauper_F","C.psittacula_P"),
  c("C.parvulus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.acutirostris_G","G.magnirostris_G"),
  c("C.parvulus_Z","P.crassirostris_Z"))

hmga2_species <- unique(unlist(hmga2_vectors))
popmap <- data.frame(rbindlist(lapply(hmga2_species,function(species){
  tmp <- read.table(paste0("data/finch_species_",species,".popmap"))
  tmp$species <- species
  return(tmp)
})))

# Input VCF path
vcf_path <- "data/111_sample_maf0.05.vcf.gz"

# Plot the hmga2 gene region

# Genotype Plot and Cluster
hmga2_genotypes <- genotype_plot(vcf=vcf_path,
                                 chr="JH739900",
                                 start=7007107,
                                 end=7120688,
                                 cluster=T,
                                 invariant_filter = T,
                                 popmap=popmap,
                                 snp_label_size = 20000)

# Dendrogram
dendrogram_metadata <- data.frame(ind = hmga2_genotypes$dendro_labels)
dendrogram_metadata$tip <- 1:nrow(dendrogram_metadata)

# Build metadata
dendrogram_metadata$species <- NA
dendrogram_metadata$beak <- "Small"
species <- unique(popmap$species)
for(i in 1:length(species)){
  dendrogram_metadata[dendrogram_metadata$ind %in% popmap[popmap$species == species[i],1],"species"] <- species[i]
}
large_species <- unlist(lapply(hmga2_vectors,'[[',2))
dendrogram_metadata[dendrogram_metadata$species %in% large_species,"beak"] <- "Large"

# Add tip labels
hmga2_genotypes$dendrogram <- hmga2_genotypes$dendrogram + 
  geom_jitter(aes(y=-2.5,x=dendrogram_metadata$tip,colour=dendrogram_metadata$species,shape=dendrogram_metadata$beak),
              size = 2.2,height = 2.2,width=0,alpha=0.75)+
  scale_shape_manual(breaks=c("Small","Large"),
                     values=c(17,19))+
  theme(legend.position="left",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  labs(colour="Species",shape="Beak Shape")+
  guides(colour = guide_legend(override.aes = list(size=8)),
         shape = guide_legend(override.aes = list(size=8)))

# Test expansion
# Combine and save
pdf("figs/FigureSX_hmga2_gene_finches.pdf",width = 10,height=8)
plot_grid(combine_genotype_plot(hmga2_genotypes),ggplot()+theme_void(),rel_widths = c(10,1))
dev.off()


