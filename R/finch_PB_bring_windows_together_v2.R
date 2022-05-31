################################################ 
# For finch P-B data, evaluate window sizes and bring data together + analyse
################################################ 

lib <- c("ggtree","pbmcapply","cowplot","RColorBrewer","tidyverse","ggplot2","parallel","regioneR","bedr","pegas","GenotypePlot")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")
source("/Volumes/jimwhiting_external/Exeter/five_aside/five_aside_functions.R")

# Read in data for all windows
window_snps <- c(50,100,200,500)

# Set comparison
comparison <- "alx1_finch_final"

# Set chrs
finch_fai <- read.table("data/finch_genome.fai")
large_chrs <- finch_fai[finch_fai$V2 > 5000000,1]

# Read in the rds results if needs be...
window_outliers <- pbmclapply(window_snps,function(wind){
  print(wind)
  chr_res <- lapply(large_chrs,function(chr){
    #  print(chr)
    return(readRDS(paste0("outputs/",comparison,"_v2_chr",chr,"_windowsize_",wind,".rds")))
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
                                       cutoffs = c(0.95,0.99,0.999,0.9999,1))
  
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
  
  # And fetch all of the pvalues
  eig_pval_res <- eigen_pvals(eigen_res_list_test,all_nulls)
  
  return(list(med_wind_size,signif_windows[[1]],signif_windows[[2]],eigenvec1_window_summaries,eigenvec2_window_summaries,eig_pval_res))
},mc.cores = 4)

# Plot Finch tree with comparisons ----------------------------------------
# Base tree is a simple newick tree capturing the topology as shown in Han et al. 2017
finch_tree <- read.tree("data/galapagos_finch_tree.nwk")

# Fetch the vectors from afvaper res
vectors <- list(# Allopatry (Cocos)
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.parvulus_Z","G.magnirostris_M"),
  c("C.pallidus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.propinqua_G","G.magnirostris_G"),
  c("G.acutirostris_G","G.magnirostris_G"))
names(vectors) <- sapply(vectors,function(x) paste0(x[1],"-",x[2]))

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
phenos <- c("Pointed","Blunt")
phenotypes <- matrix(ncol=1,nrow=length(pruned.tree$tip.label))
colnames(phenotypes) <- "phenotype"
rownames(phenotypes) <- pruned.tree$tip.label
for(i in 1:length(phenos)){
  phenotypes[rownames(phenotypes) %in% sapply(vectors,'[[',i),"phenotype"] <- phenos[i]
}

# Add to tree
final_pheno_tree <- gheatmap(final_tree, phenotypes,offset=2, width=0.2,colnames=FALSE, legend_title="genotype")+
  scale_fill_manual(breaks = phenos,
                    values = brewer.pal(n = 8, name = "Dark2")[2:3])+
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

# Get outlier count for each window size...
for(i in 1:length(window_outliers)){
  print(length(eig1_outliers[[i]]))
}

# Convert regions to bed
regions2bed <- function(regions){
  if(length(regions) > 0){
    out <- matrix(nrow=length(regions),ncol=3)
    split1 <- strsplit(regions,":")
    out[,1] <- sapply(split1,"[[",1)
    split2 <-  strsplit(sapply(split1,"[[",2),"-")
    out[,2] <- sapply(split2,"[[",1)
    out[,3] <- sapply(split2,"[[",2)
    colnames(out) <- c("chr","start","end")
    out <- data.frame(out)
    class(out$start) <- "integer"
    class(out$end) <- "integer"
    return(out)
    
  } else {
    return(NULL)
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
  
  # Clean
  winds_to_keep <- window_snps[!(sapply(outlier_regions, anyNA))]
  outlier_regions <- outlier_regions[!sapply(outlier_regions, anyNA)]
  
  # Get Beds
  eig_bed <- lapply(outlier_regions,regions2bed)
  names(eig_bed) <- winds_to_keep
  
  # Filter null
  eig_bed[sapply(eig_bed, is.null)] <- NULL
  comps_to_make <- combn(1:length(eig_bed),2)
  
  # sort lexographically
  eig_bed_sort <- lapply(eig_bed,function(bed){
    bedr.sort.region(bed,check.chr = FALSE)
  })
  
  # If we actually have outliers
  if(length(eig_bed_sort) != 0){
    names(eig_bed_sort) <- names(eig_bed)
    
    # Get the intersection
    outlier_intersect <- bedr.join.multiple.region(x = eig_bed_sort,check.chr = F)
    
    # Keep overlapping
    eig_overlaps <- outlier_intersect[order(-as.integer(outlier_intersect$n.overlaps)),]
    eig_overlaps <- eig_overlaps[as.integer(eig_overlaps$n.overlaps) > 1,]
    
    # Return
    return(eig_overlaps)
  } else {
    return(NULL)
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
# These are based on the largest window when looking at overlapping windows...
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

# Fetch all of the pvals
for(i in 1:nrow(eig1_supp)){
  window_tmp <- eig1_supp$window_id[i]
  snp_tmp <- eig1_supp$window_size_snps[i]
  res_tmp <- which(window_snps == as.integer(snp_tmp))
  eig1_supp$empPvalue[i] <- window_outliers[[res_tmp]][[6]][window_tmp,1]
}

eig1_supp <- eig1_supp[,c("window_id","eigenvector","eigenvalue","total_variance","empPvalue","parallel_lineages","parallel_pops","antiparallel_pops","window_size_snps")]

# Write these
write.table(eig1_supp,
            "tables/TableSX_PB_finch_Eig1_outliers.txt",
            row.names = F,quote = F,sep="\t")

# # Eig2
# eig2_supp <- unique(data.frame(rbindlist(lapply(rownames(overlapping_outliers[[2]][overlapping_outliers[[2]]$n.overlaps == 3,]),function(window_id){
#   
#   # Isolate the biggest window size
#   biggest_window <- tail(strsplit(overlapping_outliers[[2]][window_id,"names"],",")[[1]],1)
#   window_indices <- 1:4
#   names(window_indices) <- window_snps
#   
#   # Fetch the relevant windows from the window SNP res
#   window_res_tmp <- window_outliers[[window_indices[biggest_window]]][[5]]
#   window_res_tmp <- window_res_tmp[window_res_tmp$window_id %in% subset_by_regions(window_res_tmp$window_id,window_id),]
#   
#   # Add in the window size
#   window_res_tmp$window_size_snps <- biggest_window
#   return(window_res_tmp)
# }))))
# eig2_supp <- eig2_supp[order(-eig2_supp$eigenvalue),]
# eig2_supp$total_variance <- round((eig2_supp$eigenvalue/6)*100,2)
# eig2_supp$eigenvalue <- round(eig2_supp$eigenvalue,3)
# 
# # Fetch all of the pvals
# for(i in 1:nrow(eig2_supp)){
#   window_tmp <- eig2_supp$window_id[i]
#   snp_tmp <- eig2_supp$window_size_snps[i]
#   res_tmp <- which(window_snps == as.integer(snp_tmp))
#   eig2_supp$empPvalue[i] <- window_outliers[[res_tmp]][[6]][window_tmp,1]
# }
# 
# eig2_supp <- eig2_supp[,c("window_id","eigenvector","eigenvalue","total_variance","empPvalue","parallel_lineages","parallel_pops","antiparallel_pops","window_size_snps")]
# 
# # Write these
# write.table(eig2_supp,
#             "tables/TableSX_PB_finch_eig2_outliers.txt",
#             row.names = F,quote = F,sep="\t")

################################################################################  
##### Plot the ALX1 locus for 200 SNP windows... #####
wind200_res <- pbmclapply(large_chrs,function(chr){
  readRDS(paste0("outputs/",comparison,"_chr",chr,"_windowsize_200.rds"))
},mc.cores=4)

# Pull all the eigen_res together
wind200_eigen_res <- lapply(wind200_res,"[[",3)
allchr_eigen_res <- merge_eigen_res(wind200_eigen_res)

# Fetch nulls
wind200_null_input_list <- lapply(wind200_res,"[[",2)
wind200_all_nulls <- unlist(wind200_null_input_list, recursive=FALSE)
wind200_all_null_cutoffs <- find_null_cutoff(null_res = wind200_all_nulls,
                                             cutoffs = c(0.95,0.99,0.999,0.9999))

# Plot full genome
full_genome <- eigenval_plot(allchr_eigen_res,
                             #cutoffs = wind200_all_null_cutoffs[,3],
                             null_vectors = wind200_all_nulls,
                             plot.pvalues=T,keep_chr_scale = F)[[1]]
full_genome

# Plot together
alx_fig <- eigenval_plot(allchr_eigen_res[grep("JH739921",names(allchr_eigen_res))],
                         cutoffs = wind200_all_null_cutoffs[,3],
                         null_vectors = wind200_all_nulls,
                         plot.pvalues=T)[[1]]
alx1_region <- c(411351,428776)
alx_fig <- alx_fig +
  annotate("rect", xmin = alx1_region[1], xmax = alx1_region[2], ymin = -Inf, ymax = Inf,alpha = .5,fill="red2") +
  annotate("text", x = 1e6, y = 3.5, label = "italic(ALX1)~gene",parse = TRUE)+
  xlab("Contig JH739921 Position (Mb)")+
  theme(title = element_blank())

# Full summary figure for main manuscript...
pdf("figs/FigureX_PB_finch_results.pdf",width=16,height=6)
plot_grid(final_pheno_tree,
          plot_grid(full_genome+
                      theme(panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            panel.grid.minor.y = element_blank())+
                      xlab("Contig"),
                    alx_fig,
                    ncol=1,nrow=2,align = "v",axis = "tblr",labels=c("E","F"),label_size=30,rel_heights = c(1.8,1.3)),
          ncol=2,nrow=1,rel_widths = c(3,5),labels = c("D","",""),label_size = 30)
dev.off()

################################################################################  
# What is happening at the ALX1 gene? - genotypeplot
# ALX1 haplotypes, pointed > blunt
alx1_vectors <- list(# Allopatry (Cocos)
  c("P.inornata_C","G.magnirostris_M"),
  # Allopatry
  c("C.parvulus_Z","G.magnirostris_M"),
  c("C.pallidus_Z","G.magnirostris_M"),
  # Sympatry
  c("G.propinqua_G","G.magnirostris_G"),
  c("G.acutirostris_G","G.magnirostris_G"))

alx1_species <- unique(unlist(alx1_vectors))
popmap <- data.frame(rbindlist(lapply(alx1_species,function(species){
  tmp <- read.table(paste0("data/finch_species_",species,".popmap"))
  tmp$species <- species
  return(tmp)
})))

# Input VCF path
vcf_path <- "data/111_sample_maf0.05.vcf.gz"

# Plot the alx1 gene region

# Genotype Plot and Cluster
alx1_genotypes <- genotype_plot(vcf=vcf_path,
                                chr="JH739921",
                                start=411351,
                                end=428776,
                                cluster=T,
                                invariant_filter = T,
                                popmap=popmap)

# Dendrogram
dendrogram_metadata <- data.frame(ind = alx1_genotypes$dendro_labels)
dendrogram_metadata$tip <- 1:nrow(dendrogram_metadata)

# Build metadata
dendrogram_metadata$species <- NA
dendrogram_metadata$beak <- "Pointed"
species <- unique(popmap$species)
for(i in 1:length(species)){
  dendrogram_metadata[dendrogram_metadata$ind %in% popmap[popmap$species == species[i],1],"species"] <- species[i]
}
dendrogram_metadata[grep("magnirostris",dendrogram_metadata$species),"beak"] <- "Blunt"

# Add tip labels
alx1_genotypes$dendrogram <- alx1_genotypes$dendrogram + 
  geom_jitter(aes(y=-2.5,x=dendrogram_metadata$tip,colour=dendrogram_metadata$species,shape=dendrogram_metadata$beak),
              size = 2.2,height = 2.2,width=0,alpha=0.75)+
  scale_shape_manual(breaks=c("Pointed","Blunt"),
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
pdf("figs/FigureSX_alx1_gene_finches.pdf",width = 8,height=6)
combine_genotype_plot(alx1_genotypes)
dev.off()

# Plot a comparable interesting candidate with full parallelism and interesting candidate genes...
# JH739896:10888810-10913977
RSPO3_genotypes <- genotype_plot(vcf=vcf_path,
                                chr="JH739896",
                                start=10888810,
                                end=10913977,
                                cluster=F,
                                plot_allele_frequency = T,
                                invariant_filter = T,
                                popmap=popmap)

# Dendrogram
dendrogram_metadata <- data.frame(ind = RSPO3_genotypes$dendro_labels)
dendrogram_metadata$tip <- 1:nrow(dendrogram_metadata)

# Build metadata
dendrogram_metadata$species <- NA
dendrogram_metadata$beak <- "Pointed"
species <- unique(popmap$species)
for(i in 1:length(species)){
  dendrogram_metadata[dendrogram_metadata$ind %in% popmap[popmap$species == species[i],1],"species"] <- species[i]
}
dendrogram_metadata[grep("magnirostris",dendrogram_metadata$species),"beak"] <- "Blunt"

# Add tip labels
RSPO3_genotypes$dendrogram <- RSPO3_genotypes$dendrogram + 
  geom_jitter(aes(y=-2.5,x=dendrogram_metadata$tip,colour=dendrogram_metadata$species,shape=dendrogram_metadata$beak),
              size = 2.2,height = 2.2,width=0,alpha=0.75)+
  scale_shape_manual(breaks=c("Pointed","Blunt"),
                     values=c(17,19))+
  theme(legend.position="left",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  labs(colour="Species",shape="Beak Shape")+
  guides(colour = guide_legend(override.aes = list(size=8)),
         shape = guide_legend(override.aes = list(size=8)))

# Merge
combine_genotype_plot(RSPO3_genotypes)

################################################################################  
# Find the genes in all outlier regions
# Here we map the geoFor genome onto the zebrafinch genome and use the mapped regions to call genes from BioMart
#### Get genes from ensembl ####
library(biomaRt)

# Set up Ensembl
ensembl <- useMart("ensembl")
zebrafinch <- useDataset("tguttata_gene_ensembl",mart=ensembl)

# Fetch outliers from ensembl
outlier_res <- lapply(eig1_supp$window_id,function(window_id){
  aln_block_length=10000
  aln_qual=60
  
  # Expand by 50kb up/downstream
  pos <- strsplit(window_id,":")[[1]][2]
  start <- as.integer(strsplit(pos,"-")[[1]][1]) - 50000
  end <- as.integer(strsplit(pos,"-")[[1]][2]) + 50000
  if(start < 0){
    start <- 0
  }
  window_id2 <- paste0(strsplit(window_id,":")[[1]][1],":",start,"-",end)
  
  # Firstly, take all the outlier regions and make a multi-fasta
  system(paste0("samtools faidx ~/Exeter/Genomes/geoFor1.fa.gz ",window_id2," > outputs/tmp_multi.fa"))
  
  # Align to old genome
  system2('/Users/jimwhiting/bin/minimap2',
          args=c(paste0("-t ",6),"~/Exeter/Genomes/Taeniopygia_guttata.bTaeGut1_v1.p.dna_sm.toplevel.fa","outputs/tmp_multi.fa"),
          stdout ="outputs/tmp_multi.paf",wait=T)
  
  # Fetch regions
  aln<-read.table("outputs/tmp_multi.paf",fill=T)
  
  # Keep tidy!
  system(paste0("rm -f outputs/tmp_multi*"))
  
  # Only keep "high support alignment regions"
  aln<-aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]
  regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))
  
  # Pull uniprot genes from biomaRt for each region
  tmp_biomart<-getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","external_gene_name"),
                     filters= "chromosomal_region",
                     values=regions,
                     mart=zebrafinch)
  if(nrow(tmp_biomart)==0){
    tmp_biomart2 <- data.frame(matrix(NA,ncol=ncol(tmp_biomart),nrow=1))
    colnames(tmp_biomart2) <- colnames(tmp_biomart)
  } else {
    tmp_biomart2 <- tmp_biomart
  }
  
  cbind(outlier_id = rep(window_id,nrow(tmp_biomart2)),
        zebrafinch_region=paste(regions,collapse = ", "),
        tmp_biomart2)
})
names(outlier_res) <- eig1_supp$window_id
alx1_outlier_candidates <- data.frame(rbindlist(outlier_res))

# Format...

# Write to a table...
write.table(alx1_outlier_candidates,
            "tables/TableSX_PB_finch_candidate_genes.tsv",
            row.names = F,quote = F,sep = "\t")

