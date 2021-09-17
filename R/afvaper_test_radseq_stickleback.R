################################################ 
# Genome-wide Test of AF-vapeR analysis using RAD-seq data from stickleback
################################################ 

devtools::load_all("~/Exeter/afvaper/")
library(tidyverse)
library(ggplot2)
library(cowplot)
library(rphylopic)

# Test inputs
stickleback_rad_test <- "~/Nottingham/NatEvol_Ecol_Genomic_Parallelism/data/MxF_VCFs/Bc.vcf.gz"
wind_sizes <- c(10,25,50)
CORES <- 6
permutations <- 10000

# Make popmap
vcf_inds <- system(paste0("bcftools query -l ",stickleback_rad_test),intern = T)
popmap <- data.frame(ind=vcf_inds,
                          pop=gsub('[0-9]+', '', vcf_inds))
popmap <- popmap[!(popmap$pop %in% c("OBSM","NYPS","MUD")),]

# Make the vectors
fresh_pops <- unique(popmap[,2])[!(unique(popmap[,2]) %in% c("OBSM","NYPS","LICA","MUD"))]
test_vectors <- lapply(fresh_pops,function(fresh){return(c("LICA",fresh))})
names(test_vectors) <- fresh_pops

# Just read in the VCF once...
stickleback_vcf <- read.vcfR(stickleback_rad_test)

# Run over the following chr - omitting the sex chromosome groupXIX. Also include scaffold_324 which includes pitx1 (it maps to the end of groupVII)
chrs <- c(paste0("group",as.roman(c(1:18,20:21))),"scaffold_324")
chr_sizes <- sapply(chrs,function(chr){
  return(max(as.integer(stickleback_vcf@fix[stickleback_vcf@fix[,1] == chr ,2])))
})
chr_props <- chr_sizes/sum(chr_sizes)
chr_perms <- ceiling(chr_props * permutations)
names(chr_perms) <- chrs

# Loop over all chromosomes and write an output for each...
for(chr in chr_res){
  
  message(paste0("STARTING CHR ",which(chrs == chr)," OF ",length(chrs)))
  
  # Subset the VCF
  vcf_in <- stickleback_vcf[stickleback_vcf@fix[,1] == chr,]
  
  # Loop over window sizes
  for(window_snps in wind_sizes){
    
    # Calculate allele frequency matrices
    AF_input <- calc_AF_vectors(vcf = vcf_in,
                                window_size = window_snps,
                                popmap = popmap,
                                vectors = test_vectors,
                                n_cores = CORES)
    
    # # Calculate null frequencies
    null_input <- calc_AF_vectors(vcf = vcf_in,
                                  window_size = window_snps,
                                  popmap = popmap,
                                  vectors = test_vectors,
                                  n_cores = CORES,
                                  null_perms = chr_perms[chr])
    
    # Calculate Eigen statistics
    eigen_res <- lapply(AF_input,eigen_analyse_vectors)
    
    # Save both
    saveRDS(list(AF_input,null_input,eigen_res),
            paste0("outputs/stickleback_rad_test_",chr,"_AF_eigen_res_windsize_",window_snps,".rds"))
  }
}
