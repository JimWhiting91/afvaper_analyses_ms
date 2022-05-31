################################################ 
# Genome-wide Test of AF-vapeR analysis for Experimental Drosophila...
################################################ 

lib <- c("tidyverse","ggplot2","parallel","data.table")
lapply(lib,library,character.only=T)
devtools::load_all("~/Exeter/afvaper/")

# Test inputs
drosophila_test <- "data/kelly_hughes_supp_data.xlsx"
CORES <- 6
permutations <- 10000
wind_sizes <- c(25,50,100,200)

# Fetch SNPs
snps_input <- data.table(na.omit(readxl::read_xlsx(drosophila_test,skip = 1)))
colnames(snps_input) <- c("chr","pos","ref","alt",
                          "A0_reads","A0_pR",
                          "B0_reads","B0_pR",
                          "C0_reads","C0_pR",
                          "A7_reads","A7_pR",
                          "B7_reads","B7_pR",
                          "C7_reads","C7_pR",
                          "LRT-parallel")

# Fetch pop names and all that
pops <- grep("reads",colnames(snps_input),value=T)
pops <- gsub("_reads","",pops)

# Use Pops to calculate allele frequency input...
AF_input <- snps_input[,.(chr,pos,A0_pR,B0_pR,C0_pR,A7_pR,B7_pR,C7_pR)]
AF_input$chr <- gsub("Scf_","",AF_input$chr)

# Make the vectors
test_vectors <- list(c("A0_pR","A7_pR"),
                     c("B0_pR","B7_pR"),
                     c("C0_pR","C7_pR"))
names(test_vectors) <- c("A","B","C")

# Get per chrom perms for a total of N permutations
chr_sizes <- AF_input[,.(size=max(.SD$pos)),by=chr]
chr_props <- chr_sizes$size/sum(chr_sizes$size)
chr_perms <- ceiling(chr_props * permutations)

# Run over the following chr
chrs <- chr_sizes$chr
names(chr_perms) <- chrs 

#chr_res <- lapply(chrs,function(chr){
for(chr in chrs){
  
  message(paste0("STARTING CHR ",which(chrs == chr)," of ",length(chrs)))
  
  
  # Loop over window sizes
  for(window_snps in wind_sizes){
    
    # Calculate allele frequency matrices
    AFV <- calc_AF_vectors(vcf = data.frame(AF_input)[AF_input$chr == chr,],
                           data_type = "freq",
                           window_size = window_snps,
                           popmap = NULL,
                           vectors = test_vectors,
                           n_cores = CORES)
    
    # Calculate null frequencies
    null_AFV <- calc_AF_vectors(vcf = data.frame(AF_input)[AF_input$chr == chr,],
                                data_type = "freq",
                                window_size = window_snps,
                                popmap = NULL,
                                vectors = test_vectors,
                                n_cores = CORES,
                                null_perms = chr_perms[chr])
    
    # Calculate Eigen statistics
    eigen_res <- lapply(AFV,eigen_analyse_vectors)
    
    # Get null matrix
    null_cutoffs <- find_null_cutoff(null_res = null_AFV,
                                     cutoffs = c(0.95,0.99,0.999))
    
    # Save both
    saveRDS(list(AFV,null_AFV,eigen_res,null_cutoffs),
            paste0("outputs/drosophila_test_",chr,"_AF_eigen_res_windsize_",window_snps,".rds"))
  }
}
