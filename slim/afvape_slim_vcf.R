##########################################################
# Analyse slim output VCFs...
lib <- c("ggplot2","vcfR","data.table","parallel")
sapply(lib,library,character.only=T)
devtools::load_all("~/software/afvaper")

# Read our VCF from the command lines
args <- commandArgs(TRUE)
vcf <- as.character(args[1])
vcf_in <- read.vcfR(vcf)

# Now eigen variables
eigen_window =as.integer(args[2])
slim_vectors = list(c("pop1","pop2"),
                    c("pop1","pop3"),
                    c("pop1","pop4"),
                    c("pop1","pop5"))
names(slim_vectors) <- paste0("Daughter_",1:4)
CORES = as.integer(args[3])
permutations = as.integer(args[4])


# Make popmap
popmap <- data.frame(ind=colnames(vcf_in@gt)[2:ncol(vcf_in@gt)],
                     pop=rep(paste0("pop",1:5),each=ncol(vcf_in@gt)/5))

# Make our AF_matrices
slim_AF <- calc_AF_vectors(vcf = vcf_in,
                           window_size = eigen_window,
                           popmap = popmap,
                           vectors = slim_vectors,
                           n_cores = CORES)

# Calculate null frequencies
null_input <- calc_AF_vectors(vcf = vcf_in,
                              window_size = eigen_window,
                              popmap = popmap,
                              vectors = slim_vectors,
                              null_perms = permutations)

# Calculate Eigen statistics 
eigen_res <- lapply(slim_AF,eigen_analyse_vectors)

# Export the results
saveRDS(list(slim_AF,eigen_res,null_input),
        paste0(vcf,"_wind",eigen_window,".EIGEN_RESULTS.rds"))
