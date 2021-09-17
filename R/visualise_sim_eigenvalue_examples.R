##################################################
# Visualise eigenvalues using simulated results
##################################################

lib <- c("ggplot2","data.table","pbmcapply","dplyr","cowplot")
lapply(lib,library,character.only=T)

# Make a plotting function
plot_sim_eigenvals <- function(results_dir,eigenvector=1,window_size=10000){
  
  # Fetch all of the afvaper results
  res_files <- list.files(results_dir,pattern = "wind200.EIGEN_RESULTS.rds")
  
  # Loop over these and read them in
  all_res <- data.frame(rbindlist(pbmclapply(res_files,function(x){
    tmp <- readRDS(paste0(results_dir,"/",x))
    tmp_eigenvals <- lapply(tmp[[2]],'[[',1)
    
    # Sum the eigenvals as we would in a normal analysis
    tmp_eigenvals <- lapply(tmp_eigenvals,cumsum)
    
    # Fetch all positions
    all_pos <- gsub("chr1:","",names(tmp_eigenvals))
    pos_mat <- matrix(nrow=length(tmp_eigenvals),ncol=3)
    pos_mat[,1] <- as.integer(sapply(strsplit(all_pos,"-"),'[[',1))
    pos_mat[,2] <- as.integer(sapply(strsplit(all_pos,"-"),'[[',2))
    pos_mat[,3] <- rowMeans(pos_mat[,1:2])
    
    # Loop over all the eigenvectors we are interested in and tabulate these...
    all_eigs <- data.frame(rbindlist(lapply(1:eigenvector,function(eig){
      out_mat <- matrix(ncol=2,nrow=length(tmp_eigenvals))
      out_mat[,1] <- x
      out_mat[,2] <- sapply(tmp_eigenvals, function(wind) return(wind[eig]))
      
      # Merge with pos
      out_mat <- data.frame(cbind(out_mat,pos_mat))
      colnames(out_mat) <- c("filepath","eigenvalue","start","end","mid")
      out_mat$eigenvalue <- as.numeric(out_mat$eigenvalue)
      out_mat$start <- as.numeric(out_mat$start)
      out_mat$end <- as.numeric(out_mat$end)
      out_mat$mid <- as.numeric(out_mat$mid)
      
      # Add eigenvector ID
      out_mat$eigenvector <- eig
      
      # Now we can just simplify these to consistent windows...
      # Windowise
      start <- seq(0,max(out_mat$end),window_size)
      end <-   start+window_size
      
      # Fetch the avg
      window_calcs <- sapply(start,function(pos){
        mean(out_mat[out_mat$end >= pos &
                       out_mat$start <= pos+window_size,"eigenvalue"])
      })
      
      # Set up a new output
      wind_out <- data.frame(start=start,
                             end=end,
                             mid=start + window_size/2,
                             eigenvalue=window_calcs,
                             eigenvector=eig,
                             filepath=x)
      
      return(wind_out)
    })))
    
    # Return this...
    return(all_eigs)
  },mc.cores=6)))
  
  # Now we can average these again within windows across all iterations with sd
  all_iter_avg <- all_res %>% group_by(mid,eigenvector) %>% summarise(eig_mean = mean(eigenvalue),
                                                          eig_sd = sd(eigenvalue))
  
  # And make our plot
  if(eigenvector == 1){
    ggplot(all_iter_avg,aes(mid,eig_mean))+
      geom_line()+
      geom_ribbon(aes(ymin=eig_mean-eig_sd,ymax=eig_mean+eig_sd),alpha=0.5)+
      theme_minimal()+
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=16),
            strip.background = element_rect(colour="black"))+
      scale_x_continuous(breaks = c(0,10e6,20e6),
                         labels = c(0,10,20))+
      labs(x="Sim Chr Position (Mb)",y=expression(Parallelism~(sum(lambda[italic(i)], i=1, m))))
      
      
  } else {
    all_iter_avg$eigenvector <- paste0("Eig ",all_iter_avg$eigenvector)
    ggplot(all_iter_avg,aes(mid,eig_mean))+
      geom_line()+
      geom_ribbon(aes(ymin=eig_mean-eig_sd,ymax=eig_mean+eig_sd),alpha=0.5)+
      theme_minimal()+
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=16),
            strip.background = element_rect(colour="black"),
            strip.text = element_text(size=14))+
      scale_x_continuous(breaks = c(0,10e6,20e6),
                         labels = c(0,10,20))+
      labs(x="Sim Chr Position (Mb)",y=expression(Parallelism~(sum(lambda[italic(i)], i=1, m))))+
    facet_wrap(~eigenvector,ncol=1,strip.position = "right",scales = "free_y")
  }
}

# Where are our results?
full_parallel_dir <- "outputs/slim/21,08,18_Rerunning_full_parallel_with_200_evolving_gens_15/final_results/"
full_parallel_plots <- plot_sim_eigenvals(full_parallel_dir,eigenvector = 1,window_size = 10000)

two_parallel_dir <- "outputs/slim/21,08,20,Rerunning_two_parallel_with_200_evolving_gens_15/final_results/"
two_parallel_plots <- plot_sim_eigenvals(two_parallel_dir,eigenvector = 2,window_size = 10000)

divergent_dir <- "outputs/slim/21,09,07,Running_a_full_divergent_example_15/final_results/"
divergent_plots <- plot_sim_eigenvals(divergent_dir,eigenvector = 3,window_size = 10000)

# Combine
pdf("figs/FigureX_eigenvalue_peaks_simulated.pdf",width=12,height=4)
plot_grid(full_parallel_plots+ggtitle("Full Parallel")+theme(title=element_text(size=20))+geom_hline(yintercept = 4,linetype="dotted"),
          two_parallel_plots+ggtitle("Multi-Parallel")+theme(title=element_text(size=20))+geom_hline(yintercept = 4,linetype="dotted"),
          divergent_plots+ggtitle("Divergent")+theme(title=element_text(size=20))+geom_hline(yintercept = 4,linetype="dotted"),
          nrow=1,
          labels="AUTO",label_size=32)
dev.off()
