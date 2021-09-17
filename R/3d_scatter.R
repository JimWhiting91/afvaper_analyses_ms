# Script for scatter
pdf("~/Exeter/angle_grinder/figs/empty_3D.pdf",width=5,height=5)
scatter3D(`SNP1-AF`, `SNP2-AF`, `SNP3-AF`,alpha=0,phi=0,bty = "g",
          xlab = "SNP-1 AF", ylab = "SNP-2 AF", zlab = "SNP-3 AF",
          colkey = NULL)
dev.off()

# Make some dummy eigenvalue examples...
library(ggplot2)
full_parallel <- data.frame(x=1:4,
                            y=c(3.5,3.8,3.9,4.0),
                            y2=c(3.5,0.3,0.1,0.01))

divergent_parallel <- data.frame(x=1:4,
                                 y=c(2.0,3.8,3.9,4.0),
                                 y2=c(2.0,1.8,0.1,0.01))
full_divergent <- data.frame(x=1:4,
                             y=c(1.5,2.5,3.4,4.0),
                             y2=c(1.5,1.0,0.9,0.6))
fp <- ggplot(full_parallel,aes(x,y))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(sum(lambda[italic(i)], i=1, m)),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Full Parallel")

fp2 <- ggplot(full_parallel,aes(x,y2))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(lambda[italic(i)]),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Full Parallel")

dp <- ggplot(divergent_parallel,aes(x,y))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(sum(lambda[italic(i)], i=1, m)),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Divergent Parallel")

dp2 <- ggplot(divergent_parallel,aes(x,y2))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(lambda[italic(i)]),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Divergent Parallel")

fd <- ggplot(full_divergent,aes(x,y))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(sum(lambda[italic(i)], i=1, m)),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Full Divergent")

fd2 <- ggplot(full_divergent,aes(x,y2))+
  geom_line()+
  geom_point(colour="red2",shape=18,size=7)+
  labs(y=expression(lambda[italic(i)]),x=expression(Eigenvector[italic(i)]))+
  theme_minimal()+
  theme(axis.title.y=element_text(size=34),
        axis.title.x=element_text(size=28),
        axis.text = element_text(size=26),
        title = element_text(size=20))+
  ylim(0,4)
 # ggtitle("Full Divergent")

# Plot together
library(cowplot)
pdf("figs/demo_eigenvalue_plots.pdf",width=8.5,height=10)
plot_grid(plot_grid(fp2,dp2,fd2,
                    ncol=1,align = "v"),
          plot_grid(fp,dp,fd,
                    ncol=1,align = "v"),
          ncol=2,align="h")
dev.off()

