#vostok/migration rate plot
library(ggplot2);library(magrittr);library(data.table);library(cowplot);library(plyr)
setwd("~/Dropbox/speciation_periodic_migration/")

vostok <- fread("data/Vostok.txt") #V1=depth;V2=ice_age;V3=gas_age;V4=C02ppm
mig_epoch <- c();mig <- T;epoch <- 1
for(i in 1:length(vostok$V4)){
  newmig <- vostok$V4[i]>250
  if(newmig != mig) epoch <- epoch+1
  mig <- newmig
  mig_epoch[i] <- epoch
}
vostok$mig_epoch <- mig_epoch

vostokfig <- ggplot(data=vostok,aes(x=V3,y=V4))+
  theme_classic()+theme(axis.text=element_text(size=7),
                        axis.title=element_text(size=7))+
  geom_line(lwd=0.5)+
  geom_path(data=subset(vostok,V4<250),aes(group=mig_epoch),col="red")+
  annotate(geom="text",label=expression(italic(m)==0.025),col="red",x=4.8e4,y=240,size=2.5)+
  annotate(geom="text",label=expression(italic(m)==0),col="black",x=5e4,y=275,size=2.5)+
  #geom_hline(yintercept = 250)+
  xlab("Years Before Present")+
  ylab(expression(CO[2]~ppm))

pops <- data.frame(x=c(1,2),y=c(1,1))
popfig <- ggplot(data=pops,aes(x=x,y=y))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        )+
  geom_point(shape=1,size=8)+
  annotate(geom="segment",x=1.3,xend=1.7,y=1.1,yend=1.1,col="red",
           arrow = arrow(type = "open", angle = 30, length = unit(0.15, "cm")))+
  annotate(geom="segment",x=1.7,xend=1.3,y=0.9,yend=0.9,col="red",
           arrow = arrow(type ="open", angle = 30, length = unit(0.15, "cm")))+
  xlim(0.8,2.2)+ylim(0.5,1.5)
 
pdf("~/Dropbox/speciation_periodic_migration/figures/vostok_w_sim.pdf",width=3.5,height=2.5)      
ggdraw()+
  draw_plot(vostokfig,0,0,1,1)+
  draw_plot(popfig,0.135,0.82,0.3,0.2)
dev.off()
 
#weighted average migration rate assuming migration when CO2 < 250
tmp <- ddply(vostok,.(mig_epoch),summarize,nyears=max(V3)-min(V3),mig=V4[1]<250)
migyears <- sum(tmp$nyears[tmp$mig==T])
nyears <- max(vostok$V3)
1/((migyears/nyears)/0.025+((nyears-migyears)/nyears)/1e-9)
#summary stats for constant v periodic migration
ss_p <- fread("data/sumstats/periodic.txt")
ss_p$migration <- "periodic"
ss_c <- fread("data/sumstats/constant.txt")
ss_c$migration <- "continuous"
ss <- rbind(ss_p,ss_c)

# mann whitney U tests
wilcox.test(ss_p$segsites, ss_c$segsites) # p-value < 2.2e-16
wilcox.test(ss_p$pi, ss_c$pi) # p-value < 2.2e-16
wilcox.test(ss_p$thetaW, ss_c$thetaW) # p-value < 2.2e-16
wilcox.test(ss_p$tajD, ss_c$tajD) # p-value < 2.2e-16
wilcox.test(ss_p$het_o, ss_c$het_o) # p-value < 2.2e-16
wilcox.test(ss_p$dxy, ss_c$dxy)  # p-value < 2.2e-16
wilcox.test(ss_p$ibs_bw_mean, ss_c$ibs_bw_mean) # p-value < 2.2e-16
wilcox.test(ss_p$ibs_wi_mean, ss_c$ibs_wi_mean) # p-value < 2.2e-16
wilcox.test(ss_p$ibs_wi_skew, ss_c$ibs_wi_skew) # p-value < 2.2e-16
wilcox.test(ss_p$ibs_bw_skew, ss_c$ibs_bw_skew) # p-value = 0.04978
mean(ss_p$ibs_bw_skew) #221.3612
mean(ss_c$ibs_bw_skew) #413.3099
wilcox.test(ss_p$ibs_bw_blocks_over_1e5, ss_c$ibs_bw_blocks_over_1e5)
wilcox.test(ss_p$ibs_wi_blocks_over_1e5, ss_c$ibs_wi_blocks_over_1e5)


pd <- melt(ss,id.vars="migration")
pd <- subset(pd,!variable %in% c("ibs_wi_skew","ibs_bw_skew"))
p <- ggplot(data=pd,aes(x=value,fill=migration))+
  theme_classic()+
  theme(strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=6),
        strip.text=element_text(size=7))+
  facet_wrap(~variable,scales="free")+
  scale_fill_manual(values=c("dodgerblue","orange"))+
  geom_density(col=NA,alpha=0.7)

pdf("~/Dropbox/speciation_periodic_migration/figures/sumstats.pdf",width=6.5,height=3)
print(p)
dev.off()

# black and white version
pd <- melt(ss,id.vars="migration")
pd <- subset(pd,!variable %in% c("fst"))

# parsing special characters
pd$variable <- factor(pd$variable, 
                      levels = c("segsites",
                                 "pi",
                                 "thetaW",
                                 "tajD",
                                 "het_o",
                                 "fst",
                                 "dxy",
                                 "ibs_bw_mean",
                                 "ibs_bw_skew",
                                 "ibs_bw_blocks_over_1e5",
                                 "ibs_wi_mean",
                                 "ibs_wi_skew",
                                 "ibs_wi_blocks_over_1e5"),
                       ordered = TRUE, 
                      labels = c("SNPs",expression(pi), paste(expression(theta),'[W]'),'D[Tajima]','H[O]',
                                                   'F[ST]','D[XY]','Mean~IBS[bewteen]','IBS[between]~skew','Mean~IBS[between]~long~blocks','Mean~IBS[within]','IBS[within]~skew',
                                                   'Mean~IBS[within]~long~blocks'))


p <- ggplot(data=pd,aes(x=value,fill=migration))+
  theme_bw()+
  theme(panel.grid = element_blank()) +
  theme(strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text=element_text(size=4),
        strip.text=element_text(size=7))+
  facet_wrap(~variable,labeller = label_parsed,scales="free")+
  scale_fill_manual(values=c("gray90","gray30"))+
  geom_density(color="black",alpha=0.7)

pdf("~/Dropbox/speciation_periodic_migration/figures/Fig_5_sumstats.pdf",width=8,height=4)
print(p)
dev.off()

# look at vcf from sims
sim.med <- read.vcfR("data/vcf/vostok_sim_10388971.trees.vcf")
dna <- vcfR2DNAbin(sim.med, unphased_as_NA = F, consensus = T, extract.haps = F)
sim <- DNAbin2genind(dna)
samples <- as.character(read.table("samples.txt")[[1]]) 
rownames(sim@tab) <- samples
sim.scaled <- scaleGen(sim,NA.method="mean",scale=F)
pca <- prcomp(sim.scaled,center=F,scale=F)
pc <- data.frame(pca$x[,1:3])
ggplot(data=pc,aes(x=PC1,y=PC2))+geom_text(aes(label=samples)) #pop assignments were correct


# plot moments output
pm <- read.table("~/Dropbox/speciation_periodic_migration/data/PM_params.txt")
im <- read.table("~/Dropbox/speciation_periodic_migration/data/IM_params.txt")
pm.ll <- pm$V14 # grab log likelihoods
im.ll <- im$V7
pm.ll <- cbind.data.frame(pm.ll, rep("periodic",nrow(pm)))
colnames(pm.ll) <- c("log_likelihood", "model")
im.ll <- cbind.data.frame(im.ll, rep("continuous",nrow(im)))
colnames(im.ll) <- c("log_likelihood", "model")
#im.ll <- im.ll[sample(nrow(im.ll), 21), ] #subset to same length while models run
ll.df <- rbind.data.frame(pm.ll,im.ll)

q <- ggplot(ll.df, aes(x=model, y=log_likelihood)) + 
  theme_bw()+
  theme(panel.grid = element_blank()) +
  theme(axis.title=element_text(size=12)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  #geom_jitter(shape=16, position=position_jitter(0.2), color="#8A008A") +
  geom_boxplot(alpha=0.8) +
  ylab("Log-likelihood") +
  xlab("Model")

pdf("~/Dropbox/speciation_periodic_migration/figures/Fig_7_likelihoods.pdf",width=3,height=3)
print(q)
dev.off()

# test difference
pm.ll <- pm$V14 # regrab log likelihoods
im.ll <- im$V7
wilcox.test(im.ll, pm.ll) # W = 259, p-value = 0.9478

