devtools::install_github("thomasp85/patchwork")
install.packages("ggsave")
library(ggplot2);library(patchwork)

c0 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/allopatry_neutral.txt")
c1 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/flux50_neutral.txt")
c2 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/parapatry_neutral.txt")
df0 <- cbind.data.frame(c0, rep("isolation", 100))
colnames(df0) <- c("generations","model")
df1 <- cbind.data.frame(c1, rep("periodic", 100))
colnames(df1) <- c("generations","model")
df2 <- cbind.data.frame(c2, rep("continuous", 100))
colnames(df2) <- c("generations","model")
df <- rbind.data.frame(df0,df1,df2)

gav <- 2/1e-4
gavp <- 0.025/((1e-4)^2)
gav50 <- 1/((0.5*(1e-4/2))+(0.5*(((1e-4)^2)/0.025)))
a <- ggplot(df, aes(x=model, y=generations)) + 
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(aes(yintercept=gav),linetype="dashed",lwd=0.75)+
  geom_hline(aes(yintercept=gavp),linetype="solid",lwd=0.75)+
  geom_hline(aes(yintercept=gav50),linetype="dotted",lwd=0.75)+
  #geom_jitter(shape=16, position=position_jitter(0.2), color="#8A008A") +
  geom_boxplot(alpha=0.8) +
  scale_y_continuous(trans='log10') +
  ylab("Time to speciation") +
  xlab("Model") +
  labs(subtitle = expression(alpha~"= 0.50"))

c0 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/allopatry_neutral.txt")
c1 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/flux10_neutral.txt")
c2 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/parapatry_neutral.txt")
df0 <- cbind.data.frame(c0, rep("isolation", 50))
colnames(df0) <- c("generations","model")
df1 <- cbind.data.frame(c1, rep("periodic", 50))
colnames(df1) <- c("generations","model")
df2 <- cbind.data.frame(c2, rep("continuous", 50))
colnames(df2) <- c("generations","model")
df <- rbind.data.frame(df0,df1,df2)

gav <- 2/1e-4
gavp <- 0.025/((1e-4)^2)
gav10 <- 1/((0.1*(1e-4/2))+(0.9*(((1e-4)^2)/0.025)))
b <- ggplot(df, aes(x=model, y=generations)) + 
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(aes(yintercept=gav),linetype="dashed",lwd=0.75)+
  geom_hline(aes(yintercept=gavp),linetype="solid",lwd=0.75)+
  geom_hline(aes(yintercept=gav10),linetype="dotted",lwd=0.75)+
  #geom_jitter(shape=16, position=position_jitter(0.2), color="#8A008A") +
  geom_boxplot(alpha=0.8) +
  scale_y_continuous(trans='log10') +
  ylab("Time to speciation") +
  xlab("Model") +
  labs(subtitle = expression(alpha~"= 0.1"))

c0 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/allopatry_neutral.txt")
c1 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/flux01_neutral.txt")
c2 <- read.table("~/Dropbox/speciation_cyclical_migration/simple_model_output/parapatry_neutral.txt")
df0 <- cbind.data.frame(c0, rep("isolation", 100))
colnames(df0) <- c("generations","model")
df1 <- cbind.data.frame(c1, rep("periodic", 100))
colnames(df1) <- c("generations","model")
df2 <- cbind.data.frame(c2, rep("continuous", 100))
colnames(df2) <- c("generations","model")
df <- rbind.data.frame(df0,df1,df2)

gav <- 2/1e-4
gavp <- 0.025/((1e-4)^2)
gav01 <- 1/((0.01*(1e-4/2))+(0.99*(((1e-4)^2)/0.025)))
c <- ggplot(df, aes(x=model, y=generations)) + 
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(aes(yintercept=gav),linetype="dashed",lwd=0.75)+
  geom_hline(aes(yintercept=gavp),linetype="solid",lwd=0.75)+
  geom_hline(aes(yintercept=gav01),linetype="dotted",lwd=0.75)+
  #geom_jitter(shape=16, position=position_jitter(0.2), color="#8A008A") +
  geom_boxplot(alpha=0.8) +
  scale_y_continuous(trans='log10') +
  ylab("Time to speciation") +
  xlab("Model") +
  labs(subtitle = expression(alpha~"= 0.01"))

fig <- a + b + c
ggsave("~/Dropbox/speciation_cyclical_migration/figures/simulation_waiting_time.pdf", plot=fig, width = 10, height = 4, units = "in")

  

