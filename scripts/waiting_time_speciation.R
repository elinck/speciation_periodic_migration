### models of expected waiting time to speciation

library(ggplot2);library(patchwork);library(ggsave)

# assign parameters for plots
# m_x = migration rate for plot x
# a_x = proportion of time in allopatry
# p_x = proportion of time exchanging genes

m_1 <- 0.25
a_1 <- 0.5
p_1 <- 1-a_1

fun.1 <- function(x) 2/x
fun.2 <- function(x) m_1/x^2
fun.3 <- function(x) 1/((a_1*(x/2))+(p_1*(x^2/m_1)))

p1 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  theme_bw() +
  stat_function(fun = fun.1, aes(linetype="isolation")) +
  stat_function(fun = fun.2, aes(linetype="continuous")) +
  stat_function(fun = fun.3, aes(linetype="periodic")) +
  scale_linetype_discrete(guide=FALSE) +  
  scale_x_continuous(trans='log10',limits=c(1e-10,1e-5)) +
  scale_y_continuous(trans='log10',limits=c(1e+05,1e+20)) +
  theme(panel.grid = element_blank()) +
  xlab(expression(mu)) +
  labs(subtitle = expression(alpha~"= 0.50")) + 
  ylab("Time to speciation")

m_2 <- 0.25
a_2 <- 0.1
p_2 <- 1-a_2

fun.4 <- function(x) 2/x
fun.5 <- function(x) m_2/x^2
fun.6 <- function(x) 1/((a_2*(x/2))+(p_2*(x^2/m_2)))

p2 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  theme_bw() +
  stat_function(fun = fun.4, aes(linetype="isolation")) +
  stat_function(fun = fun.5, aes(linetype="continuous")) +
  stat_function(fun = fun.6, aes(linetype="periodic")) +
  scale_linetype_discrete(guide=FALSE) +  
  scale_x_continuous(trans='log10',limits=c(1e-10,1e-5)) +
  scale_y_continuous(trans='log10',limits=c(1e+05,1e+20)) +
  theme(panel.grid = element_blank()) +
  xlab(expression(mu)) +
  ylab(element_blank()) +
  labs(subtitle = expression(alpha~"= 0.1")) + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

m_3 <- 0.25
a_3 <- 0.01
p_3 <- 1-a_3

fun.7 <- function(x) 2/x
fun.8 <- function(x) m_3/x^2
fun.9 <- function(x) 1/((a_3*(x/2))+(p_3*(x^2/m_3)))

p3 <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
  theme_bw() +
  stat_function(fun = fun.7, aes(linetype="isolation")) +
  stat_function(fun = fun.8, aes(linetype="continuous")) +
  stat_function(fun = fun.9, aes(linetype="periodic")) +
  scale_linetype_discrete(name = "gene flow") +  
  scale_x_continuous(trans='log10',limits=c(1e-10,1e-5)) +
  scale_y_continuous(trans='log10',limits=c(1e+05,1e+20)) +
  theme(panel.grid = element_blank()) +
  xlab(expression(mu)) +
  ylab(element_blank()) +
  labs(subtitle = expression(alpha~"= 0.01")) + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

fig <- p1 + p2 + p3
ggsave("~/Dropbox/speciation_cyclical_migration/figures/waiting_time.pdf", plot=fig, width = 10, height = 4, units = "in")
