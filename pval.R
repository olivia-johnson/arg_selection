library(data.table)
library(ggplot2)

setwd(setwd("~/Documents/arg_selection/")

dt=fread("windowedlr_single_0_s1.0_sT17000_sE20000_sP1_cF1.0_cFT17499_admix0.0_sSize100.txt")
dt[, pval:=pchisq(l_ratio, df = 1, lower.tail = FALSE)]


tdiff=ggplot(dt, aes(x=Pos, y=l_ratio))+
  geom_line() + 
  theme_bw()
ggsave("tdiff_0_s1.0_sT17000_sE20000_sP1_cF1.0_cFT17499_admix0.0_sSize100.pdf", tdiff)

plot=ggplot(dt, aes(x=Pos, y=pval))+
  geom_line() + scale_y_log10() +
  theme_bw()
ggsave("pval_s1.0_s1.0_sT17000_sE20000_sP1_cF1.0_cFT17499_admix0.0_sSize100.pdf", plot)
