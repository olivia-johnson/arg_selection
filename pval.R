#!/usr/bin/env Rscript
filename <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(ggplot2)

setwd("~/Documents/arg_selection/")

dt=fread(filename)
dt[, pval:=pchisq(l_ratio, df = 1, lower.tail = FALSE)]
dt[,ppval:=-log10(pval)]


tdiff=ggplot(dt, aes(x=Pos, y=l_ratio))+
  geom_line() + 
  theme_bw()
ggsave(paste0("tdiff_",filename,".pdf"), tdiff)

plot=ggplot(dt, aes(x=Pos, y=ppval))+
  geom_line() + 
  theme_bw()
ggsave(paste0("pval_",filename,".pdf"), plot)
