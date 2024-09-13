## GRoSS analysis ##

filename <- commandArgs(trailingOnly = TRUE)


library(data.table)

setwd("~/Documents/arg_selection/")

dt=fread(filename)
dt[, count:=paste0(A,",",a), by=c("pop","time","pos")]
dt[, pop_id:=paste0("p", pop,"_",time), by=c("pop", "time")]
grossinput=dcast(dt, pos~pop_id, value.var='count',fill = "0,0")
grossinput[, CHROM:=1]
grossinput[, SNPID:=pos, by="pos"]
setnames(grossinput, 'pos', 'POS')
setcolorder(grossinput, c("CHROM","POS","SNPID"))

fwrite(grossinput, filename, sep = "\t")


