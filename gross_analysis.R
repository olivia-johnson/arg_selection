## GRoSS analysis ##

library(data.table)

setwd("~/Documents/arg_selection/")

ac_files=list.files(path ="~/Documents/arg_selection/",pattern ="ac_")
dt=fread(ac_files[1])
dt[, count:=paste0(A,",",a), by=c("pop","time","pos")]
dt[, pop_id:=paste0("p", pop,"_",time), by=c("pop", "time")]
grossinput=dcast(dt, pos~pop_id, value.var='count')
grossinput[, CHROM:=1]
setnames(grossinput, 'pos', 'POS')

fwrite(grossinput, ac_files[1], sep = "\t")
