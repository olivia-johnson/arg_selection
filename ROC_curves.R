## gLike ROC

library(data.table)
library(ggplot2)

path="~/Documents/arg_selection/"
setwd(path)

## Collate files
data = NULL 
f_list  <- list.files(path =path, pattern ="twindowedlr_single")
unwanted <- list.files(path =path, pattern =".pdf")
f_list=setdiff(f_list, unwanted)
for (i in f_list){
  dd=fread(i)
  split=strsplit(i, "windowed")[[1]]
  dd[,method:=split[1]]
  split2=strsplit(split[2], "_")[[1]]
  dd[,rep:=as.numeric(split2[3])]
  dd[,selCoeff:=as.numeric(strsplit(split2[4],"s")[[1]][2])]
  dd[,selStart:=as.numeric(strsplit(split2[5],"sT")[[1]][2])]
  dd[,selEnd:=as.numeric(strsplit(split2[6],"sE")[[1]][2])]
  dd[,selPop:=as.numeric(strsplit(split2[7],"sP")[[1]][2])]
  dd[,condFreq:=as.numeric(strsplit(split2[8],"cF")[[1]][2])]
  dd[, condFreqTime:=as.numeric(strsplit(split2[9],"cFT")[[1]][2])]
  dd[, admix:=as.numeric(strsplit(split2[10],"admix")[[1]][2])]
  dd[, sampSize:=as.numeric(strsplit(strsplit(split2[11],"sSize")[[1]][2], ".txt")[[1]][1])]
  dd[, selLength:=selEnd-selStart, by=c("selStart", "selEnd")]
  dd[,lab:=paste(split2[5], split2[6],split2[7],split2[8],split2[9], paste0("sSize",sampSize), sep="_")]
  data=rbind(dd, data)
}

data[, max:=max(l_ratio), by=c("method", "selCoeff", "lab", "rep")]



## check missing replicates
for (m in unique(data$method)) {
  for (s in unique(data[selCoeff!=0.0]$selCoeff))  {
    for (l in unique(data$lab)) {
      neutral = data[selCoeff == 0.0 & method == m  & lab == l]
      selection = data[selCoeff == s & method == m & lab == l]
      if (neutral[Pos==25000,.N]!=10){print(paste("Neutral", l, "only", neutral[Pos==25000,.N], "replicates - Missing rep", (setdiff(0:9, unique(neutral$rep)))))}
      if (selection[Pos==25000,.N]!=10){print(paste("Sel", s, l, "only", selection[Pos==25000,.N], "replicates- Missing rep", (setdiff(0:9, unique(selection$rep)))))}
    }}}
      
## Calc TPR and FPR     
ROC_data = list()
for (m in unique(data$method)) {
  for (s in unique(data[selCoeff!=0.0]$selCoeff))  {
    for (l in unique(data$lab)) {
      neutral = data[selCoeff == 0.0 & method == m  & lab == l]
      selection = data[selCoeff == s & method == m & lab == l]     
      if (min(neutral$max) >= 0 & nrow(neutral) >=1 & min(selection$max) >= 0 & nrow(selection) >=1) {
        LR_cutoff = seq(0, max(selection$max), 0.1)
        FPR <- numeric(length(LR_cutoff))
        TPR <- numeric(length(LR_cutoff))
      } else { next
      }
      for (i in seq_along(LR_cutoff)) {
        cat(paste0("calculating for ", s, " ",m, l, " at cutoff ", i,"\n"))
        neutral_positive = neutral$max >= LR_cutoff[i]
        sel_positive = selection$max >= LR_cutoff[i]
        TP = sum(sel_positive)
        FP = sum(neutral_positive)
        TN = length(neutral$max) - FP
        FN = length(selection$max) - TP
        FPR[i] = ifelse(FP + TN == 0, 0, 1 - (TN / (FP + TN)))
        TPR[i] = ifelse(TP + FN == 0, 0, TP / (TP + FN))
      }
      dt = data.table(selCoeff = s, method = m, lab = l, selStart=unique(selection$selStart),selEnd=unique(selection$selEnd),cF=unique(selection$condFreq),admixture=unique(selection$admix),Threshold = LR_cutoff, FPR = FPR, TPR = TPR)
      ROC_data[[paste(s, m, l, sep = "_")]] = dt
    }
  }
}

ROC_data=rbindlist((ROC_data))
cbp1 <- c("#000000", "#E69F00", "#0072B2", "#009E73",
          "#CC79A7","#F0E442", "#56B4E9", "#D55E00", "purple", "#882255")

plot = ggplot(ROC_data[admixture==0&cF==1.0][order(FPR, TPR)], 
       aes(FPR, TPR, color = as.character(selCoeff))) + 
  scale_linetype_manual(values = c(1, 6)) +
  geom_line() + facet_grid(selStart+selEnd + method~ selCoeff, scales = "free") +
  theme_bw() + labs(col="Selection\nCoefficient")+ lims(x=c(0,1), y=c(0,1))+
   scale_colour_manual(values = cbp1) +ggtitle("ROC curve - no admixture, conditional frequency 1")
plot
ggsave(filename="ROC_noadmix.pdf", plot, width = 12, height =8)
ggsave(filename="ROC_noadmix.jpg", plot, width = 18, height =10)