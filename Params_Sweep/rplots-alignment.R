if (!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
  # display.brewer.all() To see all the available colors
  library(dplyr)
}

setwd("~/Dropbox/University life/PHD quest/Projects_Current/GATK optimization/My scripts")


parameters = data.frame(
  symbol=c('A','B','c','d','D','E','k','L','m','O','r','T','U','w','W'), 
  meaning=c('Matching score','Mismatch penalty','Max Occurance','zDropoff','Chain fraction', 
            'Gap Extention penalty', 'Minimum seed length', 'Clipping penalty','Mate SW rounds', 
            'Gap Open penalty', 'Seed split ratio', 'Min alignment quality', 'Unpair penalty',
            'Bandwidth','Seed bases count'),
  stringsAsFactors = F)

rawdata = read.table( pipe('ssh -l hpcbiointern01 biocluster.igb.illinois.edu "cat /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/results/changing_parameters.txt"'),
                   stringsAsFactors = F,header = T)

data = rawdata
data = cbind(data,fraction=data$Total_aligned/data$Total_reads)
default = data[(data$parameter=='default'),]
data = data[!(data$parameter=='default'),]

m = match(data$parameter, parameters$symbol)
data$header = parameters$meaning[m]

data = data[,c('header','parameter','value','Mean_MAPQ','fraction','Time')]

data = arrange(data,parameter,value)

attempts = split(data,data$parameter)

pdf('myplots.pdf',title = 'Effect of parameter chages in the alignment file')
par(mfrow = c(2,2))
for (i in seq(1:length(attempts))) {
  plot(attempts[[i]][,4],attempts[[i]][,5], main=attempts[[i]][1,1], xlab='Mean MAPQ',
       ylab='Fraction mapped', col=brewer.pal(nrow(attempts[[i]]),"PuOr"),pch=19)
  text(attempts[[i]][,4],attempts[[i]][,5],labels=attempts[[i]][,'value'],cex=.7,pos=3)
  points(default$Mean_MAPQ,default$fraction,pch=13,col=brewer.pal(3,"Set1"))
  
  text(default$Mean_MAPQ,default$fraction,labels='default',cex=.7,pos=1,
       col=brewer.pal(3,"Set1"))
  abline(h=default$fraction,col=brewer.pal(3,"Set1"))
  abline(v=default$Mean_MAPQ,col=brewer.pal(3,"Set1"))
  
  cx = default$Mean_MAPQ + abs((default$Mean_MAPQ - max(attempts[[i]][,4]))/2)
  cy = default$fraction + abs((default$fraction-max(attempts[[i]][,5]))/2)
  
  text(cx,1.2,labels='optimal area',cex=.7,pos=1, col=brewer.pal(3,"Set1"))
  cl = col2rgb( brewer.pal(3,"Set1"), alpha =F) [,1] /255
  cl = rgb(cl[1],cl[2],cl[3],.1)
  rect(xleft = default$Mean_MAPQ, ybottom = default$fraction,
       xright = max(attempts[[i]][,'Mean_MAPQ']), ytop = 1.4,
       col = cl, border=NA)
}
dev.off()
