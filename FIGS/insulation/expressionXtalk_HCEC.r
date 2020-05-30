library(Cairo)


CairoSVG( 'expressionXtalk_HCEC', width=6, height=2.5 )
par( mfrow=c(1,4), mai=c(1,.5,.5,0) )

pValCutoff = .05

hcec2   = read.table( 'W:/2013_06_wnt_tgfb_crosstalk/qPCR/HCEC/HCEC_2hr.csv',sep=',',header=T,stringsAsFactors=F)
hcec2$Celltype = rep( 'HCEC',dim(hcec2)[1])
hcec2$Time     = rep( 2,dim(hcec2)[1])
hcec6   = read.table( 'W:/2013_06_wnt_tgfb_crosstalk/qPCR/HCEC/HCEC_6hr.csv',sep=',',header=T,stringsAsFactors=F)
hcec6$Celltype = rep( 'HCEC',dim(hcec6)[1])
hcec6$Time     = rep( 6,dim(hcec6)[1])
data   = rbind(hcec2,hcec6)


treatments = c("ctrl+ctrl","BMP4+ctrl","ctrl+TGFB3","BMP4+TGFB3",
               "TGFB3+TGFB3","Wnt3A+ctrl","Wnt3A+TGFB3")
outNames   = c("ctrl","BMP4","TGFB3","BMP4+TGFB3",
               "TGFB3+TGFB3","Wnt3A","Wnt3A+TGFB3")
readouts   = c('AXIN2','SMAD7')

colors        = c('darkred','blue')


for( rIdx in 1:length(readouts)){
  for( timepoint in c(2,6) ){   
    readout      = readouts[rIdx]
    subdata      = data[data$Detector==readouts[rIdx] &
                        data$Time==timepoint, ]  
    subdata$RNA  = 2^(-subdata$dCT)
    ctrls        = subdata[subdata$Treatment == 'ctrl+ctrl', 'RNA']
    subdata$RNA  = subdata$RNA / mean(ctrls) 
    ctrls        = ctrls / mean(ctrls) 
  
    subdata$RNA  = log2(subdata$RNA)
    
    ylims        = range(subdata$RNA[subdata$Treatment %in% treatments])*1.1
    yvals        = (ceiling(ylims[1])):(floor(ylims[2]))
    
    plot( NA,NA, xlim=c(.5,length(treatments)+.5), xlab='',
         ylim=ylims,
         main=paste(timepoint,'hr'), yaxt='n',
         ylab=paste('log2(relative',readout,'mRNA)'),xaxt='n')
    axis( 1, at=1:length(treatments),treatments,las=3)
    axis( 2, at=yvals )
    abline(h=0,col='gray',lty=3)
    
    for( x in 1:length(treatments)){
      values = subdata$RNA[subdata$Treatment==treatments[x]]
      lines( c(x,x), c(-1,1)*sd(values)+mean(values), col=colors[rIdx] )
      points( x, mean(values), pch=16,col=colors[rIdx])
    
      # SIGNIFICANCE
      isEqCtrl = t.test(subdata$RNA[subdata$Treatment==treatments[x]],
                        subdata$RNA[subdata$Treatment=='ctrl+ctrl'])$p.value
      if( isEqCtrl < pValCutoff){
        isEqCtrl = '*'
      }else{
        isEqCtrl = ''
      }
      
      axis(1,at=x,isEqCtrl,tick=F,col.axis=colors[rIdx])
    }
  }
}
 



dev.off()



