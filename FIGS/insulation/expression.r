library(Cairo)


CairoSVG( 'expression', width=3, height=4.5 )
layout(matrix(c(1,1,1,1,2,2,2,
                3,3,3,3,4,4,4),nrow=2,byrow=T))
par(mai=c(.75,.75,.75,0))

pValCutoff = .05

hcec  = read.table('W:/2013_06_wnt_tgfb_crosstalk/qPCR/HCEC/HCEC_2hr.csv',sep=',',header=T,stringsAsFactors=F)
hcec$Celltype = rep('HCEC',dim(hcec)[1])
skmel = read.table("W:/2013_06_wnt_tgfb_crosstalk/qPCR/SKMEL2+MALME3/SKEML2+MALME3M.csv",sep=',',header=T,stringsAsFactors=F)
skmel = skmel[skmel$Celltype=='SKMEL2',]

data  = rbind(hcec,skmel)


###### PLOT SMAD7 ######
treatments = list(tgfb=c('ctrl+ctrl','ctrl+TGFB3','BMP4+ctrl'),
                  wnt =c('ctrl+ctrl','Wnt3A+ctrl'))
outNames   = list(tgfb=c('ctrl','TGFB3','BMP4'),
                  wnt =c('ctrl','Wnt3A'))
readouts   = list(tgfb='SMAD7',wnt='AXIN2')

colors        = c('black','darkred','blue','darkgreen')
names(colors) = c('ctrl+ctrl','Wnt3A+ctrl','ctrl+TGFB3','BMP4+ctrl')

for( celltype in c('HCEC','SKMEL2')){
  for( input in c('tgfb','wnt')){
    readout      = readouts[[input]]
    subdata      = data[data$Detector==readouts[input] &
                        data$Celltype==celltype, ]  
    subdata$RNA  = 2^(-subdata$dCT)
    ctrls        = subdata[subdata$Treatment == 'ctrl+ctrl', 'RNA']
    subdata$RNA  = subdata$RNA / mean(ctrls) 
    ctrls        = ctrls / mean(ctrls) 
  
    subdata$RNA  = log2(subdata$RNA)
    
    ylims        = range(subdata$RNA[subdata$Treatment %in% treatments[[input]]])*1.1
    yvals        = 0:(floor(ylims[2]))
    
    plot( NA,NA, xlim=c(.5,length(treatments[[input]])+.5), xlab='',
         ylim=ylims,
         main=celltype, yaxt='n',
         ylab=paste('log2(relative',readout,'mRNA)'),xaxt='n')
    axis( 1, at=1:length(treatments[[input]]),treatments[[input]],las=3)
    axis( 2, at=yvals )
    abline(h=0,col='gray',lty=3)
    
    for( x in 1:length(treatments[[input]])){
      values = subdata$RNA[subdata$Treatment==treatments[[input]][x]]
      lines( c(x,x), c(-1,1)*sd(values)+mean(values), col=colors[treatments[[input]][x]] )
      points( x, mean(values), pch=16, col=colors[treatments[[input]][x]])
    
      # SIGNIFICANCE
      isEqCtrl = t.test(subdata$RNA[subdata$Treatment==treatments[[input]][x]],
                        subdata$RNA[subdata$Treatment=='ctrl+ctrl'])$p.value
      if( isEqCtrl < pValCutoff){
        isEqCtrl = '*'
      }else{
        isEqCtrl = ''
      }
      
      axis(1,at=x,isEqCtrl,tick=F,col.axis=colors[treatments[[input]][x]])
    }
  }
}
 



dev.off()



