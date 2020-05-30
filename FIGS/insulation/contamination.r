library(Cairo)


CairoSVG( 'contamination', width=2.5, height=5 )
par( mfrow=c(2,1), mai=c(1.5,.75,0,0) )

data = read.table("W:/2013_06_wnt_tgfb_crosstalk/TE20/TE20_summaryStats.csv",
                  sep=',',header=T,stringsAsFactors=F)
pValCutoff = .05
treatments = c( 'ctrl','Wnt3A.CFHP','Wnt3A',
                'Wnt3A+aTGFB','Wnt5A','Wnt5A.CF' )
outNames   = c( 'ctrl','Wnt3A(CF)','Wnt3A(LP)',
                'Wnt3A+aTGFB','Wnt5A','Wnt5A(LP/CF)' )
channels   = 2:3
readouts   = c('Smad2/3','B-catenin')
names(channels) = readouts
colors     = c('blue','darkred')
names(colors) = readouts

for( readout in readouts ){
  # Normalize the data to [0,1]
  feature    = paste('sum_nucleus.median',channels[readout],sep='.')
  values     = data[,feature]
  values     = values - mean(data[data$treatment == 'ctrl',feature])
  values     = values / (mean(data[data$treatment == 'Wnt3A',feature])-
                         mean(data[data$treatment == 'ctrl',feature]))
  

  
  
  plot( NA,NA, xlim=c(.5,length(treatments)+.5), xaxt='n',
        ylim=range(values)*1.05, yaxt='n', ylab=readout, xlab='')
  axis(1,at=1:length(treatments),outNames,las=3)
  axis(2,at=c(0,1))
  abline(h=0,col='gray',lty=3)
  abline(h=1,col='gray',lty=3)
  
  for( tIdx in 1:length(treatments)){
    theseValues = values[data$treatment == treatments[tIdx]]
    lines(c(tIdx,tIdx),mean(theseValues)+c(-1,1)*sd(theseValues),col=colors[readout])
    points( tIdx, mean(theseValues), col = colors[readout], pch=16)
    
    isEqCtrl = t.test(data[data$treatment == 'ctrl',feature],
                      data[data$treatment == treatments[tIdx],feature])$p.value
    if(isEqCtrl<pValCutoff){
      isEqCtrl = '*'
    }else{isEqCtrl = ''}
    
    axis(1,at=tIdx,isEqCtrl,tick=F)
  }
}



dev.off()



