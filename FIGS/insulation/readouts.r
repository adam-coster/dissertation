library(Cairo)


CairoSVG( 'readouts', width=8, height=3 )
par(mfrow=c(1,5))

data = read.table("W:/2013_06_wnt_tgfb_crosstalk/TI08/TI08_summaryStats.csv",
                  sep=',',header=T,stringsAsFactors=F)
data = data[data$experiment %in% c('TGFB.readouts','BMP.readouts'),]
ylims      = c(-.25,1.25)
channel    = 2
feature    = 'sum_nucleus.median.2'
treatments = c('ctrl','BMP4','TGFB3')
readouts   = c('Smad2','Smad23','pSmad23',
               'Smad1','pSmad158')
antibodies = c('Smad2','Smad2/3','pSmad2/3',
               'Smad1','pSmad1/5')
names(antibodies) = readouts

colors        = c('black','darkgreen','blue')
names(colors) = treatments

for( readout in readouts ){
  # Normalize the data to [0,1]
  subdata    = data[data$marker.2==readout,]
  ctrls      = subdata$sum_nucleus.median.2[subdata$treatment=='ctrl']
  subdata[subdata$marker.2==readout,feature] = subdata[subdata$marker.2==readout,feature]-mean(ctrls)
  vals       = subdata$sum_nucleus.median.2[subdata$treatment!='ctrl']
  subdata[subdata$marker.2==readout,feature] = subdata[subdata$marker.2==readout,feature]/mean(vals)
  
  treatment = unique(subdata$treatment[subdata$treatment!='ctrl'])
  
  plot( NA,NA, xlim=c(.75,2.25), xaxt='n',
        ylim=ylims,main=readout,
        yaxt='n', ylab='',
        xlab=treatment)
  axis(1,at=1:2,c('-','+'))
  axis(2,at=c(0,1))
  abline(h=0,col='gray',lty=3)
  abline(h=1,col='gray',lty=3)
  
  ctrls = subdata$sum_nucleus.median.2[subdata$treatment=='ctrl']
  lines(c(1,1),mean(ctrls)+c(-1,1)*sd(ctrls),col='gray',lwd=2)
  points( 1, mean(ctrls), col = colors['ctrl'], pch=16)
  
  vals = subdata$sum_nucleus.median.2[subdata$treatment==treatment]
  lines(c(2,2),mean(vals)+c(-1,1)*sd(vals),col='gray',lwd=2)
  points( 2, mean(vals), col = colors[treatment], pch=16)
  
}




dev.off()



