### 20140127 by Adam D Coster
### PURPOSE
###   To analyze the TE17.1,2,5,6 data (in form of dose-response curves)

library(Cairo)

celltypes  = c('HCEC','SKMEL2')
treatments = c('TGFB3','Wnt3A','BMP4')
jitter     = c(-.01,.01)
names(jitter) = celltypes

cols = c('black','darkred')
names(cols) = celltypes


CairoSVG( 'time', width=6, height=2 )
par( mfrow=c(1,3) )


########### TIMECOURSES ############

plates     = paste('W:/2013_06_wnt_tgfb_crosstalk/TI05/TI05.',1:2,'_summaryStats.csv',sep='')

allData    = read.table( plates[1], sep=',', header=T, stringsAsFactors=F )
allData    = rbind( allData, read.table( plates[2], sep=',', header=T, stringsAsFactors=F ))


for( treatment in treatments ){
  
  concentrations   = sort( unique(allData$time.treatment[allData$treatment==treatment]))
  readouts         = unique(allData[allData$treatment==treatment,c('marker.2','marker.3')])
  markers          = names(readouts)[readouts != '']
  readout          = readouts[readouts != '']
  channel          = as.numeric(sub('marker.','',markers))
  feature          = paste('sum_nucleus.median',channel,sep='.')
  
  logUnits         = concentrations
  
  plot( NA,NA, xlab='time (min)', xaxt='n',
        ylab = readout, xlim=range(logUnits), ylim=c(-.2,1.2),yaxt='n')
  axis(2,at=c(0,.5,1))
  axis(1,at=logUnits)
  abline( h=.5, col='gray', lty=3 )
  abline( h=1, col='gray', lty=3 )
  abline( h=0, col='gray', lty=3 )
  
  for( celltype in celltypes ){
    
    dataSubset = allData[allData$celltype==celltype & allData$treatment==treatment, c(feature,'time.treatment')]
    
    dataSubset = dataSubset[order(dataSubset$time.treatment),]
    
    responses  = matrix(dataSubset[,feature],ncol=3,byrow=T)
    
    # Normalize to 0 - 1
    
    
    responses = responses - min(apply(responses,1,mean))
    responses = responses / max(apply(responses,1,mean))
    
    means     = apply(responses,1,mean)
    sds       = apply(responses,1,sd)
    
    for( idx in 1:length(logUnits)){
      lines( rep(logUnits[idx],2)+jitter[celltype], means[idx] + c(-1,1)*sds[idx], col= cols[celltype] )
    }
    
    lines( logUnits, means, lwd=2, col= cols[celltype] )
  }
  legend('topleft',celltypes,text.col=cols[celltypes],bty='n')
}  



dev.off()



