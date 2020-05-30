######
# DOSE-RESPONSE CURVES for SKMEL2 and HCEC
######

library(Cairo)

celltypes  = c('HCEC','SKMEL2')
treatments = c('TGFB3','Wnt3A','BMP4')
jitter     = c(-.01,.01)
names(jitter) = celltypes

cols = c('black','darkred')
names(cols) = celltypes


CairoSVG( 'dose', width=6, height=2 )
par( mfrow=c(1,3) )




########### TGFB3 DOSE-RESPONSE ############



plate    = 'W:/2013_06_wnt_tgfb_crosstalk/TI03/TI03_summaryStats.csv'

allData  = read.table( plate, sep=',', header=T, stringsAsFactors=F )

treatments = c('TGFB3')
channel    = 2
marker     = 'marker.2'

for( treatment in treatments ){
  
  concentrations   = sort( unique(allData$concentration[allData$treatment==treatment]))
  readouts         = unique(allData[allData$treatment==treatment,c('marker.2')])
  readout          = readouts[readouts != '']
  feature          = paste('sum_nucleus.median',channel,sep='.')
  
  logUnits         = log10(concentrations)
  logUnits[is.infinite(logUnits)]  = min(logUnits[!is.infinite(logUnits)])-1
  
  plot( NA,NA, xlab=paste('log10(',treatment,')',sep=''), xaxt='n',
        ylab = readout, xlim=range(logUnits), ylim=c(-.2,1.2),yaxt='n')
  axis(2,at=c(0,.5,1))
  axis(1,at=seq(min(round(logUnits)),max(round(logUnits)),by=1))
  abline( h=.5, col='gray', lty=3 )
  abline( h=1, col='gray', lty=3 )
  abline( h=0, col='gray', lty=3 )
  
  
  for( celltype in celltypes ){
    
    dataSubset = allData[allData$celltype==celltype & allData$treatment==treatment, c(feature,'concentration')]
    
    dataSubset = dataSubset[order(dataSubset$concentration),]
    
    responses  = matrix(dataSubset[,feature],ncol=3,byrow=T)
    
    # Normalize to 0 - 1
    responses = responses - mean(responses[1,])
    responses = responses / mean(responses[dim(responses)[1],])
    
    means     = apply(responses,1,mean)
    sds       = apply(responses,1,sd)
    
    for( idx in 1:length(logUnits)){
      lines( rep(logUnits[idx],2)+jitter[celltype], means[idx] + c(-1,1)*sds[idx], col= cols[celltype] )
    }
    
    lines( logUnits, means, lwd=2, col= cols[celltype] )
  }
  legend('topleft',celltypes,text.col=cols[celltypes],bty='n')
}  







########### Wnt3A and BMP4 DOSE-RESPONSE ############




plates     = paste('W:/2013_06_wnt_tgfb_crosstalk/TI01/TI01.',1:2,'_summaryStats.csv',sep='')

allData    = read.table( plates[1], sep=',', header=T, stringsAsFactors=F )
allData    = rbind( allData, read.table( plates[2], sep=',', header=T, stringsAsFactors=F ))
treatments = c('Wnt3A','BMP4')


for( treatment in treatments ){
  
  concentrations   = sort( unique(allData$concentration[allData$treatment==treatment]))
  readouts         = unique(allData[allData$treatment==treatment,c('marker.2','marker.3')])
  markers          = names(readouts)[readouts != '']
  readout          = readouts[readouts != '']
  channel          = as.numeric(sub('marker.','',markers))
  feature          = paste('sum_nucleus.median',channel,sep='.')
  
  logUnits         = log10(concentrations)
  logUnits[is.infinite(logUnits)]  = min(logUnits[!is.infinite(logUnits)])-1
  
  plot( NA,NA, xlab=paste('log10(',treatment,')',sep=''), xaxt='n',
        ylab = readout, xlim=range(logUnits), ylim=c(-.1,1.1),yaxt='n')
  axis(2,at=c(0,.5,1))
  axis(1,at=seq(min(round(logUnits)),max(round(logUnits)),by=1))
  abline( h=.5, col='gray', lty=3 )
  abline( h=1, col='gray', lty=3 )
  abline( h=0, col='gray', lty=3 )
  
  
  for( celltype in celltypes ){
    
    dataSubset = allData[allData$celltype==celltype & allData$treatment==treatment, c(feature,'concentration')]
    
    dataSubset = dataSubset[order(dataSubset$concentration),]
    
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

