## Need to generate a heatmap that has all 4 celltypes
## at the 1, 2,and 18hr marks. Need to generate one
## heatmap per timepoint and readout. The readouts are
## Smad23, and CTNNB1.
library(Cairo)

timepoints = c('1','2','18')
datafiles  = c('data/TI04.60_summaryStats.csv',
               'data/TI04.120_summaryStats.csv',
               'data/TI16.18_TW_summaryStats.csv')
reps       = 3
pValCutoff = 0.05
readouts   = c('Smad23','CTNNB1')
colors     = c('blue','darkred')
channels   = 2:3
celltypes  = c('HCEC','SKMEL2','MALME3M','IEC6')
allTreatments = list(Smad23=c('ctrl+ctrl','ctrl+Wnt3A',
                           'TGFB3+ctrl','TGFB3+Wnt3A'),
                     CTNNB1=c('ctrl+ctrl','TGFB3+ctrl',
                           'ctrl+Wnt3A','TGFB3+Wnt3A'))

setToOne   = c('TGFB3+ctrl','ctrl+Wnt3A')
feature    = 'sum_nucleus'
allYlims   = list(Smad23 = list(HCEC   =c(-.15, 1.3 ),
                                SKMEL2 =c(-.15, 1.3 ),
                                MALME3M=c(-.15, 1.3 ),
                                IEC6   =c(-.15, 1.3 )),
                  CTNNB1 = list(HCEC   =c(-.25, 3.5 ),
                                SKMEL2 =c(-.25, 1.6 ),
                                MALME3M=c(-.25, 1.6 ),
                                IEC6   =c(-.25, 1.6 )) )

CairoSVG( 'wntTgfbInsulation',width=8,height=4)
par(mfrow=c(2,4))

for( rIdx in 1:length(readouts) ){
  treatments = allTreatments[[readouts[rIdx]]]
  
  for( cIdx in 1:length(celltypes) ){

    readout   = readouts[ rIdx]
    channel   = channels[ rIdx]
    marker    = paste('marker',channel,sep='.')
    celltype  = celltypes[cIdx]
    ylims     = allYlims[[readouts[rIdx]]][[celltype]]
    
    plot(NA,NA,ylim=ylims, xlim=c(.75,length(treatments)*length(timepoints)+.25),
         main=celltype,ylab=readouts[rIdx],xaxt='n',xlab='',yaxt='n')
    abline(h=0,col='gray',lty=3)
    abline(h=1,col='gray',lty=3)
    axis(2, at=c(0,1))

    
    for( treatIdx in 1:length(treatments)){
      
      for( tIdx in 1:length(timepoints)){
        
        
        data = read.table(datafiles[tIdx],sep=',',header=T,
                          stringsAsFactors=F )
      
        treatment = treatments[treatIdx]
        
        subdata  = data[data$celltype==celltype &
                        data$marker.2==readouts[1]  &
                        data$marker.3==readouts[2]  &
                        data$treatment %in% treatments,]
        
        # Convert to 0-1 scale
        thisFeature = paste(feature,'median',channel,sep='.')
        ctrls   = subdata[subdata$treatment == 'ctrl+ctrl',
                          thisFeature]
        subdata[,thisFeature] = subdata[,thisFeature]-mean(ctrls)
        subdata[,thisFeature] = subdata[,thisFeature]/
          mean(subdata[subdata$treatment==setToOne[rIdx],thisFeature])
        
        # Collect statistics
        isEqCtrl = t.test(subdata[subdata$treatment == 'ctrl+ctrl',thisFeature],
                          subdata[subdata$treatment == treatment,thisFeature])$p.value
        isEqOne  = t.test(subdata[subdata$treatment == setToOne[rIdx],thisFeature],
                          subdata[subdata$treatment == treatment,thisFeature])$p.value
        if(isEqCtrl<pValCutoff){
          isEqCtrl = '*'
        }else{isEqCtrl = '  '}
        
        if(isEqOne<pValCutoff){
          isEqOne = '*'
        }else{isEqOne = ' '}
        
        
        mu = mean(subdata[subdata$treatment == treatment,thisFeature])
        s  = sd(subdata[subdata$treatment == treatment,thisFeature])
        
        # Plot with first two treatment error bars up, second two down
        x  = (treatIdx-1)*length(timepoints)+tIdx
        errorY = mu+s*c(1,-1)
        points( x, mu, pch=16, col=c('gray',colors[rIdx],'gray',colors[rIdx])[treatIdx])
        lines( c(x,x), errorY, lwd=1,
               col=c('darkgray',colors[rIdx],'darkgray',colors[rIdx])[treatIdx])
        axis(1,at=x,
             paste(isEqOne,isEqCtrl,
                   sep=''), tick=F, family='mono', las =3, cex.axis=1.5,
             col.axis=c('darkgray',colors[rIdx],'darkgray',colors[rIdx])[treatIdx])
      }
    }
  }
}

dev.off()

