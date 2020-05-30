library(Cairo)

timepoints = c('1','2','18')
datafiles  = c('data/TI04.60_summaryStats.csv',
               'data/TI04.120_summaryStats.csv',
               'data/TI16.18_TB_summaryStats.csv')
reps       = 3
pValCutoff = 0.05
readouts   = c('Smad23')
colors     = c('blue')
colnames(colors)=readouts
channels   = c(2,2)
celltypes  = c('HCEC','SKMEL2','MALME3M','IEC6')
allTreatments = list(Smad23=c('ctrl+ctrl','ctrl+BMP4',
                           'TGFB3+ctrl','TGFB3+BMP4'),
                     pSmad158=c('ctrl+ctrl','TGFB3+ctrl',
                           'ctrl+BMP4','TGFB3+BMP4'))

setToOne   = c('TGFB3+ctrl','ctrl+BMP4')
feature    = 'sum_nucleus'
allYlims   = list(Smad23  = c(-.15, 1.3 ),
                  pSmad15 = c(-.35, 3.3) )

CairoSVG( 'bmpTgfbInsulationSmad23',width=8,height=2)
par(mfrow=c(1,4))

for( rIdx in 1:length(readouts) ){
  ylims      = allYlims[[readouts[rIdx]]]
  treatments = allTreatments[[readouts[rIdx]]]
  
  for( cIdx in 1:length(celltypes) ){

    readout   = readouts[ rIdx]
    channel   = channels[ rIdx]
    marker    = paste('marker',channel,sep='.')
    celltype  = celltypes[cIdx]
    
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
                        !(data$experiment %in% c(1,4,5)) &
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

