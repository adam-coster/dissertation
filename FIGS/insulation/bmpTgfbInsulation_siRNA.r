library(Cairo)

siRNAs    = c('scramble','Smad4')
datafile  = 'data/TI10_summaryStats.csv'
data       = read.table( datafile, sep=',', header=T, stringsAsFactors=F )
reps       = 3
pValCutoff = 0.05
readouts   = c('Smad23','pSmad158')
colors     = c('blue','darkgreen')
channels   = c(2,2,3)
celltypes  = c('HCEC')
allTreatments = list( Smad23 =c('ctrl+ctrl','ctrl+BMP4','TGFB3+ctrl','TGFB3+BMP4'),
                  pSmad158=c('ctrl+ctrl','TGFB3+ctrl','ctrl+BMP4','TGFB3+BMP4'))

setToOne   = c('TGFB3+ctrl','ctrl+BMP4')
feature    = 'sum_nucleus'
allYlims   = list(Smad23   = c(-.1, 1.2 ),
                  pSmad158 = c(-.1, 2.2) )

CairoSVG( 'bmpTgfbInsulation_siRNA',width=6,height=3)
par(mfrow=c(1,2))

for( rIdx in 1:length(readouts) ){
  
  readout   = readouts[ rIdx]
  channel   = channels[ rIdx]
  marker    = paste('marker',channel,sep='.')
  
  treatments = allTreatments[[readout]]
  
  ylims      = allYlims[[readouts[rIdx]]]
  plot( NA,NA,ylim=ylims, xlim=c(.75,length(treatments)*length(siRNAs)+.25),
        main='',ylab=readouts[rIdx],xaxt='n',xlab='',yaxt='n')
  abline(h=0,col='gray',lty=3)
  abline(h=1,col='gray',lty=3)
  axis(2, at=c(0,1))
  
  for( rnaIdx in 1:length(siRNAs) ){

    siRNA     = siRNAs[ rnaIdx]

    
    for( treatIdx in 1:length(treatments)){
      
      treatment = treatments[treatIdx]
      
      subdata  = data[data$siRNA==siRNA &
                      data[,marker]==readout  &
                      data$treatment %in% treatments,]
      
      # Convert to 0-1 scale
      thisFeature = paste(feature,'median',channel,sep='.')
      ctrls   = data[data$siRNA=='scramble' &
                        data[,marker]==readout  &
                        data$treatment == 'ctrl+ctrl',
                        thisFeature]
      subdata[,thisFeature] = subdata[,thisFeature]-mean(ctrls)
      subdata[,thisFeature] = subdata[,thisFeature]/(
                mean(data[data$siRNA=='scramble' &
                            data[,marker]==readout  &
                            data$treatment == setToOne[rIdx],
                            thisFeature]) - mean(ctrls))
      
      pRef = (data[data$siRNA=='scramble' &
                  data[,marker]==readout  &
                  data$treatment == treatment,
                  thisFeature]- mean(ctrls)) /(
                    mean(data[data$siRNA=='scramble' &
                                data[,marker]==readout  &
                                data$treatment == setToOne[rIdx],
                              thisFeature]) - mean(ctrls))
      
      isEqScr  = t.test(pRef,
                        subdata[subdata$treatment == treatment,thisFeature])$p.value
      #isEqCtrl = t.test(subdata[subdata$treatment == 'ctrl+ctrl',thisFeature],
      #                  subdata[subdata$treatment == treatment,thisFeature])$p.value
      #isEqOne  = t.test(subdata[subdata$treatment == setToOne[rIdx],thisFeature],
      #                  subdata[subdata$treatment == treatment,thisFeature])$p.value
      
      
      if(isEqScr<pValCutoff){
        isEqScr = '*'
      }else{isEqScr = ' '}
      
      mu = mean(subdata[subdata$treatment == treatment,thisFeature])
      s  = sd(subdata[subdata$treatment == treatment,thisFeature])
      
      # Plot with first two treatment error bars up, second two down
      x  = (rnaIdx-1)*length(treatments)+treatIdx
      
      if( siRNA == 'scramble'){
        theCols =c('darkgray',colors[rIdx],'darkgray',colors[rIdx])[treatIdx]
      }
      else{
        theCols = colors[rIdx]
      }
      
      
      errorY = mu+s*c(1,-1)
      points( x, mu, pch=16, col=theCols)
      lines( c(x,x), errorY, lwd=1,
             col=c('darkgray',colors[rIdx],'darkgray',colors[rIdx])[treatIdx])
      axis(1,at=x,
           paste(isEqScr, sep=''), tick=F, family='mono',
           las =3, cex.axis=1.5,
           col.axis= theCols )
    }
  }
}

dev.off()

