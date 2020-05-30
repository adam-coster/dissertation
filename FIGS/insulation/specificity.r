library(Cairo)


CairoSVG( 'specificity', width=4, height=1.5 )
par( mfrow=c(1,3),mai=c(1,.75,0,0) )

ylims = c(-.3,1.3)
xlims = c(.5,3.5)
pValCutoff = 0.05

#### TGFB3 & Wnt3A specificity ####
TWdata = read.table("W:/2013_06_wnt_tgfb_crosstalk/TI16/TI16_summaryStats.csv",sep=',',header=T,stringsAsFactors=F)

#### BMP4 specificity ####
Bdata  = read.table("W:/2013_06_wnt_tgfb_crosstalk/TI08/TI08_summaryStats.csv",sep=',',header=T,stringsAsFactors=F)


#### PLOT #####

## TGFB3->SMAD2/3 ##
ctrl.TGFB3 = c('ctrl','ctrl')
plus.TGFB3 = c('TGFB3','ctrl')
blok.TGFB3 = c('TGFB3','aTGFB')
xlab.TGFB3 = c('ctrl','TGFB3','TGFB3+aTGFB')
exp.TGFB3  = 'TI16.WaT'
cell.TGFB3 = 'HCEC'
feat.TGFB  = 'sum_nucleus.median.2' # Smad2/3 is in channel 2
Tdata      = TWdata[TWdata$experiment==exp.TGFB3,]
pair.TGFB3 = apply(Tdata[,c('treatment.1','treatment.2')],
                   1,paste,collapse='+')
# Transform the data to [0,1]
Tdata[,feat.TGFB] = Tdata[,feat.TGFB] - 
                    mean(Tdata[pair.TGFB3==paste(ctrl.TGFB3,collapse='+'),feat.TGFB])
Tdata[,feat.TGFB] = Tdata[,feat.TGFB] / 
  mean(Tdata[pair.TGFB3==paste(plus.TGFB3,collapse='+'),feat.TGFB])

plot(NA,NA,xlab='',ylab='Nuclear Smad2/3 (a.u.)',
     xaxt='n', xlim=xlims, ylim=ylims, yaxt='n' )
axis(2, at=c(0,1))
axis(1, at=1:3, xlab.TGFB3, las=3, tick=F )
abline(h=1,col='gray',lty=3)
abline(h=0,col='gray',lty=3)


for( tIdx in 1:3 ){
  treatment = apply(rbind(ctrl.TGFB3,plus.TGFB3,blok.TGFB3),
                    1,paste,collapse='+')[tIdx]
  values = Tdata[pair.TGFB3==treatment,feat.TGFB]
  lines(c(tIdx,tIdx),mean(values)+c(-1,1)*sd(values),col='blue')
  points( tIdx, mean(values), pch=16, col='blue' )

  # Check against control
  isEqCtrl = t.test( Tdata[pair.TGFB3==treatment,feat.TGFB],
                     Tdata[pair.TGFB3=='ctrl+ctrl',feat.TGFB])$p.value
  if(isEqCtrl<pValCutoff){
    isEqCtrl = '*'
  }else{isEqCtrl = ''}
  
  axis(1, at = tIdx, isEqCtrl, tick = F,
       las = 3, cex.axis = 1.5)
}








## Wnt3A->CTNNB1 ##
ctrl.Wnt3A = c('ctrl','ctrl')
plus.Wnt3A = c('Wnt3A','ctrl')
blok.Wnt3A = c('Wnt3A','DKK')
xlab.Wnt3A = c('ctrl','Wnt3A','Wnt3A+DKK')
exp.Wnt3A  = 'TI16.WaW'
cell.Wnt3A = 'SKMEL2'
feat.Wnt3A = 'sum_nucleus.median.3' # CTNNB1 is in channel 3
Wdata      = TWdata[TWdata$experiment==exp.Wnt3A,]
pair.Wnt3A = apply(Wdata[,c('treatment.1','treatment.2')],
                   1,paste,collapse='+')
# Transform the data to [0,1]
Wdata[,feat.Wnt3A] = Wdata[,feat.Wnt3A] - 
  mean(Wdata[pair.Wnt3A==paste(ctrl.Wnt3A,collapse='+'),feat.Wnt3A])
Wdata[,feat.Wnt3A] = Wdata[,feat.Wnt3A] / 
  mean(Wdata[pair.Wnt3A==paste(plus.Wnt3A,collapse='+'),feat.Wnt3A])

plot(NA,NA,xlab='',ylab='Nuclear B-catenin (a.u.)',
     xaxt='n', xlim=xlims, ylim=ylims, yaxt='n' )
axis(2, at=c(0,1))
axis(1, at=1:3, xlab.Wnt3A, las=3 )
abline(h=1,col='gray',lty=3)
abline(h=0,col='gray',lty=3)

for( tIdx in 1:3 ){
  treatment = apply(rbind(ctrl.Wnt3A,plus.Wnt3A,blok.Wnt3A),
                    1,paste,collapse='+')[tIdx]
  values = Wdata[pair.Wnt3A==treatment,feat.Wnt3A]
  lines(c(tIdx,tIdx),mean(values)+c(-1,1)*sd(values),col='darkred')
  points( tIdx, mean(values), pch=16, col='darkred' )
  # Check against control
  isEqCtrl = t.test( Wdata[pair.Wnt3A==treatment,feat.Wnt3A],
                     Wdata[pair.Wnt3A=='ctrl+ctrl',feat.Wnt3A])$p.value
  if(isEqCtrl<pValCutoff){
    isEqCtrl = '*'
  }else{isEqCtrl = ''}
  
  axis(1, at = tIdx, isEqCtrl, tick = F,
       las = 3, cex.axis = 1.5)
}





## BMP4->pSmad1/5 ##
ctrl.BMP4 = c('ctrl','ctrl')
plus.BMP4 = c('ctrl','BMP4')
blok.BMP4 = c('Noggin','BMP4')
xlab.BMP4 = c('ctrl','BMP4','BMP4+Noggin')
exp.BMP4  = 'TI16.WaW'
cell.BMP4 = 'SKMEL2'
feat.BMP4 = 'sum_nucleus.median.2' # pSmad1/5 is in channel 2
Bdata      = Bdata[Bdata$experiment=='specificity',]
pair.BMP4 = apply(Bdata[,c('treatment.1','treatment.2')],
                   1,paste,collapse='+')
# Transform the data to [0,1]
Bdata[, feat.BMP4] = Bdata[,feat.BMP4] - 
  mean(Bdata[pair.BMP4==paste(ctrl.BMP4,collapse='+'),feat.BMP4])
Bdata[, feat.BMP4] = Bdata[,feat.BMP4] / 
  mean(Bdata[pair.BMP4==paste(plus.BMP4,collapse='+'),feat.BMP4])

plot(NA,NA,xlab='',ylab='Nuclear pSmad1/5 (a.u.)',
     xaxt='n', xlim=xlims, ylim=ylims, yaxt='n' )
axis(2, at=c(0,1))
axis(1, at=1:3, xlab.BMP4, las=3 )
abline(h=1,col='gray',lty=3)
abline(h=0,col='gray',lty=3)

for( tIdx in 1:3 ){
  treatment = apply(rbind(ctrl.BMP4,plus.BMP4,blok.BMP4),
                    1,paste,collapse='+')[tIdx]
  values = Bdata[pair.BMP4==treatment,feat.BMP4]
  lines(c(tIdx,tIdx),mean(values)+c(-1,1)*sd(values),col='darkgreen' )
  points( tIdx, mean(values), pch=16, col='darkgreen' )
  # Check against control
  isEqCtrl = t.test( Wdata[pair.BMP4==treatment,feat.BMP4],
                     Wdata[pair.BMP4=='ctrl+ctrl',feat.BMP4])$p.value
  if(isEqCtrl<pValCutoff){
    isEqCtrl = '*'
  }else{isEqCtrl = ''}
  
  axis(1, at = tIdx, isEqCtrl, tick = F,
       las = 3, cex.axis = 1.5)
}





dev.off()



