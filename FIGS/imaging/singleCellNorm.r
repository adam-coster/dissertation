library(Cairo)
library(Hmisc)

params   = read.table('data/cycleParams.csv',sep=',',header=T)
filtered = read.table('data/filteredSingleCellData.csv',sep=',',header=T)

attach(params)

subset = filtered$well == 'D6'
histBins = c(-Inf, seq( from=14, to=16.5, length.out=20), Inf)

CairoSVG('singleCellNorm.svg',height=3,width=12)
par(mfrow=c(1,4))


ylims = c(14,17)
xlims = c(13.2,13.7)


plot( log2(filtered$totalIntensity.1[subset]),
      log2(filtered$totalIntensity.2[subset]),
      pch=16, col=rgb(0,0,0,.2), xlim=xlims, ylim=ylims )
fit = lm( log2(filtered$totalIntensity.2[subset]) ~
          log2(filtered$totalIntensity.1[subset]) )
abline(fit,col='gray',lwd=3)
legend('topleft', paste('r =',round( 
    rcorr( log2(filtered$totalIntensity.2[subset]),
           log2(filtered$totalIntensity.1[subset]),
          type='pearson')$r[2],2)),text.col='red',bty='n')


plot( log2(filtered$totalIntensity.1[subset]),
      log2(filtered$denoisedTotalIntensityMult.2[subset]),
      pch=16, col=rgb(0,0,0,.2),main='multiplicative', xlim=xlims, ylim=ylims )
fit = lm( log2(filtered$denoisedTotalIntensityMult.2[subset]) ~
            log2(filtered$totalIntensity.1[subset]) )
abline(fit,col='red',lwd=2)
legend('topleft', paste('r =',round( 
  rcorr( log2(filtered$denoisedTotalIntensityMult.2[subset]),
         log2(filtered$totalIntensity.1[subset]),
         type='pearson')$r[2],2)),text.col='red',bty='n')


plot( log2(filtered$totalIntensity.1[subset]),
      log2(filtered$denoisedTotalIntensity.2[subset]),
      pch=16, col=rgb(0,0,0,.2),main='regression', xlim=xlims, ylim=ylims )
fit = lm( log2(filtered$denoisedTotalIntensity.2[subset]) ~
            log2(filtered$totalIntensity.1[subset]) )
abline(fit,col='dodgerblue',lwd=2)
legend('topleft', paste('r =',round( 
  rcorr( log2(filtered$denoisedTotalIntensity.2[subset]),
         log2(filtered$totalIntensity.1[subset]),
         type='pearson')$r[2],2)),text.col='red',bty='n')


plot( log2(filtered$totalIntensity.2[subset]),
      log2(filtered$denoisedTotalIntensity.2[subset]),
      pch=16, main='shuffling',col= rgb(0,0,0,.2),cex=.5,
      xlim=ylims, ylim=ylims )

dev.off()