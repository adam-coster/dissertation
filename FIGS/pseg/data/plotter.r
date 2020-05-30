# Goal: obtain plot of CFP vs RFP for mixed infected pSeg population
library(Hmisc)
library(Cairo)

dna  = read.table('dapi.csv',header=T,sep=',',stringsAsFactors=F)
cfp  = read.table('cfp.csv',header=T,sep=',',stringsAsFactors=F)
rfp  = read.table('rfp.csv',header=T,sep=',',stringsAsFactors=F)

# Quick QC of cell cycle

hist(log2(dna$sum_nucleus),100     )
maxDNA = 20
minDNA = 19.75

keep   = log2(dna$sum_nucleus) < maxDNA &
         log2(dna$sum_nucleus) > minDNA

dna  = log2(dna[keep,'sum_nucleus'])
cfp  = log2(cfp[keep,'sum_nucleus'])
rfp  = log2(rfp[keep,'sum_nucleus'])

# Correct using DNA
# RFP first
plot( NA, NA, xlim=c(19.75,20),
      ylim=c(13.5,22))

rk = kmeans( rfp, 2 )
for( m in 1:2 ){
  indices = rk$cluster == m
  points( dna[indices],rfp[indices],
          pch=16,
          col=rgb(c(0,1)[m],0,c(1,0)[m],.5))
  fit = lm(rfp[indices]~dna[indices])
  abline(fit, col='black', lwd=3)

  legend(c('topleft','bottomleft')[m],
         paste(round(rcorr(rfp[indices],dna[indices],type='pearson')$r[2],2)),
         text.col=rgb(c(0,1)[m],0,c(1,0)[m]),bty='n')
  
  rfp[indices] = median(rfp[indices])+fit$residuals
}

# CFP next
plot( NA, NA, xlim=c(19.75,20),
      ylim=c(16.8,20.5))

ck = kmeans( cfp, 2 )
for( m in 1:2 ){
  indices = ck$cluster == m
  points( dna[indices],cfp[indices],
          pch=16,
          col=rgb(c(0,1)[m],0,c(1,0)[m],.5))
  fit = lm(cfp[indices]~dna[indices])
  abline(fit, col='black', lwd=3)
  
  legend(c('topleft','bottomleft')[m],
         paste(round(rcorr(cfp[indices],dna[indices],type='pearson')$r[2],2)),
         text.col=rgb(c(0,1)[m],0,c(1,0)[m]),bty='n')
  
  cfp[indices] = median(cfp[indices])+fit$residuals
}


#######


#

ak = kmeans(cbind(cfp,rfp),4,iter.max=100,nstart=10)


CairoSVG( 'mixed.svg', width=4,height=4)


plot(NA,NA, xlim=c(min(cfp),max(cfp)),
     ylim=c(min(rfp),max(rfp)),
     xlab='CFP',ylab='RFP')

corrs  = c()
cols   = c()
counts = c()
for( m in 1:4 ){
  indices = ak$cluster == m
  col = rgb(c(0,1,0,0)[m],c(0,0,.5,0)[m],c(1,0,1,0)[m],.25)
  points( cfp[indices],rfp[indices],
          pch=16,
          col=col)
  cols = c(cols,col)
  corrs = c(corrs,
            round(rcorr(cfp[indices],rfp[indices],'pearson')$r[2],2))
  counts = c(counts,sum(indices))
}

legend( 'bottomright',
        paste(round(counts/sum(counts)*100),'%',sep=''),
        text.col=cols,bty='n')

dev.off()
