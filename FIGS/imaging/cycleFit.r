library(Cairo)

params = read.table('data/cycleParams.csv',sep=',',header=T)
data   = read.table('data/log2DNA.csv',sep=',',header=T)

attach(params)


CairoSVG('cycleFit.svg',width=6,height=3)

par(mfrow=c(1,2))


s2p   = 1/sqrt(2*pi)
minFreq = .01

xlims = c(12.8,15)
ylims = c(0,2.5)

x=data$x
y=data$y

# SHOW EACH PIECE
plot( x, y, type='l',xlim=xlims, ylim=ylims,lwd=3,
      xlab='log2(DNA)',ylab='nuclei', col='gray')

G1.y = weight.1*s2p*(1/sd.1)*exp(-(x-mean.1)^2/(2*sd.1^2))
lines( x[G1.y>minFreq], G1.y[G1.y>minFreq], col='red',lty=3,lwd=2)

G2.y = weight.2*s2p*(1/sd.2)*exp(-(x-mean.2)^2/(2*sd.1^2))
lines( x[G2.y>minFreq], G2.y[G2.y>minFreq], col='dodgerblue',lty=3,lwd=2)

S.y  = (x>mean.1 & x<(mean.2))*v
lines( x[S.y>minFreq], S.y[S.y>minFreq], col='black',lty=3,lwd=2)

abline( v=mean.1, col='red' )
abline( v=mean.2, col='dodgerblue')


# COMBINED
plot( x, y, type='l', xlim=xlims, ylim=ylims,lwd=3,
      xlab='log2(DNA)',ylab='nuclei', col='gray')
all.y = G1.y + G2.y + S.y
lines( x[all.y>minFreq],all.y[all.y>minFreq], col='darkred',lty=3,lwd=2)

dev.off()