library(Cairo)

G1Cells = function( wellPath, G1, DNAChannel=1, zCutoff=2.5 ){
  DNA    = read.table( paste( wellPath, '-', DNAChannel, '.csv',sep='' ), stringsAsFactors=F, header=T, sep=',')
  cvs    = DNA$std_nucleus / DNA$mean_nucleus
  keep   = cvs < G1$CV.max & cvs > G1$CV.min
  
  logDNA = log( DNA$sum_nucleus )
  keep   = keep & logDNA < G1$mean.1 + zCutoff*G1$sd.1 & logDNA > G1$mean.1 - zCutoff*G1$sd.1
  
  return( keep )
}


CairoSVG('readoutInformation',width=10,height=5)
par(mfrow=c(2,4))

#### Plot distributions from TI08 to show MI of TGFB3 or BMP4
#### and their readouts.

##### The following file will only exist after this script has been
##### run, followed by the mutual_information_readouts.m script
mutualInfoFile  = 'data/readoutInformation.csv'
#####

data        = read.table('W:/2013_06_wnt_tgfb_crosstalk/TI08/TI08_summaryStats.csv',header=T,sep=',',stringsAsFactors=F)
segRoot     = 'W:/2013_06_wnt_tgfb_crosstalk/plates/140404_TI08/segmentation/NUCLEI_DATA/'
G1          = read.table("W:/2013_06_wnt_tgfb_crosstalk/plates/140404_TI08/segmentation/SKMEL2_G1.csv",header=T,sep=',',stringsAsFactors=F)
experiments = c('TGFB.readouts','BMP.readouts')
incols      = c('blue','darkgreen')
names(incols) = experiments
readouts    = matrix(c('pSmad23','Smad23','pSmad158','Smad1'),nrow=2)
colnames(readouts) = experiments
features    = c('median_nucleus','sum_nucleus')
channel     = '2'
bins        = 15

for( feature in features){
  for( experiment in experiments ){
    for( readout in readouts[,experiment]){
      subdata    = data[data$experiment==experiment &
                        data[,paste('marker',channel,sep='.')]==readout,]
      treatments = unique(subdata$treatment)
      wells = c()
      xlims = c(Inf,-Inf)
      for( treatment in treatments ){
        # Find well with median well-level value
        wells  = c(wells,subdata$well[subdata$treatment==treatment])
      }
      wells = matrix(wells,nrow=3)
      colnames(wells)=treatments
      
      hists = c()
      xlims = c(Inf,-Inf)
      for( treatment in treatments ){
        allVals = c()
        for( well in wells[,treatment]){
          wellpath    = paste(segRoot,well,sep='')
          singleCells = read.table( paste( wellpath,'-',channel,'.csv',sep='' ),
                                    sep=',', header=T, stringsAsFactors=F )
          dna         = read.table( paste( wellpath,'-1.csv',sep='' ),
                                    sep=',', header=T, stringsAsFactors=F )
          keepers     = G1Cells(wellpath,G1)
          singleCells = singleCells[keepers,paste(feature)]
          dna         = dna[keepers, paste(feature)]
          # Normalize!
          fit         = lm( log2(singleCells) ~ log2(dna) )
          singleCells = 2^(median(log2(singleCells)) + fit$residuals)
          allVals     = c(allVals,singleCells)
          
        }
        allVals = log2(allVals)
        # Trim the outliers
        allVals = allVals[allVals<(median(allVals)+3*mad(allVals)) &
                          allVals>(median(allVals)-3*mad(allVals))]
        # Write to file for idiot-Matlab
        write.csv(allVals,paste('data/',readout,treatment,feature,'csv',sep='.'),
                  row.names=F)
        hists[[treatment]] = hist(allVals,breaks=bins,plot=F)
        xlims = c(min(xlims[1],hists[[treatment]]$mids),
                  max(xlims[2],hists[[treatment]]$mids))
      }
      
      # Set up a plot
      plot(NA,NA,ylim=c(0,1),xlim=2^xlims, main=feature,
           xlab=readout,
           ylab='Frequency', yaxt='n')
      cols = c('black',incols[experiment])
      names(cols) = treatments
      for( treatment in treatments ){
        lines( 2^hists[[treatment]]$mids,
               hists[[treatment]]$density/max(hists[[treatment]]$density),
               col=cols[treatment],lwd=2)
      }
      legend('topright',treatments,text.col=cols[treatments],bty='n')
      if( file.exists(mutualInfoFile) ){
        mutualInfo = read.table(mutualInfoFile,header=T,sep=',')
        miMean     = mutualInfo[paste(readout,feature,'meanMI',sep='.')]
        miSd       = mutualInfo[paste(readout,feature,'sdMI',sep='.')]
        legend('right',c(paste(round(miMean,2),'+/-',
                               round(miSd,  2)),
                         paste('n=(',sum(hists[[1]]$counts),',',
                               sum(hists[[2]]$counts),')',sep='')),
               bty='n')}
    }
  }
}
dev.off()


