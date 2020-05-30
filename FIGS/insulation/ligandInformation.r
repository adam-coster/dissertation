library(Cairo)

G1Cells = function( wellPath, G1, DNAChannel=1, zCutoff=2.5 ){
  DNA    = read.table( paste( wellPath, '-', DNAChannel, '.csv',sep='' ), stringsAsFactors=F, header=T, sep=',')
  cvs    = DNA$std_nucleus / DNA$mean_nucleus
  keep   = cvs < G1$CV.max & cvs > G1$CV.min
  
  logDNA = log( DNA$sum_nucleus )
  keep   = keep & logDNA < G1$mean.1 + zCutoff*G1$sd.1 & logDNA > G1$mean.1 - zCutoff*G1$sd.1
  
  return( keep )
}


CairoSVG('ligandInformation',width=5,height=5)
par(mfrow=c(2,2))

#### Plot distributions from TI01/3 to show MI of TGFB3 or BMP4


##### The following file will only exist after this script has been
##### run, followed by the mutual_information_readouts.m script
mutualInfoFile  = 'data/ligandInformation.csv'
#####

data        = read.table('W:/2013_06_wnt_tgfb_crosstalk/TI03/TI03_summaryStats.csv',header=T,sep=',',stringsAsFactors=F)

segRoot     = paste( 'W:/2013_06_wnt_tgfb_crosstalk/plates/',
                     '20140317_TI03',
                    '/segmentation/NUCLEI_DATA/',sep='')
G1          = read.table("W:/2013_06_wnt_tgfb_crosstalk/plates/20140317_TI03/segmentation/SKMEL2_G1.csv",header=T,sep=',',stringsAsFactors=F)
incol       = 'blue'
readout     = 'Smad23'
features    = c('sum_nucleus')
channel     = '2'
bins        = 20
concentrations = c(0,12.15)

for( feature in features){
  subdata    = data[data$celltype=='SKMEL2' &
                      data[,paste('marker',channel,sep='.')]==readout,]
  wells = c()
  xlims = c(Inf,-Inf)
  treatment = unique(subdata$treatment)
  for( concentration in concentrations ){
    # Find well with median well-level value
    wells  = c(wells,subdata$well[subdata$concentration==concentration])
  }
  wells    = matrix(wells,nrow=3)
  colnames(wells) = as.character(concentrations)
  
  hists = c()
  xlims = c(Inf,-Inf)
  for( concentration in concentrations ){
    allVals = c()
    for( well in wells[,as.character(concentration)]){
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
    write.csv(allVals,paste('data/',treatment,concentration,feature,'csv',sep='.'),
              row.names=F)
    hists[[as.character(concentration)]] = hist(allVals,breaks=bins,plot=F)
    xlims = c(min(xlims[1],hists[[as.character(concentration)]]$mids),
              max(xlims[2],hists[[as.character(concentration)]]$mids))
  }
  
  # Set up a plot
  plot(NA,NA,ylim=c(0,1),xlim=2^xlims, main=feature,
       xlab=readout,
       ylab='Frequency', yaxt='n')
  cols = c('black',incol)
  names(cols) = as.character(concentrations)
  for( concentration in concentrations ){
    lines( 2^hists[[as.character(concentration)]]$mids,
           hists[[as.character(concentration)]]$density/max(hists[[as.character(concentration)]]$density),
           col=cols[as.character(concentration)],lwd=2)
  }
  legend('topright',treatment,text.col=cols[as.character(concentration)],bty='n')
  if( file.exists(mutualInfoFile) ){
    mutualInfo = read.table(mutualInfoFile,header=T,sep=',')
    miMean     = mutualInfo[paste(treatment,feature,'meanMI',sep='.')]
    miSd       = mutualInfo[paste(treatment,feature,'sdMI',sep='.')]
    legend('right',c(paste(round(miMean,2),'+/-',
                           round(miSd,  2)),
                     paste('n=(',sum(hists[[1]]$counts),',',
                           sum(hists[[2]]$counts),')',sep='')),
           bty='n')
  }
  
}




data   = read.table('W:/2013_06_wnt_tgfb_crosstalk/TI01/TI01.1_summaryStats.csv',header=T,sep=',',stringsAsFactors=F)
data   = rbind(data,
               read.table('W:/2013_06_wnt_tgfb_crosstalk/TI01/TI01.2_summaryStats.csv',header=T,sep=',',stringsAsFactors=F))

segRoots    = paste( 'W:/2013_06_wnt_tgfb_crosstalk/plates/20140317_TI01.',
                     1:2,
                     '/segmentation/NUCLEI_DATA/',sep='')

G1          = read.table("W:/2013_06_wnt_tgfb_crosstalk/plates/20140317_TI01.1/segmentation/SKMEL2_G1.csv",
                         header=T,sep=',',stringsAsFactors=F)
readouts    = c('Smad23','pSmad158','CTNNB1')
incols      = c('blue','darkgreen','darkred')
names(incols) = readouts
channels    = c('2','2','3')
names(channels) = readouts
bins        = 20
allConcentrations = matrix(c(0,12.15,
                          0,8.1,
                          0,729 ),nrow=2)
colnames( allConcentrations ) = readouts


for( readout in readouts ){
  for( feature in features){
    channel    = channels[readout]
    subdata    = data[data$celltype=='SKMEL2' &
                        data[,paste('marker',channel,sep='.')]==readout,]
    treatment  = unique(subdata$treatment)
    wells = c()
    xlims = c(Inf,-Inf)
    concentrations = allConcentrations[,readout]
    names(segRoots) = as.character(concentrations)
    
    for( concentration in concentrations ){
      # Find well with median well-level value
      wells  = c(wells,subdata$well[subdata$concentration==concentration])
    }
    wells    = matrix(wells,nrow=3)
    colnames(wells) = as.character(concentrations)
    
    hists = c()
    xlims = c(Inf,-Inf)
    for( concentration in concentrations ){
      allVals = c()
      for( well in wells[,as.character(concentration)]){
        wellpath    = paste(segRoots[as.character(concentration)],well,sep='')
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
      write.csv(allVals,paste('data/',treatment,concentration,feature,'csv',sep='.'),
                row.names=F)
      hists[[as.character(concentration)]] = hist(allVals,breaks=bins,plot=F)
      xlims = c(min(xlims[1],hists[[as.character(concentration)]]$mids),
                max(xlims[2],hists[[as.character(concentration)]]$mids))
    }
    
    # Set up a plot
    plot(NA,NA,ylim=c(0,1),xlim=2^xlims, main=feature,
         xlab=readout,
         ylab='Frequency', yaxt='n')
    cols = c('black',incols[readout])
    names(cols) = as.character(concentrations)
    for( concentration in concentrations ){
      lines( 2^hists[[as.character(concentration)]]$mids,
             hists[[as.character(concentration)]]$density/max(hists[[as.character(concentration)]]$density),
             col=cols[as.character(concentration)],lwd=2)
    }
    legend('topright',treatment,text.col=cols[as.character(concentration)],bty='n')
    if( file.exists(mutualInfoFile) ){
      mutualInfo = read.table(mutualInfoFile,header=T,sep=',')
      miMean     = mutualInfo[paste(treatment,feature,'meanMI',sep='.')]
      miSd       = mutualInfo[paste(treatment,feature,'sdMI',sep='.')]
      legend('right',c(paste(round(miMean,2),'+/-',
                             round(miSd,  2)),
                       paste('n=(',sum(hists[[1]]$counts),',',
                             sum(hists[[2]]$counts),')',sep='')),
             bty='n')
    }
  }
}



dev.off()


