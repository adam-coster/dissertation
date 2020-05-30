## Need to generate a heatmap that has all 4 celltypes
## at the 1, 2,and 18hr marks. Need to generate one
## heatmap per timepoint and readout. The readouts are
## Smad23, and CTNNB1.

timepoints = c('1hr','2hr','18hr')
datafiles  = c('data/TI04.60_summaryStats.csv',
               'data/TI04.120_summaryStats.csv',
               'data/TI16.18_summaryStats.csv')
reps       = 1
readouts   = c('Smad23','CTNNB1')
colors     = matrix(c(rgb(0,(0:10)*.1,0),
                      rgb((0:10)*.1,0,0)),ncol=2)
colnames(colors)=readouts
channels   = 2:3
celltypes  = c('HCEC','SKMEL2','MALME3M','IEC6')
treatments = c('ctrl+ctrl','TGFB3+ctrl',
               'ctrl+Wnt3A','TGFB3+Wnt3A')
setToOne   = c('TGFB3+ctrl','ctrl+Wnt3A')
feature    = 'mean_nucleus'
minHeat    = -.2
maxHeat    = 1.5

for( rIdx in 1:length(readouts) ){
  heatmap = matrix(0,nrow=reps*length(treatments),
                   ncol=length(timepoints)*length(celltypes))
  
  for( tIdx in 1:length(timepoints)){
    data = read.table(datafiles[tIdx],sep=',',header=T,
                      stringsAsFactors=F )
    for( cIdx in 1:length(celltypes) ){
      readout  = readouts[ rIdx]
      channel  = channels[ rIdx]
      marker   = paste('marker',channel,sep='.')
      celltype = celltypes[cIdx]
      
      subdata  = data[data$celltype==celltype &
                      data$marker.2==readouts[1]  &
                      data$marker.3==readouts[2]  &
                      data$treatment %in% treatments,]
      
      # Convert to ctrl-based Z-scores
      thisFeature = paste(feature,'median',channel,sep='.')
      ctrls   = subdata[subdata$treatment == 'ctrl+ctrl',
                        thisFeature,]
      subdata[,thisFeature] = subdata[,thisFeature]-mean(ctrls)
      subdata[,thisFeature] = subdata[,thisFeature]/
        mean(subdata[subdata$treatment==setToOne[rIdx],thisFeature])
      
      # Get it into the correct matrix coords
      for( treatIdx in 1:length(treatments)){
        treatment = treatments[treatIdx]
        rows      = (treatIdx-1)*reps+1:reps
        cols      = (tIdx-1)*length(celltypes)+ cIdx
        heatmap[rows,cols] = mean(subdata[subdata$treatment==treatment,thisFeature])
      
      }
    }
  }
  tempMap = heatmap
  tempMap[tempMap < minHeat] = minHeat
  tempMap[tempMap > maxHeat] = maxHeat
  heatmap(tempMap,Rowv=NA,Colv=NA,col=colors[,readout],scale='none',
          labRow = rep(treatments,each=reps),
          labCol = paste(rep(celltypes,length(timepoints)),
                   rep(timepoints,each=length(celltypes))))
}