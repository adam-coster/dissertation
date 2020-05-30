library(Cairo)

params   = read.table('data/cycleParams.csv',sep=',',header=T)
filtered = read.table('data/filteredSingleCellData.csv',sep=',',header=T)
allData  = read.table('data/features.csv',sep=',',header=T)
thumbRoot = 'W:/Publications/2014_05_Wnt_TGFB_BMP_insulation/FIGURES/InputOutput/figure_data/140130_te20_plate_2013012037/thumbs/'
outRoot   = 'data/nuclei/'
attach(params)

numCells = 25

# Collect IDs for random cells in each of the various bins

cells.G1 = allData$cellID[abs(mean.1 - log2(allData$totalIntensity.1)) < 2*sd.1]
cells.G2 = allData$cellID[abs(mean.2 - log2(allData$totalIntensity.1)) < 2*sd.1]
cells.S  = log2(allData$totalIntensity.1) < (mean.2-3*sd.1) &
           log2(allData$totalIntensity.1) > (mean.1+3*sd.1)
cells.S  = allData$cellID[cells.S]




for( type in c('G1','G2','S')){
  cells = sample(eval(parse(text=paste('cells',type,sep='.'))))
  count = 0
  while(T){
    inName = paste(thumbRoot,cells[count+1],'-1.png',sep='')
    if (file.exists( inName )){
      file.copy(from=inName,
                to=paste(outRoot,type,'-',count+1,'-1.png',sep='') )
      count = count + 1
    }
    if( count == numCells ){
      break
    }
  }
}

