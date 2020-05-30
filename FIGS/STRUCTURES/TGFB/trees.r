library(Cairo)
library(ape)

treefiles = dir()
treefiles = treefiles[grep('(phylip|nexus)',treefiles)]

CairoSVG( 'trees.svg', width=3*length(treefiles),height=4)
par(mfrow=c(1,length(treefiles)))
for( treefile in treefiles ){
  outname = gsub('(phylip|nexus)','svg',treefile)
  tree    = read.tree( file=treefile )
  
  plot( tree )
  
  
}
dev.off()
