args = commandArgs(trailingOnly=TRUE)
library('gsubfn')

#args[1] is filepath to output folder; subsequent argument(s) are distance score columns e.g. DistanceScore_d9

myreport<-read.table(gsubfn('%1', list('%1'=args[1]),'%1/output/distancestats.tsv'), header = TRUE, sep='\t')

#order data by sample name
myreport<-myreport[order(myreport$Sample1,myreport$Sample2),]

#edit sample names
#myreport$Sample1<-sapply(as.list(as.vector(myreport$Sample1)), function(x) substr(x, 1, nchar(x)-37))
#myreport$Sample2<-sapply(as.list(as.vector(myreport$Sample2)), function(x) substr(x, 1, nchar(x)-37)) #OLD:function(x) paste(substr(x, 1, nchar(x)-39),substr(x,nchar(x),nchar(x)),sep='_')

distargs=args[2:length(args)]

for (distarg in distargs) {
  ###make distance matrix
  samples1<-as.character(sort(unique(myreport$Sample1)))
  samples<-sort(unique(c(as.character(unique(myreport$Sample1)),as.character(unique(myreport$Sample2)))))
  myreportmatrix<-matrix(NA, ncol=length(samples), nrow=length(samples))
  myreportlist<-list()

  #first make list of samples1:samples2/scores
  colnamesmyreport<-colnames(myreport)
  scorecols<-c(which(colnamesmyreport=='Sample2'),which(colnamesmyreport==distarg))
  for (i in 1:length(samples1))  {
    myreportlist[[i]]<-myreport[myreport$Sample1==samples1[i],scorecols]
  }
  names(myreportlist)<-samples1


  #make matrix

  #in M-k 2014 when there are no hsps, distance is set as twice maximum observed distance (rather than just setting to 1)
  #head(rev(sort(myreport$DistanceScore_d6)),100); many distances are close to 1, so using twice max helps to distinguish samples with no similarity to others, from low similarity samples

  #maxdist<-2*(max(myreport$DistanceScore_d9))
  maxdist<-2*(max(myreport[,distarg]))

  for(i in 1:length(samples)) {
    sample<-samples[i]
    if (sample %in% samples1) {
      for(j in 1:length(samples)) {
	if (j<=i) { #upper-right triangle or diagonal (automatically set to NA and 0 respectively)
	  next
	}
	#print(c(i,j))
	value<-myreportlist[[sample]][which(myreportlist[[sample]][,1]==samples[j]),2]
	#print(value)
	if (length(value)==0) {
	  value=maxdist #set distance to maximum (no blast hits) for given pairwise sample comparison
	}
	myreportmatrix[j,i]<-value #this assigns values to the lower triangle of the matrix
      }
    }
    else {
      for(j in 1:length(samples)) {
	if (j<=i) {
	  next
	}
	#print(c(i,j))
	value=maxdist #set distance to maximum (no blast hits)
	myreportmatrix[j,i]<-value #this assigns values to the lower triangle of the matrix
      }
    }
  }


  colnames(myreportmatrix)<-samples
  rownames(myreportmatrix)<-samples
  myreportmatrix[upper.tri(myreportmatrix, diag=F)]<-NA #assign upper triangle values (excluding diagonals) to NA
  diag(myreportmatrix)<-0 #assign diagnonals to 0 distance
  d<-as.dist(myreportmatrix) #convert distance matrix to dist structure that can be used as input for hclust


  #make ggplot dendrogram / save dendrogram to file
  clust = hclust(d) #complete linkage used by default; not using FastME which is used by M-k.  
  clustdendro<-as.dendrogram(clust)
  dendfilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/dendrogram_%2.rd')
  saveRDS(clustdendro,dendfilepath)

  #plot

  numsamples<-length(samples1)
  print(numsamples)
  print('number of samples')
  height=(numsamples%/%5)*5
  width=(numsamples%/%5)*3
  if (height<15) {
     height=15
  }
  if (width<10) {
     width=10
  }

  writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/dendrogram_%2.pdf')
  pdf(writefilepath,width,height)
  par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
  plot(clustdendro,horiz=T)
  dev.off()
}