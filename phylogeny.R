args = commandArgs(trailingOnly=TRUE)
library('gsubfn')
library('ape')
library('foreach')
library('doParallel')

#args[1] is filepath to output folder; args[2] is threads; args[3] is bootstrap; subsequent argument(s) are distance score columns e.g. DistanceScore_d9

threads=as.integer(args[2])
boot=as.integer(args[3])
distargs=args[4:length(args)]


getphy<-function(myreport,distarg,getnumsamples=TRUE) {
  samples1<-as.character(sort(unique(myreport$Sample1)))
  samples<-sort(unique(c(as.vector(myreport$Sample1),as.vector(myreport$Sample2))))
  numsamples<-length(samples)
  myreportmatrix<-matrix(NA, ncol=length(samples), nrow=length(samples))
  myreportlist<-list()
  
  #first make list of samples1:samples2/scores
  colnamesmyreport<-colnames(myreport)
  scorecols<-c(which(colnamesmyreport=='Sample2'),which(colnamesmyreport==distarg))
  for (mysample in samples1)  {
    myreportlist[[mysample]]<-myreport[myreport$Sample1==mysample,scorecols]
  }
  
  #make matrix
  maxdist<-2*(max(myreport[,distarg]))
  
  for(i in 1:length(samples)) {
    sample<-samples[i]
    if (sample %in% samples1) {
      for(j in 1:length(samples)) {
        if (j<=i) { #upper-right triangle or diagonal (automatically set to NA and 0 respectively)
          next
        }
        value<-myreportlist[[sample]][which(myreportlist[[sample]][,1]==samples[j]),2]
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
  
  #make phylogeny
  myphy<-fastme.bal(d, nni = TRUE, spr = TRUE, tbr = TRUE)
  if (getnumsamples==TRUE) {
    return(list(myphy,numsamples))
  } else {
    return(myphy)
  }
}


#if no bootstrapping is specified, plot tree using distancestats.tsv data, otherwise use bootstrapped distance scores in order to plot distancestats.tsv tree + bootstrap confidence values

#read distancestats and make sure data is ordered data by sample name
myreport<-read.table(gsubfn('%1', list('%1'=args[1]),'%1/output/distancestats.tsv'), header = TRUE, sep='\t')
myreport<-myreport[order(myreport$Sample1,myreport$Sample2),]


for (distarg in distargs) {
  #get original phylogeny
  phyout<-getphy(myreport,distarg)
  originalphy<-phyout[[1]]
  originalphy$edge.length[originalphy$edge.length<0]<-0  
  numsamples<-phyout[[2]]
  height=(numsamples%/%5)*5
  width=(numsamples%/%5)*3
  if (height<15) {
    height=15
  }
  if (width<10) {
    width=10
  }
  
  if (boot==0) {
    writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/tree_%2.pdf')
    pdf(writefilepath,width,height)
    par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
    plot(originalphy)
    dev.off()
    writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/tree_%2.rds')
    saveRDS(originalphy,writefilepath)
  } else {
    #get bootstrapped phylogenies
    myreportboot<-read.table(gsubfn('%1', list('%1'=args[1]),'%1/output/distancestats_bootstrapped.tsv'), header = TRUE, sep='\t')
    myreportbootsplit<-split(myreportboot, myreportboot$bootstrap)
    myreportbootsplit<-mclapply(myreportbootsplit, function(x) x<-x[order(x$Sample1,x$Sample2),], mc.cores=threads)
  
    cl<-makeCluster(as.integer(threads))
    registerDoParallel(cl)

    bootphys<-foreach(i=1:length(names(myreportbootsplit)), .packages=c('ape')) %dopar% {
      bootname<-names(myreportbootsplit)[i]
      phyout<-getphy(myreportbootsplit[[bootname]],distarg,getnumsamples=FALSE)
      print(list(phyout))
    }

    stopCluster(cl)
    names(bootphys)<-names(myreportbootsplit)
    bootphys<-lapply(bootphys, function(l) l[[1]])

    boot.clades<-prop.clades(originalphy, bootphys) #if a clade of the original tree is not represented in any of the bootstrap trees, the node support value will be 'NA'; convert to 0; express values as percentage; then re-convert 0s + convert other low values to NA
    boot.clades[is.na(boot.clades)]<-0
    boot.clades<-round(boot.clades/boot,2)*100 #express as percentage
    boot.clades[boot.clades<50]<-NA
    
    writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/tree_%2_bootstrapped.pdf')
    pdf(writefilepath,width,height)
    par(mar = c(5,2,2,10)) 
    plot(originalphy)
    nodelabels(boot.clades,bg="white",frame="none",cex=0.8,adj=c(1.1,-0.4))
    dev.off()
    writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/tree_%2.rds')
    saveRDS(originalphy,writefilepath)
    writefilepath=gsubfn('%1|%2', list('%1'=args[1],'%2'=distarg), '%1/output/tree_%2_bootstrapped.rds')
    saveRDS(boot.clades,writefilepath)
  }
}
