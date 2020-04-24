args<-commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressMessages(library('gsubfn')))
suppressWarnings(suppressMessages(library('ape')))
suppressWarnings(suppressMessages(library('foreach')))
suppressWarnings(suppressMessages(library('doParallel')))
#library('phangorn')

#args[1] is filepath to output folder; args[2] is threads; args[3] is bootstrap; args[4]/args[5] are dendrogram/phylogeny booleans; args[6] is a '|'joined string of distance arguments; args[7] is a '|'joined string of matrix arguments

outputpath<-as.character(args[1])
threads<-as.integer(args[2])
boot<-as.integer(args[3])
#network<-as.character(args[4]) #true or false
dendrogram<-as.character(args[4]) #true or false
phylogeny<-as.character(args[5]) #true or false
distargs<-as.character(args[6])
distargs<-unique(sort(unlist(strsplit(distargs,'|',fixed=TRUE))))
matrixargs<-as.character(args[7])
matrixargs<-unique(sort(unlist(strsplit(matrixargs,'|',fixed=TRUE))))
combinedargs<-unique(sort(c(distargs,matrixargs)))


getmatrix<-function(myreport,arg,bootstraprun=FALSE) {
  NAinmatrix<-FALSE
  samples1<-as.character(sort(unique(myreport$Genome1)))
  samples<-sort(unique(c(as.vector(myreport$Genome1),as.vector(myreport$Genome2))))
  numsamples<-length(samples)
  myreportmatrix<-matrix(NA, ncol=length(samples), nrow=length(samples))
  myreportlist<-list()
  #first make list of samples1:samples2/scores
  colnamesmyreport<-colnames(myreport)
  scorecols<-c(which(colnamesmyreport=='Genome2'),which(colnamesmyreport==arg))
  for (mysample in samples1)  {
    myreportlist[[mysample]]<-myreport[myreport$Genome1==mysample,scorecols]
  }
  #make ANI or distance matrix
  for(i in 1:length(samples)) {
    sample<-samples[i]
    if (sample %in% samples1) {
      for(j in 1:length(samples)) {
        if (j<=i) { #upper-right triangle or diagonal (automatically set to NA and 0 respectively)
          next
        }
        value<-myreportlist[[sample]][which(myreportlist[[sample]][,1]==samples[j]),2]
        if (length(value)==0) {
          value<-NA
          NAinmatrix<-TRUE #if a (distance) score is missing (i.e. no blast hits for comparison), don't plot tree 
        }
        myreportmatrix[j,i]<-value #this assigns values to the lower left triangle of the matrix
        myreportmatrix[i,j]<-value #this assigns values to the upper right triangle of the matrix
      }
    }
    else {
      for(j in 1:length(samples)) {
        if (j<=i) {
          next
        }
        value<-NA
        NAinmatrix<-TRUE
        myreportmatrix[j,i]<-value #this assigns values to the lower left triangle of the matrix
        myreportmatrix[i,j]<-value #this assigns values to the upper right triangle of the matrix
      }
    }
  }
  colnames(myreportmatrix)<-samples
  rownames(myreportmatrix)<-samples
  #myreportmatrix[upper.tri(myreportmatrix, diag=F)]<-NA #assign upper triangle values (excluding diagonals) to NA #!update: now creating symmetric matrix
  if (arg!='Average_nucleotide_identity') { #assign diagnonals to distance=0 or ANI=1
    diag(myreportmatrix)<-0 
  } else {
    diag(myreportmatrix)<-1
  }
  if (bootstraprun==FALSE) {
    return(list(myreportmatrix,NAinmatrix,numsamples))
  } else { ###this is for bootstrapped trees - no need to check NAinmatrix or get numsamples, just return matrix as dist object
    d<-as.dist(myreportmatrix)
    return(d)
  }

}

getdend<-function(d) {
  mydend<-as.phylo(hclust(d)) #complete linkage used by default
  return(mydend)
}

getphy<-function(d) {
  myphylo<-fastme.bal(d, nni = TRUE, spr = TRUE, tbr = TRUE)
  myphylo$edge.length[myphylo$edge.length<0]<-0
  return(myphylo)
}
#N.B for trees: if no bootstrapping is specified, plot tree using distancestats.tsv data, otherwise use bootstrapped distance scores in order to plot distancestats.tsv tree + bootstrap confidence values


#read distancestats
myreport<-read.table(gsubfn('%1', list('%1'=outputpath),'%1/output/distancestats.tsv'), header = TRUE, sep='\t')
#myreport<-myreport[order(myreport$Genome1,myreport$Genome2),] #already ordered in granges.R
#N.B the distancestats.tsv report produced by granges.R is already ordered as follows: Genome1/Genome2 pairs are sorted alphabetically; rows are sorted according to the Genome1 and Genome 2 columns.



for (arg in combinedargs) {
  #get ANI/distance matrix; output as matrix if arg is specified in matrix args; output as tree if arg is specified in distargs (and no NAs present in matrix i.e. no missing comparison data)
  matrixout<-getmatrix(myreport,arg)
  matrixobj<-matrixout[[1]]
  NAinmatrix<-matrixout[[2]]
  numsamples<-matrixout[[3]]
  if (NAinmatrix==TRUE) {
    dendrogram<-FALSE
    phylogeny<-FALSE
  }
  if (arg %in% matrixargs) {
    writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=arg), '%1/output/matrix_%2.tsv')
    matrixdf<-as.data.frame(matrixobj)
    write.table(matrixdf,writefilepath,sep='\t',col.names = NA,row.names = TRUE,quote=FALSE)
  }

  if (dendrogram=='True' || phylogeny=='True') {
    if (arg %in% distargs) {
      distarg<-arg
      distobj<-as.dist(matrixobj)
      #set plotting dimensions based on number of samples
      height<-(numsamples%/%5)*5
      width<-(numsamples%/%5)*3
      if (height<15) {
        height<-15
      }
      if (width<10) {
        width<-10
      }
      #get original dendrogram/phylogeny
      if (dendrogram=='True') {
        originaldend<-getdend(distobj)
      }
      
      if (phylogeny=='True') {
        originalphy<-getphy(distobj)
      }
      ##plot trees...
      if (boot==0) {
        #no bootstrapping, plot trees
        if (dendrogram=='True') {
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/dend_%2.pdf')
          pdf(writefilepath,width,height)
          par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
          plot(originaldend)
          dev.off()
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/dend_%2.rds')
          saveRDS(originaldend,writefilepath)
        }
        
        if (phylogeny=='True') {
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/phylo_%2.pdf')
          pdf(writefilepath,width,height)
          par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
          plot(originalphy)
          dev.off()
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/phylo_%2.rds')
          saveRDS(originalphy,writefilepath)
        }
        
      } else {
        #get bootstrapped phylogenies
        myreportboot<-read.table(gsubfn('%1', list('%1'=outputpath),'%1/output/distancestats_bootstrapped.tsv'), header = TRUE, sep='\t')
        myreportbootsplit<-split(myreportboot, myreportboot$bootstrap)
        myreportbootsplit<-mclapply(myreportbootsplit, function(x) x<-x[order(x$Genome1,x$Genome2),], mc.cores=threads)
        
        cl<-makeCluster(as.integer(threads))
        registerDoParallel(cl)
        
        bootdists<-foreach(i=1:length(names(myreportbootsplit)), .packages=c('ape')) %dopar% {
          bootname<-names(myreportbootsplit)[i]
          distout<-getmatrix(myreportbootsplit[[bootname]],distarg,bootstraprun=TRUE)
          print(distout)
        }
        
        stopCluster(cl)
        
        names(bootdists)<-names(myreportbootsplit)
        
        if (dendrogram=='True') {
          bootdends<-lapply(bootdists, getdend)
          bootdendspart<-prop.part(bootdends,check.labels = TRUE)
          boot.clades<-prop.clades(originaldend, part=bootdendspart) #if a clade of the original tree is not represented in any of the bootstrap trees, the node support value will be 'NA'; convert to 0; express values as percentage; then re-convert 0s + convert other low values to NA
          boot.clades[is.na(boot.clades)]<-0
          boot.clades<-round(boot.clades/boot,2)*100 #express as percentage
          #boot.clades[boot.clades<50]<-NA #don't label clade support bootstrap values <50%
          
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/dend_%2_bootstrapped.pdf')
          pdf(writefilepath,width,height)
          par(mar = c(5,2,2,10)) 
          plot(originaldend)
          nodelabels(boot.clades,bg="white",frame="none",cex=0.8,adj=c(1.1,-0.4))
          dev.off()
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/dend_%2.rds')
          saveRDS(originaldend,writefilepath)
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/dend_%2_bootstrapped.rds')
          saveRDS(boot.clades,writefilepath)
        }
        
        if (phylogeny=='True') {
          bootphys<-lapply(bootdists, getphy)
          boot.clades<-prop.clades(originalphy, bootphys) #if a clade of the original tree is not represented in any of the bootstrap trees, the node support value will be 'NA'; convert to 0; express values as percentage; then re-convert 0s + convert other low values to NA
          boot.clades[is.na(boot.clades)]<-0
          boot.clades<-round(boot.clades/boot,2)*100 #express as percentage
          #boot.clades[boot.clades<50]<-NA #don't label clade support bootstrap values <50%
          
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/phylo_%2_bootstrapped.pdf')
          pdf(writefilepath,width,height)
          par(mar = c(5,2,2,10)) 
          plot(originalphy)
          nodelabels(boot.clades,bg="white",frame="none",cex=0.8,adj=c(1.1,-0.4))
          dev.off()
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/phylo_%2.rds')
          saveRDS(originalphy,writefilepath)
          writefilepath<-gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/phylo_%2_bootstrapped.rds')
          saveRDS(boot.clades,writefilepath)
        }
      }
    }
  }
}


if (dendrogram==FALSE && phylogeny==FALSE) {
  print('notree')
}




#OLD CODE - making 2d network

# for (distarg in distargs) {
#   #get original phylogeny
#   if (network=='True') {
#     phyout<-getphy(myreport,distarg,getmatrix=TRUE)
#     originalphy<-phyout[[1]]
#     originalphy$edge.length[originalphy$edge.length<0]<-0
#     mymatrix=phyout[[2]]
#     nnetout<-neighborNet(mymatrix)
#     numsamples<-phyout[[3]]
#   } else {
#     phyout<-getphy(myreport,distarg)
#     originalphy<-phyout[[1]]
#     originalphy$edge.length[originalphy$edge.length<0]<-0  
#     numsamples<-phyout[[2]]
#   }

#   height=(numsamples%/%5)*5
#   width=(numsamples%/%5)*3
#   if (height<15) {
#     height=15
#   }
#   if (width<10) {
#     width=10
#   }
  
#   if (boot==0) {
#     writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/tree_%2.pdf')
#     pdf(writefilepath,width,height)
#     par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
#     plot(originalphy)
#     dev.off()
#     writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/tree_%2.rds')
#     saveRDS(originalphy,writefilepath)

#     if (network=='True') {
#        writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/network_%2.pdf')
#        pdf(writefilepath,width,height)
#        par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
#        plot(nnetout)
#        dev.off()
#        writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/network_%2.rds')
#        saveRDS(nnetout,writefilepath)
#     }
#   } else {
#     #get bootstrapped phylogenies
#     myreportboot<-read.table(gsubfn('%1', list('%1'=outputpath),'%1/output/distancestats_bootstrapped.tsv'), header = TRUE, sep='\t')
#     myreportbootsplit<-split(myreportboot, myreportboot$bootstrap)
#     myreportbootsplit<-mclapply(myreportbootsplit, function(x) x<-x[order(x$Sample1,x$Sample2),], mc.cores=threads)
  
#     cl<-makeCluster(as.integer(threads))
#     registerDoParallel(cl)

#     bootphys<-foreach(i=1:length(names(myreportbootsplit)), .packages=c('ape')) %dopar% {
#       bootname<-names(myreportbootsplit)[i]
#       phyout<-getphy(myreportbootsplit[[bootname]],distarg,bootstraprun=TRUE)
#       print(list(phyout))
#     }

#     stopCluster(cl)
#     names(bootphys)<-names(myreportbootsplit)
#     bootphys<-lapply(bootphys, function(l) l[[1]])

#     boot.clades<-prop.clades(originalphy, bootphys) #if a clade of the original tree is not represented in any of the bootstrap trees, the node support value will be 'NA'; convert to 0; express values as percentage; then re-convert 0s + convert other low values to NA
#     boot.clades[is.na(boot.clades)]<-0
#     boot.clades<-round(boot.clades/boot,2)*100 #express as percentage
#     boot.clades[boot.clades<50]<-NA
    
#     writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/tree_%2_bootstrapped.pdf')
#     pdf(writefilepath,width,height)
#     par(mar = c(5,2,2,10)) 
#     plot(originalphy)
#     nodelabels(boot.clades,bg="white",frame="none",cex=0.8,adj=c(1.1,-0.4))
#     dev.off()
#     writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/tree_%2.rds')
#     saveRDS(originalphy,writefilepath)
#     writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/tree_%2_bootstrapped.rds')
#     saveRDS(boot.clades,writefilepath)

#     if (network=='True') {
#        writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/network_%2.pdf')
#        pdf(writefilepath,width,height)
#        par(mar = c(5,2,2,10)) #bottom left top right N.B this must come below pdf command 
#        plot(nnetout)
#        dev.off()
#        writefilepath=gsubfn('%1|%2', list('%1'=outputpath,'%2'=distarg), '%1/output/network_%2.rds')
#        saveRDS(nnetout,writefilepath)
#     }

#   }
# }
