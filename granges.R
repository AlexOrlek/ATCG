args = commandArgs(trailingOnly=TRUE)
library('gsubfn')
library('GenomicRanges')
library('purrr')
reduce<-GenomicRanges::reduce
simplify<-IRanges::simplify
library(data.table)
transpose<-purrr::transpose

#args[1] is filepath to pipeline output folder; args[2] is threads; args[3] is breakpoint stats; args[4] is bootstrap number

cores=as.integer(args[2])
breakpoint=as.character(args[3])
alnlenstats=as.character(args[4])
boot=as.integer(args[5])

###define functions

#stats functions
gethspidpositions<-function(x) {
  #returns the number of identical nucleotide alignment positions for a given set of alignments ('HSPs'); to avoid bias, it is applied after overlaps have been excluded from the alignment set
  return(sum(round(width(x)*(x$pid/100))))
}

gethsplength<-function(x) {
  #returns the total coverage length of a set of alignments ('HSPs'); to avoid bias, it is applied after overlaps have been excluded from the alignment set
  return(sum(width(x)))
}

getstats<-function(x) {
  #returns raw statistics                                                        
  hspidpositions<-gethspidpositions(x)
  hsplength<-gethsplength(x)
  return(c(hspidpositions,hsplength))                                                                     
}


#trimming functions
reducefunction<-function(x) {
  #reduces (joins contiguous) disjoint ranges, after splitting by input hsp   
  myinputhsp<-unique(mcols(x)$inputhsp)
  output<-reduce(x)
  mcols(output)$inputhsp<-myinputhsp
  return(output)
}

splitreducecombine<-function(x) {
  #takes list of ranges for a given paired sample; splits into needs reducing / doesn't need reducing; reduces; combines 
  isduplicate<-duplicated(mcols(x)$inputhsp) | duplicated(mcols(x)$inputhsp, fromLast=TRUE)
  numdup<-sum(isduplicate)
  if(numdup>0) {
    if(numdup==length(isduplicate)) {
      duplicates<-x
      splitduplicates<-split(duplicates,mcols(duplicates)$inputhsp)
      reducedduplicates<-lapply(splitduplicates, function(x) x=reducefunction(x))
      combinedreduced<-do.call(getMethod(c, "GenomicRanges"), reducedduplicates)
      final<-combinedreduced
    } else {
      duplicates<-x[isduplicate]
      nonduplicates<-x[!isduplicate]
      splitduplicates<-split(duplicates,mcols(duplicates)$inputhsp)
      reducedduplicates<-lapply(splitduplicates, function(x) x=reducefunction(x))
      combinedreduced<-do.call(getMethod(c, "GenomicRanges"), reducedduplicates)
      final<-append(nonduplicates, combinedreduced)
    }
  } else {
    final<-x
  }
  return(final)
}


addcols<-function(x) {
  #add pid and strand to reduced output (apply to each list element i.e. each qname split)     
  myinputhsps<-mcols(x)$inputhsp
  mcols(x)$pid<-pid[myinputhsps]
  mcols(x)$mystrand<-mystrand[myinputhsps]
  mcols(x)$snames<-snames[myinputhsps] #!!!ADDED
  return(x)
}


trimalignments<-function(qfinal,sfinal,qtrimonly=FALSE) {
  #filter alignments of qfinal based on sfinal and vice-versa
  finalhsps<-sort(intersect(mcols(qfinal)$inputhsp,mcols(sfinal)$inputhsp))
  qfinal<-qfinal[mcols(qfinal)$inputhsp %in% finalhsps]
  sfinal<-sfinal[mcols(sfinal)$inputhsp %in% finalhsps]
  finalalignments<-report[finalhsps,]
  #order alignments
  qfinal<-qfinal[order(mcols(qfinal)$inputhsp)]
  sfinal<-sfinal[order(mcols(sfinal)$inputhsp)]
  #trim QUERY alignments
  copyqfinal<-qfinal
  addstart<-start(sfinal)-finalalignments$sstart
  minusend<-finalalignments$send-end(sfinal)
  mydiff<-addstart+minusend
  for (i in 1:length(mydiff)) {
    if (mydiff[i]>0) {
      #calculate new start/end positions
      if (finalalignments$strand[i]=='+') {
        newstart<-start(qfinal[i])+addstart[i]
        newend<-end(qfinal[i])-minusend[i]
      } else {
        newstart<-start(qfinal[i])+minusend[i]
        newend<-end(qfinal[i])-addstart[i]
      }
      #apply new start/end positions in order to trim alignments
      if (newstart>=end(qfinal[i])) {
        start(qfinal[i])<-end(qfinal[i])
      } else {
        start(qfinal[i])<-newstart
      }
      if (newend<=start(qfinal[i])) {
        end(qfinal[i])<-start(qfinal[i])
      } else {
        end(qfinal[i])<-newend
      }
    }
  }
  if (qtrimonly==FALSE) {
  #trim SUBJECT alignments
  addstart<-start(copyqfinal)-finalalignments$qstart
  minusend<-finalalignments$qend-end(copyqfinal)
  mydiff<-addstart+minusend
  for (i in 1:length(mydiff)) {
    if (mydiff[i]>0) {
      #calculate new start/end positions
      if (finalalignments$strand[i]=='+') {
        newstart<-start(sfinal[i])+addstart[i]
        newend<-end(sfinal[i])-minusend[i]
      } else {
        newstart<-start(sfinal[i])+minusend[i]
        newend<-end(sfinal[i])-addstart[i]
      }
      #apply new start/end positions in order to trim alignments
      if (newstart>=end(sfinal[i])) {
        start(sfinal[i])<-end(sfinal[i])
      } else {
        start(sfinal[i])<-newstart
      }
      if (newend<=start(sfinal[i])) {
        end(sfinal[i])<-start(sfinal[i])
      } else {
        end(sfinal[i])<-newend
      }
    }
  }
  mylist<-list("qfinal"=qfinal,"sfinal"=sfinal)
  return(mylist)
  } else {
  return(qfinal)
  }
}


#breakpoint caluclation functions
makepairs <- function(x) paste(head(x, -1), tail(x, -1)) 
BP <- function(x, y) makepairs(x)[!(makepairs(x) %in% makepairs(y) | makepairs(x) %in% makepairs(rev(y*-1)))]

BPfunc<-function(x,y) {  #apply to granges object
  stopifnot(nrow(x)==nrow(y))
  alncount<-length(mcols(x)$inputhsp)
  bpcount<-length(BP(mcols(x)$inputhsp,mcols(y)$inputhsp))
  return(list("bpcount"=bpcount,"alncount"=alncount))
}

breakpointcalc<-function(qtrimmed,strimmed,mydf) {
  qcontigsplit<-lapply(qtrimmed, function(x) split(x,as.vector(mcols(x)$sname))) 
  scontigsplit<-lapply(strimmed, function(x) split(x,as.vector(mcols(x)$sname))) #nested lists- need to split by subject contigs before calculating breakpoints; convert mcols()$sname to vector otherwise fator levels are used
  bpout<-list()
  for (qname in names(qcontigsplit)) {
    #sort alignments by position
    qcontigsplit[[qname]]<-lapply(qcontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    scontigsplit[[qname]]<-lapply(scontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    out<-mapply(BPfunc,qcontigsplit[[qname]],scontigsplit[[qname]]) #apply breakpoint function to sublists
    bpcount<-sum(unlist((out["bpcount",]))) #sum across subject contig splits to get count per query
    alncount<-sum(unlist((out["alncount",])))
    bpout[[qname]]<-c(bpcount,alncount)
  }
  mydfbp<-as.data.frame(do.call(rbind, bpout)) 
  mydfbp<-cbind(querysample=rownames(mydfbp),subjectsample=rep(sample,nrow(mydfbp)),mydfbp)
  #merge dataframes
  myfinaldf<-merge(mydf,mydfbp,by=c("querysample","subjectsample"))
  return(myfinaldf)
}


#bootstrapping functions

resample<-function(x) {
  indices<-sample(1:length(x),replace=T)
  return(x[indices])
}


combinebootfunc<-function(x) {
  bootlist<-lapply(x, function(l) l[[2]])
  bootlist<-lapply(bootlist, function(x) as.data.frame(rbindlist(x,idcol="index")))
  return(as.data.frame(do.call(rbind,bootlist)))
}


#alignment length distribution function

getalnlenstats<-function(x) {
  alnlens<-rev(sort(width(x)))
  totalalnlen<-sum(alnlens)
  nxvector<-vector()
  lxvector<-vector()
  quartiles<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  for (i in 1:length(quartiles)) {
    quartile<-quartiles[i]
    lx<-which(cumsum(alnlens)>(totalalnlen*quartile))[1]
    nx<-alnlens[lx]
    lxvector[i]<-lx
    nxvector[i]<-nx
  }
  return(list(lxvector,nxvector))
}

#final ditance stats calculation function
statsfunc<-function(stats, breakpoint,mygenomelen,mymingenomelen,alnlenstats='False') {
  lxcols<-c('l10','l20','l30','l40','l50','l60','l70','l80','l90')
  nxcols<-c('n10','n20','n30','n40','n50','n60','n70','n80','n90')
  d0<-as.numeric(1-(stats["hsplength"]/mygenomelen))
  d1<-as.numeric(1-(stats["hsplength"]/mymingenomelen))
  d2<-as.numeric(-log(stats["hsplength"]/mygenomelen))
  d3<-as.numeric(-log(stats["hsplength"]/mymingenomelen))
  d4<-as.numeric(1-(stats["hspidpositions"]/stats["hsplength"]))
  d5<-as.numeric(-log(stats["hspidpositions"]/stats["hsplength"]))
  d6<-as.numeric(1-(stats["hspidpositions"]/mygenomelen))
  d7<-as.numeric(1-(stats["hspidpositions"]/mymingenomelen))
  d8<-as.numeric(-log(stats["hspidpositions"]/mygenomelen))
  d9<-as.numeric(-log(stats["hspidpositions"]/mymingenomelen))
  percentid<-as.numeric(stats["hspidpositions"]/stats["hsplength"])
  covbreadthmin<-as.numeric(stats["hsplength"]/mymingenomelen)
  if (breakpoint=='True' && alnlenstats=='True') {
    bpdist<-as.numeric(stats["breakpoints"]/stats["alignments"])
    breakpoints<-as.numeric(stats["breakpoints"])
    alignments<-as.numeric(stats["alignments"])
    lxstats<-as.integer(stats[lxcols])
    nxstats<-as.integer(stats[nxcols])
    return(c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadthmin,bpdist,breakpoints,alignments,lxstats,nxstats))
  } else if (breakpoint=='True') {
    bpdist<-as.numeric(stats["breakpoints"]/stats["alignments"])
    breakpoints<-as.numeric(stats["breakpoints"])
    alignments<-as.numeric(stats["alignments"])
    return(c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadthmin,bpdist,breakpoints,alignments))
  } else if (alnlenstats=='True') {
    lxstats<-as.integer(stats[lxcols])
    nxstats<-as.integer(stats[nxcols])
    return(c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadthmin,lxstats,nxstats))
  } else {
    return(c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadthmin))
  }
}


###iterate through samples, to get initial raw statistics for sample-pairs

#read seqlength file
seqlenreport<-read.table(gsubfn('%1',list('%1'=args[1]),'%1/seqlengths.tsv'),sep='\t',header=FALSE)
colnames(seqlenreport)<-c('sequence','length')
seqlenreport<-seqlenreport[order(seqlenreport$sequence),]

#read samples file
library('foreach')
library('doParallel')

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)

samples<-read.table(gsubfn('%1',list('%1'=args[1]),'%1/included.txt'),sep='\t',header=FALSE)
samples<-as.vector(samples[,1])
#samples<-samples[1:6]
allsampledflist<-list()

allsampledflist<-foreach(i=1:length(samples), .packages = c('gsubfn','GenomicRanges','purrr')) %dopar% {
    #read alignmnents file for given sample
    sample<-samples[i]
    report<-read.table(gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/blast/%2/alignments.tsv'),sep='\t',header=FALSE) #same subject, different queries
    colnames(report)<-c('qname','sname','pid','alnlen','mismatches','gapopens','qstart','qend','sstart','send','evalue','bitscore','qcov','qcovhsp','qlength','slength','strand')
    report$qname<-sapply(strsplit(as.vector(report$qname),"|",fixed=T),function(x) x=x[1]) #!!!ADDED
    report$sname<-sapply(strsplit(as.vector(report$sname),"|",fixed=T),function(x) x=x[1])   
    ###disjoin method trimming
    pid<-report$pid
    snames<-report$sname #!!!ADDDED
    #alnlen<-report$alnlen
    mystrand<-report$strand  #strand info is lost after disjoin - need to retain strand info for toppid
    qgr<-GRanges(seqnames = report$qname, ranges = IRanges(start=(report$qstart), end = (report$qend)), strand=(report$strand))
    sgr<-GRanges(seqnames = report$qname, ranges = IRanges(start=(report$sstart), end= (report$send)), strand = (report$strand))
    qalnlen<-width(qgr)
    salnlen<-width(sgr)
    #shift subject ranges where there are multiple contigs
    samplecontigs<-levels(as.factor(report$sname))
    samplecontiglens<-seqlenreport[seqlenreport$sequence %in% samplecontigs,2]
    sampleindices<-as.numeric(as.factor(report$sname))
    for (j in 1:length(samplecontiglens)) {
      if (j==1) {
        next
      }
      myindices<-which(sampleindices==j)
      addlen<-sum(samplecontiglens[c(1:(j-1))])
      sgr[myindices]<-shift(sgr[myindices],addlen)
    }
    #query disjoin
    gr2<-disjoin(qgr,with.revmap=TRUE,ignore.strand=TRUE)
    revmap<-gr2$revmap
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(qalnlen[x])]))
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample                                            
    qreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x))
    #subject disjoin
    gr2<-disjoin(sgr,with.revmap=TRUE,ignore.strand=TRUE)
    revmap<-gr2$revmap
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(salnlen[x])]))
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample
    sreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x))
    #add pid and strand
    qfinal<-lapply(qreducedoutput, function(x) x=addcols(x))
    sfinal<-lapply(sreducedoutput, function(x) x=addcols(x))
    #trim query alignments
    if (breakpoint=='True') {
      trimmedalignments<-transpose(mapply(trimalignments,qfinal,sfinal,SIMPLIFY = F)) #!ADDED
      qtrimmed<-trimmedalignments$qfinal
      strimmed<-trimmedalignments$sfinal
    } else {
      qtrimmed<-mapply(trimalignments,qfinal,sfinal,qtrimonly=TRUE)
    }
    #trimmedalignments<-mapply(trimalignments,qfinal,sfinal)
    #trimmedalignments<-lapply(split(1:nrow(trimmedalignments), rownames(trimmedalignments)), function(i) trimmedalignments[i,])  #!!!ADDED
    #get hsp id stats
    mystats<-lapply(qtrimmed,getstats) #!!!CHANGED
    mydf<-as.data.frame(do.call(rbind, mystats)) #convert list of vectors to dataframe                                               
    mydf<-cbind(querysample=rownames(mydf),subjectsample=rep(sample,nrow(mydf)),mydf)
    #get breakpoint stats
    if (breakpoint=='True') {
      myfinaldf<-breakpointcalc(qtrimmed,strimmed,mydf)
    } else {
      myfinaldf<-mydf
    }
    #get alignment length distribution stats
    if (alnlenstats=='True') {
       alnlenstatslist<-lapply(qtrimmed, getalnlenstats)
       alnlenstatsdf<-cbind(as.data.frame(do.call(rbind, lapply(alnlenstatslist, function(l) l[[1]]))), as.data.frame(do.call(rbind, lapply(alnlenstatslist, function(l) l[[2]]))))
       myfinaldf<-cbind(myfinaldf,alnlenstatsdf)
    }
    #IF NO BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS
    if (boot==0) {
      print(myfinaldf)
    } else {
      #IF BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS + BOOTSTRAPPED STATS
      myfinaldfbootlist<-list()
      for (i in 1:boot) {
        qtrimmedboot<-lapply(qtrimmed,resample)
        strimmedboot<-lapply(strimmed,resample)
        mystatsboot<-lapply(qtrimmedboot,getstats) #!!!CHANGED
        mydfboot<-as.data.frame(do.call(rbind, mystatsboot)) #convert list of vectors to dataframe                                               
        mydfboot<-cbind(querysample=rownames(mydfboot),subjectsample=rep(sample,nrow(mydfboot)),mydfboot)
        #get breakpoint stats
        if (breakpoint=='True') {
          myfinaldfboot<-breakpointcalc(qtrimmedboot,strimmedboot,mydfboot)
        } else {
          myfinaldfboot<-mydfboot
        }
        myfinaldfbootlist[[i]]<-myfinaldfboot
      }
      print(list(myfinaldf, myfinaldfbootlist))
    }
}

stopCluster(cl)


lxcols<-c('l10','l20','l30','l40','l50','l60','l70','l80','l90')
nxcols<-c('n10','n20','n30','n40','n50','n60','n70','n80','n90')
if (boot==0) {
  allsampledf<-as.data.frame(do.call(rbind, allsampledflist))
  if (breakpoint=='True' && alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments',lxcols, nxcols)
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments",lxcols,nxcols)
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments')
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments")
  } else if (alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength',lxcols,nxcols)
    statscols<-c("hspidpositions","hsplength",lxcols,nxcols)
  } else {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength')
    statscols<-c("hspidpositions","hsplength")
  }
} else {
  allsampledf<-as.data.frame(do.call(rbind, lapply(allsampledflist, function(l) l[[1]])))
  allsampledfboot<-combinebootfunc(allsampledflist)
  if (breakpoint=='True' && alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments',lxcols,nxcols)
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments",lxcols,nxcols)
    statscolsboot<-c("hspidpositions","hsplength","breakpoints","alignments")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments')
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments')
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments")
    statscolsboot<-c("hspidpositions","hsplength","breakpoints","alignments")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments')
  } else if (alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength',lxcols,nxcols)
    statscols<-c("hspidpositions","hsplength",lxcols,nxcols)
    statscolsboot<-c("hspidpositions","hsplength")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  } else {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength')
    statscols<-c("hspidpositions","hsplength")
    statscolsboot<-c("hspidpositions","hsplength")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  }
}

###get final stats (distance scores) for each pairwise sample combination

#aggregate seqlen repot at sample level
seqlenreport$sequence<-sapply(strsplit(as.vector(seqlenreport$sequence),"|",fixed=T),function(x) x=x[1])
sampleseqlen<-aggregate(seqlenreport$length, by=list(seqlenreport$sequence), FUN=sum)
colnames(sampleseqlen)<-c('sample','length')

#get final stats
samples2<-samples
allsampledflist<-list()
counter1=0
for (sample in samples) {
  print(sample)
  sampleaseqlen<-sampleseqlen[which(sampleseqlen[,"sample"]==sample),"length"]
  counter1=counter1+1
  counter2=0
  output<-list()
  for (sampleb in samples2) {
    if (sample==sampleb) {
      #counter2=counter2+1 #if you want self-self distances to be included, unhash these two lines
      #output[[counter2]]<-c(sample,sampleb,0,0,0)
      next
    }
    counter2=counter2+1
    #print(c(sample, sampleb))
    stats1row<-which(allsampledf[,"querysample"]==sample & allsampledf[,"subjectsample"]==sampleb)
    stats2row<-which(allsampledf[,"querysample"]==sampleb & allsampledf[,"subjectsample"]==sample)
    if (length(stats1row)==0 && length(stats2row)==0) {
      print(c(sampleb,'no pairwise matches found for this sample'))
      next
    }
    samplebseqlen<-sampleseqlen[which(sampleseqlen[,"sample"]==sampleb),"length"]
    mygenomelen<-sampleaseqlen+samplebseqlen
    mymingenomelen<-min(sampleaseqlen,samplebseqlen)*2
    if (length(stats1row)==0) {
      stats<-data.frame(rbind(colSums(allsampledf[stats2row,statscols])),row.names = NULL)
      if (stats["hsplength"]>mymingenomelen) { #this shouldn't be necessary given that only one genome has alignments
        stats["hsplength"]<-mymingenomelen
      }
      if (stats["hspidpositions"]>mymingenomelen) {
        stats["hspidpositions"]<-mymingenomelen
      }
      statsout<-statsfunc(stats,breakpoint,mygenomelen,mymingenomelen,alnlenstats)

    } else if (length(stats2row)==0) {
      stats<-data.frame(rbind(colSums(allsampledf[stats1row,statscols])),row.names = NULL)
      if (stats["hsplength"]>mymingenomelen) { #this shouldn't be necessary given that only one genome has alignments
        stats["hsplength"]<-mymingenomelen
      }
      if (stats["hspidpositions"]>mymingenomelen) {
        stats["hspidpositions"]<-mymingenomelen
      }
      statsout<-statsfunc(stats,breakpoint,mygenomelen,mymingenomelen,alnlenstats)

    } else { ##combine both
      stats1<-data.frame(rbind(colSums(allsampledf[stats1row,statscols])),row.names = NULL)
      stats2<-data.frame(rbind(colSums(allsampledf[stats2row,statscols])),row.names = NULL)
      mergedstats<-stats1+stats2
      if (mergedstats["hsplength"]>mymingenomelen) {
        mergedstats["hsplength"]<-mymingenomelen
      }
      if (mergedstats["hspidpositions"]>mymingenomelen) {
        mergedstats["hspidpositions"]<-mymingenomelen
      }
      statsout<-statsfunc(mergedstats,breakpoint,mygenomelen,mymingenomelen,alnlenstats)

    }
    output[[counter2]]<-c(sample,sampleb,sampleaseqlen,samplebseqlen,statsout)

  }
  finaldf<-as.data.frame(do.call(rbind,output))
  allsampledflist[[counter1]]<-finaldf
  samples2<-setdiff(samples2,sample) #this means non-symmetric matrix of distance scores - this is what I want because by using merged stats, I'm averaging over both directions
}

finaldf<-as.data.frame(do.call(rbind,allsampledflist))

if (breakpoint=='True' && alnlenstats=='True') {
  colnames(finaldf)<-c('Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome','Breakpoint_distance','Breakpoints','Alignments',lxcols,nxcols)
} else if (breakpoint=='True') {
  colnames(finaldf)<-c('Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome','Breakpoint_distance','Breakpoints','Alignments')
} else if (alnlenstats=='True') {
  colnames(finaldf)<-c('Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome',lxcols,nxcols)
} else {
  colnames(finaldf)<-c('Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome')
}

write.table(finaldf, file=gsubfn('%1',list('%1'=args[1]),'%1/output/distancestats.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)


###if there is bootstrapping, need to produce additional stats table - split by bootstrap, apply above code to get stats to each boostrap, then combine

allsampledfbootsplit<-split(allsampledfboot, allsampledfboot$bootstrap)
finaldflist<-list()
statscols<-statscolsboot
alnlenstats='False'

for (i in names(allsampledfbootsplit)){
  print(i)
  allsampledf<-allsampledfbootsplit[[i]]
  samples2<-samples
  allsampledflist<-list()
  counter1=0
  for (sample in samples) {
    print(c(sample,'sample'))
    sampleaseqlen<-sampleseqlen[which(sampleseqlen[,"sample"]==sample),"length"]
    counter1=counter1+1
    counter2=0
    output<-list()
    for (sampleb in samples2) {
      if (sample==sampleb) {
        #counter2=counter2+1 #if you want self-self distances to be included, unhash these two lines
        #output[[counter2]]<-c(sample,sampleb,0,0,0)
        next
      }
      counter2=counter2+1
      stats1row<-which(allsampledf[,"querysample"]==sample & allsampledf[,"subjectsample"]==sampleb)
      stats2row<-which(allsampledf[,"querysample"]==sampleb & allsampledf[,"subjectsample"]==sample)
      if (length(stats1row)==0 && length(stats2row)==0) {
        print(c(sampleb,'no pairwise matches found for this sample'))
        next
      }
      samplebseqlen<-sampleseqlen[which(sampleseqlen[,"sample"]==sampleb),"length"]
      mygenomelen<-sampleaseqlen+samplebseqlen
      mymingenomelen<-min(sampleaseqlen,samplebseqlen)*2
      if (length(stats1row)==0) {
        stats<-data.frame(rbind(colSums(allsampledf[stats2row,statscols])),row.names = NULL)
        if (stats["hsplength"]>mymingenomelen) { #this shouldn't be necessary given that only one genome has alignments
          stats["hsplength"]<-mymingenomelen
        }
        if (stats["hspidpositions"]>mymingenomelen) {
          stats["hspidpositions"]<-mymingenomelen
        }
        statsout<-statsfunc(stats,breakpoint,mygenomelen,mymingenomelen)
        
      } else if (length(stats2row)==0) {
        stats<-data.frame(rbind(colSums(allsampledf[stats1row,statscols])),row.names = NULL)
        if (stats["hsplength"]>mymingenomelen) { #this shouldn't be necessary given that only one genome has alignments
          stats["hsplength"]<-mymingenomelen
        }
        if (stats["hspidpositions"]>mymingenomelen) {
          stats["hspidpositions"]<-mymingenomelen
        }
        statsout<-statsfunc(stats,breakpoint,mygenomelen,mymingenomelen)
        
      } else { ##combine both
        stats1<-data.frame(rbind(colSums(allsampledf[stats1row,statscols])),row.names = NULL)
        stats2<-data.frame(rbind(colSums(allsampledf[stats2row,statscols])),row.names = NULL)
        mergedstats<-stats1+stats2
        if (mergedstats["hsplength"]>mymingenomelen) {
          mergedstats["hsplength"]<-mymingenomelen
        }
        if (mergedstats["hspidpositions"]>mymingenomelen) {
          mergedstats["hspidpositions"]<-mymingenomelen
        }
        statsout<-statsfunc(mergedstats,breakpoint,mygenomelen,mymingenomelen)
        
      }
      output[[counter2]]<-c(sample,sampleb,sampleaseqlen,samplebseqlen,statsout)
      
    }
    finaldf<-as.data.frame(do.call(rbind,output))
    allsampledflist[[counter1]]<-finaldf
    samples2<-setdiff(samples2,sample) #this means non-symmetric matrix of distance scores - this is what I want because by using merged stats, I'm averaging over both directions
  }
  finaldflist[[i]]<-allsampledflist
}



finaldflist2<-lapply(finaldflist, function(x) as.data.frame(do.call(rbind,x)))
finaldfboot<-rbindlist(finaldflist2,idcol = "index")

if (breakpoint=='True' && alnlenstats=='True') {
  colnames(finaldfboot)<-c('bootstrap','Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome','Breakpoint_distance','Breakpoints','Alignments',lxcols,nxcols)
} else if (breakpoint=='True') {
  colnames(finaldfboot)<-c('bootstrap','Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome','Breakpoint_distance','Breakpoints','Alignments')
} else if (alnlenstats=='True') {
  colnames(finaldfboot)<-c('bootstrap','Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome',lxcols,nxcols)  
} else {
  colnames(finaldfboot)<-c('bootstrap','Sample1','Sample2','Sample1_length','Sample2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth_mingenome')
}

write.table(finaldfboot, file=gsubfn('%1',list('%1'=args[1]),'%1/output/distancestats_bootstrapped.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

