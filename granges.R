args = commandArgs(trailingOnly=TRUE)
library('gsubfn')
library('GenomicRanges')
library('purrr')
reduce<-GenomicRanges::reduce
shift<-GenomicRanges::shift
library(data.table)
rbindlist<-data.table::rbindlist
transpose<-purrr::transpose

#args[1] is filepath to pipeline output folder; args[2] is threads; args[3] is breakpoint stats; args[4] is bootstrap number

cores=as.integer(args[2])
breakpoint=as.character(args[3])
alnlenstats=as.character(args[4])
boot=as.integer(args[5])

###define functions

#stats functions
gethspidpositions<-function(x) {
  #returns the number of identical nucleotide alignment positions for a given set of alignments ('HSPs'); it is applied after overlaps have been excluded from the alignment set
  return(sum(round(width(x)*(x$pid/100))))
}

gethsplength<-function(x) {
  #returns the total coverage length of a set of alignments ('HSPs'); it is applied after overlaps have been excluded from the alignment set
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
  #add pid, strand, and subject names to reduced output (apply to each list element i.e. each qname split)     
  myinputhsps<-mcols(x)$inputhsp
  mcols(x)$pid<-pid[myinputhsps]
  mcols(x)$mystrand<-mystrand[myinputhsps]
  mcols(x)$snames<-snames[myinputhsps] #!!!ADDED
  return(x)
}


pastefunction<-function(x) {
  if (length(x)==1) {
    x=x[1]
  } else {
    x=paste(x[1],x[2],sep='|')
  }
}



canceldoubletrim<-function(x,y) {
  #x is query, y is subject; return x,y after cancelling double trim
  xcancel<-y-x
  ycancel<-x-y
  xcancel[xcancel<0]<-0
  ycancel[ycancel<0]<-0
  return(list(xcancel,ycancel))
}

canceldoubletrimwrapper<-function(qaddstart,qminusend,saddstart,sminusend,alldiffindices,posindices,negindices) {
  cancelout<-canceldoubletrim(qaddstart[alldiffindices][posindices],saddstart[alldiffindices][posindices])
  qaddstart[alldiffindices][posindices]<-cancelout[[1]]
  saddstart[alldiffindices][posindices]<-cancelout[[2]]
  
  cancelout<-canceldoubletrim(qminusend[alldiffindices][posindices],sminusend[alldiffindices][posindices])
  qminusend[alldiffindices][posindices]<-cancelout[[1]]
  sminusend[alldiffindices][posindices]<-cancelout[[2]]
  
  cancelout<-canceldoubletrim(qaddstart[alldiffindices][negindices],sminusend[alldiffindices][negindices])
  qaddstart[alldiffindices][negindices]<-cancelout[[1]]
  sminusend[alldiffindices][negindices]<-cancelout[[2]]
  
  cancelout<-canceldoubletrim(qminusend[alldiffindices][negindices],saddstart[alldiffindices][negindices])
  qminusend[alldiffindices][negindices]<-cancelout[[1]]
  saddstart[alldiffindices][negindices]<-cancelout[[2]]
  
  return(list(qaddstart,qminusend,saddstart,sminusend))
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
  #get addstart/minusend vectors for query and subject
  qaddstart<-start(sfinal)-finalalignments$sstart #based on subject trimming, how much should be trimmed from query 
  qminusend<-finalalignments$send-end(sfinal)
  saddstart<-start(copyqfinal)-finalalignments$qstart #based on query trimming, how much should be trimmed from subject - need to know this at query trimming stage to avoid double trimming
  sminusend<-finalalignments$qend-end(copyqfinal)
  #get indices that require position reassignment
  qmydiffindices<-which((qaddstart > 0 | qminusend > 0)==TRUE)  
  smydiffindices<-which((saddstart > 0 | sminusend > 0)==TRUE)
  alldiffindices<-intersect(qmydiffindices,smydiffindices) #indices where there is trimming on both query and subject - this is used to prevent double trimming
  posindices<-which(mcols(qfinal[alldiffindices])$mystrand=='+') #of the alldiffindices, which are positive strand
  negindices<-which(mcols(qfinal[alldiffindices])$mystrand=='-')
  cancelout<-canceldoubletrimwrapper(qaddstart,qminusend,saddstart,sminusend,alldiffindices,posindices,negindices)
  qaddstart<-cancelout[[1]];qminusend<-cancelout[[2]];saddstart<-cancelout[[3]];sminusend<-cancelout[[4]]
  #calculate new start/end positions
  lenmydiff<-length(qmydiffindices)
  newstart<-numeric(lenmydiff)
  newend<-numeric(lenmydiff)
  strandpos<-which(finalalignments$strand[qmydiffindices]=='+') #of the mydiffindices, which are positive strand; strandpos is for subsetting mydiffindices/indexing newstart/end vectors; mydiffindices is for subsetting alignments
  lenstrandpos<-length(strandpos)
  if (lenstrandpos==0) { #all negative
    newstart<-start(qfinal[qmydiffindices])+qminusend[qmydiffindices]
    newend<-end(qfinal[qmydiffindices])-qaddstart[qmydiffindices]
  } else if (lenstrandpos==lenmydiff) { #all positive
    newstart<-start(qfinal[qmydiffindices])+qaddstart[qmydiffindices]
    newend<-end(qfinal[qmydiffindices])-qminusend[qmydiffindices]
  } else {
    mydiffindicespos<-qmydiffindices[strandpos]
    mydiffindicesneg<-qmydiffindices[-strandpos] #error - now corrected !  instead of -
    newstart[strandpos]<-start(qfinal[mydiffindicespos])+qaddstart[mydiffindicespos]
    newend[strandpos]<-end(qfinal[mydiffindicespos])-qminusend[mydiffindicespos]
    newstart[-strandpos]<-start(qfinal[mydiffindicesneg])+qminusend[mydiffindicesneg]
    newend[-strandpos]<-end(qfinal[mydiffindicesneg])-qaddstart[mydiffindicesneg]
  }
  #set bounds for new start/ends
  boundnewstartindices<-which(newstart>end(qfinal[qmydiffindices]))
  if (length(boundnewstartindices)>0) {
    newstart[boundnewstartindices]<-end(qfinal[qmydiffindices][boundnewstartindices])
  }
  start(qfinal)[qmydiffindices]<-newstart #reassign
  boundnewendindices<-which(newend<start(qfinal[qmydiffindices]))
  if (length(boundnewendindices)>0) {
    newend[boundnewendindices]<-start(qfinal[qmydiffindices][boundnewendindices])
  }
  end(qfinal)[qmydiffindices]<-newend #reassign
  
  if (qtrimonly==FALSE) {
    #trim SUBJECT alignments
    #calculate new start/end positions
    lenmydiff<-length(smydiffindices)
    newstart<-numeric(lenmydiff)
    newend<-numeric(lenmydiff)
    strandpos<-which(finalalignments$strand[smydiffindices]=='+') #strandpos is for subsetting mydiffindices/indexing newstart/end vectors; mydiffindices is for subsetting alignments
    lenstrandpos<-length(strandpos)
    if (lenstrandpos==0) { #all negative
      newstart<-start(sfinal[smydiffindices])+sminusend[smydiffindices]
      newend<-end(sfinal[smydiffindices])-saddstart[smydiffindices]
    } else if (lenstrandpos==lenmydiff) { #all positive
      newstart<-start(sfinal[smydiffindices])+saddstart[smydiffindices]
      newend<-end(sfinal[smydiffindices])-sminusend[smydiffindices]
    } else {
      mydiffindicespos<-smydiffindices[strandpos]
      mydiffindicesneg<-smydiffindices[-strandpos] #error - now corrected !  instead of -
      newstart[strandpos]<-start(sfinal[mydiffindicespos])+saddstart[mydiffindicespos]
      newend[strandpos]<-end(sfinal[mydiffindicespos])-sminusend[mydiffindicespos]
      newstart[-strandpos]<-start(sfinal[mydiffindicesneg])+sminusend[mydiffindicesneg]
      newend[-strandpos]<-end(sfinal[mydiffindicesneg])-saddstart[mydiffindicesneg]
    }
    #set bounds for new start/ends
    boundnewstartindices<-which(newstart>end(sfinal[smydiffindices]))
    if (length(boundnewstartindices)>0) {
      newstart[boundnewstartindices]<-end(sfinal[smydiffindices][boundnewstartindices])
    }
    start(sfinal)[smydiffindices]<-newstart #reassign
    boundnewendindices<-which(newend<start(sfinal[smydiffindices]))
    if (length(boundnewendindices)>0) {
      newend[boundnewendindices]<-start(sfinal[smydiffindices][boundnewendindices])
    }
    end(sfinal)[smydiffindices]<-newend #reassign
    
    mylist<-list("qfinal"=qfinal,"sfinal"=sfinal)
    return(mylist)
  } else {
    return(qfinal)
  }
}



#breakpoint caluclation functions
makepairs<-function(x) mapply(c, head(x,-1), tail(x,-1), SIMPLIFY = FALSE)
BP<-function(x,y) { #this function works on signed permuations (numeric vectors with +/- indicated)
  out1<-makepairs(x)[!(makepairs(x) %in% makepairs(y) | makepairs(x) %in% makepairs(rev(y*-1)))]
  out2<-makepairs(x)[unlist(lapply(makepairs(x), function(z) length(unique(sign(z)))))>1]
  all<-c(out1,out2)
  return(all[!duplicated(all)]) #deduplicate to aovid double counting breakpoints that occur in out1 and out2
}

BPfunc<-function(x,y) {  #apply to granges objects
  stopifnot(nrow(x)==nrow(y))
  alncount<-length(mcols(x)$inputhsp)
  xinputhsps<-mcols(x)$inputhsp
  yinputhsps<-mcols(y)$inputhsp
  xstrand<-as.vector(mcols(x)$mystrand)
  ystrand<-as.vector(mcols(y)$mystrand)
  xstrand[xstrand=='+']<-1  #convert inputhsp vectors to signed permuations by multipling by +1/-1 strand vector
  xstrand[xstrand=='-']<--1
  ystrand[ystrand=='+']<-1
  ystrand[ystrand=='-']<--1
  xinputhsps<-xinputhsps*as.numeric(xstrand)
  yinputhsps<-yinputhsps*as.numeric(ystrand)
  bpcount<-length(BP(xinputhsps,yinputhsps))
  return(list("bpcount"=bpcount,"alncount"=alncount))
}


breakpointcalc<-function(qtrimmed,strimmed,mydf) {
  qcontigsplit<-lapply(qtrimmed, function(x) split(x,bpsplit[mcols(x)$inputhsp])) 
  scontigsplit<-lapply(strimmed, function(x) split(x,bpsplit[mcols(x)$inputhsp])) #nested lists- need to split by query and subject contigs before calculating breakpoints
  bpout<-list()
  for (qname in names(qcontigsplit)) {
    #sort alignments by position
    qcontigsplit[[qname]]<-lapply(qcontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    scontigsplit[[qname]]<-lapply(scontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    numsplits<-length(qcontigsplit[[qname]])
    out<-mapply(BPfunc,qcontigsplit[[qname]],scontigsplit[[qname]]) #apply breakpoint function to sublists
    bpcount<-sum(unlist((out["bpcount",]))) #sum across subject contig splits to get count per query
    alncount<-sum(unlist((out["alncount",])))
    bpout[[qname]]<-c(bpcount,alncount,numsplits)
  }
  mydfbp<-as.data.frame(do.call(rbind, bpout))
  colnames(mydfbp)<-c('breakpoints','alignments','numsplits')
  mydfbp<-cbind(querysample=rownames(mydfbp),subjectsample=rep(sample,nrow(mydfbp)),mydfbp)
  #merge dataframes
  myfinaldf<-merge(mydf,mydfbp,by=c("querysample","subjectsample"))
  return(myfinaldf)
}


bpdistcalc<-function(bps,alns,numsplits) {
  if ((alns-numsplits)==0) {
    return(as.numeric(0))
  } else {
    return(as.numeric(bps/(alns-numsplits)))
  }
}


#bootstrapping functions

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

#final ditance stats calculation functions
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
    bpdist<-bpdistcalc(stats["breakpoints"],stats["alignments"],stats["numsplits"])
    breakpoints<-as.numeric(stats["breakpoints"])
    alignments<-as.numeric(stats["alignments"])
    lxstats<-as.integer(stats[lxcols])
    nxstats<-as.integer(stats[nxcols])
    return(c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadthmin,bpdist,breakpoints,alignments,lxstats,nxstats))
  } else if (breakpoint=='True') {
    bpdist<-bpdistcalc(stats["breakpoints"],stats["alignments"],stats["numsplits"])
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


mergestats<-function(x,y) {
  merge<-data.frame()
  for (i in 1:length(colnames(x))) {
    mycol<-colnames(x)[i]
    if (mycol=="hspidpositions" || mycol=="hsplength") {
      merge[1,i]<-as.integer(x[i])+as.integer(y[i])
    } else {
      merge[1,i]<-mean(as.integer(x[i]),as.integer(y[i]))
    }
  }
  colnames(merge)<-colnames(x)
  return(merge)
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

lxcols<-c('l10','l20','l30','l40','l50','l60','l70','l80','l90')
nxcols<-c('n10','n20','n30','n40','n50','n60','n70','n80','n90')

allsampledflist<-foreach(i=1:length(samples), .packages = c('gsubfn','GenomicRanges','purrr')) %dopar% {
    #read alignmnents file for given sample
    sample<-samples[i]
    report<-read.table(gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/blast/%2/alignments.tsv'),sep='\t',header=FALSE) #same subject, different queries
    colnames(report)<-c('qname','sname','pid','alnlen','mismatches','gapopens','qstart','qend','sstart','send','evalue','bitscore','qcov','qcovhsp','qlength','slength','strand')
    #get information for shifting subject ranges where there are multiple contigs (below)
    seqlenreportseqs<-sapply(strsplit(as.vector(seqlenreport$sequence),'|',fixed=T),pastefunction)
    reformattedsname<-as.factor(sapply(strsplit(as.vector(report$sname),"|",fixed=T),pastefunction))
    samplescontigs<-levels(reformattedsname)
    samplescontiglens<-seqlenreport[seqlenreportseqs %in% samplescontigs,2]
    samplesindices<-as.numeric(reformattedsname)
    #get information for shifting query ranges where there are multiple contigs (below)
    reformattedqname<-as.factor(sapply(strsplit(as.vector(report$qname),"|",fixed=T),pastefunction))
    sampleqcontigs<-levels(reformattedqname)
    sampleqcontiglens<-seqlenreport[seqlenreportseqs %in% sampleqcontigs,2]
    sampleqindices<-as.numeric(reformattedqname)
    #
    if (breakpoint=='True') {
        bpsplit<-paste(reformattedqname,reformattedsname,sep='_') #required for breakpoint calculation
    }
    #remove any contig information
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
    for (j in 1:length(samplescontiglens)) {
      if (j==1) {
        next
      }
      myindices<-which(samplesindices==j)
      addlen<-sum(samplescontiglens[c(1:(j-1))])
      sgr[myindices]<-shift(sgr[myindices],addlen)
    }
    #shift query ranges where there are multiple contigs
    for (j in 1:length(sampleqcontiglens)) {
      if (j==1) {
        next
      }
      myindices<-which(sampleqindices==j)
      addlen<-sum(sampleqcontiglens[c(1:(j-1))])
      qgr[myindices]<-shift(qgr[myindices],addlen)
    }
    #overwrite report with shifted forward ranges
    if (length(samplescontiglens)>1) {
       report$sstart<-start(sgr)
       report$send<-end(sgr)
    }
    if (length(sampleqcontiglens)>1) {
       report$qstart<-start(qgr)
       report$qend<-end(qgr)
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
    #add pid, strand, and subject name
    qfinal<-lapply(qreducedoutput, function(x) x=addcols(x))
    sfinal<-lapply(sreducedoutput, function(x) x=addcols(x))
    #trim alignments
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
    colnames(mydf)<-c('hspidpositions','hsplength')
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
       alnlenstatsdf<-cbind(rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[1]]),as.data.frame.list),idcol="querysample"),rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[2]]),as.data.frame.list)))
       colnames(alnlenstatsdf)<-c('querysample',lxcols,nxcols)
       myfinaldf<-merge(myfinaldf,alnlenstatsdf,by="querysample")
    }
    #IF NO BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS
    if (boot==0) {
      print(myfinaldf)
    } else {
      #IF BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS + BOOTSTRAPPED STATS
      myfinaldfbootlist<-list()
      for (z in 1:boot) {
        #stopifnot(sapply(qtrimmed, function(x) length(x))==sapply(strimmed, function(x) length(x)))
	indices<-lapply(qtrimmed, function(x) sample(1:length(x), replace=T)) #get indices for resampling qtrimmed/strimmed
	qtrimmedboot<-lapply(1:length(qtrimmed), FUN=function(x, list1, list2) list1[[x]][list2[[x]]] , list1=qtrimmed, list2=indices) #resample qtrimmed using indices
	names(qtrimmedboot)<-names(qtrimmed)
        mystatsboot<-lapply(qtrimmedboot,getstats) #!!!CHANGED
	mydfboot<-rbindlist(lapply(mystatsboot,as.data.frame.list),idcol="querysample") ##convert list of vectors to dataframe
        mydfboot<-cbind(subjectsample=rep(sample,nrow(mydfboot)),mydfboot)
	colnames(mydfboot)<-c('subjectsample','querysample','hspidpositions','hsplength')
        #get breakpoint stats
        #if (breakpoint=='True') {
	#  strimmedboot<-lapply(1:length(strimmed), FUN=function(x, list1, list2) list1[[x]][list2[[x]]] , list1=strimmed, list2=indices) #resample strimmed using indices
	#  names(strimmedboot)<-names(strimmed)
        #  myfinaldfboot<-breakpointcalc(qtrimmedboot,strimmedboot,mydfboot)
        #} else {
        #  myfinaldfboot<-mydfboot
        #}
	myfinaldfboot<-mydfboot #NO LONGER RESAMPLING BREAKPOINT STATS
        myfinaldfbootlist[[z]]<-myfinaldfboot
      }
      print(list(myfinaldf, myfinaldfbootlist))
    }
}

stopCluster(cl)


if (boot==0) {
  allsampledf<-as.data.frame(do.call(rbind, allsampledflist))
  if (breakpoint=='True' && alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits',lxcols, nxcols)
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits",lxcols,nxcols)
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits')
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits")
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
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits',lxcols,nxcols)
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits",lxcols,nxcols)
    #statscolsboot<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits") #NO LONGER BOOTSTRAPPING BREAKPOINTS
    #colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits')
    statscolsboot<-c("hspidpositions","hsplength")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits')
    statscols<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits")
    #statscolsboot<-c("hspidpositions","hsplength","breakpoints","alignments","numsplits")
    #colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','numsplits')
    statscolsboot<-c("hspidpositions","hsplength")
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
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
  #print(sample)
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
      #print(c(sampleb,'no pairwise matches found for this sample'))
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
      mergedstats<-mergestats(stats1,stats2)
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

if (boot!=0) {

  allsampledfbootsplit<-split(allsampledfboot, allsampledfboot$bootstrap)
  finaldflist<-list()
  statscols<-statscolsboot
  alnlenstats='False'
  breakpoint='False' #NO LONGER DOING BOOTSTRAPPING OF BREAKPOINTS

  for (i in names(allsampledfbootsplit)){
    #print(i)
    allsampledf<-allsampledfbootsplit[[i]]
    samples2<-samples
    allsampledflist<-list()
    counter1=0
    for (sample in samples) {
      #print(c(sample,'sample'))
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
	  #print(c(sampleb,'no pairwise matches found for this sample'))
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
	  mergedstats<-mergestats(stats1,stats2)
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

}