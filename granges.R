args = commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressMessages(library('gsubfn',quietly=TRUE)))
suppressWarnings(suppressMessages(library('GenomicRanges',quietly=TRUE)))  #gives masking messages even when warn.conflict=FALSE
suppressWarnings(suppressMessages(library('purrr',quietly=TRUE)))
reduce<-GenomicRanges::reduce
shift<-GenomicRanges::shift
suppressWarnings(suppressMessages(library('data.table',quietly=TRUE)))
rbindlist<-data.table::rbindlist
transpose<-purrr::transpose
map2<-purrr::map2

#args[1] is filepath to pipeline output folder; args[2] is threads; args[3] is breakpoint stats; args[4] is alignment length stats; arg[5] is bootstrap number; arg[6] is best alignment selection criterion; arg[7] is length threshold for filtering alignments prior to breakpoint calculation; arg[8] is trimmed alignment output

cores=as.integer(args[2])
breakpoint=as.character(args[3])
alnlenstats=as.character(args[4])
boot=as.integer(args[5])
alnrankmethod=as.character(args[6])
lengthfilter=as.integer(args[7])
outputbestblastalignments=as.character(args[8])
outputnonoverlappingalignments=as.character(args[9])
outputtrimmedalignments=as.character(args[10])
bidirectionalblast=as.character(args[11])
statsfromtrimmed=as.character(args[12]) #False by default
keepbisectedrangesarg=as.character(args[13])
alnlenstatsquantiles=as.character(args[14])
alnlenstatsquantiles<-unique(sort(unlist(strsplit(alnlenstatsquantiles,'|',fixed=TRUE))))

lxcols<-vector()
nxcols<-vector()
for (i in 1:length(alnlenstatsquantiles)) {
  lxcols[i]<-paste('L',alnlenstatsquantiles[i],sep='')
  nxcols[i]<-paste('N',alnlenstatsquantiles[i],sep='')
}

###define functions (N.B see bottom of script for deprecated functions)

#stats functions
gethspidpositions<-function(x) {
  #returns the number of identical nucleotide alignment positions for a given set of alignments ('HSPs'); it is applied after overlaps have been excluded from the alignment set
  #return(sum(round(width(x)*(x$pid/100))))
  return(round(sum(width(x)*(x$pid/100))))
}

gethsplength<-function(x) {
  #returns the total coverage length of a set of alignments ('HSPs'); it is applied after overlaps have been excluded from the alignment set
  return(sum(width(x)))
}

getstatsfromdisjoint<-function(x,y) {
  #returns raw statistics from disjoint alignments (qtrimmed,strimmed)
  widthx<-width(x)
  widthy<-width(y)
  mywidths<-ifelse(widthx<=widthy,widthx,widthy) #for each hsp, get width of shortest range from qtrimmed/strimmed
  mypids<-x$pid/100
  hspidpositions<-round(sum(mywidths*mypids))
  hsplength<-sum(mywidths)
  return(c(hspidpositions,hsplength))
}

getstatsfromtrimmed<-function(x) {
  #returns raw statistics from trimmed alignments (for simplicity, using qtrimmed only)                 
  hspidpositions<-gethspidpositions(x)
  hsplength<-gethsplength(x)
  return(c(hspidpositions,hsplength))                                                                     
}

#range shifting function
rangeshifting<-function(reformattedname, seqlenreport, seqlenreportseqs) { #reformatted name is seqnames in sample|contig format; seqlenreport is original seqlen report; seqlenreportseqs are the sequences in sample|contig format from the seqlenreport
  genomes<-as.vector(sapply(strsplit(as.vector(reformattedname),'|',fixed=T), function(x) x[1]))
  genomes<-factor(genomes,levels=unique(genomes)) #query genomes not necessarily ordered alphabeticlally so need to specify factor level ordering to match actual order before using factor levels to split
  genomesplit<-split(reformattedname,genomes)
  #get contiglens and indices
  samplecontiglenslist<-list()
  sampleindiceslist<-list()
  for (a in 1:length(genomesplit)) {
    samplename<-as.factor(as.vector(genomesplit[[a]]))  #!need to convert to vector first
    samplecontigs<-levels(samplename) #levels (corresponding to contigs of a genome) are ordered alphabetically; rangeshifting will follow this alphabetical order through the contigs
    sampleindices<-as.numeric(samplename)
    samplecontiglens<-numeric(length(samplecontigs))
    for (b in 1:length(samplecontigs)) {
      mycontig<-samplecontigs[b]
      contiglen<-as.vector(unlist(seqlenreport[which(seqlenreportseqs==mycontig),2]))
      samplecontiglens[b]<-contiglen
    }
    samplecontiglenslist[[a]]<-samplecontiglens
    sampleindiceslist[[a]]<-sampleindices
  }
  #get vector of lengths to add from nested per genome lists of contig lengths and indices 
  alladdlens<-numeric()
  for (x in 1:length(samplecontiglenslist)) {
    samplecontiglens<-samplecontiglenslist[[x]]
    sampleindices<-sampleindiceslist[[x]]
    addlens<-numeric(length(sampleindices))
    for (y in 1:length(samplecontiglens)) {
      if (y==1) {
        addlen<-0
        myindices<-which(sampleindices==y)
        addlens[myindices]<-addlen
        next
      }  
      myindices<-which(sampleindices==y)
      #addlen<-sum(samplecontiglens[c(1:(y-1))])+1 #!NO this is incorrect
      addlen<-sum(samplecontiglens[c(1:(y-1))])
      addlens[myindices]<-addlen
    }
    alladdlens<-c(alladdlens,addlens)  #append addlens from each genome
  }
  return(alladdlens)
}

#function which filters original blast report to include only best alignments (retained in qfinal and/or sfinal after selecting best hsps from disjoins)
getbestblasthits<-function(qfinal,sfinal) {
  besthitindices<-unique(qfinal$inputhsp,sfinal$inputhsp)
  besthitreport<-originalreport[besthitindices,]
  return(besthitreport)
}

#function which filters disjoint alignments to include only those present in both qfinal and sfinal; then sorts alignments by inputhsp
intersectsortalignments<-function(qfinal,sfinal) {
  #filter qfinal/sfinal alignments based on intersect of qfinal and sfinal hsps
  finalhsps<-sort(intersect(mcols(qfinal)$inputhsp,mcols(sfinal)$inputhsp))
  qfinal<-qfinal[mcols(qfinal)$inputhsp %in% finalhsps]
  sfinal<-sfinal[mcols(sfinal)$inputhsp %in% finalhsps]
  finalalignments<-report[finalhsps,]
  #sort alignments
  qfinal<-qfinal[order(mcols(qfinal)$inputhsp)]
  sfinal<-sfinal[order(mcols(sfinal)$inputhsp)]
  mylist<-list('qfinal'=qfinal,'sfinal'=sfinal,'finalalignments'=finalalignments)
  return(mylist)
}

#function to convert qfinal/sfinal and qtrimmed/strimmed lists to list of single dataframes which can then be written to file
grtodf<-function(qgr,sgr) {  #same but made generic
  myqdf<-as.data.frame(qgr)
  mysdf<-as.data.frame(sgr)
  colnames(myqdf)<-c('seqnames','qstart','qend', 'qwidth', 'strand',names(mcols(qgr)))
  colnames(mysdf)<-c('seqnames','sstart','send', 'swidth', 'strand',names(mcols(sgr)))
  combineddf<-cbind(myqdf[c("seqnames","qcontignames","snames","scontignames","inputhsp","alnlen","pid","bitscore","mystrand","qstart","qend","qwidth")],mysdf[c("sstart","send","swidth")])
  colnames(combineddf)<-c("qname","qcontig","sname","scontig","originalhspindex","originalalnlen","pid","bitscore","strand","qstart","qend","qwidth","sstart","send","swidth")
  combineddf[c("qname","qcontig","sname","scontig")]<-sapply(combineddf[c("qname","qcontig","sname","scontig")],as.vector)
  return(combineddf)
}


#trimming functions
reducefunction<-function(x,keepbisectedranges) {
  #reduces (joins contiguous) disjoint ranges, after splitting by input hsp; if there are discontiguous ranges for a given hsp (due to alignment being split in two - if it's longer on query/subject but is considered suboptimal e.g. lower bitscore, if using bitscore as best hit selection criterion) 1) selects longest range of hsp ranges if keepbisectedranges=True 2) excludes hsp ranges if keepbisectedranges=False. N.B there should be at least one best alignment retained for every query-sample pair, even with keepbisectedranges=False, so empty granges object bug shouldn't occur
  myinputhsp<-unique(mcols(x)$inputhsp)
  output<-reduce(x)
  if (length(output)>1) { #discontiguous (bisected) range
    if (keepbisectedranges=='False') {
      return(NULL)
    }
    output<-output[which.max(width(output)),]
  }
  mcols(output)$inputhsp<-myinputhsp
  return(output)
}


splitreducecombine<-function(x,keepbisectedranges) {  
  #takes list of ranges for a given paired sample; splits into needs reducing / doesn't need reducing; reduces; combines 
  isduplicate<-duplicated(mcols(x)$inputhsp) | duplicated(mcols(x)$inputhsp, fromLast=TRUE)
  numdup<-sum(isduplicate)
  if(numdup>0) {
    if(numdup==length(isduplicate)) {
      duplicates<-x
      splitduplicates<-split(duplicates,mcols(duplicates)$inputhsp)
      reducedduplicates<-lapply(splitduplicates, function(x) x=reducefunction(x,keepbisectedranges))
      includedalignmentindices<-lapply(reducedduplicates,length)>0
      reducedduplicates<-reducedduplicates[includedalignmentindices] #at least one best alignment should be remain
      combinedreduced<-do.call(getMethod(c, "GenomicRanges"), reducedduplicates)
      final<-combinedreduced
    } else {
      duplicates<-x[isduplicate]
      nonduplicates<-x[!isduplicate]
      splitduplicates<-split(duplicates,mcols(duplicates)$inputhsp)
      reducedduplicates<-lapply(splitduplicates, function(x) x=reducefunction(x,keepbisectedranges))
      includedalignmentindices<-lapply(reducedduplicates,length)>0
      if (all(includedalignmentindices==FALSE)) {
        final<-nonduplicates
      } else {
        reducedduplicates<-reducedduplicates[includedalignmentindices]
        combinedreduced<-do.call(getMethod(c, "GenomicRanges"), reducedduplicates)
        final<-append(nonduplicates, combinedreduced)
      }
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
  mcols(x)$bitscore<-bitscore[myinputhsps]
  mcols(x)$mystrand<-mystrand[myinputhsps]
  mcols(x)$snames<-snames[myinputhsps]
  mcols(x)$alnlen<-alnlen[myinputhsps] #added so that short alignments can be filtered prior to breakpoint calculation
  mcols(x)$qcontignames<-qcontignames[myinputhsps]
  mcols(x)$scontignames<-scontignames[myinputhsps]
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
  xcancel<-x-y
  ycancel<-y-x
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


trimalignments<-function(qfinal,sfinal,finalalignments,qtrimonly=FALSE) {
  #function trims disjoint alignments (after applying intersectsortalignments); trims query alignments based on subject trimming and vice-versa
  #trim QUERY alignments
  #get addstart/minusend vectors for query and subject
  qaddstart<-start(sfinal)-finalalignments$sstart #based on subject trimming, how much should be trimmed from query 
  qminusend<-finalalignments$send-end(sfinal)
  saddstart<-start(qfinal)-finalalignments$qstart #based on query trimming, how much should be trimmed from subject - need to know this at query trimming stage to avoid double trimming
  sminusend<-finalalignments$qend-end(qfinal)
  #get indices of alignments that require position reassignment
  qmydiffindices<-which((qaddstart > 0 | qminusend > 0)==TRUE)  
  smydiffindices<-which((saddstart > 0 | sminusend > 0)==TRUE)
  alldiffindices<-intersect(qmydiffindices,smydiffindices) #indices where there is trimming on both query and subject - this is used to prevent double trimming
  posindices<-which(mcols(qfinal[alldiffindices])$mystrand=='+') #of the alldiffindices, which are positive strand
  negindices<-which(mcols(qfinal[alldiffindices])$mystrand=='-')
  #edit addstart/minusend vectors to cancel out double trimming - where an alignment is trimmed from same flank on both query and subject genome - in this case just want to trim to meet the maximal trim
  cancelout<-canceldoubletrimwrapper(qaddstart,qminusend,saddstart,sminusend,alldiffindices,posindices,negindices)
  qaddstart<-cancelout[[1]];qminusend<-cancelout[[2]];saddstart<-cancelout[[3]];sminusend<-cancelout[[4]]
  #calculate new start/end positions
  lenmydiff<-length(qmydiffindices)
  if (lenmydiff>0) {
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
      mydiffindicesneg<-qmydiffindices[-strandpos]
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
  }
  if (qtrimonly==FALSE) {
    #trim SUBJECT alignments
    #calculate new start/end positions
    lenmydiff<-length(smydiffindices)
    if (lenmydiff>0) {
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
        mydiffindicesneg<-smydiffindices[-strandpos]
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
    }
    
    #after any necessary filtering/trimming, return qfinal (and sfinal if qtrimonly=FALSE)
    mylist<-list("qfinal"=qfinal,"sfinal"=sfinal)
    return(mylist)
  } else {
    return(qfinal)
  }
}


reformatcombineddf<-function(combineddf) {  #disjointdf/trimmeddf has separate genome/contig columns - merge these into a single column
  combineddf$qname<-ifelse(is.na(combineddf$qcontig), combineddf$qname, paste(combineddf$qname,combineddf$qcontig,sep='|'))
  combineddf$sname<-ifelse(is.na(combineddf$scontig), combineddf$sname, paste(combineddf$sname,combineddf$scontig,sep='|'))
  combineddf<-combineddf[c("qname","sname","originalhspindex","originalalnlen","pid","strand","qstart","qend","qwidth","sstart","send","swidth")]
  return(combineddf)
}


#breakpoint caluclation functions
makepairs<-function(x) mapply(c, head(x,-1), tail(x,-1), SIMPLIFY = FALSE)

BP<-function(x,y,contigsbool=FALSE,contigs=NULL,samecontigadjacencies=NULL,numpairs=NULL) { #this function works on signed permuations (numeric vectors with +/- indicated)
  if (contigsbool==FALSE) {
    makepairsx<-makepairs(x)
    out<-makepairsx[!(makepairsx %in% makepairs(y) & unlist(lapply(makepairsx, function(x) sum(sign(x))))==2 | makepairsx %in% makepairs(rev(y)) & unlist(lapply(makepairsx, function(x) sum(sign(x))))==-2)]
    bps<-length(out)
    numpairs<-length(makepairsx) #==alignments-1
  } else {
    makepairsx<-makepairs(x)
    out<-makepairsx[!(makepairsx %in% makepairs(y) & unlist(lapply(makepairsx, function(x) sum(sign(x))))==2 | makepairsx %in% makepairs(rev(y)) & unlist(lapply(makepairsx, function(x) sum(sign(x))))==-2) & samecontigadjacencies]
    bps<-length(out)
  }
  return(list('bps'=bps,'numpairs'=numpairs))
}

BPfunc<-function(x,y) {  #apply to granges objects
  stopifnot(nrow(x)==nrow(y))
  qcontigs<-mcols(x)$qcontignames
  scontigs<-mcols(y)$scontignames
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
  if (all(is.na(c(qcontigs,scontigs)))) {  #no contigs
    out<-BP(xinputhsps,yinputhsps,contigsbool=FALSE)
  } else {
    qsamecontigadjacencies<-unlist(lapply(makepairs(qcontigs), function(x) length(unique(x))==1))
    qnumpairs<-sum(qsamecontigadjacencies)
    ssamecontigadjacencies<-unlist(lapply(makepairs(scontigs), function(x) length(unique(x))==1))
    snumpairs<-sum(ssamecontigadjacencies)
    minpairsindex<-which.min(c(qnumpairs,snumpairs))
    if (minpairsindex==1) {  #calculate stats based on samples with fewest adjacencies (qcontigs)
      out<-BP(xinputhsps,yinputhsps,contigsbool=TRUE,contigs=qcontigs,samecontigadjacencies=qsamecontigadjacencies,numpairs=qnumpairs)
    } else {
      out<-BP(yinputhsps,xinputhsps,contigsbool=TRUE,contigs=scontigs,samecontigadjacencies=ssamecontigadjacencies,numpairs=snumpairs)
    }
  }
  bpcount<-out$bps
  pairscount<-out$numpairs
  return(c("bpcount"=bpcount,"alncount"=alncount,"pairscount"=pairscount))
}

breakpointcalc<-function(qtrimmed,strimmed,mydf) {
  #first filter short alignments (based on post-trimming alignment length)
  includedindices<-mapply(filtershortwidthalignments,qtrimmed,strimmed,widthfilter=lengthfilter,SIMPLIFY=FALSE)
  qtrimmed<-map2(qtrimmed,includedindices, `[`) #filter list of alignments using list of included indices
  strimmed<-map2(strimmed,includedindices, `[`)
  includedindices<-lapply(qtrimmed,length)>0 #use included indices to remove list elements with no alignments after mapping includedindices list
  if (all(includedindices==FALSE)) { #if no alignments remain across all queries after applying length filter, fill dataframe with NAs
    mydfbp<-cbind(querysample=rownames(mydf),subjectsample=rep(sample,nrow(mydf)),breakpoints=rep(NA,nrow(mydf)),alignments=rep(NA,nrow(mydf)),pairs=rep(NA,nrow(mydf)))
    myfinaldf<-merge(mydf,mydfbp,by=c("querysample","subjectsample"),all=TRUE) #all=TRUE is redundant here
    return(myfinaldf)
  }
  qtrimmed<-qtrimmed[includedindices]
  strimmed<-strimmed[includedindices]
  #sort by inputhsp
  qtrimmedsorted<-lapply(qtrimmed, function(x) sort(x,ignore.strand=T))
  strimmedsorted<-lapply(strimmed, function(x) sort(x,ignore.strand=T))
  #calculate breakpoint stats
  bpout<-mapply(BPfunc,qtrimmedsorted,strimmedsorted,SIMPLIFY=TRUE)
  mydfbp<-as.data.frame(t(bpout))
  colnames(mydfbp)<-c('breakpoints','alignments','pairs')
  mydfbp<-cbind(querysample=rownames(mydfbp),subjectsample=rep(sample,nrow(mydfbp)),mydfbp)
  #merge dataframes; replace missing breakpoint data with NA where necessary (if all alignments have been filtered due to filtered) - use all=TRUE argument to achieve NA missing cell replacement
  myfinaldf<-merge(mydf,mydfbp,by=c("querysample","subjectsample"),all=TRUE) #all=TRUE means keep all and fill with NAs
  return(myfinaldf)
}


filtershortwidthalignments<-function(qtrimmed,strimmed,widthfilter) {
  qwidth<-width(qtrimmed)
  swidth<-width(strimmed)
  includedindices<-qwidth > widthfilter & swidth > widthfilter
  return(includedindices)
}


bpdistcalc<-function(bps,pairs) {
  if (pairs==0) {
    return(as.numeric(0))
  } else {
    return(as.numeric(bps/pairs))
  }
}


#bootstrapping functions

combinebootfunc<-function(x) {
  bootlist<-lapply(x, function(l) l[[2]])
  bootlist<-lapply(bootlist, function(x) as.data.frame(rbindlist(x,idcol="index")))
  return(as.data.frame(do.call(rbind,bootlist)))
}


#alignment length distribution function

getalnlenstats<-function(x,widthfilter,quantiles) {
  alnlens<-rev(sort(width(x)))
  alnlenbool<-alnlens>widthfilter
  if(all(alnlenbool==FALSE)) {
    lxvector<-rep(NA,length(lxcols))
    nxvector<-rep(NA,length(nxcols))
    return(list(lxvector,nxvector))
  }
  alnlens<-alnlens[alnlenbool] #calculations will be based on length filtered alignments
  totalalnlen<-sum(alnlens)
  numquant<-length(quantiles)
  nxvector<-integer(numquant)
  lxvector<-integer(numquant)
  for (i in seq_len(numquant)) {
    lx<-which(cumsum(alnlens)>(totalalnlen*quantiles[i]))[1]
    nx<-alnlens[lx]
    lxvector[i]<-lx
    nxvector[i]<-nx
  }
  return(list(lxvector,nxvector))
}


#final ditance stats calculation functions
statsfunc<-function(stats, breakpoint,mygenomelenvector,mygenomelen,mymingenomelen,bootstrap,bidirectionalblast,alnlenstats='False') {
  hsplength<-as.numeric(stats$hsplength)
  hspidpositions<-as.numeric(stats$hspidpositions)
  covbreadth<-as.numeric(hsplength/mygenomelen)
  covbreadthmin<-as.numeric(hsplength/mymingenomelen)
  d0<-as.numeric(1-covbreadth)
  d1<-as.numeric(1-covbreadthmin)
  d2<-as.numeric(-log(covbreadth))
  d3<-as.numeric(-log(covbreadthmin))
  percentid<-as.numeric(hspidpositions/hsplength)
  d4<-as.numeric(1-percentid)
  d5<-as.numeric(-log(percentid))
  hspidgenlen<-as.numeric(hspidpositions/mygenomelen)
  hspidmingenlen<-as.numeric(hspidpositions/mymingenomelen)
  d6<-as.numeric(1-hspidgenlen)
  d7<-as.numeric(1-hspidmingenlen)
  d8<-as.numeric(-log(hspidgenlen))
  d9<-as.numeric(-log(hspidmingenlen))
  returnvector<-c(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,percentid,covbreadth,covbreadthmin)
  if (bootstrap=='True') {
    return(returnvector)
  } else {
    sample1len<-mygenomelenvector[1]
    sample2len<-mygenomelenvector[2]
    sample1hsplengthpretrim<-as.numeric(stats$qhsplengthpretrim)
    sample2hsplengthpretrim<-as.numeric(stats$shsplengthpretrim)
    sample1hspidpositionspretrim<-as.numeric(stats$qhspidpositionspretrim)
    sample2hspidpositionspretrim<-as.numeric(stats$shspidpositionspretrim)
    sample1percentid<-as.numeric(sample1hspidpositionspretrim/sample1hsplengthpretrim)
    sample2percentid<-as.numeric(sample2hspidpositionspretrim/sample2hsplengthpretrim)
    if (bidirectionalblast=='False') {
      sample1covbreadth<-as.numeric(sample1hsplengthpretrim/sample1len)
      sample2covbreadth<-as.numeric(sample2hsplengthpretrim/sample2len)
      bpdist2divisor=1000
    } else {
      sample1covbreadth<-as.numeric(sample1hsplengthpretrim/(2*sample1len))
      sample2covbreadth<-as.numeric(sample2hsplengthpretrim/(2*sample2len))
      bpdist2divisor=2000
    }
    returnvector<-c(returnvector,sample1covbreadth,sample2covbreadth,sample1percentid,sample2percentid)
    if (breakpoint=='True') {
      breakpoints<-as.numeric(stats$breakpoints)
      alignments<-as.numeric(stats$alignments)
      pairs<-as.numeric(stats$pairs)
      if (is.na(breakpoints) || is.na(pairs)) {
        bpdist=NA
        bpdist2=NA
      } else {
         bpdist<-bpdistcalc(breakpoints,pairs)
         bpdist2<-as.numeric(breakpoints/as.numeric((hsplength/bpdist2divisor)))#express per kb
      }
      returnvector<-c(returnvector,bpdist,bpdist2,breakpoints,alignments,pairs)
    }
    if (alnlenstats=='True') {
      lxstats<-as.integer(stats[lxcols])
      nxstats<-as.integer(stats[nxcols])
      returnvector<-c(returnvector,lxstats,nxstats)
    }
    return(returnvector)
  }
}


mergestats<-function(x,y) {
  merge<-data.frame()
  for (i in seq_along(colnames(x))) {
    mycol<-colnames(x)[i]
    if (mycol=="hspidpositions" || mycol=="hsplength" || mycol=="qhsplengthpretrim" || mycol=="shsplengthpretrim" || mycol=="qhspidpositionspretrim" || mycol=="shspidpositionspretrim") {
      merge[1,i]<-as.integer(x[i])+as.integer(y[i])
    } else {
      #merge[1,i]<-mean(as.integer(x[i]),as.integer(y[i]))
      merge[1,i]<-sum(as.integer(x[i]),as.integer(y[i]))/2
    }
  }
  colnames(merge)<-colnames(x)
  return(merge)
}


#functions for converting allsampledflist to finaldf

reorderallsampledf<-function(df) { #sort query/subject samples alphabetically and sort pretrim hsplengths/hpsidpositions accordingly; after reordering, the terms 'query'/'sample' lose meaning
  ordervec<-order(df[c(1,2)])
  df<-df[c(ordervec,ordervec+4,ordervec+6)]
  return(df)
}

getseqlens<-function(sample1,sample2) {
  #provide genome names; returns genome length stats from sampleseqlen dataframe
  seqlens<-c(sampleseqlen[which(sampleseqlen$sample==sample1),"length"],sampleseqlen[which(sampleseqlen$sample==sample2),"length"])
  return(seqlens)
}

applystatscalc<-function(mystats,mystatscols,bootstrap='False') {
  #wrapper function for statsfunc(); apply to a list of dataframes, split by genome names; each dataframe ('mystats') will have 1 or 2 rows (usually 2 rows when there is data for both blast directions)
  sample1<-mystats$querysample[1];sample2<-mystats$subjectsample[1]
  mygenomelenvector<-getseqlens(sample1,sample2) #mygenomelenvector is the sample1,sample2 length vector
  mygenomelen<-sum(mygenomelenvector);mymingenomelen<-min(mygenomelenvector)*2
  if(nrow(mystats)==1) {
    stats<-mystats[,mystatscols]
    if (bidirectionalblast=='False') {
      mygenomelen<-mygenomelen/2;mymingenomelen<-mymingenomelen/2
    }
  } else {
    stats<-mergestats(mystats[1,mystatscols],mystats[2,mystatscols])
  }
  if (stats$hsplength>mymingenomelen) {
    stats$hsplength<-mymingenomelen
  }
  if (stats$hspidpositions>mymingenomelen) {
    stats$hspidpositions<-mymingenomelen
  }
  statsout<-statsfunc(stats,breakpoint,mygenomelenvector,mygenomelen,mymingenomelen,bootstrap,bidirectionalblast,alnlenstats)
  return(c(sample1,sample2,mygenomelenvector,statsout))
}



###iterate through samples, to get initial raw statistics for sample-pairs

#read seqlength file
seqlenreport<-fread(gsubfn('%1',list('%1'=args[1]),'%1/seqlengths.tsv'),sep='\t',colClasses=c('character','integer'))
colnames(seqlenreport)<-c('sequence','length')
seqlenreport<-seqlenreport[order(seqlenreport$sequence),]

#read samples file
samples<-read.table(gsubfn('%1',list('%1'=args[1]),'%1/includedsubjects.txt'),sep='\t',header=FALSE)
samples<-as.character(samples[,1])
#samples<-samples[1:6]


suppressWarnings(suppressMessages(library('foreach',quietly=TRUE)))
suppressWarnings(suppressMessages(library('doParallel',quietly=TRUE)))

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)

allsampledflist<-list()
allsampledflist<-foreach(i=1:length(samples), .packages = c('gsubfn','GenomicRanges','purrr','data.table')) %dopar% {
    #read alignmnents file for given sample
    sample<-samples[i]
    report<-fread(gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/blast/%2/alignments.tsv'),select=c('qname','sname','pid','alnlen','mismatches','gapopens','qstart','qend','sstart','send','evalue','bitscore','strand'),colClasses = list('character'=c('qname','sname','strand'), 'numeric'=c('pid','evalue','bitscore'),'integer'=c('alnlen','mismatches','gapopens','qstart','qend','sstart','send')),header=TRUE,sep='\t')
    originalreport<-report
    originalreport$originalhspindex<-rownames(originalreport)
    #get information for shifting query and subject ranges where there are multiple contigs)
    seqlenreportseqs<-sapply(strsplit(as.vector(seqlenreport$sequence),'|',fixed=T),pastefunction)
    reformattedqname<-as.factor(sapply(strsplit(as.vector(report$qname),"|",fixed=T),pastefunction))
    reformattedsname<-as.factor(sapply(strsplit(as.vector(report$sname),"|",fixed=T),pastefunction))
    addqlens<-rangeshifting(reformattedqname, seqlenreport, seqlenreportseqs)
    addslens<-rangeshifting(reformattedsname, seqlenreport, seqlenreportseqs)
    #remove any contig information; save contig info to variables first
    qcontignames<-sapply(strsplit(as.vector(report$qname),"|",fixed=T),function(x) x=x[2])
    scontignames<-sapply(strsplit(as.vector(report$sname),"|",fixed=T),function(x) x=x[2])
    report$qname<-sapply(strsplit(as.vector(report$qname),"|",fixed=T),function(x) x=x[1])
    report$sname<-sapply(strsplit(as.vector(report$sname),"|",fixed=T),function(x) x=x[1])   
    ###disjoin method trimming
    pid<-report$pid
    bitscore<-report$bitscore
    snames<-report$sname
    alnlen<-report$alnlen
    mystrand<-report$strand  #strand info is lost after disjoin - need to retain strand info for toppid
    qgr<-GRanges(seqnames = report$qname, ranges = IRanges(start=(report$qstart), end = (report$qend)), strand=(report$strand))
    sgr<-GRanges(seqnames = report$qname, ranges = IRanges(start=(report$sstart), end= (report$send)), strand = (report$strand))
    #qalnlen<-width(qgr)
    #salnlen<-width(sgr)
    #shift query and subject ranges where there are multiple contigs
    qgr<-shift(qgr,addqlens)
    sgr<-shift(sgr,addslens)
    #overwrite report with shifted forward ranges
    if (sum(addslens)>0) {
       report$sstart<-start(sgr)
       report$send<-end(sgr)
    }
    if (sum(addqlens)>0) {
       report$qstart<-start(qgr)
       report$qend<-end(qgr)
    }
    #query disjoin
    gr2<-disjoin(qgr,with.revmap=TRUE,ignore.strand=TRUE)
    revmap<-gr2$revmap
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(report[,get(alnrankmethod)][x])]))
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample                                            
    qreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x,keepbisectedrangesarg))
    #subject disjoin
    gr2<-disjoin(sgr,with.revmap=TRUE,ignore.strand=TRUE)
    revmap<-gr2$revmap
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(report[,get(alnrankmethod)][x])]))
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample
    sreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x,keepbisectedrangesarg))
    #add pid, strand, subject name, and alignment length
    qfinal<-lapply(qreducedoutput, function(x) x=addcols(x))
    sfinal<-lapply(sreducedoutput, function(x) x=addcols(x))
    #write best blast alignments to file
    if (outputbestblastalignments=='True') {
      bestblasthitslist<-mapply(getbestblasthits,qfinal,sfinal,SIMPLIFY = FALSE)
      bestblasthitsdf<-do.call(rbind,bestblasthitslist)
      #adhere to blast outfmt 6 (-ve strand indicated implicitly by flipping sstart/send)
      strandispositive<-bestblasthitsdf$strand=='+'
      blastsstart<-ifelse(strandispositive,bestblasthitsdf$sstart,bestblasthitsdf$send)
      blastsend<-ifelse(strandispositive,bestblasthitsdf$send,bestblasthitsdf$sstart)
      bestblasthitsdf$sstart<-blastsstart
      bestblasthitsdf$send<-blastsend
      bestblasthitsdf$strand<-NULL
      #
      bestblasthitsdf<-bestblasthitsdf[with(bestblasthitsdf,order(bestblasthitsdf$qname,bestblasthitsdf$sname,-bestblasthitsdf$bitscore)),]
      write.table(bestblasthitsdf, file=gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/output/bestblastalignments_%2.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
    }
    #get pre-trimmmed hsplength and hspidpositions
    qhsplenpretrim<-as.data.frame(do.call(rbind,lapply(qfinal,gethsplength)))
    shsplenpretrim<-as.data.frame(do.call(rbind,lapply(sfinal,gethsplength)))
    colnames(qhsplenpretrim)<-'qhsplenpretrim'
    colnames(shsplenpretrim)<-'shsplenpretrim'
    qhspidpositionspretrim<-as.data.frame(do.call(rbind,lapply(qfinal,gethspidpositions)))
    shspidpositionspretrim<-as.data.frame(do.call(rbind,lapply(sfinal,gethspidpositions)))
    colnames(qhspidpositionspretrim)<-'qhspidpositionspretrim'
    colnames(shspidpositionspretrim)<-'shspidpositionspretrim'
    #filter disjoint alignments - remove alignments not present in both qfinal and sfinal; then sort by inputhsp
    intersectsortoutput<-transpose(mapply(intersectsortalignments,qfinal,sfinal,SIMPLIFY = FALSE))
    qfinal<-intersectsortoutput$qfinal
    sfinal<-intersectsortoutput$sfinal
    finalalignments<-intersectsortoutput$finalalignments
    #write disjoint (non-overlapping) alignments to file
    if (outputnonoverlappingalignments=='True') {
      disjointdflist<-mapply(grtodf,qfinal,sfinal,SIMPLIFY=FALSE)
      disjointdf<-do.call(rbind,disjointdflist)
      minusqlens<-addqlens[disjointdf$originalhspindex]
      minusslens<-addslens[disjointdf$originalhspindex]
      disjointdf$qstart<-disjointdf$qstart-minusqlens
      disjointdf$qend<-disjointdf$qend-minusqlens
      disjointdf$sstart<-disjointdf$sstart-minusslens
      disjointdf$send<-disjointdf$send-minusslens
      disjointdf<-reformatcombineddf(disjointdf)
      write.table(disjointdf, file=gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/output/nonoverlappingalignments_%2.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
    }
    if (statsfromtrimmed=='False') {
      #get hsp id/len stats from disjoint alignments
      mystats<-mapply(getstatsfromdisjoint, qfinal,sfinal,SIMPLIFY = FALSE)
      mydf<-as.data.frame(do.call(rbind, mystats)) #convert list of vectors to dataframe
      mydf<-cbind(querysample=rownames(mydf),subjectsample=rep(sample,nrow(mydf)),mydf,qhsplenpretrim,shsplenpretrim,qhspidpositionspretrim,shspidpositionspretrim)
      #get breakpoint stats
      if (breakpoint=='True') {
        myfinaldf<-breakpointcalc(qfinal,sfinal,mydf)
      } else {
        myfinaldf<-mydf
      }
      #get alignment length distribution stats
      if (alnlenstats=='True') {
        myquantiles=as.numeric(alnlenstatsquantiles)/100
        alnlenstatslist<-lapply(qfinal, getalnlenstats, widthfilter=lengthfilter,quantiles=myquantiles)
        alnlenstatsdf<-cbind(rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[1]]),as.data.frame.list),idcol="querysample",use.names=FALSE),rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[2]]),as.data.frame.list),use.names=FALSE))
        colnames(alnlenstatsdf)<-c('querysample',lxcols,nxcols)
        myfinaldf<-merge(myfinaldf,alnlenstatsdf,by="querysample")
      }
    }
    #trim alignments
    if (statsfromtrimmed=='True' || outputtrimmedalignments=='True') {
      trimmedalignments<-transpose(mapply(trimalignments,qfinal,sfinal,finalalignments,SIMPLIFY = FALSE))
      qtrimmed<-trimmedalignments$qfinal
      strimmed<-trimmedalignments$sfinal
      includedalignmentindices<-lapply(qtrimmed,length)>0 & lapply(strimmed,length)>0 #safeguards against bug due to empty list element (probably unecessary)
      qtrimmed<-qtrimmed[includedalignmentindices]
      strimmed<-strimmed[includedalignmentindices]
    } else {
      qtrimmed<-mapply(trimalignments,qfinal,sfinal,finalalignments,qtrimonly=TRUE,SIMPLIFY=FALSE)
      qtrimmed<-qtrimmed[lapply(qtrimmed,length)>0]
    }
    #write trimmed alignments to file
    if (outputtrimmedalignments=='True') {
      trimmeddflist<-mapply(grtodf, qtrimmed, strimmed, SIMPLIFY = FALSE)
      trimmeddf<-do.call(rbind,trimmeddflist)
      #shift back query/subject positions and reformat df so genome/contig columns are merged to single column
      minusqlens<-addqlens[trimmeddf$originalhspindex]
      minusslens<-addslens[trimmeddf$originalhspindex]
      trimmeddf$qstart<-trimmeddf$qstart-minusqlens
      trimmeddf$qend<-trimmeddf$qend-minusqlens
      trimmeddf$sstart<-trimmeddf$sstart-minusslens
      trimmeddf$send<-trimmeddf$send-minusslens
      trimmeddf<-reformatcombineddf(trimmeddf)
      write.table(trimmeddf, file=gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/output/trimmedalignments_%2.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
    }
    if (statsfromtrimmed=='True') {
      #get hsp id stats
      mystats<-lapply(qtrimmed,getstatsfromtrimmed)
      mydf<-as.data.frame(do.call(rbind, mystats)) #convert list of vectors to dataframe
      mydf<-cbind(querysample=rownames(mydf),subjectsample=rep(sample,nrow(mydf)),mydf,qhsplenpretrim,shsplenpretrim,qhspidpositionspretrim,shspidpositionspretrim)
      #get breakpoint stats
      if (breakpoint=='True') {
        myfinaldf<-breakpointcalc(qtrimmed,strimmed,mydf)
      } else {
        myfinaldf<-mydf
      }
      #get alignment length distribution stats
      if (alnlenstats=='True') {
        myquantiles=as.numeric(alnlenstatsquantiles)/100
        alnlenstatslist<-lapply(qtrimmed, getalnlenstats, widthfilter=lengthfilter,quantiles=myquantiles)
        alnlenstatsdf<-cbind(rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[1]]),as.data.frame.list),idcol="querysample"),rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[2]]),as.data.frame.list)))
        colnames(alnlenstatsdf)<-c('querysample',lxcols,nxcols)
        myfinaldf<-merge(myfinaldf,alnlenstatsdf,by="querysample")
      }
    }
    #IF NO BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS
    if (boot==0) {
      print(myfinaldf)
    } else {
      #IF BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS + BOOTSTRAPPED STATS; N.B NO LONGER RESAMPLING BREAKPOINT/ALNLEN STATS
      myfinaldfbootlist<-vector("list",boot)
      if (statsfromtrimmed=='False') {
        for (z in seq_len(boot)) {
          indices<-lapply(qfinal, function(x) sample(1:length(x), replace=T)) #get indices for resampling qfinal/sfinal
          qfinalboot<-map2(qfinal,indices, `[`)
          sfinalboot<-map2(sfinal,indices, `[`)
          mystatsboot<-mapply(getstatsfromdisjoint, qfinalboot,sfinalboot,SIMPLIFY = FALSE)
          mydfboot<-as.data.frame(do.call(rbind,mystatsboot))
          mydfboot<-cbind(rownames(mydfboot),rep(sample,nrow(mydfboot)),mydfboot)
          myfinaldfboot<-mydfboot
          myfinaldfbootlist[[z]]<-myfinaldfboot
        }
      } else {
        for (z in seq_len(boot)) {
          indices<-lapply(qtrimmed, function(x) sample(1:length(x), replace=T)) #get indices for resampling qtrimmed/strimmed
          qtrimmedboot<-map2(qtrimmed,indices, `[`)
          names(qtrimmedboot)<-names(qtrimmed)
          mystatsboot<-lapply(qtrimmedboot,getstatsfromtrimmed)
          mydfboot<-as.data.frame(do.call(rbind,mystatsboot))
          mydfboot<-cbind(rownames(mydfboot),rep(sample,nrow(mydfboot)),mydfboot)
          myfinaldfboot<-mydfboot 
          myfinaldfbootlist[[z]]<-myfinaldfboot
        }
      }
      print(list(myfinaldf, myfinaldfbootlist))
    }
}

stopCluster(cl)


#create allsampledf and allsampledfboot colname vectors + statscols and statscolsboot vectors
allsampledfcolnames<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','qhspidpositionspretrim','shspidpositionspretrim')
statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','qhspidpositionspretrim','shspidpositionspretrim')
allsampledfbootcolnames<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
statscolsboot<-c('hspidpositions','hsplength')
if (breakpoint=='True') {
  allsampledfcolnames<-c(allsampledfcolnames,'breakpoints','alignments','pairs')
  statscols<-c(statscols,'breakpoints','alignments','pairs')
}
if (alnlenstats=='True') {
  allsampledfcolnames<-c(allsampledfcolnames,lxcols,nxcols)
  statscols<-c(statscols,lxcols,nxcols)
}

#convert allsampledf lists to data frames and use colname vectors to label columns
if (boot==0) {
  allsampledf<-as.data.frame(do.call(rbind, allsampledflist))
  colnames(allsampledf)<-allsampledfcolnames
} else { #N.B for bootstrapped data, there is only hspid/hsplength based on bootstrapped trimmmed alignments; no longer bootstrapping breakpoints
  allsampledf<-as.data.frame(do.call(rbind, lapply(allsampledflist, function(l) l[[1]])))
  colnames(allsampledf)<-allsampledfcolnames
  allsampledfboot<-combinebootfunc(allsampledflist)
  colnames(allsampledfboot)<-allsampledfbootcolnames
}

###get final stats (distance scores) for each pairwise sample combination

#aggregate seqlen repot at sample level
seqlenreport$sequence<-sapply(strsplit(as.vector(seqlenreport$sequence),"|",fixed=T),function(x) x=x[1])
sampleseqlen<-aggregate(seqlenreport$length, by=list(seqlenreport$sequence), FUN=sum)
colnames(sampleseqlen)<-c('sample','length')

#get final stats
allsampledf[c(1,2,5,6,7,8)]<-t(apply(allsampledf,1,reorderallsampledf)) #reorder query/subject sample alphabetically; reorder qhsplenpretrim/shsplenpretrim and qhspidpositionspretrim/shspidpositionspretrim accoringly (N.B. the terms query/subject now lose their meaning)
splitsampledf<-split(allsampledf, list(allsampledf$querysample,allsampledf$subjectsample),drop=TRUE) #split by pairwise combination

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)
clusterExport(cl,c("applystatscalc","getseqlens","sampleseqlen","mergestats","statsfunc","breakpoint","alnlenstats","bpdistcalc","lxcols","nxcols","bidirectionalblast"))

finaldf<-parLapply(cl, splitsampledf,applystatscalc, mystatscols=statscols)

stopCluster(cl)


finaldf<-as.data.frame(do.call(rbind,finaldf))

#N.B 'Genome1'/'Genome2' are used for column names of comparisonstats.tsv/comparisonstats_bootstrapped.tsv files, but in code in this script the term 'sample' instead of 'genome' is used

#create finaldfcolnames vector
finaldfcolnames<-c('Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Average_nucleotide_identity','Coverage_breadth_trimmed','Coverage_breadth_mingenome_trimmed','Coverage_breadth_Genome1','Coverage_breadth_Genome2','Percent_identity_Genome1','Percent_identity_Genome2')
if (breakpoint=='True') {
  finaldfcolnames<-c(finaldfcolnames,'Breakpoint_distance_d0','Breakpoint_distance_d1','Breakpoints','Alignments','Alignment_pairs')
}
if (alnlenstats=='True') {
  finaldfcolnames<-c(finaldfcolnames,lxcols,nxcols)
}

#use finaldfcolnames vector to label columns
colnames(finaldf)<-finaldfcolnames

#order by genome names and write to file
finaldf<-finaldf[with(finaldf,order(finaldf$Genome1,finaldf$Genome2)),]
write.table(finaldf, file=gsubfn('%1',list('%1'=args[1]),'%1/output/comparisonstats.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)


###if there is bootstrapping, need to produce additional stats table - split by bootstrap, apply same code as above to get stats for each boostrap, then combine

if (boot!=0) {

  allsampledfbootsplit<-split(allsampledfboot, allsampledfboot$bootstrap)
  finaldflist<-vector("list", as.integer(boot))
  alnlenstats='False'
  breakpoint='False' #NO LONGER DOING BOOTSTRAPPING OF BREAKPOINTS

  cl<-makeCluster(as.integer(args[2]))
  registerDoParallel(cl)
  
  finaldflist<-foreach(i=names(allsampledfbootsplit)) %dopar% {
    allsampledf<-allsampledfbootsplit[[i]]
    allsampledf[,c(2,3)]<-t(apply(allsampledf[,c(2,3)],1,sort)) #sort samples alphabetically
    splitsampledf<-split(allsampledf, list(allsampledf$querysample,allsampledf$subjectsample),drop=TRUE) #split by pairwise combination 
    finaldfboot<-lapply(splitsampledf,applystatscalc,mystatscols=statscolsboot,bootstrap='True')
    finaldfboot<-as.data.frame(do.call(rbind,finaldfboot))
    print(finaldfboot)
  }
  
  stopCluster(cl)
  
  finaldfboot<-rbindlist(finaldflist,idcol = "index")

  #create finaldfbootcolnames vector
  finaldfbootcolnames<-c('bootstrap','Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Average_nucleotide_identity','Coverage_breadth_trimmed','Coverage_breadth_mingenome_trimmed')

  colnames(finaldfboot)<-finaldfbootcolnames

  finaldfboot<-finaldfboot[with(finaldfboot,order(as.numeric(finaldfboot$bootstrap),finaldfboot$Genome1,finaldfboot$Genome2)),]
  write.table(finaldfboot, file=gsubfn('%1',list('%1'=args[1]),'%1/output/comparisonstats_bootstrapped.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

}




###OLD CODE


###depreacated code

#code for getting bootstrap resampled breakpoint/alnlenstats statistics:

          #colnames(mydfboot)<-c('querysample','subjectsample','hspidpositions','hsplength')
          #get breakpoint stats
          #if (breakpoint=='True') {
          #  strimmedboot<-lapply(1:length(strimmed), FUN=function(x, list1, list2) list1[[x]][list2[[x]]] , list1=strimmed, list2=indices) #resample strimmed using indices
          #  names(strimmedboot)<-names(strimmed)
          #  myfinaldfboot<-breakpointcalc(qtrimmedboot,strimmedboot,mydfboot)
          #} else {
          #  myfinaldfboot<-mydfboot
          #}


###deprecated fuctions

# #breakpoint calculation functions
# filtershortalignments<-function(x,lengthfilter) { #this function is deprecated - now using filtershortwidthalignments instead - using post-trimmed alignment length rather than pre-trimmed alignment length
#   includedindices<-mcols(x)$alnlen>lengthfilter
#   return(x[includedindices])
# }




