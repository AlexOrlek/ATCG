args = commandArgs(trailingOnly=TRUE)
library('gsubfn')
library('GenomicRanges')
library('purrr')
reduce<-GenomicRanges::reduce
shift<-GenomicRanges::shift
library(data.table)
rbindlist<-data.table::rbindlist
transpose<-purrr::transpose

#args[1] is filepath to pipeline output folder; args[2] is threads; args[3] is breakpoint stats; args[4] is alignment length stats; arg[5] is bootstrap number; arg[6] is best alignment selection criterion; arg[7] is length threshold for filtering alignments prior to breakpoint calculation; arg[8] is trimmed alignment output

cores=as.integer(args[2])
breakpoint=as.character(args[3])
alnlenstats=as.character(args[4])
boot=as.integer(args[5])
alnrankmethod=as.character(args[6])
lengthfilter=as.integer(args[7])
outputtrimmedalignments=as.character(args[8])

###define functions

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

getstats<-function(x) {
  #returns raw statistics                                                        
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




#function to convert qtrimmed/strimmed lists to list of single dataframes which can then be written to file

grtodf<-function(qtrimmed,strimmed) {
  myqdf<-as.data.frame(qtrimmed)
  mysdf<-as.data.frame(strimmed)
  colnames(myqdf)<-c('seqnames','qstart','qend', 'qwidth', 'strand',names(mcols(qtrimmed)))
  colnames(mysdf)<-c('seqnames','sstart','send', 'swidth', 'strand',names(mcols(strimmed)))
  trimmeddf<-cbind(myqdf[c("seqnames","qcontignames","snames","scontignames","inputhsp","alnlen","pid","bitscore","mystrand","qstart","qend","qwidth")],mysdf[c("sstart","send","swidth")])
  colnames(trimmeddf)<-c("qname","qcontig","sname","scontig","originalhspindex","originalalnlen","pid","bitscore","strand","qstart","qend","qwidth","sstart","send","swidth")
  trimmeddf[c("qname","qcontig","sname","scontig")]<-sapply(trimmeddf[c("qname","qcontig","sname","scontig")],as.vector)
  return(trimmeddf)
}



#trimming functions
reducefunction<-function(x) {
  #reduces (joins contiguous) disjoint ranges, after splitting by input hsp; if there are discontiguous ranges (due to alignment being split in two - if it's longer on query/subject but is considered suboptimal according to blast table alignment length) selects longest range
  myinputhsp<-unique(mcols(x)$inputhsp)
  output<-reduce(x)
  if (length(output)>1) {
    output<-output[which.max(width(output)),]
  }
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
  mcols(x)$bitscore<-bitscore[myinputhsps]
  mcols(x)$mystrand<-mystrand[myinputhsps]
  mcols(x)$snames<-snames[myinputhsps] #!!!ADDED
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


trimalignments<-function(qfinal,sfinal,qtrimonly=FALSE) {
  #function first filters alignments to include only those present in both qfinal and sfinal; then trims alignments - trims query alignments based on subject trimming and vice-versa
  #filter alignments of qfinal based on sfinal and vice-versa
  finalhsps<-sort(intersect(mcols(qfinal)$inputhsp,mcols(sfinal)$inputhsp))
  qfinal<-qfinal[mcols(qfinal)$inputhsp %in% finalhsps]
  sfinal<-sfinal[mcols(sfinal)$inputhsp %in% finalhsps]
  finalalignments<-report[finalhsps,]
  #order alignments
  qfinal<-qfinal[order(mcols(qfinal)$inputhsp)]
  sfinal<-sfinal[order(mcols(sfinal)$inputhsp)]
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


reformattrimmeddf<-function(trimmeddf) {  #trimmeddf has separate genome/contig columns - merge these into a single column
  trimmeddf$qname<-ifelse(is.na(trimmeddf$qcontig), trimmeddf$qname, paste(trimmeddf$qname,trimmeddf$qcontig,sep='|'))
  trimmeddf$sname<-ifelse(is.na(trimmeddf$scontig), trimmeddf$sname, paste(trimmeddf$sname,trimmeddf$scontig,sep='|'))
  trimmeddf<-trimmeddf[c("qname","sname","originalhspindex","originalalnlen","pid","strand","qstart","qend","qwidth","sstart","send","swidth")]
  return(trimmeddf)
}


#breakpoint caluclation functions
makepairs<-function(x) mapply(c, head(x,-1), tail(x,-1), SIMPLIFY = FALSE)
BP<-function(x,y) { #this function works on signed permuations (numeric vectors with +/- indicated)
  out<-makepairs(x)[!(makepairs(x) %in% makepairs(y) & unlist(lapply(makepairs(x), function(x) sum(sign(x))))==2 | makepairs(x) %in% makepairs(rev(y)) & unlist(lapply(makepairs(x), function(x) sum(sign(x))))==-2)]
  return(out)
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
  #split by contig
  qcontigsplit<-lapply(qtrimmed, function(x) split(x,bpsplit[mcols(x)$inputhsp])) 
  scontigsplit<-lapply(strimmed, function(x) split(x,bpsplit[mcols(x)$inputhsp])) #nested lists- need to split by query and subject contigs before calculating breakpoints
  bpout<-list()
  for (qname in names(qcontigsplit)) {
    #sort alignments by position
    qcontigsplit[[qname]]<-lapply(qcontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    scontigsplit[[qname]]<-lapply(scontigsplit[[qname]], function(x) sort(x,ignore.strand=T))
    numsplits<-length(qcontigsplit[[qname]])
    out<-mapply(BPfunc,qcontigsplit[[qname]],scontigsplit[[qname]],SIMPLIFY=TRUE) #apply breakpoint function to sublists
    bpcount<-sum(unlist((out["bpcount",]))) #sum across subject contig splits to get count per query
    alncount<-sum(unlist((out["alncount",])))
    pairscount<-alncount-numsplits #ADDED
    bpout[[qname]]<-c(bpcount,alncount,pairscount)
  }
  mydfbp<-as.data.frame(do.call(rbind, bpout))
  colnames(mydfbp)<-c('breakpoints','alignments','pairs')
  mydfbp<-cbind(querysample=rownames(mydfbp),subjectsample=rep(sample,nrow(mydfbp)),mydfbp)
  #merge dataframes; replace missing breakpoint data with NA where necessary (if all alignments have been filtered due to filtered) - use all=TRUE argument to achieve NA missing cell replacement
  myfinaldf<-merge(mydf,mydfbp,by=c("querysample","subjectsample"),all=TRUE) #all=TRUE means keep all and fill with NAs
  return(myfinaldf)
}



filtershortalignments<-function(x,lengthfilter) { #this function is deprecated - now using filtershortwidthalignments instead - using post-trimmed alignment length rather than pre-trimmed alignment length
  includedindices<-mcols(x)$alnlen>lengthfilter
  return(x[includedindices])
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

getalnlenstats<-function(x,widthfilter) {
  alnlens<-rev(sort(width(x)))
  alnlenbool<-alnlens>widthfilter
  if(all(alnlenbool==FALSE)) {
    lxvector<-rep(NA,length(lxcols))
    nxvector<-rep(NA,length(nxcols))
    return(list(lxvector,nxvector))
  }
  alnlens<-alnlens[alnlenbool] #calculations will be based on length filtered alignments
  totalalnlen<-sum(alnlens)
  quartiles<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  numquart<-length(quartiles)
  nxvector<-integer(numquart)
  lxvector<-integer(numquart)
  for (i in seq_len(numquart)) {
    lx<-which(cumsum(alnlens)>(totalalnlen*quartiles[i]))[1]
    nx<-alnlens[lx]
    lxvector[i]<-lx
    nxvector[i]<-nx
  }
  return(list(lxvector,nxvector))
}


#final ditance stats calculation functions
statsfunc<-function(stats, breakpoint,mygenomelenvector,mygenomelen,mymingenomelen,bootstrap,alnlenstats='False') {
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
    sample1covbreadth<-as.numeric(sample1hsplengthpretrim/(2*sample1len))
    sample2covbreadth<-as.numeric(sample2hsplengthpretrim/(2*sample2len))
    returnvector<-c(returnvector,sample1covbreadth,sample2covbreadth)
    if (breakpoint=='True') {
      breakpoints<-as.numeric(stats$breakpoints)
      alignments<-as.numeric(stats$alignments)
      pairs<-as.numeric(stats$pairs)
      if (is.na(breakpoints) || is.na(pairs)) {
        bpdist=NA
        bpdist2=NA
      } else {
         bpdist<-bpdistcalc(breakpoints,pairs)
         bpdist2<-as.numeric(breakpoints/as.numeric((hsplength/2000)))#express per kb
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
    if (mycol=="hspidpositions" || mycol=="hsplength" || mycol=="qhsplengthpretrim" || mycol=="shsplengthpretrim") {
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

reorderallsampledf<-function(df) { #sort query/subject samples alphabetically and shift pretrim hsplengths accordingly; after reordering, the terms 'query'/'sample' lose meaning
  ordervec<-order(df[c(1,2)])
  df<-df[c(ordervec,ordervec+4)]
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
  } else {
    stats<-mergestats(mystats[1,mystatscols],mystats[2,mystatscols])
  }
  if (stats$hsplength>mymingenomelen) {
    stats$hsplength<-mymingenomelen
  }
  if (stats$hspidpositions>mymingenomelen) {
    stats$hspidpositions<-mymingenomelen
  }
  statsout<-statsfunc(stats,breakpoint,mygenomelenvector,mygenomelen,mymingenomelen,bootstrap,alnlenstats)
  return(c(sample1,sample2,mygenomelenvector,statsout))
}



###iterate through samples, to get initial raw statistics for sample-pairs

#read seqlength file
seqlenreport<-fread(gsubfn('%1',list('%1'=args[1]),'%1/seqlengths.tsv'),sep='\t',colClasses=c('character','integer'))
colnames(seqlenreport)<-c('sequence','length')
seqlenreport<-seqlenreport[order(seqlenreport$sequence),]

#read samples file
library('foreach')
library('doParallel')

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)

samples<-read.table(gsubfn('%1',list('%1'=args[1]),'%1/included.txt'),sep='\t',header=FALSE)
samples<-as.character(samples[,1])
#samples<-samples[1:6]
allsampledflist<-list()

lxcols<-c('l10','l20','l30','l40','l50','l60','l70','l80','l90')
nxcols<-c('n10','n20','n30','n40','n50','n60','n70','n80','n90')


allsampledflist<-foreach(i=1:length(samples), .packages = c('gsubfn','GenomicRanges','purrr','data.table')) %dopar% {
    #read alignmnents file for given sample
    sample<-samples[i]
    #report<-fread(gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/blast/%2/alignments.tsv'),select=c(1,2,3,4,7,8,9,10,12,17),sep='\t') #same subject, different queries
    report<-fread(gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/blast/%2/alignments.tsv'),select=c('qname','sname','pid','alnlen','qstart','qend','sstart','send','bitscore','strand'),colClasses = list('character'=c('qname','sname','strand'), 'numeric'=c('pid','bitscore'),'integer'=c('alnlen','qstart','qend','sstart','send')),header=TRUE,sep='\t')
    #colnames(report)<-c('qname','sname','pid','alnlen','mismatches','gapopens','qstart','qend','sstart','send','evalue','bitscore','qcov','qcovhsp','qlength','slength','strand')
    #colnames(report)<-c('qname','sname','pid','alnlen','qstart','qend','sstart','send','bitscore','strand')
    #get information for shifting query and subject ranges where there are multiple contigs)
    seqlenreportseqs<-sapply(strsplit(as.vector(seqlenreport$sequence),'|',fixed=T),pastefunction)
    reformattedqname<-as.factor(sapply(strsplit(as.vector(report$qname),"|",fixed=T),pastefunction))
    reformattedsname<-as.factor(sapply(strsplit(as.vector(report$sname),"|",fixed=T),pastefunction))
    addqlens<-rangeshifting(reformattedqname, seqlenreport, seqlenreportseqs)
    addslens<-rangeshifting(reformattedsname, seqlenreport, seqlenreportseqs)
    #
    if (breakpoint=='True') {
        bpsplit<-paste(reformattedqname,reformattedsname,sep='_') #required for breakpoint calculation
    }
    #remove any contig information; save contig info to variables first
    qcontignames<-sapply(strsplit(as.vector(report$qname),"|",fixed=T),function(x) x=x[2])
    scontignames<-sapply(strsplit(as.vector(report$sname),"|",fixed=T),function(x) x=x[2])
    report$qname<-sapply(strsplit(as.vector(report$qname),"|",fixed=T),function(x) x=x[1]) #!!!ADDED
    report$sname<-sapply(strsplit(as.vector(report$sname),"|",fixed=T),function(x) x=x[1])   
    ###disjoin method trimming
    pid<-report$pid
    bitscore<-report$bitscore
    snames<-report$sname #!!!ADDDED
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
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(report[,get(alnrankmethod)][x])])) ##!changed alnlen[x] to alnrankmethod for flexibility; !changed qalnlen[x] to alnlen[x] - improves agreement between query/subject trimming in terms of which alignments are retained, avoiding excessive alignment filtering when finalhsps are selected as an intersection of query/subject hsps.
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample                                            
    qreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x))
    #subject disjoin
    gr2<-disjoin(sgr,with.revmap=TRUE,ignore.strand=TRUE)
    revmap<-gr2$revmap
    tophsp<-unlist(lapply(revmap, function(x) x[which.max(report[,get(alnrankmethod)][x])]))
    mcols(gr2)$inputhsp<-tophsp
    mcols(gr2)$revmap<-NULL
    grsplit<-split(gr2,seqnames(gr2)) #split by seqnames i.e. one list per paired sample
    sreducedoutput<-lapply(grsplit, function(x) x=splitreducecombine(x))
    #add pid, strand, subject name, and alignment length
    qfinal<-lapply(qreducedoutput, function(x) x=addcols(x))
    sfinal<-lapply(sreducedoutput, function(x) x=addcols(x))
    #get pre-trimmmed hsplength
    qhsplenpretrim<-as.data.frame(do.call(rbind,lapply(qfinal,gethsplength)))
    shsplenpretrim<-as.data.frame(do.call(rbind,lapply(sfinal,gethsplength)))
    colnames(qhsplenpretrim)<-'qhsplenpretrim'
    colnames(shsplenpretrim)<-'shsplenpretrim'
    #trim alignments
    if (breakpoint=='True' || outputtrimmedalignments=='True') {
      trimmedalignments<-transpose(mapply(trimalignments,qfinal,sfinal,SIMPLIFY = FALSE)) #!ADDED
      qtrimmed<-trimmedalignments$qfinal
      strimmed<-trimmedalignments$sfinal
      includedalignmentindices<-lapply(qtrimmed,length)>0 & lapply(strimmed,length)>0 #safeguards against bug due to empty list element (probably unecessary)
      qtrimmed<-qtrimmed[includedalignmentindices]
      strimmed<-strimmed[includedalignmentindices]
    } else {
      qtrimmed<-mapply(trimalignments,qfinal,sfinal,qtrimonly=TRUE,SIMPLIFY=FALSE)
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
      trimmeddf<-reformattrimmeddf(trimmeddf)
      write.table(trimmeddf, file=gsubfn('%1|%2',list('%1'=args[1],'%2'=sample),'%1/output/trimmedalignments_%2.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)
    }
    #get hsp id stats
    mystats<-lapply(qtrimmed,getstats) #!!!CHANGED
    mydf<-as.data.frame(do.call(rbind, mystats)) #convert list of vectors to dataframe
    mydf<-cbind(querysample=rownames(mydf),subjectsample=rep(sample,nrow(mydf)),mydf,qhsplenpretrim,shsplenpretrim)
    #get breakpoint stats
    if (breakpoint=='True') {
      myfinaldf<-breakpointcalc(qtrimmed,strimmed,mydf)
    } else {
      myfinaldf<-mydf
    }
    #get alignment length distribution stats
    if (alnlenstats=='True') {
       alnlenstatslist<-lapply(qtrimmed, getalnlenstats, widthfilter=lengthfilter)
       alnlenstatsdf<-cbind(rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[1]]),as.data.frame.list),idcol="querysample"),rbindlist(lapply(lapply(alnlenstatslist, function(l) l[[2]]),as.data.frame.list)))
       colnames(alnlenstatsdf)<-c('querysample',lxcols,nxcols)
       myfinaldf<-merge(myfinaldf,alnlenstatsdf,by="querysample")
    }
    #IF NO BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS
    if (boot==0) {
      print(myfinaldf)
    } else {
      #IF BOOTSTRAPPING, SAVE ALL ALIGNMENT STATS + BOOTSTRAPPED STATS
      myfinaldfbootlist<-vector("list",boot)
      for (z in seq_len(boot)) {
        #stopifnot(sapply(qtrimmed, function(x) length(x))==sapply(strimmed, function(x) length(x)))
        indices<-lapply(qtrimmed, function(x) sample(1:length(x), replace=T)) #get indices for resampling qtrimmed/strimmed
        qtrimmedboot<-lapply(1:length(qtrimmed), FUN=function(x, list1, list2) list1[[x]][list2[[x]]] , list1=qtrimmed, list2=indices) #resample qtrimmed using indices
        names(qtrimmedboot)<-names(qtrimmed)
        mystatsboot<-lapply(qtrimmedboot,getstats) #!!!CHANGED
        mydfboot<-as.data.frame(do.call(rbind,mystatsboot))
        mydfboot<-cbind(rownames(mydfboot),rep(sample,nrow(mydfboot)),mydfboot)
        #colnames(mydfboot)<-c('querysample','subjectsample','hspidpositions','hsplength')
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
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs',lxcols, nxcols)
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs',lxcols,nxcols)
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs')
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs')
  } else if (alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim',lxcols,nxcols)
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim',lxcols,nxcols)
  } else {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim')
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim')
  }
} else { #N.B for bootstrapped data, there is only hspid/hsplength based on bootstrapped trimmmed alignments
  allsampledf<-as.data.frame(do.call(rbind, lapply(allsampledflist, function(l) l[[1]])))
  allsampledfboot<-combinebootfunc(allsampledflist)
  if (breakpoint=='True' && alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs',lxcols,nxcols)
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs',lxcols,nxcols)
    #statscolsboot<-c('hspidpositions','hsplength','breakpoints','alignments','pairs') #NO LONGER BOOTSTRAPPING BREAKPOINTS
    #colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','pairs')
    statscolsboot<-c('hspidpositions','hsplength')
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  } else if (breakpoint=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs')
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim','breakpoints','alignments','pairs')
    #statscolsboot<-c('hspidpositions','hsplength','breakpoints','alignments','pairs')
    #colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength','breakpoints','alignments','pairs')
    statscolsboot<-c('hspidpositions','hsplength')
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  } else if (alnlenstats=='True') {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim',lxcols,nxcols)
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim',lxcols,nxcols)
    statscolsboot<-c('hspidpositions','hsplength')
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  } else {
    colnames(allsampledf)<-c('querysample','subjectsample','hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim')
    statscols<-c('hspidpositions','hsplength','qhsplengthpretrim','shsplengthpretrim')
    statscolsboot<-c('hspidpositions','hsplength')
    colnames(allsampledfboot)<-c('bootstrap','querysample','subjectsample','hspidpositions','hsplength')
  }
}

###get final stats (distance scores) for each pairwise sample combination

#aggregate seqlen repot at sample level
seqlenreport$sequence<-sapply(strsplit(as.vector(seqlenreport$sequence),"|",fixed=T),function(x) x=x[1])
sampleseqlen<-aggregate(seqlenreport$length, by=list(seqlenreport$sequence), FUN=sum)
colnames(sampleseqlen)<-c('sample','length')

#get final stats
allsampledf[c(1,2,5,6)]<-t(apply(allsampledf,1,reorderallsampledf)) #reorder query/subject sample alphabetically; reorder qhsplenpretrim/shsplenpretrim accoringly (N.B. the terms query/subject now lose their meaning)
splitsampledf<-split(allsampledf, list(allsampledf$querysample,allsampledf$subjectsample),drop=TRUE) #split by pairwise combination

cl<-makeCluster(as.integer(args[2]))
registerDoParallel(cl)
clusterExport(cl,c("applystatscalc","getseqlens","sampleseqlen","mergestats","statsfunc","breakpoint","alnlenstats","bpdistcalc","lxcols","nxcols"))

finaldf<-parLapply(cl, splitsampledf,applystatscalc, mystatscols=statscols)

stopCluster(cl)


finaldf<-as.data.frame(do.call(rbind,finaldf))

#N.B 'Genome1'/'Genome2' are used for column names of distancestats.tsv/distancestats_bootstrapped.tsv files, but in code in this script the term 'sample' instead of 'genome' is used
if (breakpoint=='True' && alnlenstats=='True') {
  colnames(finaldf)<-c('Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Coverage_breadth_Genome1','Coverage_breadth_Genome2','Breakpoint_distance_d0','Breakpoint_distance_d1','Breakpoints','Alignments','Alignment_pairs',lxcols,nxcols)
} else if (breakpoint=='True') {
  colnames(finaldf)<-c('Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Coverage_breadth_Genome1','Coverage_breadth_Genome2','Breakpoint_distance_d0','Breakpoint_distance_d1','Breakpoints','Alignments','Alignment_pairs')
} else if (alnlenstats=='True') {
  colnames(finaldf)<-c('Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Coverage_breadth_Genome1','Coverage_breadth_Genome2',lxcols,nxcols)
} else {
  colnames(finaldf)<-c('Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Coverage_breadth_Genome1','Coverage_breadth_Genome2')
}

finaldf<-finaldf[with(finaldf,order(finaldf$Genome1,finaldf$Genome2)),]
write.table(finaldf, file=gsubfn('%1',list('%1'=args[1]),'%1/output/distancestats.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)


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

  if (breakpoint=='True' && alnlenstats=='True') {
    colnames(finaldfboot)<-c('bootstrap','Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Breakpoint_distance_d0','Breakpoint_distance_d1','Breakpoints','Alignments','Alignment_pairs',lxcols,nxcols)
  } else if (breakpoint=='True') {
    colnames(finaldfboot)<-c('bootstrap','Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome','Breakpoint_distance_d0','Breakpoint_distance_d1','Breakpoints','Alignments','Alignment_pairs')
  } else if (alnlenstats=='True') {
    colnames(finaldfboot)<-c('bootstrap','Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome',lxcols,nxcols)  
  } else {
    colnames(finaldfboot)<-c('bootstrap','Genome1','Genome2','Genome1_length','Genome2_length','DistanceScore_d0','DistanceScore_d1','DistanceScore_d2','DistanceScore_d3','DistanceScore_d4','DistanceScore_d5','DistanceScore_d6','DistanceScore_d7','DistanceScore_d8','DistanceScore_d9','Percent_identity','Coverage_breadth','Coverage_breadth_mingenome')
  }

  finaldfboot<-finaldfboot[with(finaldfboot,order(as.numeric(finaldfboot$bootstrap),finaldfboot$Genome1,finaldfboot$Genome2)),]
  write.table(finaldfboot, file=gsubfn('%1',list('%1'=args[1]),'%1/output/distancestats_bootstrapped.tsv'), sep='\t', quote=F, col.names=TRUE, row.names=FALSE)

}




#OLD CODE

#incorrect breakpoint calculation function
# BP<-function(x,y) { #this function works on signed permuations (numeric vectors with +/- indicated)
#   out1<-makepairs(x)[!(makepairs(x) %in% makepairs(y) | makepairs(x) %in% makepairs(rev(y*-1)))]
#   out2<-makepairs(x)[unlist(lapply(makepairs(x), function(z) length(unique(sign(z)))))>1]
#   all<-c(out1,out2)
#   return(all[!duplicated(all)]) #deduplicate to aovid double counting breakpoints that occur in out1 and out2
# }
