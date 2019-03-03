args = commandArgs(trailingOnly=TRUE)
library('gsubfn')
library('genoPlotR')
library('dendextend') #prune function 
library('ape') #handling trees
library('tools') #file_ext function  #tools is a base R package
library('ggplot2')
library('cowplot')
inputdir=as.character(args[1]) #alignments/tree input dir
filenamesyntax=as.character(args[2]) #syntax for alignment files and (optionally) for feature annotation files e.g. trimmedalignments_GENOMENAME.tsv,GENOMENAME.gb; GENOMENAME will be replaced
seqlengths=as.character(args[3]) #seqlength file
comparisons=as.character(args[4]) #comparisons file
outdir=as.character(args[5]) #output dir
comparisontype=as.character(args[6]) #'chain' or 'singlereference'
title=as.character(args[7]) #'my title' or NULL ##                                                 
titlepos=as.character(args[8]) #centre   
rightmargin=as.numeric(args[9]) #0.05
sequencecols=as.character(args[10]) #"light yellow,light cyan"
dnaseglabels=as.character(args[11])  #genome1,genome2 OR NULL (default) ##
dnaseglabelcex=as.numeric(args[12]) #1
dnaseglabelcol=as.character(args[13]) #black or black,black ##
dnasegline=as.character(args[14]) #false / black or false,flase / black,black ##
mingapsize=as.numeric(args[15]) #0.02
mainscale=as.character(args[16]) #true/false ##
dnasegscale=as.character(args[17]) #true/false ##
dnasegscalecex=as.numeric(args[18]) #1                 
outputheight=as.character(args[19])
outputwidth=as.character(args[20])
legendorientation=as.character(args[21])
poscolvec=as.character(args[22])
negcolvec=as.character(args[23])
featurespresent=as.character(args[24])
if (featurespresent=='featurespresent') { #there is a feature input directory provided
  featuresdir=as.character(args[25]) #features input dir
  annotationtxtname=as.character(args[26]) #auto product gene
  annotationtxtheight=as.character(args[27]) #auto or numeric
  annotationtxtrot=as.integer(args[28])
  exclusionpresent=as.character(args[29]) #exclusionabsent commandline filepath
  exclusionarg=as.character(args[30])  #placeholder, comma-separated string, path to file
  inclusionpresent=as.character(args[31])
  inclusionarg=as.character(args[32])
  casesensitive=as.character(args[33])
  annotationgenetype=as.character(args[34]) #default: side_bars
  annotationoutlinecol=as.character(args[35]) #default: black
  annotationfillcol=as.character(args[36]) #default: black 
  annotationlty=as.integer(args[37]) #default: 1
  annotationlwd=as.integer(args[38]) #default: 1
  #annotationcex=as.character(args[19]) #default: auto or numeric
} else {
  annotationtxtheight=0
}


annotationcex<-1 #not sure what these do; set to defaults
annotationpch<-8

#handle title
if (title=='NULL') {
  title<-NULL
}

#handle col vecs ###!!!added
poscolvec<-unlist(strsplit(poscolvec,',',fixed=T))
negcolvec<-unlist(strsplit(negcolvec,',',fixed=T))

#handle sequencecols, syntax
sequencecols<-unlist(strsplit(sequencecols,',',fixed=T))

syntaxvec<-unlist(strsplit(filenamesyntax,',',fixed=T))
alignmentsyntax<-syntaxvec[1]
if (length(syntaxvec)>1) {
  annotationsyntax<-syntaxvec[2]
}

#handle dnaseglabels
if (dnaseglabels=='NULL') {
  dnaseglabels<-NULL
} else {
  dnaseglabels<-unlist(strsplit(dnaseglabels,',',fixed=T))
}

#handle dnaseglabelcol
dnaseglabelcol<-unlist(strsplit(dnaseglabelcol,',',fixed=T))

#handle dnasegline
dnasegline<-unlist(strsplit(dnasegline,',',fixed=T))
for (i in 1:length(dnasegline)) {
  myline<-dnasegline[i]
  if (myline=='false') {
    dnasegline[i]<-FALSE
  }
  if (myline=='true') {
    dnasegline[i]<-'black'
  }
}

#handle mainscsale, dnasegscale
if (mainscale=='true') {
  mainscale<-TRUE
}
if (mainscale=='false') {
  mainscale<-FALSE
}
if (dnasegscale=='true') {
  dnasegscale<-TRUE
}
if (dnasegscale=='false') {
  dnasegscale<-FALSE
}

#handle inclusion/exclusion criteria
if (featurespresent=='featurespresent') {
  if (exclusionpresent=='commandline') {
    exclusionarg<-unlist(strsplit(exclusionarg,',',fixed=T))
  }
  if (exclusionpresent=='filepath') {
    exclusionfile<-read.table(exclusionarg,sep='\t',header=FALSE)
    exclusionarg<-as.vector(exclusionfile[,1])
  }

  if (inclusionpresent=='commandline') {
    inclusionarg<-unlist(strsplit(inclusionarg,',',fixed=T))
  }
  if (inclusionpresent=='filepath') {
    inclusionfile<-read.table(inclusionarg,sep='\t',header=FALSE)
    inclusionarg<-as.vector(inclusionfile[,1])
  }
}

#handle annotation settings
if (annotationgenetype=='arrows') {
  annotationgenetype<-'arrowsgrob'
}

##function definitions

getxlims<-function(x) {
  xlims<-vector()
  for (i in 1:length(x)) {
    if (i==1) {
      start<-1
      end<-x[1]
      xlims<-c(start,end)
    } else {
      start<-sum(x[c(1:i-1)])+1
      end<-sum(x[c(1:i)])
      xlims<-c(xlims,start,end)
    }
  }
  return(xlims)
}


pastefunction<-function(x) {
  if (length(x)==1) {
    x=x[1]
  } else {
    x=paste(x[1],x[2],sep='|')
  }
}


cellpresent<-function(x) {
  if (is.null(x) || is.na(x)) {
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}


makepairs<-function(x) mapply(c, head(x,-1), tail(x,-1), SIMPLIFY = FALSE)


getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


applyexclincl<-function(annotvec,terms,casesen,exclincl) {
  newannotvec<-annotvec
  for (i in 1:length(annotvec)) {
    annot<-annotvec[i]
    if (casesen!='True') {
      annot<-tolower(annot)
    }
    for (term in terms) {
      if (casesen!='True') {
        term<-tolower(term)
      }
      if (exclincl=='exclusion') {
        if (grepl(term,annot,fixed = TRUE)==TRUE) {
          newannotvec[i]<-''
	  break
        }
      } else {
        if (grepl(term,annot,fixed = TRUE)==FALSE) {
          newannotvec[i]<-''
        } else {
	  newannotvec[i]<-annotvec[i]
          break
        }
      }
    }
  }
  return(newannotvec)
}


reformatgff<-function(gff,annotationtagname=annotationtxtname,defaultoutlinecol=annotationoutlinecol,defaultgenetype=annotationgenetype,defaultfillcol=annotationfillcol,defaultlty=annotationlty,defaultlwd=annotationlwd,defaultpch=annotationpch,defaultcex=annotationcex,exclusioncriteria=exclusionarg,inclusioncriteria=inclusionarg,casesen=casesensitive) {
  #get annotation text
  if (annotationtagname=='auto') {
    name1<-getAttributeField(gff$attributes,'gene')
    name2<-getAttributeField(gff$attributes,'product')
    name<-name1
    name[is.na(name1)]<-name2[is.na(name1)]
  } else {
    name<-getAttributeField(gff$attributes,annotationtagname)    
  }
  name[is.na(name)]<-''
  name<-as.character(name)
  #apply inclusion/exlusion criteria
  if (exclusioncriteria[1]!='placeholder') {
    name<-applyexclincl(name,exclusioncriteria,casesen,'exclusion')
  }
  if (inclusioncriteria[1]!='placeholder') {
    name<-applyexclincl(name,inclusioncriteria,casesen,'inclusion')
  }
  start<-as.integer(gff$start)
  end<-as.integer(gff$end)
  strand<-as.character(gff$strand)
  strand[strand=='+']<-'1'
  strand[strand=='-']<-'-1'
  strand<-as.integer(strand)
  col<-getAttributeField(gff$attributes,'outline')
  col[is.na(col)]<-defaultoutlinecol
  col<-as.character(col)
  gene_type<-getAttributeField(gff$attributes,'gene_type')
  gene_type[is.na(gene_type)]<-defaultgenetype
  gene_type[gene_type=='arrows']<-'arrowsgrob'
  gene_type<-as.character(gene_type)
  fill<-getAttributeField(gff$attributes,'fill')
  fill[is.na(fill)]<-defaultfillcol
  fill<-as.character(fill)
  lty<-getAttributeField(gff$attributes,'lty')
  lty[is.na(lty)]<-defaultlty
  lwd<-getAttributeField(gff$attributes,'lwd')
  lwd[is.na(lwd)]<-defaultlwd
  pch<-getAttributeField(gff$attributes,'pch')
  pch[is.na(pch)]<-defaultpch
  cex<-getAttributeField(gff$attributes,'cex')
  cex[is.na(cex)]<-defaultcex
  #make genoplotr style annotation dataframe
  reformatteddf<-data.frame(name,start,end,strand,col,gene_type,fill,lty,lwd,pch,cex,stringsAsFactors = F)
  reformatteddf$lty<-as.integer(reformatteddf$lty)
  reformatteddf$lwd<-as.integer(reformatteddf$lwd)
  reformatteddf$pch<-as.integer(reformatteddf$pch)
  reformatteddf$cex<-as.integer(reformatteddf$cex)
  return(reformatteddf)
}



removeparentannotations<-function(gff) {
  childindices<-which(is.na(getAttributeField(gff$attributes,'Parent'))==F && is.na(getAttributeField(gff$attributes,'ID'))==T) #rows with parent attribute but no id (genbank location parts)
  if (length(childindices)>0) {
    ids<-getAttributeField(gff$attributes,'ID')
    childids<-getAttributeField(gff[childindices,]$attributes,'Parent')

    indicestoremove<-vector()
    for (i in 1:length(childids)) {
      childid<-childids[i]
      parentindex<-which(ids==childid)
      indicestoremove<-c(indicestoremove,parentindex)
      parentattribute<-gff[parentindex,]$attribute
      gff[childindices[i],]$attributes<-parentattribute #replace parent flag with parent attribute
    }
    gff<-gff[-unique(indicestoremove),]
    return(gff)
  } else {
    return(gff)
  }
}



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


###colour scale functions ###!!ADDED
library(ggplot2)
library(cowplot)

getbreaks<-function(pidvalues, lengthout) {
  lowpid<-floor(min(pidvalues))
  highpid<-ceiling(max(pidvalues))
  breaks<-seq(lowpid,highpid,length.out = lengthout)
  diff<-breaks[2]-breaks[1]
  if (diff>1.5) {
    breaks2<-round(breaks)
  } else {
    breaks2<-round(breaks,2)
  }
  return(breaks2)
}

getlegendplot<-function(comparison,indices,mypids,colvec,legendorientation,legendtitle) {
  dummydf<-data.frame(1:length(mypids),mypids)
  colnames(dummydf)<-c('indices','PercentIdentity')
  p<-ggplot(dummydf,aes(x=indices,y=PercentIdentity,colour=PercentIdentity)) + geom_point()
  if (length(colvec)>2) {
    p2<-p + scale_colour_gradientn(colours=colvec,limits = c(floor(min(comparison$pid)),ceiling(max(comparison$pid))), 
                                   breaks = getbreaks(comparison$pid,10)) + labs(colour="Percent\nIdentity") + theme(legend.text.align = 1, legend.text = element_text(size=rel(0.7)),legend.title = element_text(size=rel(0.7)),legend.position = c(0.5,0.5),legend.key.size = unit(30,"pt")) #legend.margin = margin(0,0,0,0),legend.box.margin = margin(0,0,0,0)  #legend.position = legendpos
  } else {
    p2<-p + scale_colour_gradient(low=colvec[1],high=colvec[2],limits = c(floor(min(comparison$pid)),ceiling(max(comparison$pid))), 
                                  breaks = getbreaks(comparison$pid,10)) + labs(colour="Percent\nIdentity") + theme(legend.text.align = 1,legend.text = element_text(size=rel(0.7)),legend.title=element_text(size=rel(0.7)),legend.position = c(0.5,0.5),legend.key.size = unit(30,"pt"))
  }
  if (legendorientation=='horizontal') {
    p2<-p2+theme(legend.direction = "horizontal", legend.text = element_text(angle = 90))
  }
  return(p2)
}



getlegendcols<-function(comparison,legendorientation,poscolvec,negcolvec,colonly=FALSE) {
  posindices<-which(comparison$strand==1)
  negindices<-which(comparison$strand==-1)
  finalcols<-numeric(length(comparison$pid))
  poslegend<-NULL
  neglegend<-NULL
  #positive alignment colours and legend
  if (length(posindices)>0) {
    mypids<-comparison$pid[posindices]
    poslegendplot<-getlegendplot(comparison,posindices,mypids,poscolvec,legendorientation)
    poscols<-ggplot_build(poslegendplot)$data[[1]]$colour
    finalcols[posindices]<-poscols
    if (colonly==FALSE) {
      poslegendplot<-poslegendplot + labs(colour="Percent identity\n(+ve orientation)")
      poslegend<-get_legend(poslegendplot)
    }
  }
  if (length(negindices)>0) {
    mypids<-comparison$pid[negindices]
    neglegendplot<-getlegendplot(comparison,negindices,mypids,negcolvec,legendorientation)
    negcols<-ggplot_build(neglegendplot)$data[[1]]$colour
    finalcols[negindices]<-negcols
    if (colonly==FALSE) {
      neglegendplot<-neglegendplot + labs(colour="Percent identity\n(-ve orientation)")
      neglegend<-get_legend(neglegendplot)
    }
  }
  return(list(poslegend,neglegend,finalcols))
}


makecoltransparent<-function(colvec,alphatransparency) {
  rgbcollist<-lapply(colvec,col2rgb)
  transparentcolvec<-unlist(lapply(rgbcollist,function(x) rgb(x[1],x[2],x[3],alpha = alphatransparency,maxColorValue = 255)))
  return(transparentcolvec)
}

poscolvec<-makecoltransparent(poscolvec,alphatransparency = 130)
negcolvec<-makecoltransparent(negcolvec,alphatransparency = 130)



assignInNamespace("colour_ramp", function(colors, na.color = NA, alpha = TRUE){
  if (length(colors) == 0) {
    stop("Must provide at least one color to create a color ramp")
  }
  colorMatrix <- grDevices::col2rgb(colors, alpha = alpha)
  structure(function(x) {
    scales:::doColorRamp(colorMatrix, x, alpha, ifelse(is.na(na.color), 
                                              "", na.color))
  }, safe_palette_func = TRUE)
}, "scales")



###annotation legend functions

genetypetoshape<-function(genetype) {
  blockindices<-which(genetype=='blocks')
  arrowindices<-which(genetype=='arrowsgrob')
  genetype[blockindices]<-22
  genetype[arrowindices]<-24
  if (length(c(blockindices,arrowindices))==0) {
    genetype[1:length(genetype)]<-21
  } else {
    genetype[-c(blockindices,arrowindices)]<-21
  }
  return(as.numeric(genetype))
}


getallannotations<-function(annotations) {  #gets annotations to be plotted on legend (i.e. annottions that have customised outline/fill/gene_type); need to deduplicate by annotation name prior to plotting
  annotationtext<-annotations$name
  annotationindices<-annotations$col!=annotationoutlinecol | annotations$gene_type!=annotationgenetype | annotations$fill!=annotationfillcol
  annotations<-annotations[annotationindices,]
  annotationtext<-annotationtext[annotationindices]
  return(list(annotations,annotationtext))
}

getannotationlegendplot<-function(annotation,annotationtext) {
  dummydf<-data.frame(column1=rep(1,length(annotationtext)),column2=rep(1,length(annotationtext)),gene_type=annotation$gene_type,fill=annotation$fill,outline=annotation$col,annotationtext,stringsAsFactors = F)
  p<-ggplot(dummydf,aes(x=column1,y=column2,fill=as.character(1:nrow(dummydf)),colour=as.character(1:nrow(dummydf)))) + geom_point() +  
    scale_colour_manual(values=dummydf$outline,labels=annotationtext) + 
    scale_fill_manual(values=dummydf$fill,labels=annotationtext) + 
    guides(colour=guide_legend(override.aes=list(shape=genetypetoshape(dummydf$gene_type)))) +  #doesn't matter whether colour or fill legend is overridden 
    labs(fill="Gene annotations",colour="Gene annotations") + 
    theme(legend.text.align = NULL,legend.text = element_text(size=rel(0.7)),legend.title=element_text(size=rel(0.7)),legend.position = c(0.5,0.5),legend.key.size = unit(30,"pt"))
  annotationlegend<-get_legend(p)
  return(annotationlegend)
}


###



#custom grobs

arrow_coord <- function(x1, x2, y=0.5, strand=NULL, width=1, head_len=100){
  # take care of strand, to get x1 as bottom and x2 as tip of arrow
  if (!is.null(strand) && strand == -1){
    x_temp <- x2
    x2 <- x1
    x1 <- x_temp
  }
  w2 <- width/4
  # if the head of the arrow is larger than half of the gene, reduce to half
  #if (head_len > abs(x1-x2)/2){
  #  head_len <- abs(x1-x2)/2
  #}
  head_len<-abs(x1-x2)  #this means arrows will be directed triangles
  # calculate xi, x "internal"
  if (x2 > x1){
    xi <- x2-head_len
  } else {
    xi <- x2+head_len
  }
  list(x=c(x1,   xi,   xi,     x2, xi,     xi,   x1),
       y=c(y-w2, y-w2, y-w2*2, y,  y+w2*2, y+w2, y+w2)
  )
}

arrowsgrob <- function(gene, head_len=Inf, i=0, grob_cex=1, ...){
  if (!is.dna_seg(gene)) stop("A dna_seg object is required")
  if (nrow(gene) > 1) stop ("gene must be single-row")
  mid <- (gene$start + gene$end)/2
  name <- paste("seg.", i, ".", gene$name, sep="")
  color <- gene$col
  fill <- gene$fill
  if (is.null(fill)) fill <- color
  if (is.null(color)) color <- fill
  arrow <- arrow_coord(x1=gene$start, x2=gene$end,y=0.5, strand=gene$strand, head_len=head_len)
  scalevec<-((arrow$y-0.5)*grob_cex)-(arrow$y-0.5) #grob cex gives height of arrow; the alignments don't shift based on annotation grob height so not recommended to use cex >1
  polygonGrob(arrow$x, arrow$y+scalevec, name=name,gp=gpar(col=color, fill=fill, lty=gene$lty,lwd=gene$lwd),default.units="native")
}



##

###read sequence lengths and make sure they are ordered alphabetically to match rangeshifting order (because for each genome, contigs were shifted in alphabetical order)
seqlengths<-read.table(seqlengths,sep='\t',header=FALSE)
colnames(seqlengths)<-c('names','lengths')
seqlengths$names<-sapply(strsplit(as.vector(seqlengths$names),'|',fixed=T),pastefunction)
factorlevels<-as.vector(levels(as.factor(seqlengths$names)))
seqlengths<-seqlengths[with(seqlengths,order(seqlengths$names)),]
seqlengthsreport<-seqlengths #copying for rangeshifting function
seqlengths<-cbind(as.data.frame(t(sapply(strsplit(as.vector(seqlengths$names),'|',fixed=T), function(x) c(x[1],x[2])))),seqlengths$lengths)
colnames(seqlengths)<-c('names','contigs','lengths')



###read a file with each row containing a column of query(queries) and a subject column
comparisonfile=read.table(comparisons,sep='\t',header=FALSE)
numcols<-ncol(comparisonfile)
cols=c('subject','queries','outputname','subjectannotation','queryannotation','subjecttxtannotation','querytxtannotation','treename')
cols<-cols[1:numcols]
colnames(comparisonfile)<-cols




for (i in 1:nrow(comparisonfile)) {
  comparison=comparisonfile[i,]
  outputname<-as.character(comparison$outputname)
  subject<-as.character(comparison$subject)
  queries<-as.character(comparison$queries)
  queryvec<-unlist(strsplit(queries,',',fixed=T))
  #initialise annotation lists for legend
  annotlgndcounter=0
  alllegendannotationtexts<-list()
  alllegendannotations<-list()
  ###get xlims
  xlim_subject<-getxlims(seqlengths[seqlengths$names==subject,3])  #seqlengths must be ordered by contig in the same way as in rangeshifting
  xlim_queries<-vector("list",length(queryvec))
  for (j in 1:length(queryvec)) {
    query<-queryvec[j]
    xlim_query<-getxlims(seqlengths[seqlengths$names==query,3]) #if there are multiple contigs, multiple corresponding ranges will be produced
    xlim_queries[[j]]<-xlim_query  
  }
  xlims<-c(list(xlim_subject),xlim_queries) #combine subject with query list
  
  
  ###get dna_segs based on xlims; loop through query vec to get list of query dna_segs; if annotations are available, append to the dna_seg
  
  #get subject dna_seg
  subject_matrix<-matrix(xlim_subject,ncol=2,byrow=T)
  nseq<-nrow(subject_matrix)
  myfill<-rep(sequencecols,nseq)[1:nseq]
  subject_seq<-data.frame(name=seqlengths[seqlengths$names==subject,2], #using contigs for names
                          start=subject_matrix[,1],
                          end=subject_matrix[,2],
                          strand=rep(1,nseq),
                          col=rep("black",nseq),
                          gene_type=rep("blocks",nseq),
                          #fill=rep("light gray",nseq),  #!!!this should be alternating to distinguish adjacent contigs
                          fill=myfill,
                          lty=rep(1,nseq),
                          lwd=rep(1,nseq),
                          pch=rep(8,nseq),
                          cex=rep(1,nseq))
  cellpresentbool<-cellpresent(comparison[i,4]) #subject annotations provided?
  if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
    annotationsfile<-as.character(comparison$subjectannotation)
    annotationsfile<-gsub("GENOMENAME", subject, annotationsyntax) #replaces GENOMENAME placeholder with subject genome name
    #subject_annotations<-read_dna_seg_from_file(gsubfn('%1|%2',list('%1'=featuresdir,'%2'=annotationsfile),'%1/%2'),gene_type = annotationgenetype)
    #subject_annotations<-as.data.frame(subject_annotations)
    #subject_annotations<-subject_annotations[c(annotationtxtname,"start","end","strand","col","gene_type","fill","lty","lwd","pch","cex")]
    #colnames(subject_annotations)<-c("name","start","end","strand","col","gene_type","fill","lty","lwd","pch","cex")
    #subject_annotations<-subject_annotations[subject_annotations$gene_type!='introns',]
    sbj_gff<-read.gff(gsubfn('%1|%2',list('%1'=featuresdir,'%2'=annotationsfile),'%1/%2'))
    sbj_gff<-sbj_gff[sbj_gff$type=='CDS',]
    sbj_gff<-removeparentannotations(sbj_gff)
    ##shift ranges if there are contigs
    reformattednames<-sapply(strsplit(as.vector(sbj_gff$seqid),'|',fixed=T),pastefunction)
    addlens<-rangeshifting(reformattednames,seqlengthsreport,seqlengthsreport$names)
    sbj_gff$start<-sbj_gff$start+addlens
    sbj_gff$end<-sbj_gff$end+addlens
    ##
    subject_annotations<-reformatgff(sbj_gff)
    subject_seg<-rbind(subject_seq,subject_annotations)
    subject_seg<-dna_seg(subject_seg)
    #get annotation text - for legend
    annotlgndcounter<-annotlgndcounter+1
    legendannotationlist<-getallannotations(subject_annotations)
    legendannotation<-legendannotationlist[[1]];legendtext<-legendannotationlist[[2]]
    alllegendannotations[[annotlgndcounter]]<-legendannotation;alllegendannotationtexts[[annotlgndcounter]]<-legendtext
  } else {
    subject_seg<-dna_seg(subject_seq)
  }


  #get annotation text - for genoplotr plot; if there are no annotations, annotation will be blank ("")
  cellpresentbool<-cellpresent(comparison[i,6]) #plot subject text annotations?
  if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
    mid_pos<-middle(subject_seg)[1:nrow(subject_seg)]
    annottext<-subject_seg$name #need to overwrite [1] with blank
    annottext[1:nrow(subject_seq)]<-""
    subjectannot<-annotation(x1=mid_pos,text=annottext,rot = annotationtxtrot)
  } else {
    subjectannot<-annotation(x1=1,text="") 
  }
  
  #get query dna_seg(s)
  query_segs<-vector("list",length(queryvec))
  query_annotation_texts<-vector("list",length(queryvec))
  
  #get queryannotation/querytxtannotation genome names (indicating which should be plotted)
  cellpresentbool<-cellpresent(comparison[i,5]) #query annotations provided?
  if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
    queryannotationfiles<-as.character(comparison$queryannotation)
    queryannotationfiles<-unlist(strsplit(queryannotationfiles,',',fixed=T))
    queryannotationfiles[queryannotationfiles=='-']<-''
  }
  cellpresentbool<-cellpresent(comparison[i,7]) #plot query text annotations?
  if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
    querytxtannotationfiles<-as.character(comparison$querytxtannotation)
    querytxtannotationfiles<-unlist(strsplit(querytxtannotationfiles,',',fixed=T))
    querytxtannotationfiles[querytxtannotationfiles=='-']<-''
  }
  #
  for (j in 1:length(queryvec)) {
    query<-queryvec[j]
    #query_matrix<-matrix(xlim_query,ncol=2,byrow=T) #!bug
    query_matrix<-matrix(xlim_queries[[j]],ncol = 2,byrow = T)
    nseq<-nrow(query_matrix)
    myfill<-rep(sequencecols,nseq)[1:nseq]
    query_seq<-data.frame(name=seqlengths[seqlengths$names==query,2], #using contigs for names
                          start=query_matrix[,1],
                          end=query_matrix[,2],
                          strand=rep(1,nseq),
                          col=rep("black",nseq),
                          gene_type=rep("blocks",nseq),
                          #fill=rep("light gray",nseq), #!!!this should be alternating to distinguish adjacent contigs
                          fill=myfill,
                          lty=rep(1,nseq),
                          lwd=rep(1,nseq),
                          pch=rep(8,nseq),
                          cex=rep(1,nseq))
    cellpresentbool<-cellpresent(comparison[i,5]) #query annotations provided?
    if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
      if (nchar(queryannotationfiles[j])>0) { #annotations to plot are provided in a comma-sep string, so a blank (nchar=0) means don't plot annotations for this query
        annotationsfile<-queryannotationfiles[j]
        annotationsfile<-gsub("GENOMENAME", query, annotationsyntax) #replaces GENOMENAME placeholder with query genome name
        #query_annotations<-read_dna_seg_from_file(gsubfn('%1|%2',list('%1'=featuresdir,'%2'=annotationsfile),'%1/%2'),gene_type = annotationgenetype)
        #query_annotations<-as.data.frame(query_annotations)
        #query_annotations<-query_annotations[c(annotationtxtname,"start","end","strand","col","gene_type","fill","lty","lwd","pch","cex")]
        #colnames(query_annotations)<-c("name","start","end","strand","col","gene_type","fill","lty","lwd","pch","cex")
        #query_annotations<-query_annotations[query_annotations$gene_type!='introns',]
        qry_gff<-read.gff(gsubfn('%1|%2',list('%1'=featuresdir,'%2'=annotationsfile),'%1/%2'))
        qry_gff<-qry_gff[qry_gff$type=='CDS',]
        qry_gff<-removeparentannotations(qry_gff)
        ##shift ranges if there are contigs
        reformattednames<-sapply(strsplit(as.vector(qry_gff$seqid),'|',fixed=T),pastefunction)
        addlens<-rangeshifting(reformattednames,seqlengthsreport,seqlengthsreport$names)
        qry_gff$start<-qry_gff$start+addlens
        qry_gff$end<-qry_gff$end+addlens
        ##
        query_annotations<-reformatgff(qry_gff)
        query_seg<-rbind(query_seq,query_annotations)
        query_seg<-dna_seg(query_seg)
        query_segs[[j]]<-query_seg
	#get annotation text - for legend
        annotlgndcounter<-annotlgndcounter+1
        legendannotationlist<-getallannotations(query_annotations)
        legendannotation<-legendannotationlist[[1]];legendtext<-legendannotationlist[[2]]
        alllegendannotations[[annotlgndcounter]]<-legendannotation;alllegendannotationtexts[[annotlgndcounter]]<-legendtext
      } else {
        query_seg<-dna_seg(query_seq)
        query_segs[[j]]<-query_seg	
      }
    } else {
      query_seg<-dna_seg(query_seq)
      query_segs[[j]]<-query_seg	
    }
    
    #get annotation text - for genoplotr plot; if there are no annotations this will be blank ("")
    cellpresentbool<-cellpresent(comparison[i,7]) #plot query text annotations?
    if (cellpresentbool==TRUE && featurespresent=='featurespresent') {
      if (nchar(querytxtannotationfiles[j])>0) { #annotations to plot are provided in a comma-sep string, so a blank (nchar=0) means don't plot annotations for this query
        mid_pos<-middle(query_seg)[1:nrow(query_seg)]  ##!N.B you can't plot annotation text if you don't have annotation (text will just be the blank(s) corresponding to contigs)
        annottext<-query_seg$name #need to overwrite [1] with blank
        annottext[1:nrow(query_seq)]<-""
        queryannot<-annotation(x1=mid_pos,text=annottext,rot = annotationtxtrot)
        query_annotation_texts[[j]]<-queryannot
      } else {
        queryannot<-annotation(x1=1,text="")
        query_annotation_texts[[j]]<-queryannot
      }
    } else {
      queryannot<-annotation(x1=1,text="")
      query_annotation_texts[[j]]<-queryannot
    }
  }
  dna_segs<-c(list(subject_seg),query_segs)
  names(dna_segs)<-c(subject,queryvec)
  annotationtexts<-c(list(subjectannot),query_annotation_texts)


  ###plot annotation legend
  alllegendannotations<-do.call(rbind,alllegendannotations)
  alllegendannotationtexts<-do.call(c,alllegendannotationtexts) 
  if (nrow(alllegendannotations)>0) {
    #deduplicate annotation data - currently data include all annotations with customised outline/fill/gene_type
    uniqueannotationindices<-!duplicated(alllegendannotationtexts)
    alllegendannotationtexts<-alllegendannotationtexts[uniqueannotationindices]
    alllegendannotations<-alllegendannotations[uniqueannotationindices,]
    #remove annotations that have blank text i.e. removed by inclusion/exclusion criteria
    nonblankindices<-!nchar(alllegendannotationtexts)==0
    alllegendannotationtexts<-alllegendannotationtexts[nonblankindices]
    alllegendannotations<-alllegendannotations[nonblankindices,]
    if (nrow(alllegendannotations)>0) {
      #plot
      writefilepath=gsubfn('%1|%2', list('%1'=outdir,'%2'=outputname), '%1/%2_annotationlegend.pdf')
      pdf(writefilepath)
      mylegend<-getannotationlegendplot(alllegendannotations,alllegendannotationtexts)
      plot(plot_grid(mylegend))
      dev.off()
    }
  }

  ###get tree; input can be a nexus or newick tree or a phylo object saved as an rds file
  treepresentbool<-cellpresent(comparison[i,8])
  if (treepresentbool==TRUE) {
    treename<-as.character(comparison[i,8])
    treepath<-gsubfn('%1|%2',list('%1'=inputdir,'%2'=treename),'%1/%2')
    extension<-file_ext(treename)
    if (extension=='RDS' || extension=='rds') { #read phylo object saved as rds
      mytree<-readRDS(treepath)
      treeclass<-class(mytree)
      if (treeclass!="phylo") {
        errormessage<-gsubfn('%1',list('%1'=treeclass),'Error: tree object is of class %1; class phylo is required')
        stop(errormessage)
      }
    }
    else { #read nexus or newick
      mytree<-try(read.tree(treepath))
      if (class(mytree)=="try-error") {
        mytree<-try(read.nexus(treepath))
        if (class(mytree)=="try-error") {
          stop('Error: unrecognised tree format')
        }
      }
    }
    ##if necessary, prune tree
    genomestoprune<-mytree$tip.label[!mytree$tip.label %in% c(subject,queryvec)]
    if (length(genomestoprune)>0) {
      mytree<-prune(mytree,genomestoprune)
    }
    
    #convert phylo tree to newick; then convert newick to phylog
    mynewick<-write.tree(mytree, file = "")
    mytree<-newick2phylog(mynewick)
  }
  
  
  #N.B until now, haven't used alignment data - need to subet according to query !also if doing chain of comparisons (sample A-B-C...) need to open different alignment data files
  ###read comparison alignment data
  if (comparisontype=='chain') {
    #comparison - chained reference - loop through pairs of genomes, starting with subject vs query1 pair
    #make pairs of genomes to loop through
    genomepairs<-makepairs(c(subject,queryvec)) #list of pairs
    comparisons<-vector("list",length(genomepairs))
    for (z in 1:length(genomepairs)) {
      pair<-genomepairs[[z]]
      sbj<-pair[1]
      qry<-pair[2]
      filename<-gsub("GENOMENAME", sbj, alignmentsyntax) #replaces GENOMENAME placeholder with subject genome name
      report<-read.table(gsubfn('%1|%2',list('%1'=inputdir,'%2'=filename),'%1/%2'),sep='\t',header=TRUE) #read subject alignment file
      report<-report[report$qname==qry,] #filter subject alignment file by query
      mystrand<-as.vector(report$strand)
      mystrand[mystrand=='+']<-1
      mystrand[mystrand=='-']<--1
      mystrand<-as.numeric(mystrand)
      compdf<-report[,c("sstart","send","qstart","qend","pid")]
      compdf<-data.frame(compdf,as.data.frame(mystrand))
      colnames(compdf)<-c('start1', 'end1', 'start2', 'end2','pid','strand')
      comparison<-comparison(compdf)
      #comparison$col<-apply_color_scheme(x=rep(0.5,length(comparison$strand)),direction=comparison$strand,"red_blue")
      #comparison$col<-apply_color_scheme(x=comparison$pid,direction=comparison$strand,"red_blue",rng=c(60,100))
      legendcols<-getlegendcols(comparison,legendorientation,poscolvec,negcolvec,colonly = TRUE)
      #comparison$col<-makecoltransparent(colvec=legendcols[[3]],alphatransparency=130)
      comparison$col<-legendcols[[3]]
      comparisons[[z]]<-comparison
    }
    #add right hand side margin
    xlimmaxindex<-which.max(sapply(xlims, function(x) max(x)))
    xlimmaxindexlast<-length(xlims[[xlimmaxindex]])
    newxlimvalue<-xlims[[xlimmaxindex]][xlimmaxindexlast]+(xlims[[xlimmaxindex]][xlimmaxindexlast]*rightmargin)
    xlims[[xlimmaxindex]][xlimmaxindexlast]<-newxlimvalue
    #calculate output dimensions, and from this calculate annotation cex
    if (annotationtxtheight=='auto') {
      annotationtxtheight<-max(nchar(do.call(rbind,annotationtexts)$text))/3.5
    } else {
      annotationtxtheight<-as.numeric(annotationtxtheight)
    }
    
    if (outputwidth=='auto') {
      outputwidth<-((max(unlist(xlims))/5000)/3)/3
      if (outputwidth<10) {
        outputwidth<-10
      }
    } else {
      outputwidth<-as.numeric(outputwidth)
    }
    if (outputheight=='auto') {
      outputheight<-(((annotationtxtheight/3)*5)+(length(dna_segs)*3))/3
      if (outputheight<10) {
        outputheight<-10
      }
    } else {
      outputheight<-as.numeric(outputheight)
    }
    
    print('dimensions: height, width')
    print(outputheight)
    print(outputwidth)
    
    #
    writefilepath=gsubfn('%1|%2', list('%1'=outdir,'%2'=outputname), '%1/%2.pdf')
    pdf(writefilepath,outputwidth,outputheight)
    
    if (treepresentbool==TRUE) {
      plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,xlims = xlims,tree=mytree,annotations = annotationtexts,annotation_height = annotationtxtheight/1.75,annotation_cex = 0.5,main=title,main_pos=titlepos,dna_seg_labels=dnaseglabels,dna_seg_label_cex=dnaseglabelcex,dna_seg_label_col=dnaseglabelcol,dna_seg_line=dnasegline,min_gap_size=mingapsize,scale=mainscale,dna_seg_scale=dnasegscale,scale_cex=dnasegscalecex)
    } else {
      plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,xlims = xlims,annotations = annotationtexts,annotation_height = annotationtxtheight/1.75,annotation_cex = 0.5,main=title,main_pos=titlepos,dna_seg_labels=dnaseglabels,dna_seg_label_cex=dnaseglabelcex,dna_seg_label_col=dnaseglabelcol,dna_seg_line=dnasegline,min_gap_size=mingapsize,scale=mainscale,dna_seg_scale=dnasegscale,scale_cex=dnasegscalecex)
    }
    dev.off()
    
    #plot legends
    writefilepath=gsubfn('%1|%2', list('%1'=outdir,'%2'=outputname), '%1/%2_alignmentlegend.pdf')
    pdf(writefilepath)
    allcomparisons<-do.call(rbind,comparisons)
    mylegends<-getlegendcols(allcomparisons,legendorientation,poscolvec,negcolvec,colonly=FALSE)
    poslegend<-mylegends[[1]];neglegend<-mylegends[[2]]
    if (legendorientation=='vertical') {
      print(plot_grid(poslegend,neglegend,nrow=2))
    } else {
      print(plot_grid(poslegend,neglegend,ncol=2))
    }
    dev.off()
  }
  
  
  if (comparisontype=='singlereference') {
    #comparison - single reference - need to overwrite xlims so that they are set to limits of the longest genome in the comparison
    #reset xlims
    maxxlim<-1
    for (xlim in xlims) {
      maxxlim<-max(maxxlim,xlim)
    }
    
    xlims2<-vector("list",length(xlims)) #problem - by using maxxlim it means genomes can't be shifted to the middle for optimal display - this is a limitation I'll accept since other options don't seem to work
    for (z in 1:length(xlims)) {
      xlims2[[z]]<-c(1,maxxlim)
    }
    
    filename<-gsub("GENOMENAME", subject, alignmentsyntax) #replaces GENOMENAME placeholder with subject genome name
    report<-read.table(gsubfn('%1|%2',list('%1'=inputdir,'%2'=filename),'%1/%2'),sep='\t',header=TRUE) #read subject alignment file
    comparisons<-vector("list",length(queryvec))
    for (z in 1:length(queryvec)) {
      query<-queryvec[z]
      qreport<-report[report$qname==query,]
      mystrand<-as.vector(qreport$strand)
      mystrand[mystrand=='+']<-1
      mystrand[mystrand=='-']<--1
      mystrand<-as.numeric(mystrand)
      compdf<-qreport[,c("sstart","send","qstart","qend","pid")]
      compdf<-data.frame(compdf,as.data.frame(mystrand))
      colnames(compdf)<-c('start1', 'end1', 'start2', 'end2','pid','strand')
      comparison<-comparison(compdf)
      #comparison$col<-apply_color_scheme(x=rep(0.5,length(comparison$strand)),direction=comparison$strand,"red_blue")
      #comparison$col<-apply_color_scheme(x=comparison$pid,direction=comparison$strand,"red_blue",rng=c(60,100))
      legendcols<-getlegendcols(comparison,legendorientation,poscolvec,negcolvec,colonly = TRUE)
      #comparison$col<-makecoltransparent(colvec=legendcols[[3]],alphatransparency=130)
      comparison$col<-legendcols[[3]]
      comparisons[[z]]<-comparison
    }
    
    #instead of visualising contigs using xlims, visualise contigs using dna_segs (need to change feature type from line) so contigs can be viewed
    #dna_segs2<-dna_segs
    #for (z in 1:length(dna_segs2)) {
    #  dna_segs2[[z]]$gene_type<-"arrows"
    #}
    #sequencecols<-c("light gray","light blue")
    #sequencecols<-rep(sequencecols,length(dna_segs2))[1:length(dna_segs2)]
    #for (i in 1:length(dna_segs2)) {
    #  dna_segs2[[z]][1,]$gene_type<-"blocks"
    #  dna_segs2[[z]][1,]$fill<-sequencecols[z]
    #}
    
    #use offsets to ensure correct layout
    offsets<-rep(0,length(dna_segs))
    
    #add right hand side margin
    xlimmaxindex<-which.max(sapply(xlims2, function(x) max(x)))
    xlimmaxindexlast<-length(xlims2[[xlimmaxindex]])
    newxlimvalue<-xlims2[[xlimmaxindex]][xlimmaxindexlast]+(xlims2[[xlimmaxindex]][xlimmaxindexlast]*rightmargin)
    xlims2[[xlimmaxindex]][xlimmaxindexlast]<-newxlimvalue
    
    #calculate output dimensions, and from this calculate annotation cex
    if (annotationtxtheight=='auto') {
      annotationtxtheight<-max(nchar(do.call(rbind,annotationtexts)$text))/3.5
    } else {
      annotationtxtheight<-as.numeric(annotationtxtheight)
    }
    
    if (outputwidth=='auto') {
      outputwidth<-((max(unlist(xlims))/5000)/3)/3
      if (outputwidth<10) {
        outputwidth<-10
      }
    } else {
      outputwidth<-as.numeric(outputwidth)
    }
    if (outputheight=='auto') {
      outputheight<-(((annotationtxtheight/3)*5)+(length(dna_segs)*3))/3
      if (outputheight<10) {
        outputheight<-10
      }
    } else {
      outputheight<-as.numeric(outputheight)
    }
    
    print('dimensions: height, width')
    print(outputheight)
    print(outputwidth)
    #
    writefilepath=gsubfn('%1|%2', list('%1'=outdir,'%2'=outputname), '%1/%2.pdf')
    pdf(writefilepath,outputwidth,outputheight)
    
    if (treepresentbool==TRUE) {
      #plot_gene_map(dna_segs = dna_segs2, comparisons = comparisons,xlims = xlims2,tree=mytree,annotations = annotationtexts,annotation_height = max(nchar(do.call(rbind,annotationtexts)$text))/4,annotation_cex = 0.5,plot_new = T)
      plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,xlims = xlims2,tree=mytree,annotations = annotationtexts,annotation_height = annotationtxtheight/1.75,annotation_cex = 0.5,offsets=offsets,main=title,main_pos=titlepos,dna_seg_labels=dnaseglabels,dna_seg_label_cex=dnaseglabelcex,dna_seg_label_col=dnaseglabelcol,dna_seg_line=dnasegline,min_gap_size=mingapsize,scale=mainscale,dna_seg_scale=dnasegscale,scale_cex=dnasegscalecex)
    } else {
      #plot_gene_map(dna_segs = dna_segs2, comparisons = comparisons,xlims = xlims2,annotations = annotationtexts,annotation_height = max(nchar(do.call(rbind,annotationtexts)$text))/4,annotation_cex = 0.5,plot_new = T)
      plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,xlims = xlims2,annotations = annotationtexts,annotation_height = annotationtxtheight/1.75,annotation_cex = 0.5,offsets=offsets,main=title,main_pos=titlepos,dna_seg_labels=dnaseglabels,dna_seg_label_cex=dnaseglabelcex,dna_seg_label_col=dnaseglabelcol,dna_seg_line=dnasegline,min_gap_size=mingapsize,scale=mainscale,dna_seg_scale=dnasegscale,scale_cex=dnasegscalecex)
    }
    dev.off()
    
    #plot legends
    writefilepath=gsubfn('%1|%2', list('%1'=outdir,'%2'=outputname), '%1/%2_alignmentlegend.pdf')
    pdf(writefilepath)
    allcomparisons<-do.call(rbind,comparisons)
    mylegends<-getlegendcols(allcomparisons,legendorientation,poscolvec,negcolvec,colonly=FALSE)
    poslegend<-mylegends[[1]];neglegend<-mylegends[[2]]
    if (legendorientation=='vertical') {
      print(plot_grid(poslegend,neglegend,nrow=2))
    } else {
      print(plot_grid(poslegend,neglegend,ncol=2))
    }
    dev.off()
  }
}

