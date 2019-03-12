################################
# Annotation class and methods
################################
# annotation is a set of text and one or two positions for each
# if one, other must be NA
annotation <- function(x1, x2=NA, text, rot=0, col="black"){
  if (missing(x1) | missing(text)) stop("Args x1 and text must be provided")
  if (!is.numeric(x1)) stop("x1 must be numeric")
  if (!is.na(x2) && !is.numeric(x2)) stop("x2 must be numeric")
  if (!is.character(text)) stop("text must be character")
  as.annotation(data.frame(x1=x1, x2=x2, text=text, stringsAsFactors=FALSE),
                rot=rot, col=col)
}
as.annotation <- function(df, x2=NA, rot=0, col="black"){
  if (is.annotation(df)) return(df)
  if (!all(c("x1", "text") %in% names(df)))
    stop("Data frame should have at least a x1 and text column")
  # attributes x2, col and arg to all rows if not defined
  if (is.null(df$x2)) df$x2 <- x2
  if (is.null(df$color)) df$color <- col
  if (is.null(df$rot)) df$rot <- rot
  class(df) <- c("annotation", "data.frame")
  df
}
is.annotation <- function(annotation){
  inherits(annotation, "annotation")  
}
range.annotation <- function(x, ...){
  annotation <- x
  range(annotation$x1, annotation$x2, na.rm=TRUE)
}
trim.annotation <- function(x, xlim=NULL, ...){
  annotation <- x
  xlim <- as.numeric(xlim)
  if (!is.null(xlim)){
    if (!is.numeric(xlim)) stop("xlim must be numeric")
    if (length(xlim) != 2) stop("xlim must be length 2")
    # to be accepted, x1 > xlim1 and, if x2=NA, xlim1 also < xlim1 or,
    # x2 < xlim2
    annotation <- annotation[annotation$x1 >= xlim[1] &
                               ((is.na(annotation$x2) &
                                   annotation$x1 <= xlim[2]) |
                                  (!is.na(annotation$x2) &
                                     annotation$x2 <= xlim[2])),]
  }
  annotation
}


#####
# create gene grobs
gene_grob <- function(gene, head_len=200, i=0, ...){
  if (!is.dna_seg(gene)) stop("A dna_seg object is required")
  if (nrow(gene) > 1) stop ("gene must be single-row")
  mid <- (gene$start + gene$end)/2
  name <- paste("seg.", i, ".", gene$name, sep="")
  color <- gene$col
  fill <- gene$fill
  if (is.null(fill)) fill <- color
  if (is.null(color)) color <- fill
  # arrows
  if (gene$gene_type == "arrows" || gene$gene_type == "headless_arrows"){
    if (gene$gene_type == "arrows"){
      arrow <- arrow_coord(x1=gene$start, x2=gene$end,
                           y=0.5, strand=gene$strand, head_len=head_len)
    }
    else {
      arrow <- block_coord(gene$start, gene$end, strand=1, y=0.25)
    }
    grob <- polygonGrob(arrow$x, arrow$y, name=name,
                        gp=gpar(col=color, fill=fill, lty=gene$lty,
                                lwd=gene$lwd),
                        default.units="native")
  }
  # blocks
  else if (gene$gene_type == "blocks" || gene$gene_type == "side_blocks"){
    if (gene$gene_type == "side_blocks"){
      block <- block_coord(gene$start, gene$end, strand=gene$strand, y=0.5)
    }
    else {
      block <- block_coord(gene$start, gene$end, strand=2, y=0)
    }
    grob <- polygonGrob(block$x, block$y, name=name,
                        gp=gpar(col=color, fill=fill, lty=gene$lty,
                                lwd=gene$lwd),
                        default.units="native")  
  }
  # lines
  else if (gene$gene_type == "lines" || gene$gene_type == "side_lines"){
    x <- c(gene$start, gene$end)
    if (gene$gene_type == "side_lines") {
      y <- gene$strand/4 + 0.5
    }
    else {
      y <- 0.5
    }
    grob <- segmentsGrob(x0=gene$start, y0=y, x1=gene$end, y1=y, name=name,
                         gp=gpar(col=color, fill=fill, lwd=gene$lwd,
                                 lty=gene$lty),
                         default.units="native")
  }
  # exons
  else if (gene$gene_type == "exons" || gene$gene_type == "side_exons"){
    if (gene$gene_type == "side_exons"){
      block <- exon_coord(gene$start, gene$end, gene$strand)
    }
    else {
      block <- exon_coord(gene$start, gene$end, 0)
    }
    grob <- polygonGrob(block$x, block$y, name=name,
                        gp=gpar(fill=fill, col=color, lty=gene$lty,
                                lwd=gene$lwd),
                        default.units="native")  
  }
  
  # bars
  else if (gene$gene_type == "bars" || gene$gene_type == "side_bars") {
    if (gene$gene_type == "side_bars"){
      y0 <- 0.5; y1 <- 0.5+gene$strand/2
    }
    else {
      y0 <- 0; y1 <- 1
    }
    grob <- segmentsGrob(x0=mid, y0=y0, x1=mid, y1=y1, name=name,
                         gp=gpar(col=color, fill=fill,lwd=gene$lwd,
                                 lty=gene$lty),
                         default.units="native")
  }
  
  # introns
  else if (gene$gene_type == "introns") {
    intron <- list(x0=c(gene$start, mid), x1=c(mid, gene$end),
                   y0=c(0.3, 0.5)*gene$strand + 0.5,
                   y1=c(0.5, 0.3)*gene$strand + 0.5)
    grob <- segmentsGrob(x0=intron$x0, y0=intron$y0,
                         x1=intron$x1, y1=intron$y1, name=name,
                         gp=gpar(lty=gene$lty, lwd=gene$lwd),
                         default.units="native")
  }
  # points
  else if (gene$gene_type == "points" || gene$gene_type == "side_points") {
    if (gene$gene_type == "side_points"){
      y <- 0.5+gene$strand/4
    }
    else {
      y <- 0.5
    }
    ## For pch with filled symbols (21-25), "col" refers to fill
    ## whereas col becomes
    grob <- pointsGrob(x=mid, y=y, name=name,
                       pch=gene$pch, size=unit(gene$cex/2, "char"),
                       gp=gpar(col=color, fill=fill),
                       default.units="native")
  }
  # text
  else if (gene$gene_type == "text" || gene$gene_type == "side_text") {
    if (gene$gene_type == "side_text"){
      just <- c("centre", c("top", "bottom")[gene$strand/2 + 1.5])
    }
    else {
      just <- c("left", "left") #!CHANGED from left left (2)
    }
    grob <- textGrob(label=gene$name, x=mid, y=0.5, name=name,
                     just=just, gp=gpar(col=color, cex=gene$cex),
                     default.units="native")
  }
  else {
    grob <- try(do.call(gene$gene_type, list(gene, ...)), silent=FALSE)
    if (!(is.grob(grob) || all(sapply(grob, is.grob))))
      stop(paste(gene$gene_type, "is an invalid gene_type or",
                 "does not return a grob"))
  }
  grob
}
# create dna_seg grobs
dna_seg_grob <- function(dna_seg, ...){
  if(!is.dna_seg(dna_seg)) stop("A dna_seg object is required")
  grob_list <- gList()
  if (nrow(dna_seg) < 1) return(grob_list)
  for (i in 1:nrow(dna_seg)){
    gene <- as.dna_seg(dna_seg[i,])
    grob_list[[i]] <- gene_grob(gene, ...)
  }
  grob_list
}
# create similarity grobs
similarity_grob <- function(similarity, i){
  if (!is.comparison(similarity)) stop("A comparison object is required")  
  if (nrow(similarity) > 1) stop ("gene must be single-row")
  if (is.null(similarity$col)) similarity$col <- grey(0.5)
  x1 <- c(similarity$start1, similarity$end1)
  x2 <- c(similarity$end2, similarity$start2)
  polygonGrob(x=c(x1, x2), y=c(1, 1, 0, 0),
              name=paste("comp.", i, ".",
                         similarity$start1, "-", similarity$end1, "_",
                         similarity$start2, "-", similarity$end2, sep=""),
              gp=gpar(fill=similarity$col, col=similarity$col, lwd=0.1),
              default.units="native")
}
# create comparisons grobs
comparison_grob <- function(comparison, ...){
  if (!is.comparison(comparison)) stop("A comparison object is required")
  grob_list <- gList()
  if (nrow(comparison) < 1) return(grob_list)
  # go in the reverse order to plot strongest comparison last
  for (i in 1:nrow(comparison)){
    grob_list[[i]] <- similarity_grob(comparison[i,], ...)
  }
  grob_list
}
# create annot grob
label_grob <- function(label, cex=0.8){
  y <- 0
  w <- 0.1
  if (!is.annotation(label)) stop("An annotation object is required")
  if (nrow(label) > 1) stop("A single-line annotation is required")
  grob_list <- gList()
  # range
  if (!is.na(label$x2)){
    bracket_coord <- bracket_coord(label$x1, label$x2, y=y, w=w)
    grob_list[[2]] <- linesGrob(x=bracket_coord$x, y=bracket_coord$y,
                                name=paste("annot", "line",
                                           gsub(" ", "_", label$text), sep="."),
                                default.units="native",
                                gp=gpar(col=label$col))
    x <- mean(c(label$x1, label$x2))
    w <- w*2
  } else {
    x <- label$x1
  }
  if (label$rot == 0){
    just <- c(0.5, 0)  #could try "centre"?
  } else {
    just <- "left"   #c(-0.1, 0.5)
  }
  grob_list[[1]] <- textGrob(label$text, x=x, y=y+w, just=just,
                             name=paste("annot", "label",
                                        gsub(" ", "_", label$text), sep="."),
                             rot=label$rot,
                             default.units="native",
                             gp=gpar(col=label$col, cex=cex))
  grob_list
}
# create annotation grob
annotation_grob <- function(annotation, ...){
  if (!is.annotation(annotation)) stop("An annotation object is required")
  grob_list <- gList()
  if (nrow(annotation) < 1) return(grob_list)
  for (i in 1:nrow(annotation)){
    label <- as.annotation(annotation[i,])
    grob_list[[i]] <- label_grob(label, ...)
  }
  grob_list
}


# create tree grob
dna_seg_label_grob <- function(labels, cex, col){
  n_label <- length(labels)
  y <- seq(1, 0, len=n_label)
  labelGrobs <- gList()
  for (i in 1:n_label) {
    labelGrobs[[i]] <-
      textGrob(x=0, y=y[i], name=paste("label", i, sep="."),
               label=labels[i], just="left", gp=gpar(cex=cex, col=col[i]),
               default.units="native")
  }
  width <- unit(1, "grobwidth", labelGrobs[[which.max(nchar(labels))]])
  labelTree <- gTree(children=labelGrobs,
                     vp=viewport(xscale=c(0, 1), yscale=c(0, 1),
                                 width=width, name="labels"), name="labelsTree")
  list(grob=labelTree, width=width)
}
# create scale grob
scale_grob <- function(max_length){
  rng <- diff(pretty(c(0, max_length), n=8))[1]
  gList(segmentsGrob(x0=max_length, y0=0, x1=max_length-rng, y1=0,
                     name="scale.lines", default.units="native"),
        textGrob(label=human_nt(rng)$text, x=max_length-rng/2, y=0.5,
                 name="scale.text", default.units="native")
  )
}
# create dna_seg scale grob
dna_seg_scale_grob <- function(range, cex=0.6, unit, i, j){
  range <- as.numeric(range)
  if (length(range) != 2 && !is.numeric(range))
    stop("range must be numeric and length 2")
  x0 <- ceiling(range[1]/unit)*unit
  x1 <- floor(range[2]/unit)*unit
  if (x1 < x0) x1 <- x0
  ticks <- seq(x0, x1, by=unit)
  labels <- human_nt(ticks)
  gList(segmentsGrob(x0=ticks, x1=ticks, y0=0, y1=1,
                     gp=gpar(col=grey(0.3)),
                     name=paste("dna_seg_scale", i, j, "lines", sep="."),
                     default.units="native"),
        textGrob(labels$text, x=ticks, y=0.5, hjust=-0.05,
                 gp=gpar(col=grey(0.3), cex=cex),
                 name=paste("dna_seg_scale", i, j, "labels", sep="."),
                 default.units="native")
  )
}
# create gap grob
gap_grob <- function(w, m, i, j){
  segmentsGrob(x0=c(m-w/4, m-w/8),
               x1=c(m+w/8, m+w/4),
               y0=c(0.2,0.2),
               y1=c(0.8,0.8),
               gp=gpar(col=grey(0.3)),
               default.units="native",
               name=paste("gap", i, j, sep="."))
}
## yaxis grob
yaxis_grob <- function(ylim=c(0, 1), cex=0.6, n=3, i){
  at <- pretty(ylim, n=n)
  at <- at[at >= ylim[1] & at <= ylim[2]]
  coords <- yaxis_coords(at, x0=0, x1=0.5*cex)
  gList(segmentsGrob(x0=unit(coords$x0, "lines"), x1=unit(coords$x1, "lines"),
                     y0=unit(coords$y0, "native"), y1=unit(coords$y1, "native"),
                     name=paste("yaxis.segments", i, sep="."),
                     default.units="native"),
        textGrob(at, x=unit(1, "lines"), y=unit(at, "native"),
                 gp=gpar(cex=cex), just=c("left", "left"), #!CHANGED from left centre (3)
                 name=paste("yaxis.labels", i, sep="."),
                 default.units="native")
  )
}


########

# calculate arrow coordinates from gene coordinates
arrow_coord <- function(x1, x2, y=0.5, strand=NULL, width=1, head_len=100){
  # take care of strand, to get x1 as bottom and x2 as tip of arrow
  if (!is.null(strand) && strand == -1){
    x_temp <- x2
    x2 <- x1
    x1 <- x_temp
  }
  w2 <- width/4
  # if the head of the arrow is larger than half of the gene, reduce to half
  if (head_len > abs(x1-x2)/2){
    head_len <- abs(x1-x2)/2
  }
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
# coords for a block
block_coord <- function(start, end, strand, y=0.5){
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, y + strand/2, y + strand/2, y)
  list(x=x, y=y)
}

# exon coord
exon_coord <- function(start, end, strand){
  x <- c(rep(start, 2), rep(end, 2))
  if (strand == 0 ){ y <- c(0.2, 0.8, 0.8, 0.2) }
  if (strand == 1 ){ y <- c(0.5, 0.8, 0.8, 0.5) }
  if (strand == -1 ){ y <- c(0.2, 0.5, 0.5, 0.2) }
  list(x=x, y=y)
}

# coords for a zone annotation
bracket_coord <- function(start, end, y=0, w=0.1){
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, rep(y+w, 2), y)
  list(x=x, y=y)
}
# axis coords
yaxis_coords <- function(at, x0=0, x1=0.5){
  n <- length(at)
  list(x0 = c(rep(x0, n), x1),
       x1 = rep(x1, n+1),
       y0 = c(at, at[1]),
       y1 = c(at, at[n]))
}
# human readable coordinates
human_nt <- function(nt, signif=FALSE){
  tag <- "nt"
  mult <- 1
  med <- median(nt)
  if (med >= 1e9){
    nt <- nt/1e9
    tag <- "Gb"
    mult <- 1e9
  } else if (med >= 1e6){
    nt <- nt/1e6
    tag <- "Mb"
    mult <- 1e6
  } else if (med >= 1e3){
    nt <- nt/1e3
    tag <- "kb"
    mult <- 1e3
  }
  if (signif) nt <- signif(nt, signif)
  list(n=nt, tag=tag, mult=mult, text=paste(nt, tag))
}
# calculate comparison coordinates
calc_comp_coor <- function(gap, xlim, comp, side){
  if (length(gap) != nrow(xlim))
    stop("gap should have the same length as xlim")
  if (side < 1 || side > 2) stop("side should be 1 or 2")
  # x is the moving cursor
  x <- 0
  old_start <- if (side==1) comp$start1 else comp$start2 
  old_end <- if (side==1) comp$end1 else comp$end2
  start <- old_start
  end <- old_end
  for (i in 1:nrow(xlim)){
    # increment by the gap length
    x <- x + gap[i]
    # select comps
    idx <- old_start >= xlim$x0[i] & old_end <= xlim$x1[i]
    # re-number by substracting the xlim and adding x
    if (xlim$strand[i] == 1){
      start[idx] <- old_start[idx] - xlim$x0[i] + x
      end[idx] <- old_end[idx] - xlim$x0[i] + x
    } else {
      start[idx] <- xlim$x1[i] - old_start[idx] + x
      end[idx] <- xlim$x1[i] - old_end[idx] + x
    }
    # increment x by the length of the segment
    x <- x + xlim$length[i]
  }
  # reattribute start and stop
  if (side==1) comp$start1 <- start else comp$start2 <- start
  if (side==1) comp$end1 <- end else comp$end2 <- end
  # return the modified comp
  comp
}
middle <- function(dna_seg){
  if (!is.dna_seg(dna_seg)) stop("argument should be a dna_seg object")
  apply(dna_seg[,c("start", "end")], 1, mean)
}
## Emulate artemis colors
## 0  white          (RGB values: 255 255 255)
## 1  dark grey      (RGB values: 100 100 100)
## 2  red            (RGB values: 255   0   0)
## 3  green          (RGB values:   0 255   0)
## 4  blue           (RGB values:   0   0 255)
## 5  cyan           (RGB values:   0 255 255)
## 6  magenta        (RGB values: 255   0 255)
## 7  yellow         (RGB values: 255 255   0)
## 8  pale green     (RGB values: 152 251 152)
## 9  light sky blue (RGB values: 135 206 250)
## 10 orange         (RGB values: 255 165   0)
## 11 brown          (RGB values: 200 150 100)
## 12 pale pink      (RGB values: 255 200 200)
## 13 light grey     (RGB values: 170 170 170)
## 14 black          (RGB values:   0   0   0)
## 15 mid red:       (RGB values: 255  63  63)
## 16 light red      (RGB values: 255 127 127)
## 17 pink           (RGB values: 255 191 191)
artemisColors <- function(){
  names <- c("white", "dark grey", "red", "green",
             "blue", "cyan", "magenta", "yellow", "pale green",
             "light sky blue", "orange", "brown", "pale pink",
             "light grey", "black", "mid red", "light red", "pink")
  numbers <- 0:(length(names)-1)
  r <- c(255, 100, 255, 0, 0, 0, 255, 255, 152, 135, 255, 200, 255,
         170, 0, 255, 255, 255)
  g <- c(255, 100, 0, 255, 0, 255, 0, 255, 251, 206, 165, 150, 200,
         170, 0, 63, 127, 191)
  b <- c(255, 100, 0, 0, 255, 255, 255, 0, 152, 250, 0, 100, 200,
         170, 0, 63, 127, 191)
  colors <- rgb(r, g, b, maxColorValue=255)
  data.frame(n=numbers, names=names, colors=colors, r=r, g=g, b=b,
             stringsAsFactors=FALSE)
}


# Minimize comparison sizes
################################################################################
# from a set of comparisons and lengths, determine the best possible
# arrangement of maps
minimize_comps <- function(comparisons, xlims, lengths, prel_offsets,
                           fixed_gap_length=FALSE){
  # function to minimize. Calculates the mean of the absolute differences
  # between starts and ends for direct hits only
  mean_w_fixed_gaps <- function(first_offset, fixed_offsets, ...){
    mean_w_gaps(c(first_offset, fixed_offsets), ...)
  }
  mean_w_gaps <- function(offsets, offsets_ref, xlim, xlim_ref, comp, side_ref){
    if (nrow(comp) == 0) return(0)
    side_test <- if (side_ref == 2) 1 else 2
    # recalc for ref side
    comp <- calc_comp_coor(gap=offsets_ref, xlim=xlim_ref, comp=comp,
                           side=side_ref)
    # test
    comp <- calc_comp_coor(gap=offsets, xlim=xlim, comp=comp, side=side_test)
    direction <- sign(comp$end1-comp$start1)*sign(comp$end2-comp$start2)
    dists <- c(abs(comp$start1-comp$start2),
               abs(comp$end1-comp$end2))#[direction > 0]
    lengths <- abs(comp$end1-comp$start1) + abs(comp$end2-comp$start2)
    lengths[lengths==0] <- 1
    weighted.mean(dists, c(lengths, lengths))
  }
  
  n_org <- length(lengths)
  offsets <- prel_offsets
  if (length(comparisons) < 1) return(offsets)
  idx_ref <- which.max(lengths)
  max_len <- max(lengths)
  
  # go up from ref
  if (idx_ref > 1){
    # comp i is between org i and i+1
    for (i in (idx_ref-1):1){
      # optimise
      # if fixed_gap_length, optimse only on the first gap
      if (fixed_gap_length && length(offsets[[i]]) > 1){
        fixed_offsets <- offsets[[i]][2:length(offsets[[i]])]
        opt <- optim(par=offsets[[i]][1], fn=mean_w_fixed_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i]][1],
                     fixed_offsets=fixed_offsets, offsets_ref=offsets[[i+1]],
                     xlim=xlims[[i]], xlim_ref=xlims[[i+1]],
                     comp=comparisons[[i]], side_ref=2)
        offsets[[i]][1] <- opt$par
      }
      # else optimise on all offsets
      else {
        opt <- optim(par=offsets[[i]], fn=mean_w_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i]],
                     offsets_ref=offsets[[i+1]],
                     xlim=xlims[[i]], xlim_ref=xlims[[i+1]],
                     comp=comparisons[[i]], side_ref=2)
        offsets[[i]] <- opt$par
      }
    }
  }
  # go down
  if (idx_ref < n_org){
    for (i in idx_ref:(n_org-1)){
      # optimise
      # if fixed_gap_length, optimse only on the first gap
      if (fixed_gap_length && length(offsets[[i]]) > 1){
        fixed_offsets <- offsets[[i]][2:length(offsets[[i]])]
        opt <- optim(par=offsets[[i+1]][1], fn=mean_w_fixed_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i+1]][1],
                     fixed_offsets=fixed_offsets, offsets_ref=offsets[[i]],
                     xlim=xlims[[i+1]], xlim_ref=xlims[[i]],
                     comp=comparisons[[i]], side_ref=1)
        offsets[[i+1]][1] <- opt$par
      }
      # else optimise on all offsets
      else {
        opt <- optim(par=offsets[[i+1]], fn=mean_w_gaps,
                     method="L-BFGS-B",
                     lower=offsets[[i+1]],
                     offsets_ref=offsets[[i]],
                     xlim=xlims[[i+1]], xlim_ref=xlims[[i]],
                     comp=comparisons[[i]], side_ref=1)
        offsets[[i+1]] <- opt$par
      }
    }
  }
  offsets
}

################################
# dna_seg class and methods
################################
dna_seg <- function(x, ...){
  # support for list ?
  if (is.data.frame(x)){
    return(as.dna_seg(x, ...))
  } else if (is.list(x)) {
    dna_seg(as.data.frame(x, stringsAsFactors=FALSE))
  } else {
    stop(paste("Cannot coerce class", class(x), "to dna_seg"))
  }
}
# convert to dna_seg format. 
as.dna_seg <- function(df, col="blue", fill="blue", lty=1, lwd=1, pch=8,
                       cex=1, gene_type="arrows"){
  # check for class dna_seg, list, df
  if (is.dna_seg(dna_seg)) return(df)
  #if (is.list(df) && !is.data.frame(df)) df <- dna_seg(df)
  if (is.data.frame(df)) {
    # check that it has rows
    if (nrow(df) < 1){
      stop("Number of rows is 0, check data input")
    }
    # check for columns
    names <- c("name", "start", "end", "strand")
    if (!identical(names(df)[1:4], names))
      stop("Col names should start with name, start, end, strand")
    if (!all(sapply(df, function(x) all(!is.null(x)))))
      stop("NULL values not allowed in data.frame")
    if (!is.numeric(df$start) | !is.numeric(df$end)){
      stop("Start and end must be numeric")
    }
    if (is.factor(df$name)) df$name <- as.character(df$name)
    if (is.factor(df$strand)) df$strand <- as.character(df$strand)
    if (is.factor(df$col)) df$col <- as.character(df$col)
    if (is.factor(df$fill)) df$fill <- as.character(df$fill)
    if (is.factor(df$gene_type)) df$gene_type <- as.character(df$gene_type)
    # care for strand
    if (is.character(df$strand)) {
      df$strand[df$strand=="+"] <- 1
      df$strand[df$strand=="-"] <- -1
      df$strand <- as.numeric(df$strand)
    }
    if (!all(df$strand %in% c(-1, 1)))
      stop("Strand vector must be composed of 1, -1, - and +, uniquely")
    # col
    if (is.null(df$col)) df$col <- col
    if (is.null(df$fill)) df$fill <- fill
    # lwd & lty
    if (is.null(df$lty)) df$lty <- lty
    if (is.null(df$lwd)) df$lwd <- lwd
    if (is.null(df$pch)) df$pch <- pch
    if (is.null(df$cex)) df$cex <- cex
    # gene_type: not given in gene
    if (is.null(df$gene_type)){
      if (gene_type == "auto") gene_type <- auto_gene_type(nrow(df))
      df$gene_type <- gene_type
    }
    if (is.null(df$gene_type)) df$gene_type <- gene_type
    # check for correct argument types
    if (!is.character(df$name)) stop("Non-character name")
    if (!(is.numeric(df$start) && is.numeric(df$end)))
      stop("Non-numeric start or end")
    if (!all(is.numeric(c(df$lwd, df$lty, df$pch, df$cex))))
      stop("lwd, lty, pch and cex must be numeric")
    if (!is.character(df$gene_type))
      stop(paste("gene_type must be a character vector, made of:",
                 paste(gene_types(), collapse=", "), "or a function name"))
  }
  else {
    stop("Unable to handle this format")
  }
  class(df) <- c("dna_seg", "data.frame")
  return(df)
}
is.dna_seg <- function(dna_seg){
  inherits(dna_seg, "dna_seg")
}
# calculates dna_seg range
range.dna_seg <- function(x, ...){
  dna_seg <- x
  range(dna_seg$start, dna_seg$end, na.rm=FALSE)
}
#range.data.frame <- range.dna_seg
# trim dna_seg given x limits
trim.dna_seg <- function(x, xlim=NULL, ...){
  dna_seg <- x
  if (!is.null(xlim)){
    if (!is.numeric(xlim)) stop("xlim must be numeric")
    if (length(xlim) != 2) stop("xlim must be length 2")
    dna_seg <- dna_seg[dna_seg$start >= xlim[1] & dna_seg$end <= xlim[2],]
  }
  dna_seg
}
# reverse dna_seg
reverse.dna_seg <- function(x, ...){
  dna_seg <- x
  start <- -dna_seg$end
  dna_seg$end <- -dna_seg$start
  dna_seg$start <- start
  dna_seg$strand <- -dna_seg$strand
  dna_seg
}
# merges two dna_segs. So far returns the minimal common set of columns, maybe
# changed in future
c.dna_seg <- function(...){
  # a little helper function to grab colnames
  #fold <- function(f, x, L) (for(e in L) x <- f(x, e))
  fold <- function(x, fun) {
    if (length(x) == 1) return(fun(x))
    accumulator <- fun(x[[1]], x[[2]])
    if (length(x) == 2) return(accumulator)
    for(i in 3:length(x)) {
      accumulator <- fun(accumulator, x[[i]])
    }
    accumulator
  }
  # parse args
  x <- list(...)
  n <- length(x)
  if (!all(sapply(x, is.dna_seg))) 
    stop("All elements must be of class dna_seg")
  if (n == 1) return(x)
  cols <- lapply(x, names)
  com_cols <- fold(cols, intersect)
  dna_seg <- x[[1]][com_cols]
  for (i in 2:n){
    dna_seg <- rbind(dna_seg, x[[i]][com_cols])
  }
  dna_seg
}


################################
# Comparison class and methods
################################
# A comparison is mostly a bunch of similarities
comparison <- function(x){
  # if data.frame
  if (is.data.frame(x)){
    return(as.comparison(x))
  } else if (is.list(x)){
    return(as.comparison(as.data.frame(x, stringsAsFactors=FALSE)))
  } else {
    stop(paste("Cannot coerce class", class(x), "to comparison"))
  }
}
as.comparison <- function(df){
  # check for class comparison, list, df
  if (is.comparison(comparison)) return(df)
  if (is.list(df) && !is.data.frame(df)) df <- comparison(df)
  if (is.data.frame(df)) {
    # check for columns
    names <- c("start1", "end1", "start2", "end2")
    if (!identical(names(df)[1:4], names))
      stop("Col names should start with start1, end1, start2 and end2")
    if (!all(sapply(df, function(x) all(!is.null(x)))))
      stop("NULL values not allowed in data.frame")
    if (!is.numeric(df$start1) | !is.numeric(df$end1) |
        !is.numeric(df$start2) | !is.numeric(df$end2)){
      stop("Starts and ends must be numeric")
    }
    # add direction column
    df$direction <- ifelse(sign(df$start1-df$end1)
                           * sign(df$start2-df$end2) > 0, 1, -1)
    # check color (character and default)
    if (is.factor(df$col)) df$col <- as.character(df$col)
  }
  else {
    stop("Unable to handle this format")
  }
  class(df) <- c("comparison", "data.frame")
  df
}
is.comparison <- function(comparison){
  inherits(comparison, "comparison")
}
# gives the range of a comparison
range.comparison <- function(x, overall=TRUE, ...){
  comparison <- x
  if (overall){
    range <- range(comparison$start1, comparison$end1,
                   comparison$start2, comparison$end2, na.rm=FALSE)
  } else {
    xlim1 <- range(comparison$start1, comparison$end1, na.rm=FALSE)
    xlim2 <- range(comparison$start2, comparison$end2, na.rm=FALSE)
    range <- data.frame(xlim1=xlim1, xlim2=xlim2)
  }
  range
}
# trim comparison given x limits
trim.comparison <- function(x, xlim1=c(-Inf, Inf), xlim2=c(-Inf, Inf), ...){
  comparison <- x
  if (!is.null(xlim1) && !is.null(xlim2)){
    if (!is.numeric(xlim1) || !is.numeric(xlim2)) stop("xlims must be numeric")
    if (length(xlim1) != 2 || length(xlim2) != 2) stop("xlims must be length 2")
    ## testing to include overlapping comps
    ## direction 1
    comparison$start1[comparison$start1 < xlim1[1] &
                        comparison$end1 > xlim1[1]] <- xlim1[1]
    comparison$end1[comparison$start1 < xlim1[2] &
                      comparison$end1 > xlim1[2]] <- xlim1[2]
    comparison$start2[comparison$start2 < xlim2[1] &
                        comparison$end2 > xlim2[1]] <- xlim2[1]
    comparison$end2[comparison$start2 < xlim2[2] &
                      comparison$end2 > xlim2[2]] <- xlim2[2]
    ## direction -1
    comparison$start1[comparison$start1 > xlim1[2] &
                        comparison$end1 < xlim1[2]] <- xlim1[2]
    comparison$end1[comparison$start1 > xlim1[1] &
                      comparison$end1 < xlim1[1]] <- xlim1[1]
    comparison$start2[comparison$start2 > xlim2[2] &
                        comparison$end2 < xlim2[2]] <- xlim2[2]
    comparison$end2[comparison$start2 > xlim2[1] &
                      comparison$end2 < xlim2[1]] <- xlim2[1]
    comparison <-
      comparison[comparison$start1 >= xlim1[1] & comparison$end1 <= xlim1[2] &
                   comparison$start2 >= xlim2[1] & comparison$end2 <= xlim2[2],]
  }
  comparison
}
# reverses a comparison. side <1 for no sides, 1 for first,
# 2 for second, >2 for bth
reverse.comparison <- function(x, side=0, ...){
  comparison <- x
  if (side > 1){
    comparison$start2 <- -comparison$start2
    comparison$end2 <- -comparison$end2
  }
  if (side == 1 || side >2){
    comparison$start1 <- -comparison$start1
    comparison$end1 <- -comparison$end1    
  }
  # if only one side changed, change direction
  if (side == 1 || side == 2){
    comparison$direction <- -comparison$direction
  }
  comparison
}



#############



plot_gene_map2<-function(dna_segs, comparisons = NULL, tree = NULL, tree_width = NULL, 
    tree_branch_labels_cex = NULL, tree_scale = FALSE, legend = NULL, 
    annotations = NULL, annotation_height = 1, annotation_cex = 0.8, 
    seg_plots = NULL, seg_plot_height = 3, seg_plot_height_unit = "lines", 
    seg_plot_yaxis = 3, seg_plot_yaxis_cex = scale_cex, xlims = NULL, 
    offsets = NULL, minimum_gap_size = 0.05, fixed_gap_length = FALSE, 
    limit_to_longest_dna_seg = TRUE, main = NULL, main_pos = "centre", 
    dna_seg_labels = NULL, dna_seg_label_cex = 1, dna_seg_label_col = "black", 
    gene_type = NULL, arrow_head_len = 200, dna_seg_line = TRUE, 
    scale = TRUE, dna_seg_scale = !scale, n_scale_ticks = 7, 
    scale_cex = 0.6, global_color_scheme = c("auto", "auto", 
        "blue_red", 0.5), override_color_schemes = FALSE, plot_new = TRUE, 
    debug = 0, ...) 
{
    if (missing(dna_segs)) 
        stop("Argument dna_segs must be provided")
    if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) 
        stop("Argument dna_segs must be a list of dna_seg objects")
    n_dna_segs <- length(dna_segs)
    n_rows <- 3 * n_dna_segs - 1
    n_comparisons <- length(comparisons)
    if (n_comparisons > 1) {
        if (!is.list(comparisons) || !all(sapply(comparisons, 
            is.comparison))) 
            stop("Argument comparisons must be a list of comparison objects")
    }
    if (n_comparisons > 0 && !(n_dna_segs - n_comparisons == 
        1)) 
        stop("Number of comparisons not correct")
    if (is.null(dna_seg_labels) && !is.null(dna_segs)) {
        dna_seg_labels <- names(dna_segs)
    }
    if (!is.null(dna_seg_labels) && !(length(dna_seg_labels) == 
        n_dna_segs)) 
        stop("Argument dna_seg_labels doesn't have the same length as dna_segs")
    if (length(dna_seg_label_col) == 1) {
        dna_seg_label_col <- rep(dna_seg_label_col, n_dna_segs)
    }
    else if (!length(dna_seg_label_col) == n_dna_segs) {
        stop("Length of argument dna_seg_label_col must be 1 or as dna_segs")
    }
    if (!is.null(tree)) {
        if (!inherits(tree, "phylog")) 
            stop("Argument tree should be of class phylog (ade4)")
        if (is.null(dna_seg_labels)) 
            stop("If tree is given, label names should be provided via named list dna_segs or dna_seg_labels")
        if (length(tree$leaves) != n_dna_segs) 
            stop("Number of leaves in the tree not equal to number of dna segs")
        if (!all(dna_seg_labels %in% names(tree$leaves))) 
            stop("Tree leaves not corresponding to dna_seg labels")
        if (is.null(tree_branch_labels_cex)) {
            if (length(grep("^[^I]", names(tree$nodes)))) {
                tree_branch_labels_cex <- 0.8
            }
            else {
                tree_branch_labels_cex <- 0
            }
        }
    }
    if (!is.null(seg_plots)) {
        if (is.seg_plot(seg_plots)) {
            s_plot <- seg_plots
            seg_plots <- c(list(s_plot), rep(list(NULL), n_dna_segs - 
                1))
        }
        else if (length(seg_plots) == n_dna_segs) {
            if (!all(sapply(seg_plots, function(x) is.seg_plot(x) || 
                is.null(x)))) 
                stop("All elements of seg_plots should be NULL or seg_plots objects")
        }
        else stop("seg_plots must be of same length as dna_segs")
        seg_plot_h <- ifelse(sapply(seg_plots, is.null), 0, seg_plot_height)
    }
    else {
        seg_plot_h <- rep(0, n_dna_segs)
    }
    if (is.null(seg_plot_yaxis) || !is.numeric(seg_plot_yaxis) || 
        is.null(seg_plots)) 
        seg_plot_yaxis <- 0
    if (!is.null(annotations)) {
        if (is.annotation(annotations)) {
            annot <- annotations
            annotations <- c(list(annot), rep(list(NULL), n_dna_segs - 
                1))
        }
        else if (length(annotations) == n_dna_segs) {
            if (!all(sapply(annotations, function(x) is.annotation(x) || 
                is.null(x)))) 
                stop("All elements of annotations should be NULL or annotation objects")
        }
        else stop("annotation must be of same length as dna_segs")
        annot_h <- ifelse(sapply(annotations, is.null), 0, annotation_height)
    }
    else {
        annot_h <- rep(0, n_dna_segs)
    }
    if (!is.null(xlims)) {
        if (!is.list(xlims)) 
            stop("xlims must be a list")
        if (length(xlims) != n_dna_segs) 
            stop("xlims must be of the same length as dna_segs")
        if (!all(sapply(xlims, function(x) (length(x)%%2) == 
            0))) 
            stop("All elements of xlims should have an even number of elements")
        for (i in 1:length(xlims)) {
            xlim <- xlims[[i]]
            rng <- range(dna_segs[[i]])
            min <- rng[1] - 0.02 * diff(rng)
            max <- rng[2] + 0.02 * diff(rng)
            if (is.null(xlim)) 
                xlim <- c(min, max)
            xlim[xlim == -Inf] <- min
            xlim[xlim == Inf] <- max
            if (!is.numeric(xlim)) 
                stop("All elements of xlims should be numeric")
            xlim <- data.frame(matrix(xlim, ncol = 2, byrow = TRUE))
            names(xlim) <- c("x0", "x1")
            xlim$strand <- ifelse(xlim$x0 < xlim$x1, 1, -1)
            for (j in 1:nrow(xlim)) xlim[j, 1:2] <- sort(xlim[j, 
                1:2])
            xlims[[i]] <- xlim
        }
    }
    else {
        xlims <- list()
        for (i in 1:n_dna_segs) {
            rng <- range(dna_segs[[i]])
            xlims[[i]] <- data.frame(x0 = rng[1] - 0.02 * diff(rng), 
                x1 = rng[2] + 0.02 * diff(rng), strand = 1)
        }
    }
    if (!is.null(offsets) && length(offsets) != n_dna_segs) {
        stop("Length of offsets not equal to number of dna_segs")
    }
    if (main_pos == "centre") {
        main_x <- 0.5
        main_just <- "centre"
    }
    else if (main_pos == "left") {
        main_x <- 0
        main_just <- "left"
    }
    else if (main_pos == "right") {
        main_x <- 1
        main_just <- "right"
    }
    else {
        stop("main_pos should be one of centre, left, right")
    }
    if (is.logical(dna_seg_line)) {
        dna_seg_line <- as.character(dna_seg_line)
        dna_seg_line[dna_seg_line == "TRUE"] <- "black"
    }
    if (!is.character(dna_seg_line)) 
        stop("dna_seg_line should be eiher a logical or character giving color")
    if (length(dna_seg_line) == 1) {
        dna_seg_line <- rep(dna_seg_line, n_dna_segs)
    }
    else if (length(dna_seg_line) != n_dna_segs) {
        stop("dna_seg_line should be of length 1 or same length as dna_segs")
    }
    if (!is.null(gene_type) && !(gene_type %in% gene_types())) 
        stop(paste("gene_type muste be one of:", paste(gene_types(), 
            collapse = ", ")))
    if (is.logical(dna_seg_scale)) {
        if (length(dna_seg_scale) == 1) {
            dna_seg_scale <- rep(dna_seg_scale, n_dna_segs)
        }
        else if (length(dna_seg_scale) != n_dna_segs) {
            stop("dna_seg_scale must be the same length dna_segs")
        }
    }
    else {
        stop("dna_seg_scale must be logical")
    }
    if (length(global_color_scheme) != 4) 
        stop("global_color_scheme should be length 4")
    glob_col_sch_2_vals <- c("increasing", "decreasing", "auto")
    if (length(grep(global_color_scheme[2], glob_col_sch_2_vals)) != 
        1) {
        stop(paste("Second argument to global_color_scheme should be one of", 
            paste(glob_col_sch_2_vals, collapse = ", ")))
    }
    else {
        global_color_scheme[2] <- grep(global_color_scheme[2], 
            glob_col_sch_2_vals, value = TRUE)
    }
    h <- rep(1, n_rows)
    h[seq(2, n_rows, by = 3)] <- 1 + scale_cex * dna_seg_scale + 
        annot_h
    h[seq(1, n_rows, by = 3)] <- seg_plot_h
    dna_seg_heights <- unit(h, c(rep(c(seg_plot_height_unit, 
        "lines", "null"), n_dna_segs), "lines"))
    if (!is.null(gene_type)) {
        if (gene_type == "auto") {
            n_genes <- sapply(dna_segs, nrow)
            gene_type <- auto_gene_type(n_genes)
        }
        for (i in 1:n_dna_segs) {
            dna_segs[[i]]$gene_type <- gene_type
        }
    }
    if ((!any("col" %in% unlist(lapply(comparisons, names))) || 
        override_color_schemes) && !is.null(comparisons)) {
        num_cols <- lapply(comparisons, function(x) names(x)[sapply(x, 
            is.numeric)])
        shared_num_cols <- names(which(table(unlist(num_cols)) == 
            length(num_cols)))
        shared_num_cols <- shared_num_cols[!shared_num_cols %in% 
            c("start1", "start2", "end1", "end2")]
        if (global_color_scheme[1] == "auto") {
            names_comp_1 <- names(comparisons[[1]])
            global_color_scheme[1] <- if ("per_id" %in% shared_num_cols) {
                "per_id"
            }
            else if ("e_value" %in% shared_num_cols) {
                "e_value"
            }
            else {
                names_comp_1[names_comp_1 %in% shared_num_cols][1]
            }
        }
        else if (!global_color_scheme[1] %in% shared_num_cols) {
            stop("One or all columns don't have the indicated column for global color scheme")
        }
        if (global_color_scheme[2] == "auto") {
            global_color_scheme[2] <- if (global_color_scheme[1] %in% 
                c("mism", "gaps", "e_value")) 
                TRUE
            else FALSE
        }
        else if (global_color_scheme[2] == "decreasing") {
            global_color_scheme[2] <- TRUE
        }
        else if (global_color_scheme[2] == "increasing") {
            global_color_scheme[2] <- FALSE
        }
        else {
            stop("Invalid value for global_color_scheme[2]")
        }
        range_col_from <- c(Inf, -Inf)
        for (i in 1:n_comparisons) {
            if (nrow(comparisons[[i]]) > 0) {
                range_col_from[1] <- min(c(range_col_from[1], 
                  comparisons[[i]][[global_color_scheme[1]]]))
                range_col_from[2] <- max(c(range_col_from[2], 
                  comparisons[[i]][[global_color_scheme[1]]]))
            }
        }
        for (i in 1:n_comparisons) {
            comparisons[[i]]$col <- apply_color_scheme(x = comparisons[[i]][[global_color_scheme[1]]], 
                direction = comparisons[[i]]$direction, color_scheme = global_color_scheme[3], 
                decreasing = global_color_scheme[2], rng = range_col_from, 
                transparency = as.numeric(global_color_scheme[4]))
        }
    }
    for (i in 1:n_dna_segs) {
        xlims[[i]]$length <- xlims[[i]]$x1 - xlims[[i]]$x0
    }
    def_gap_length <- max(sapply(xlims, function(x) sum(x$length))) * 
        minimum_gap_size
    unpadded_lengths <- sapply(xlims, function(x) sum(x$length) + 
        (nrow(x) - 1) * def_gap_length)
    longest_seg <- which.max(unpadded_lengths)
    max_length <- unpadded_lengths[longest_seg]
    scale_unit <- diff(pretty(c(0, max_length), n = n_scale_ticks + 
        2)[1:2])
    seg_subplots <- list()
    dna_subsegs <- list()
    for (i in 1:n_dna_segs) {
        n_subsegs <- nrow(xlims[[i]])
        dna_subsegs[[i]] <- list()
        seg_subplots[[i]] <- list()
        for (j in 1:n_subsegs) {
            dna_subsegs[[i]][[j]] <- trim.dna_seg(dna_segs[[i]], 
                c(xlims[[i]]$x0[j], xlims[[i]]$x1[j]))
            if (seg_plot_h[[i]] > 0) {
                seg_subplots[[i]][[j]] <- trim.seg_plot(seg_plots[[i]], 
                  c(xlims[[i]]$x0[j], xlims[[i]]$x1[j]))
            }
        }
    }
    if (n_comparisons > 0) {
        for (i in 1:n_comparisons) {
            comp1 <- comparisons[[i]][0, ]
            for (j in 1:nrow(xlims[[i]])) {
                for (k in 1:nrow(xlims[[i + 1]])) {
                  comp1 <- rbind(comp1, trim.comparison(comparisons[[i]], 
                    xlim1 = as.numeric(xlims[[i]][j, c("x0", 
                      "x1")]), xlim2 = as.numeric(xlims[[i + 
                      1]][k, c("x0", "x1")])))
                }
            }
            comparisons[[i]] <- comp1
        }
    }
    if (is.null(offsets)) {
        prel_offsets <- lapply(xlims, function(x) c(0, rep(def_gap_length, 
            nrow(x) - 1)))
        offsets <- minimize_comps(comparisons, xlims, unpadded_lengths, 
            prel_offsets, fixed_gap_length)
    }
    else {
        offsets <- as.list(offsets)
        for (i in 1:n_dna_segs) {
            if (length(offsets[[i]]) == 1) {
                offsets[[i]] <- c(offsets[[i]], rep(def_gap_length, 
                  nrow(xlims[[i]]) - 1))
            }
            else if (length(offsets[[i]]) != nrow(xlims[[i]])) {
                stop("The length of each element of offsets should be either one or equal to the number of subsegments in the corresponding segment.")
            }
        }
    }
    if (limit_to_longest_dna_seg) {
        for (i in 1:n_dna_segs) {
            tot_length <- sum(c(xlims[[i]]$length, offsets[[i]]))
            if (tot_length > max_length) {
                excess <- tot_length - max_length
                for (j in length(offsets[[i]]):1) {
                  if ((offsets[[i]][j] - excess) < def_gap_length) {
                    excess <- excess - offsets[[i]][j] + def_gap_length
                    offsets[[i]][j] <- def_gap_length
                  }
                  else {
                    offsets[[i]][j] <- offsets[[i]][j] - excess
                    break
                  }
                }
            }
        }
    }
    else {
        padded_lengths <- sapply(xlims, function(x) sum(x$length)) + 
            sapply(offsets, sum)
        max_length <- max(padded_lengths)
    }
    if (n_comparisons > 0) {
        for (i in 1:n_comparisons) {
            comparisons[[i]] <- calc_comp_coor(offsets[[i]], 
                xlims[[i]], comparisons[[i]], side = 1)
            comparisons[[i]] <- calc_comp_coor(offsets[[i + 1]], 
                xlims[[i + 1]], comparisons[[i]], side = 2)
        }
    }
    padded_lengths <- sapply(xlims, function(x) sum(x$length)) + 
        sapply(offsets, sum)
    max_length <- max(padded_lengths)
    longest_segment <- which.max(padded_lengths)
    dna_seg_grobs <- list()
    dna_seg_scale_grobs <- list()
    for (i in 1:n_dna_segs) {
        dna_seg_grobs[[i]] <- list()
        dna_seg_scale_grobs[[i]] <- list()
        for (j in 1:length(dna_subsegs[[i]])) {
            if (debug > 0 && debug < nrow(dna_subsegs[[i]][[j]])) 
                dna_subsegs[[i]][[j]] <- dna_subsegs[[i]][[j]][1:debug, 
                  ]
            dna_seg_grobs[[i]][[j]] <- dna_seg_grob(dna_subsegs[[i]][[j]], 
                arrow_head_len, i, ...)
            dna_seg_scale_grobs[[i]][[j]] <- if (dna_seg_scale[[i]]) 
                dna_seg_scale_grob(range = xlims[[i]][j, c("x0", 
                  "x1")], cex = scale_cex, unit = scale_unit, 
                  i = i, j = j)
            else NULL
        }
    }
    seg_plot_grobs <- list()
    seg_plot_ylims <- list()
    for (i in 1:n_dna_segs) {
        seg_plot_grobs[[i]] <- list()
        xl_sg <- c(Inf, -Inf)
        for (j in 1:length(dna_subsegs[[i]])) {
            if (length(seg_plots[[i]]) > 0) {
                grb <- do.call(seg_subplots[[i]][[j]]$func, seg_subplots[[i]][[j]]$args)
                rng <- nice_ylim.seg_plot(seg_subplots[[i]][[j]])
                xl_sg[1] <- min(xl_sg[1], rng[1])
                xl_sg[2] <- max(xl_sg[2], rng[2])
                seg_plot_grobs[[i]][[j]] <- grb
            }
            else {
                seg_plot_grobs[[i]] <- NULL
            }
        }
        seg_plot_ylims[[i]] <- if (is.null(seg_plots[[i]]$ylim)) 
            xl_sg
        else seg_plots[[i]]$ylim
    }
    seg_plot_yaxis_grobs <- list()
    for (i in 1:n_dna_segs) {
        seg_plot_yaxis_grobs[[i]] <- if (length(seg_plots[[i]]) > 
            0 && seg_plot_yaxis > 0) 
            yaxis_grob(seg_plot_ylims[[i]], cex = seg_plot_yaxis_cex, 
                n = seg_plot_yaxis, i)
        else NULL
    }
    comparison_grobs <- list()
    if (n_comparisons > 0) {
        for (i in 1:n_comparisons) {
            if (debug > 0 && debug < nrow(comparisons[[i]])) 
                comparisons[[i]] <- comparisons[[i]][1:debug, 
                  ]
            comparison_grobs[[i]] <- comparison_grob(comparisons[[i]], 
                i)
        }
    }
    if (!is.null(annotations)) {
        annotation_grobs <- list()
        for (i in 1:n_dna_segs) {
            if (is.null(annotations[[i]])) {
                annotation_grobs[[i]] <- list(NULL)
            }
            else {
                annotation_grobs[[i]] <- list()
                for (j in 1:length(dna_subsegs[[i]])) {
                  annot <- trim.annotation(annotations[[i]], 
                    xlims[[i]][j, c("x0", "x1")])
                  annotation_grobs[[i]][[j]] <- annotation_grob(annot, 
                    annotation_cex)
                }
            }
        }
    }
    if (scale) {
        scale_grob <- scale_grob(max_length)
        scale_h <- 1
    }
    else {
        scale_h <- 0
    }
    if (!is.null(main)) {
        main_grob <- textGrob(x = main_x, label = main, gp = gpar(cex = 1.2), 
            just = main_just)
        main_h <- 1.8
    }
    else {
        main_h <- 0
    }
    if (!is.null(tree)) {
        y <- permute_tree(tree, dna_seg_labels)
        tree_grob <- phylog_grob(tree, 1 - ((y - 1)/(n_dna_segs - 
            1)), clabel.leaves = dna_seg_label_cex, clabel.nodes = tree_branch_labels_cex, 
            tree.scale = tree_scale)
        tree_w <- unit(0.1, "npc") + tree_grob$width
    }
    else if (!is.null(dna_seg_labels)) {
        tree_grob <- dna_seg_label_grob(dna_seg_labels, cex = dna_seg_label_cex, 
            col = dna_seg_label_col)
        tree_w <- tree_grob$width
    }
    else {
        tree_grob <- NULL
        tree_w <- unit(0, "npc")
    }
    if (!is.null(tree_width)) 
        tree_w <- unit(tree_width, "inches")
    if (tree_scale) 
        scale_h <- 1
    if (plot_new) 
        grid.newpage()
    pushViewport(viewport(width = unit(1, "npc") - unit(1, "lines"), 
        height = unit(1, "npc") - unit(1, "lines"), name = "oma"), 
        viewport(layout = grid.layout(2, 1, heights = unit(c(main_h, 
            1), c("lines", "null"))), name = "oma_layout"))
    if (!is.null(main)) {
        pushViewport(viewport(layout.pos.row = 1, name = "main"))
        grid.draw(main_grob)
        upViewport()
    }
    seg_plot_yaxis_w <- if (seg_plot_yaxis > 0 && !is.null(seg_plots[[longest_segment]])) 
        unit(3, "grobwidth", data = seg_plot_yaxis_grobs[[longest_segment]])
    else unit(0, "null")
    pushViewport(viewport(layout.pos.row = 2, layout = grid.layout(2, 
        3, heights = unit(c(1, scale_h), c("null", "lines")), 
        widths = unit.c(tree_w, unit(1, "null"), seg_plot_yaxis_w)), 
        name = "frame"))
    if (scale) {
        pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2, 
            xscale = c(0, max_length), name = "scale"))
        grid.draw(scale_grob)
        upViewport()
    }
    if (!is.null(tree_grob)) {
        annot_margin <- unit(if (is.null(annotations[[1]])) 
            0
        else annotation_height, "lines")
        seg_plot_margin <- unit(if (is.null(seg_plot_grobs[[1]][[1]])) 
            0
        else seg_plot_height, seg_plot_height_unit)
        hli <- unit(0.5, "lines")
        dna_scale_margin <- unit(scale_cex * dna_seg_scale[[n_dna_segs]], 
            "lines")
        pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, 
            layout = grid.layout(6, 1, heights = unit.c(seg_plot_margin, 
                annot_margin, hli, unit(n_dna_segs * (1 + seg_plot_height), 
                  "null"), hli, dna_scale_margin)), name = "tree_outer"))
        pushViewport(viewport(layout.pos.row = 4, width = unit(1, 
            "npc") - unit(1, "lines"), just = c("centre", "bottom"), 
            name = "tree"))
        grid.draw(tree_grob$grob)
        upViewport(2)
    }
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, 
        name = "plotarea_outer", clip = "on"), viewport(width = unit(1, 
        "npc") - unit(1, "lines"), height = unit(1, "npc") - 
        unit(0, "lines"), name = "plotarea", clip = "off"))
    pushViewport(viewport(layout = grid.layout(n_rows, 1, heights = dna_seg_heights), 
        name = "map"))
    if (n_comparisons > 0) {
        for (i in 1:n_comparisons) {
            pushViewport(viewport(layout.pos.row = 3 * i, yscale = c(0, 
                1), xscale = c(0, max_length), name = paste("comparison", 
                i, sep = ".")))
            grid.draw(comparison_grobs[[i]])
            upViewport()
        }
    }
    for (i in 1:n_dna_segs) {
        n_dna_subsegs <- length(dna_subsegs[[i]])
        n_cols <- n_dna_subsegs * 2 + 1
        widths <- numeric(n_cols)
        widths[1:n_dna_subsegs * 2] <- xlims[[i]]$length
        widths[1:n_dna_subsegs * 2 - 1] <- offsets[[i]]
        widths_units <- unit(widths, rep("native", n_cols))
        heights <- unit(c(annot_h[i], 1, scale_cex * dna_seg_scale[i]), 
            c("lines", "null", "lines"))
        pushViewport(viewport(layout.pos.row = 3 * i - 2, layout = grid.layout(1, 
            n_cols, widths = widths_units, just = c("left", "centre")), 
            xscale = c(0, max_length), name = paste("seg_plot", 
                i, sep = ".")))
        for (j in 1:n_dna_subsegs) {
            idx <- if (xlims[[i]]$strand[j] == 1) 
                c("x0", "x1")
            else c("x1", "x0")
            xscale <- as.numeric(xlims[[i]][j, idx])
            if (!is.null(seg_plots[[i]])) {
                pushViewport(viewport(layout.pos.col = j * 2, 
                  yscale = seg_plot_ylims[[i]], xscale = xscale, 
                  just = c("left", "centre"), name = paste("seg_subplot", 
                    i, j, sep = ".")))
                grid.draw(seg_plot_grobs[[i]][[j]])
                upViewport()
            }
        }
        if (!is.null(seg_plots[[i]]) && seg_plot_yaxis > 0) {
            pushViewport(viewport(layout.pos.col = n_cols, yscale = seg_plot_ylims[[i]], 
                width = unit(1, "grobwidth", data = seg_plot_yaxis_grobs[[i]]), 
                just = c("left", "centre"), name = paste("seg_plot_yaxis", 
                  i, sep = ".")))
            grid.draw(seg_plot_yaxis_grobs[[i]])
            upViewport()
        }
        upViewport()
        pushViewport(viewport(layout.pos.row = 3 * i - 1, layout = grid.layout(3, 
            n_cols, heights = heights, widths = widths_units, 
            just = c("left", "centre")), xscale = c(0, max_length), 
            name = paste("scale_and_dna_seg", i, sep = ".")))
        for (j in 1:n_dna_subsegs) {
            idx <- if (xlims[[i]]$strand[j] == 1) 
                c("x0", "x1")
            else c("x1", "x0")
            xscale <- as.numeric(xlims[[i]][j, idx])
            if (!is.null(annotations[[i]])) {
                pushViewport(viewport(layout.pos.row = 1, layout.pos.col = j * 
                  2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                  "centre"), name = paste("annotation", i, j, 
                  sep = ".")))
                grid.draw(annotation_grobs[[i]][[j]])
                upViewport()
            }
            if (dna_seg_scale[i]) {
                pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j * 
                  2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                  "centre"), name = paste("dna_seg_scale", i, 
                  j, sep = ".")))
                grid.draw(dna_seg_scale_grobs[[i]][[j]])
                upViewport()
            }
            pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j * 
                2, yscale = c(0, 1), xscale = xscale, just = c("left", 
                "centre"), name = paste("dna_seg", i, j, sep = ".")))
            if (!dna_seg_line[i] == "FALSE") {
                grid.segments(x0 = unit(xlims[[i]]$x0[j], "native"), 
                  y0 = unit(0.5, "native"), x1 = unit(xlims[[i]]$x1[j], 
                    "native"), y1 = unit(0.5, "native"), name = paste("dna_seg_line", 
                    i, j, sep = "."), gp = gpar(col = dna_seg_line[i]))
            }
            grid.draw(dna_seg_grobs[[i]][[j]])
            upViewport()
            if (j > 1) {
                pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j * 
                  2 - 1, yscale = c(0, 1), xscale = c(0, widths[2 * 
                  j - 1]), just = c("centre", "centre"), name = paste("gap", 
                  i, j, sep = ".")))
                grid.draw(gap_grob(w = def_gap_length, m = widths[2 * 
                  j - 1]/2, i, j))
                upViewport()
            }
        }
        upViewport()
    }
    upViewport(2)
    upViewport(2)
    upViewport(2)
}