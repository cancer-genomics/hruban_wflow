# Defining a function to annotate variants with their allele frequency in gnomAD
gnomAD_annotate <- function(chrom, pos, ref, alt) {

  # Creating a data.frame to store the output data
  out.df <- data.frame(liftover.success = rep("Y", length(chrom)),
                       in.gnomad = rep("N", length(chrom)),
                       gnomad.pass = as.character(rep(NA, length(chrom))),
                       gnomad.af = as.numeric(rep(NA, length(chrom))), stringsAsFactors = FALSE)

  # Storing the variants as a GRanges 'mut.gr'
  mut.gr <- GRanges(seqnames = chrom, ranges = IRanges(start = pos, end = pos))

  # Defining the path to the hg19 to hg38 chain for lifting over variants
  hg19ToHg38ChainPath <- "../data/hg19ToHg38.over.chain"
  hg19ToHg38Chain <- rtracklayer::import.chain(hg19ToHg38ChainPath)

  # Lifting over variants from hg19 to hg38 to match gnomAD
  hg38.mut.grl <- rtracklayer::liftOver(mut.gr, hg19ToHg38Chain) # Note that some positions will not lift over (will mark as N)
  names(hg38.mut.grl) <- 1:nrow(out.df) # The name will be the index of the variant passed to this function
  n.liftover <- elementNROWS(hg38.mut.grl)

  no.liftover <- which(n.liftover == 0)
  if (length(no.liftover) > 0) {
    out.df$liftover.success[no.liftover] <- "N" # Marking variants that don't lift over to hg38
  }

  many.liftover <- which(n.liftover > 1)
  if (length(many.liftover) > 0) {
    out.df$liftover.success[many.liftover] <- "M" # Marking variants that lift over multiple times to hg38
    hg38.mut.grl <- hg38.mut.grl[-many.liftover]  # and excluding them from the gnomAD annotation stage
  }

  hg38.mut.gr <- unlist(hg38.mut.grl)
  if (length(hg38.mut.gr) == 0) return(out.df)
  rm(list = c("hg38.mut.grl", "n.liftover", "no.liftover", "many.liftover"))

  # Marking variants if they lifted over to one position but the reference base is different between the two builds
  hg38.mut.gr$REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, hg38.mut.gr))
  different.ref <- hg38.mut.gr$REF != ref[as.numeric(names(hg38.mut.gr))]
  if (any(different.ref)) {
    different.ref.ind <- which(different.ref)
    out.df$liftover.success[as.numeric(names(hg38.mut.gr))[different.ref.ind]] <- "R"
    hg38.mut.gr <- hg38.mut.gr[-different.ref.ind]
  }
  if (length(hg38.mut.gr) == 0) return(out.df)
  rm(list = c("different.ref"))

  hg38.mut.gr$ALT <- alt[as.numeric(names(hg38.mut.gr))]

  # Getting gnomAD variant postions
  gnomad <- readRDS("../data/gnomad.genomes.r3.0.sites.sbs.rds")

  olaps <- findOverlaps(hg38.mut.gr, gnomad, type = "equal")
  hg38.mut.gr <- hg38.mut.gr[queryHits(olaps)]
  tmp.gnomad <- gnomad[subjectHits(olaps)]

  is.same.alt <- hg38.mut.gr$ALT == tmp.gnomad$ALT

  if (any(is.same.alt)) {
    hits <- as.numeric(names(hg38.mut.gr))[is.same.alt]
    out.df$in.gnomad[hits] <- "Y"
    out.df$gnomad.pass[hits] <- ifelse(tmp.gnomad$FILTER[is.same.alt] == "P", "Y", "N")
    out.df$gnomad.af[hits] <- tmp.gnomad$AF[is.same.alt]
  }

  return(out.df)

} # End of gnomAD_annotate

