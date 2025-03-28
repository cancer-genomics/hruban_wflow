library(tidyverse)
library(VariantAnnotation)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(foreach)
library(doParallel)

args <- commandArgs(trailingOnly = T)
strelkaSbsVcf <- args[1]
plasmaID <- args[2]
plasmaAln <- args[3]
outDir <- args[4]
nProcesses <- args[5] %>% as.numeric

message(plasmaID)

# Setting up the parallel computing environment
cl <- makeCluster(nProcesses)
registerDoParallel(cl)

# Reading in SBS mutations
sbs.vcf <- readVcf(strelkaSbsVcf)

# Filtering for mutations that passed Strelka QC
sbs.vcf <- sbs.vcf[rowRanges(sbs.vcf)$FILTER == "PASS"]
muts <- rowRanges(sbs.vcf)

# Creating a data.frame with information about each variant
getNucCts <- function(tier, vcf) {
  stopifnot(tier %in% c(1, 2))

  tumor.A <- geno(sbs.vcf)[["AU"]][, , tier][, "TUMOR"]
  tumor.T <- geno(sbs.vcf)[["TU"]][, , tier][, "TUMOR"]
  tumor.C <- geno(sbs.vcf)[["CU"]][, , tier][, "TUMOR"]
  tumor.G <- geno(sbs.vcf)[["GU"]][, , tier][, "TUMOR"]
  normal.A <- geno(sbs.vcf)[["AU"]][, , tier][, "NORMAL"]
  normal.T <- geno(sbs.vcf)[["TU"]][, , tier][, "NORMAL"]
  normal.C <- geno(sbs.vcf)[["CU"]][, , tier][, "NORMAL"]
  normal.G <- geno(sbs.vcf)[["GU"]][, , tier][, "NORMAL"]

  ct.df <- data.frame(tumor_A = tumor.A,
                      tumor_T = tumor.T,
                      tumor_C = tumor.C,
                      tumor_G = tumor.G,
                      normal_A = normal.A,
                      normal_T = normal.T,
                      normal_C = normal.C,
                      normal_G = normal.G) 
  return(ct.df)
}

tier1_ct.df <- getNucCts(tier = 1, vcf = sbs.vcf)
tier2_ct.df <- getNucCts(tier = 2, vcf = sbs.vcf)

sbs.df <- data.frame(variant_id = names(muts),
                     chrom = seqnames(muts) %>% as.character,
                     pos = start(muts),
                     ref = muts$REF,
                     alt = unlist(muts$ALT) %>% as.character)

# Computing tier1 coverage metrics
sbs.df$tumor_coverage_tier1 = rowSums(tier1_ct.df[,c("tumor_A", "tumor_T", "tumor_C", "tumor_G")])
sbs.df$tumor_ref_ct_tier1 <- tier1_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("tumor_", sbs.df$ref), colnames(tier1_ct.df)))]
sbs.df$tumor_alt_ct_tier1 <- tier1_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("tumor_", sbs.df$alt), colnames(tier1_ct.df)))]
sbs.df$tumor_maf_tier1 <- sbs.df$tumor_alt_ct_tier1 / sbs.df$tumor_coverage_tier1

sbs.df$normal_coverage_tier1 = rowSums(tier1_ct.df[,c("normal_A", "normal_T", "normal_C", "normal_G")])
sbs.df$normal_ref_ct_tier1 <- tier1_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("normal_", sbs.df$ref), colnames(tier1_ct.df)))]
sbs.df$normal_alt_ct_tier1 <- tier1_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("normal_", sbs.df$alt), colnames(tier1_ct.df)))]
sbs.df$normal_maf_tier1 <- sbs.df$normal_alt_ct_tier1 / sbs.df$normal_coverage_tier1

# Computing tier2 coverage metrics
sbs.df$tumor_coverage_tier2 = rowSums(tier2_ct.df[,c("tumor_A", "tumor_T", "tumor_C", "tumor_G")])
sbs.df$tumor_ref_ct_tier2 <- tier2_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("tumor_", sbs.df$ref), colnames(tier2_ct.df)))]
sbs.df$tumor_alt_ct_tier2 <- tier2_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("tumor_", sbs.df$alt), colnames(tier2_ct.df)))]
sbs.df$tumor_maf_tier2 <- sbs.df$tumor_alt_ct_tier2 / sbs.df$tumor_coverage_tier2

sbs.df$normal_coverage_tier2 = rowSums(tier2_ct.df[,c("normal_A", "normal_T", "normal_C", "normal_G")])
sbs.df$normal_ref_ct_tier2 <- tier2_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("normal_", sbs.df$ref), colnames(tier2_ct.df)))]
sbs.df$normal_alt_ct_tier2 <- tier2_ct.df[cbind(1:nrow(sbs.df) ,match(paste0("normal_", sbs.df$alt), colnames(tier2_ct.df)))]
sbs.df$normal_maf_tier2 <- sbs.df$normal_alt_ct_tier2 / sbs.df$normal_coverage_tier2

# Annotating each variant observed in tissue with its allele frequency in gnomAD
source("../functions/gnomAD_annotate.R")

sbs.df$liftover_success <- as.character(NA)
sbs.df$in_gnomad <- as.character(NA)
sbs.df$gnomad_pass <- as.character(NA)
sbs.df$gnomad_af <- as.numeric(NA)

m <- gnomAD_annotate(chrom = sbs.df$chrom,
                     pos = sbs.df$pos,
                     ref = sbs.df$ref,
                     alt = sbs.df$alt)

sbs.df[, c("liftover_success", "in_gnomad", "gnomad_pass", "gnomad_af")] <- m

# Generate pileup in plasma bam file at positions harboring known somatic mutations
pyFile <- "../functions/pileup_runner.py"

regions <- paste0(seqnames(muts), ":", start(muts), "-", end(muts))
stopifnot(length(unique(regions)) == length(regions))
regions <- split(regions, cut(seq_along(regions), nProcesses, labels = FALSE))

res <- foreach(i = 1:nProcesses) %dopar% {
  pileupFile <- file.path(outDir, "pileup", paste0(plasmaID, "_process", i, ".txt"))
  regionFile <- file.path(outDir, "regions", paste0(plasmaID, "_process", i, ".txt"))
  writeLines(regions[[i]], con = regionFile)
  pyCmd <- paste("python3", pyFile, plasmaAln, regionFile, pileupFile)
  system(pyCmd)
}

# Reading in the pileup files and coverting so that each line is a unique position/fragment combination 
# with read 1 and read 2 observations on the same line
pu.df <- file.path(outDir, "pileup", paste0(plasmaID, "_process", 1:nProcesses, ".txt")) %>% 
           lapply(., function(t) read_tsv(t, col_types = cols())) %>% 
           do.call("rbind", .) %>%
           mutate(read = paste0("read", read)) %>%
           pivot_wider(names_from = read, 
                       values_from = c(strand, obs, mapq, phred),
                       names_glue = "{read}_{.value}") %>%
           dplyr::select(fragment, chrom, pos, read1_obs, read1_strand, read1_phred, read1_mapq, read2_obs, 
                         read2_strand, read2_phred, read2_mapq) %>%
           rename_at(vars(-chrom, -pos), ~ paste0("plasma_", .x))

# Annotating each variant observed in plasma with the number of 'N' bases in each read
fragmentFile <- file.path(outDir, "fragments", paste0(plasmaID, ".txt"))
plasmaFragments <- pu.df %>% 
                     filter(!is.na(plasma_fragment)) %>% 
                     .$plasma_fragment
writeLines(plasmaFragments, con = fragmentFile)

# Counting the number of 'N' bases in each read
pyFile2 <- "../functions/count_n_bases.py"
ctFile <- file.path(outDir, "n_bases", paste0(plasmaID, ".txt"))
pyCmd <- paste("python3", pyFile2, plasmaAln, fragmentFile, ctFile)
system(pyCmd)

# Reading in the counts for the number of 'N' bases in each read
ct.list <- read_tsv(ctFile) %>%
             split(.$read)

# Adding the number of 'N' bases in each read to obs.df
pu.df <- pu.df %>%
           left_join(., ct.list[["1"]], by = c("plasma_fragment" = "fragment")) %>%
           dplyr::select(-read) %>%
           relocate(plasma_read1_n_bases = n_bases, .after = "plasma_read1_mapq") %>%
           left_join(., ct.list[["2"]], by = c("plasma_fragment" = "fragment")) %>%
           dplyr::select(-read) %>%
           relocate(plasma_read2_n_bases = n_bases, .after = "plasma_read2_mapq")

# Joining sbs.df and pu.df
summary.df <- left_join(sbs.df, pu.df, by = c("chrom", "pos"))

write.table(summary.df, file = file.path(outDir, "obs", paste0(plasmaID, ".txt")),
            row.names = F, col.names = T, sep = "\t", quote = F)



