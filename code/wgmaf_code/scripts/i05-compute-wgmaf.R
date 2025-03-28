library(tidyverse)
library(readxl)
library(tidyverse)

manifest <- read_tsv("../data/manifest.txt")

# Definining the directory that holds the mutation observation files
obsDir <- "../outDir/plasma_obs/obs"
obsFiles <- list.files(obsDir, full.names = T)

# Defining a function for filtering observations
filterObs <- function(df, minPhred, minMAPQ, minCovTumorTier1, minCovNormalTier1, minAltCovTumorTier1, minMafTumorTier1, maxMafNormalTier1, maxMafNormalTier2, maxNBases,
                      tnFilters, plasmaFilters) {

  if (tnFilters == T) {
  df <- df %>% 
    filter(tumor_coverage_tier1 >= minCovTumorTier1,
           normal_coverage_tier1 >= minCovNormalTier1,
           tumor_alt_ct_tier1 >= minAltCovTumorTier1,
           tumor_maf_tier1 >= minMafTumorTier1,
           normal_maf_tier1 <= maxMafNormalTier1,
           normal_maf_tier2 <= maxMafNormalTier2,
           liftover_success == "Y",
           is.na(gnomad_pass) | gnomad_pass == "Y",
           is.na(gnomad_af) | gnomad_af <= 1/1e4)
  }

  if (plasmaFilters == T) {
    df <- df %>%
      filter(!is.na(plasma_read1_obs),
             !is.na(plasma_read2_obs),
             plasma_read1_obs == plasma_read2_obs,
             plasma_read1_phred >= minPhred,
             plasma_read2_phred >= minPhred,
             plasma_read1_mapq >= minMAPQ,
             plasma_read2_mapq >= minMAPQ,
             plasma_read1_n_bases <= maxNBases,
             plasma_read2_n_bases <= maxNBases)
  }

  return(df)

}

plasma_ids <- filter(manifest, sample_type == "plasma")$sample_id

maf.df <- NULL

for (i in 1:length(plasma_ids)) {

  obs.df <- read_tsv(file.path(obsDir, paste0(plasma_ids[i], ".txt")), col_types = cols()) %>%
              mutate(sample_id = plasma_ids[i],
                     .before = variant_id)

  # Apply filters
  obs.df <- filterObs(obs.df, minPhred = 30, minMAPQ = 42, minCovTumorTier1 = 0, minCovNormalTier1 = 30, minAltCovTumorTier1 = 3, 
                      minMafTumorTier1 = 0.1, maxMafNormalTier1 = 0, maxMafNormalTier2 = 0.05, maxNBases = 4,
                      tnFilters = T, plasmaFilters = T)

  ref.obs <- obs.df %>%
               filter(plasma_read1_obs == ref) %>%
               summarise(n_ref = n()) 

  alt.df <- obs.df %>%
              filter(plasma_read1_obs == alt)
  alt.obs <- alt.df  %>% 
               summarise(n_alt = n()) 
  
  maf.df.i <- bind_cols(ref.obs, alt.obs) %>%
                mutate(wgmaf = n_alt / (n_alt + n_ref)) %>%
                mutate(sample_id = plasma_ids[i], .before = n_ref)

  maf.df <- bind_rows(maf.df, maf.df.i)
}

# Setting values of NaN (when there are 0 plasma observations) to 0
if (any(is.nan(maf.df$wgmaf))) {
  maf.df$wgmaf[is.nan(maf.df$wgmaf)] <- 0
}

# Saving the MAFs for each sample at each timepoint
outDir <- "../outDir"
if (!dir.exists(outDir)) dir.create(outDir)
write.table(maf.df, file = file.path(outDir, "plasma_wgmafs.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)



