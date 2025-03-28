get_cnv <- function(facet_summary_filename){
  facet_summary <- read_tsv(facet_summary_filename, show_col_types=FALSE)
  facet_summary<-facet_summary %>% mutate(length=end-start)
  facet_summary<-facet_summary %>% mutate(total_score=length*cnlr.median)
  facet_summary<-facet_summary %>%
      mutate(chrarm=case_when((chrom==1 & start>=152500000) ~ "1q",
      (chrom==1 & end<=152500000) ~ "1p",
      (chrom==2 & start>=102500000) ~ "2q",
      (chrom==2 & end<=102500000) ~"2p",
      (chrom==3 & start>=97500000) ~"3q",
      (chrom==3 & end<=97500000) ~ "3p",
      (chrom==4 & start>=57500000) ~"4q",
      (chrom==4 & end<=57500000) ~ "4p",
      (chrom==5 & start>=52500000) ~"5q",
      (chrom==5 & end<=52500000) ~"5p",
      (chrom==6 & start>=67500000) ~"6q",
      (chrom==6 & end<=67500000) ~"6p",
      (chrom==7 & start>=67500000) ~ "7q",
      (chrom==7 & end<=67500000) ~ "7p",
      (chrom==8 & start>=52500000) ~"8q",
      (chrom==8 & end<=52500000) ~"8p",
      (chrom==9 & start>=77500000) ~ "9q",
      (chrom==9 & end<=77500000) ~"9p",
      (chrom==10 & start>=57500000) ~"10q",
      (chrom==10 & end<=57500000) ~"10p",
      (chrom==11 & start>=57500000) ~"11q",
      (chrom==11 & end<=57500000) ~"11p",
      (chrom==12 & start>=42500000)~"12q",
      (chrom==12 & end<=42500000)~"12p",
      (chrom==13 & start>=22500000)~"13q",
      (chrom==14 & start>=22500000)~"14q",
      (chrom==15 & start>=27500000)~"15q",
      (chrom==16 & start>=47500000)~"16q",
      (chrom==16 & end<=47500000)~"16p",
      (chrom==17 & start>=32500000)~"17q",
      (chrom==17 & end<=32500000) ~"17p",
      (chrom==18 & start>=22500000)~"18q",
      (chrom==18 & end<=22500000) ~"18p",
      (chrom==19 & start>=32500000)~"19q",
      (chrom==19 & end<=32500000)~"19p",
      (chrom==20 & start>=32500000)~"20q",
      (chrom==20 & end<=32500000)~"20p",
      (chrom==21 & start>=17500000)~"21q",
      (chrom==22 & start>=27500000)~"22q",
      ))
  facet_summary<-facet_summary%>%filter(!is.na(chrarm))
  armcalc_facet_summary <- facet_summary %>%
      group_by(chrarm) %>%
      mutate(cnv_score = (sum(total_score))/(sum(length)))
  armcalc<-distinct(armcalc_facet_summary, chrarm, .keep_all = TRUE)
  armcalc1<-armcalc[,c('chrarm','cnv_score')]
  armcalc1<-armcalc1 %>% mutate(ID=facet_summary_filename)
  armcalc1
}

### Fast PCA using QR decomposition
.svd.qr = function(A, l=2) {
    corr <- cor(A, method="pearson", use="complete.obs")
#     corr <- cor(A, method="spearman")

    center <- rowMeans2(corr, na.rm = TRUE)
    corr <- sweep(corr, 1L, center, check.margin = FALSE) 
#     scale <- rowMeans2(corr, na.rm = TRUE)
#     corr <- sweep(corr, 1L, scale, FUN="/", check.margin = FALSE) 
#     l <- 5
    n <- ncol(corr)
    m <- nrow(corr)

    ## Fast implementation (Projection onto lower dimensional subspace)
    if(TRUE) {
    G <- matrix(rnorm(n * l, 0, 1), nrow = n, ncol = l)
    h.list <- vector("list", 2L)
#     h.list <- vector("list", i + 1L)
    h.list[[1]] <- corr %*% G
    for (j in seq(2, 1 + 1L)) {
        h.list[[j]] <- corr %*% (crossprod(corr, h.list[[j - 1L]]))
    }
#     h.list[[2]] <- corr %*% (crossprod(corr, h.list[[1L]]))
    H <- do.call("cbind", h.list) # n x [(1+1)l] matrix
    # QR algorithm
    Q <- qr.Q(qr(H, 0))
    TT <- crossprod(Q, corr)
    svd <- svd(TT)
    u <- Q %*% svd$u[,1, drop = TRUE]
    }
    ## Slow
#     svd <- svd(corr)
#     u <- svd$u[, 1, drop = TRUE]
    u
}

.standardize <- function(u) {
    u <- u - median(u)
    u <- u/sqrt(sum(u^2))
    u
}

##  For smoothing out vectors for nicer plotting of AB compartments.
##  Credit goes to JP Fortin.
.meanSmoother <- function(x, k=2, iter=2, na.rm=TRUE){
    meanSmoother.internal <- function(x, k=2, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)
        
        window.mean <- function(x, j, k, na.rm=na.rm){
            if (k>=1){
                return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
            } else {
                return(x[j])
            }    
        }
        
        for (i in (k+1):(n-k)){
            y[i] <- window.mean(x,i,k, na.rm)
        }
        for (i in 1:k){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-k+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }
    
    for (i in 1:iter){
        x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
    }
    x
}


.imputeMatrix <- function(matrix) {
  matrix[is.infinite(matrix) & matrix > 0] <- max(matrix[is.finite(matrix)])
  matrix[is.infinite(matrix) & matrix < 0] <- min(matrix[is.finite(matrix)])
  
  # Imputation of the missing values:
  missing <- which(is.na(matrix), arr.ind = TRUE)
  if (length(missing) != 0) {
    for (j in seq_len(nrow(missing))) {
      mean <- mean(
        x = matrix[missing[j, 1L], ][
          is.finite(matrix[missing[j, 1L], ])],
        na.rm = TRUE)
      matrix[missing[j, 1L], missing[j, 2L]] <- mean
    }
  }
  matrix
}

chr_wrangling <- function(dat, sel_chr = "chr14", source_list = NULL,
                          label_list = NULL) {
  if (is.null(source_list)) source_list <- sort(unique(dat$source2))
  if (is.null(label_list)) label_list <- source_list
  end_track <- length(source_list)

  sc <- dat %>%
    filter(chr == sel_chr) %>%
    filter(source2 %in% source_list) %>%
    mutate(source2 = factor(source2, source_list, label_list))

  # the bins that don't agree b/w reference tracks
  mis_bins <- sc %>%
    filter(source2 == label_list[1], lymph_agree == "no") %>%
    pull(bin)

  # trying to do victor's thing of bins with a big change in magnitude
  lihc <- dat %>% filter(source2 == source_list[1])
  lymph <- dat %>% filter(source2 == source_list[end_track])
  bins <- inner_join(lymph, lihc, by = "bin")
  bins$diff <- bins$eigen.x - bins$eigen.y
  bins$z_score <- (bins$diff - mean(bins$diff)) / sd(bins$diff)
  sig_bins <- bins %>% filter(z_score > 1.96 | z_score < -1.96) %>% pull(bin)
  ##p=.05, maybe should be corrected but not really ind tests?

  d <- setdiff(sig_bins, mis_bins)
  ## adding the discordant sign bins to the track and make the informative bins stand out
  d <- append(d, mis_bins)

  track.data <- sc %>%
    filter(chr == sel_chr) %>%
    mutate(#eigen = .meanSmoother(eigen), # having to smooth here is new vs. liver
           color = ifelse(eigen < 0, "gray50", "red4"),
           transp = ifelse(bin %in% d, .95, .1))

  ##now making the other panel thing to show mag diff
  sc <- dat %>% filter(chr == sel_chr)
  ##only the plasma tracks
  sc <- sc %>%
    filter(source2 %in% source_list[3:4])
  sc$source2 <- factor(sc$source2,
                       levels = source_list[3:4],
                       labels = label_list[3:4])

  #coloring the bins as above
  tib <- sc %>%
    filter(chr == sel_chr) %>%
    mutate(#eigen = .meanSmoother(eigen),
           color = ifelse(eigen < 0, "gray50", "red4"),
           transp = ifelse(bin %in% d, .95, .1))

  ##the title thing is weird and it needs a legend also the lines
  ##are too thick, fix later
  tib$source3 <- "Plasma -- White = PAAD, Black = Non-Cancer"
  list(tib = tib, track.data = track.data)
}

#' @export
deconvolution_boxplots <- function(data){
  mis_bins <- data %>%
    filter(source2=="1. AB Compartments PAAD" & lymph_agree=="no")
  mis_bins <- mis_bins$bin
  data2 <- data %>% filter(bin %in% mis_bins)
  ##adding a bunch of columns that could be cool
  lymph_agreement <- data2 %>% group_by(source2,chr) %>%
    count(lymph_agree=="yes") %>%
    filter(`lymph_agree == "yes"`==TRUE)
  lihc_agreement <- data2 %>%
    group_by(source2,chr) %>%
    count(lihc_agree=="yes") %>%
    filter(`lihc_agree == "yes"`==TRUE)
  all_groups <- data2 %>%
    group_by(source2,chr) %>% count()
  lymph_agreement$combo <- paste(lymph_agreement$source2,
                                 lymph_agreement$chr,sep=" ")
  lihc_agreement$combo <- paste(lihc_agreement$source2,
                                lihc_agreement$chr,sep=" ")
  all_groups$combo <- paste(all_groups$source2,
                            all_groups$chr, sep=" ")
  lymph_agreement <- inner_join(lymph_agreement,
                                all_groups,by="combo")
  lymph_agreement$percentage <- lymph_agreement$n.x/lymph_agreement$n.y
  lymph_agreement$type <- "Lymphoblast Agreement"
  lymph_agreement <- lymph_agreement %>%
    select(source2.x,combo,chr.x,percentage,type)
  lihc_agreement <- inner_join(lihc_agreement, all_groups,by="combo")
  lihc_agreement$percentage <- lihc_agreement$n.x/lihc_agreement$n.y
  lihc_agreement$type <- "PAAD Agreement"
  lihc_agreement <- lihc_agreement %>%
    select(source2.x,combo,chr.x,percentage,type)
  agree <- rbind(lihc_agreement,lymph_agreement)
  lihc_agreement <- lihc_agreement %>%
    rename(Liver_Percentage=percentage)
  lymph_agreement <- lymph_agreement %>%
    rename(Lymphoblast_Percentage=percentage)
  perc <- inner_join(lihc_agreement %>%
                       select(Liver_Percentage, combo, source2.x, chr.x),
                     lymph_agreement %>%
                       select(Lymphoblast_Percentage,combo),by="combo")
  perc2 <- perc %>%
    filter(source2.x=="2. Extracted Pancreatic component - Ratios" |
             source2.x=="4. Median Healthy Plasma Ratios" |
             source2.x=="6. Median Pancreatic Plasma Ratios" )
  chrlevels <- paste0("chr", 1:22)
  perc2$chr.x <- factor(perc2$chr.x, chrlevels)
  slevels <- c("2. Extracted Pancreatic component - Ratios",
               "6. Median Pancreatic Plasma Ratios",
               "4. Median Healthy Plasma Ratios")
  perc2 <- perc2 %>%
    mutate(source2.x=factor(source2.x,
                            slevels,
                            c("Estimated Pancreatic Component",
                              "Pancreatic Cancer Plasma",
                              "Healthy Plasma")),
           odds=(Liver_Percentage)/(1-Liver_Percentage))
  perc2
}

get_coefs <- function(model) {
  orig_coefs <- coef(model$finalModel, s = model$bestTune$lambda) * (-1)
  pr <- prep(model$recipe)
  model_input <- suppressWarnings(bake(pr, new_data = model$trainingData))

  feature_means <- model_input  %>%
      select(-c(id, type)) %>%
      colMeans()
  feature_sds <- model_input %>%
      select(-c(id, type)) %>%
      as.data.frame() %>%
      summarise_all(sd)
  feature_coefs <- data.frame(features = names(feature_sds),
                            sd = as.numeric(feature_sds))
  feature_coefs <- merge(feature_coefs,
                     data.frame(features = rownames(orig_coefs),
                                orig_coefs = as.numeric(orig_coefs)),
                     by = 'features', all.x = TRUE)
  feature_coefs$scaled_coefs <- feature_coefs$orig_coefs * feature_coefs$sd
  coefs<-feature_coefs %>% filter(scaled_coefs != 0)
  coefs
}


