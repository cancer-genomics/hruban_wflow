library(tidyverse)
library(patchwork)

# find the data type of variable(s) used to train the model
get_data_type <- function(train_data, coxph_mf, sel_col_name=NULL) {
    col_names <- names(train_data)
    data_classes <- attr(terms(coxph_mf), "dataClasses")
    # which variables in the data were used to train the model?
    col_hits <- match(col_names, names(data_classes))
    col_hits <- col_hits[!is.na(col_hits)]
    # reverse of the operation above - to create the HR table with the same sequence of variables as in the cox model
    ord_col_hits <- match(names(data_classes), col_names)
    ord_col_names <- names(data_classes)[!is.na(ord_col_hits)]
    # return the data type or recursively find the data type of all the variables used in the model
    if(!is.null(sel_col_name)) {
        if(sel_col_name %in% names(data_classes[col_hits])) return(data_classes[[sel_col_name]]) else stop("Selected column name is not present in the model object!")
    } else {
        sapply(ord_col_names, function(x) get_data_type(train_data, coxph_mf, x))
    }
}

# find the number of levels for data type of variable
# <!> the number of levels are based on the data and not the model, so ensure that you are using the data that is used to train the model 
get_nlevels <- function(train_data, coxph_mf, sel_col_name=NULL) {
    if(!is.null(sel_col_name)) {
        col_data_type <- get_data_type(train_data, coxph_mf, sel_col_name) 
        # assign/compute the number of levels depending on the data type
        switch(col_data_type,
               "numeric" = 1,
               "logical" = 2,
               "character" = length(unique(train_data[[sel_col_name]])),
               "factor" = length(levels(train_data[[sel_col_name]])))
    } else {
        # this step is wasteful, but meh!
        all_col_data_type <- get_data_type(train_data, coxph_mf)
        # find the number of levels for all the variables used in the model
        sapply(names(all_col_data_type), function(x) get_nlevels(train_data, coxph_mf, x)) 
    }
}

# extract the HR, conf intervals and p-value for variable(s)
get_hr_tbl <- function(train_data, coxph_mf, sel_col_name=NULL, status="death") {

    # cox ph object coefficients and confidence intervals
    summary_tbl <- summary(coxph_mf)$conf.int
    coef_tbl <- summary(coxph_mf)$coefficients

    # if only extracting the information for the selected variable
    # else extract the information for all the variables
    if(!is.null(sel_col_name)) {

        # counting events and totals
        col_data_type <- get_data_type(train_data, coxph_mf, sel_col_name)
        col_nlevels <- get_nlevels(train_data, coxph_mf, sel_col_name)
        if(col_nlevels > 1) {
            # for logical, character or factor data type
            nevent_tbl <- table(train_data[train_data[[status]] == 1, sel_col_name])
            ntotal_tbl <- table(train_data[[sel_col_name]])
        } else {
            # for numeric data type 
            nevent_tbl <- table(ifelse(train_data[[status]] == 1, sel_col_name, NA))
            ntotal_tbl <- table(ifelse(!is.na(train_data[[status]]), sel_col_name, NA))
        }

        # initializing the HR, confidence intervals, pvalues
        HRvec <- rep(1,col_nlevels)
        HRlowvec <- rep(NA,col_nlevels)
        HRhigvec <- rep(NA,col_nlevels)
        pvaluevec <- rep(NA,col_nlevels)
        neventvec <- rep(NA,col_nlevels)
        ntotalvec <- rep(NA,col_nlevels)

        # suffices for levels
        suffix_nlevels <- switch(col_data_type,
                                 "numeric" = "",
                                 "logical" = c("FALSE", "TRUE"),
                                 "character" = unique(train_data[[sel_col_name]]),
                                 "factor" = levels(train_data[[sel_col_name]]))

        # extracting values from the cox ph object
        for (i_num in seq_along(suffix_nlevels)) {
            i <- suffix_nlevels[i_num]
            i_full <- paste0(sel_col_name, i) # rownames used in the summary table of cox model
            # for each level, fill the HR, interval and p-value
            if(paste0(sel_col_name, i) %in% rownames(summary_tbl)) {
                HRvec[i_num] <- summary_tbl[i_full,1]
                HRlowvec[i_num] <- summary_tbl[i_full,3]
                HRhigvec[i_num] <- summary_tbl[i_full,4]
                pvaluevec[i_num] <- coef_tbl[i_full,5]
            }
            # if logical, character or factor, use the suffix
            # else use the variable name
            # or keep it as 0
            if(i %in% names(ntotal_tbl)) {
                neventvec[i_num] <- nevent_tbl[i]
                ntotalvec[i_num] <- ntotal_tbl[i]
            } else if(i_full %in% names(ntotal_tbl)) { 
                neventvec[i_num] <- nevent_tbl[i_full]
                ntotalvec[i_num] <- ntotal_tbl[i_full]
            } else {
                neventvec[i_num] <- 0 
                ntotalvec[i_num] <- 0 
            }
        }

        # adding a header row
        HRvec <- c(NA, HRvec)
        HRlowvec <- c(NA, HRlowvec)
        HRhigvec <- c(NA, HRhigvec)
        pvaluevec <- c(NA, pvaluevec)
        neventvec <- c(NA, neventvec)
        ntotalvec <- c(NA, ntotalvec)
        suffix_nlevels <- c("", suffix_nlevels)

        # createing a table for selected variable
        tbl_row <- tibble(var = paste0(sel_col_name, suffix_nlevels),
                          HR = HRvec,
                          HRlow = HRlowvec,
                          HRhig = HRhigvec,
                          pvalue = pvaluevec,
                          nevent = neventvec,
                          ntotal = ntotalvec)
    } else {
        all_col_nlevels <- get_nlevels(train_data, coxph_mf)
        lapply(names(all_col_nlevels), function(x) get_hr_tbl(train_data, coxph_mf, x, status)) %>% bind_rows()
    }
}

# add 0 to the right
r_pad <- function(x, n) sprintf(paste0("%.", n, "f"), round(x, n))
    
# make the forest plot
plot_hr <- function(hr_tbl, file_name = "coxph_table.png", var_labels=NULL, xlim=NULL, xticks_rotate=FALSE, scientific=FALSE, fraction=FALSE, fig_height = 7) {

    # find the title rows and fix shape of the point and interval
    title_row <- which(is.na(hr_tbl$ntotal))
    hr_tbl$shape <- 15
    hr_tbl$se <- 1
    hr_tbl$shape[title_row] <- NA
    hr_tbl$se[title_row] <- 0

    if(is.null(var_labels)) var_labels <- hr_tbl$var # keep the var names as is if nothing is provided

    # adding a heading row to the table - this table will be used to make the intervals plot and put "in" between the actual table
    hr_tbl_in <- hr_tbl %>%
        add_row(var = NA, HR = NA, HRlow = NA, HRhig = NA, pvalue = NA, nevent = NA, ntotal = NA, shape = NA, se = NA, .before=1)

    # clean up the headings and reorganize the table for printing
    hr_tbl_disp <- hr_tbl %>%
        mutate(var = var_labels,
               #events_total = paste(nevent, ntotal, sep="/"),
               events_total = ntotal, 
               hru = ifelse(is.na(HR), "", ifelse(is.na(HRlow), HR, paste0(r_pad(HR,2), " (", r_pad(HRlow,2),"-",r_pad(HRhig,2),")"))),
               pvalue = sapply(pvalue, function(x) ifelse(is.na(x), "", pvalue_print(x))),
               events_total = ifelse(row_number() %in% title_row, "", events_total)) %>% 
    select(var, events_total, hru, pvalue) %>%
    rename(Subgroup = var,
           #`Events/Patients` = events_total,
           `Patients` = events_total,
           `Hazard Ratio (95% CI)` = hru,
           `P value` = pvalue)

    # adding a heading row to the table for display
    hr_tbl_lr <- hr_tbl_disp %>%
        add_row(Subgroup = "Subgroup",
#                `Events/Patients` = "Events/Patients",
                `Patients` = "Patients",
                `Hazard Ratio (95% CI)` = "Hazard Ratio (95% CI)",
                `P value` = "P value", .before= 1) %>% 
    mutate(row_order = row_number(),
           row_order = factor(row_order, levels = rev(row_order)))

    # defining margins - no method here - just trial and error
    pos1 <- 0; mar1 <- pos1-0.1
    pos2 <- 4; mar2 <- pos2+0.5
    pos3 <- 2; mar3 <- pos3-3.0
    pos4 <- 4; mar4 <- pos4+1.0

    # overwrite title_row as a new row has been added
    title_row <- which(is.na(hr_tbl_in$ntotal))

    # keep variable names and events / patients on the left side
    fig_left <- ggplot(data = hr_tbl_lr)+
        geom_text(aes(x = pos1, y = row_order, label = Subgroup, 
                      fontface = ifelse(row_order %in% title_row, "bold", "plain")), hjust=0)+
    #geom_text(aes(x = pos2, y = row_order, label = `Events/Patients`,
    geom_text(aes(x = pos2, y = row_order, label = `Patients`,
                  fontface = ifelse(row_order %in% title_row, "bold", "plain")), hjust=1)+
    theme_void()+
    coord_cartesian(xlim = c(mar1, mar2))

    # keep hazard ratios and p-value on the right side
    fig_right <- ggplot(data = hr_tbl_lr)+
        geom_text(aes(x = pos3, y = row_order, label = `Hazard Ratio (95% CI)`,
                      fontface = ifelse(row_order %in% title_row, "bold", "plain")), hjust=1)+
        geom_text(aes(x = pos4, y = row_order, label = `P value`,
                      fontface = ifelse(row_order %in% title_row, "bold", "plain")), hjust=1)+
        theme_void()+
        coord_cartesian(xlim = c(mar3, mar4))
    
    # keeping the HR table in the middle
    if(is.null(xlim)) {
        xmin <- round(log(min(hr_tbl_in$HRlow, na.rm=TRUE))/log(2))
        xmax <- round(log(max(hr_tbl_in$HRhig, na.rm=TRUE))/log(2))
        xlim <- c(xmin, xmax)
        xrange <- pretty(c(xmin, xmax))
    } else {
        xrange <- pretty(c(xlim[1], xlim[2]))
    } 
    
    # forest plot
    fig <- hr_tbl_in %>%
        mutate(var_num = row_number(), 
               var_num = factor(var_num, levels = rev(var_num)),
               se = as.character(se),
               log_HR = log(HR)/log(2),
               log_HRlow = log(HRlow)/log(2),
               log_HRhig = log(HRhig)/log(2)) %>%
        replace_na(list(HR=0, HRlow=0, HRhig=0)) %>%
        ggplot(aes(x = log_HR, y = var_num))+
        geom_point(aes(shape=shape))+
        geom_errorbar(aes(xmin = log_HRlow, xmax = log_HRhig, linetype=se, width=0.2))+
        scale_linetype_manual(breaks = c("0","1"), values = c("blank", "solid"))+
        geom_vline(xintercept = 0, linetype="dotted")+
        scale_shape_identity()

    # if the range is too long, use the scientific notation - works well with text rotation
    if(!scientific | fraction) {
        if(fraction) {
            labels_xrange <- c()
            for(xunit in xrange) {
                if(xunit < 0) {
                    labels_xrange <- c(labels_xrange, frac_str(1, 2^(-xunit)))
                } else if(xunit == 0) {
                    labels_xrange <- c(labels_xrange, 1)
                } else {
                    labels_xrange <- c(labels_xrange, 2^xunit)
                }
            }
        } else {
            labels_xrange <- 2^xrange
        } 

        fig <- fig+
            scale_x_continuous(breaks = xrange, labels = labels_xrange)
    } else {
        fig <- fig+
            scale_x_continuous(breaks = xrange, labels = format(2^xrange, scientific=TRUE, digits=2))
    }

    # remove legend and apply x-axis limits
    fig <- fig+
        coord_cartesian(xlim = xlim)+
        labs(x="", y="")+
        theme_classic()+
        theme(legend.position="none",
              axis.line.y= element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    
    if(xticks_rotate) fig <- fig+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # arrange the different panels
    layout <- c(
                area(t=0, b=7, l=0, r=3),  # left of forest plot
                area(t=0, b=7, l=4, r=8), # forest plot
                area(t=0, b=7, l=9, r=12), # right of forest plot
                area(t=0, b=7, l=0, r=12)  # to make a line after the heading
    )
    
    # to plot a line under the title row
    line_pos <- mean(as.numeric(hr_tbl_in %>% 
                     mutate(var_num = row_number(), 
                            var_num = factor(var_num, levels = rev(var_num))) %>%
                     select(var_num) %>% 
                     slice_head(n=2) %>% 
                     pull(var_num)))

    # plotting the title row
    fig_title <- hr_tbl_in %>%
        mutate(var_num = row_number(), 
               var_num = factor(var_num, levels = rev(var_num))) %>%
        select(var_num) %>% 
        mutate(shape = NA) %>%
        ggplot()+
        geom_point(aes(x = shape, y = var_num, shape = shape))+
        scale_shape_identity()+
        geom_hline(yintercept = line_pos)+
        theme_void()

    # final plot arrangement
    intg_fig <- fig_left+fig+fig_right+fig_title+plot_layout(design = layout)
    
#    if(dirname(file_name) == ".") {
#        dir.create("results")
#        ggsave(file.path("results", paste0(file_name, ".png")), intg_fig, height = fig_height, width = 12)
#        # there must be a better way, but NP
#    } else {
#        ggsave(paste0(file_name, ".png"), intg_fig, height = fig_height, width = 12)
#    }
    
    invisible(intg_fig)

}

# To create a vertical fraction
frac_str <- function(x,y) {
    parse(text = paste0("frac(", x, ",", y, ")"))
}

# To print p-values upto 4 decimal places or switch to scientific format
pvalue_print <- function(x, type = c("categ", "round")) {
  type <- match.arg(type)
  if (x < 0.0001) {
    if (type == "round") {
      xf <- round(x, 5)
    } else {
      xf <- "< 0.0001"
    }
  } else if (x < 0.001) {
    if (type == "round") {
      xf <- round(x, 4)
    } else {
      xf <- "< 0.001"
    }
  } else if (x < 0.01) {
    if (type == "round") {
      xf <- round(x, 3)
    } else {
      xf <- "< 0.01"
    }
  } else {
    if (type == "round") {
      xf <- round(x, 2)
    } else {
      xf <- as.character(round(x, 2))
    }
  }
  return(xf)
}
