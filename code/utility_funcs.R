add.legend <- function(title="", x=0, y=0, xadj=0, labels, col_palette=NULL, ...) {
    dots <- list(...)
    if(!is.null(dots$border)) border <- dots$border else border <- NA
    if(is.null(col_palette)) col_palette <- scales::hue_pal()(length(labels))
    text(x+xadj,y,title,adj=0,font=2)
    legend(x,y,legend=labels,
           xjust=0, yjust=1, xpd=TRUE,
           title.adj=0, adj=0,
           bty="n", border=border,
           fill=col_palette)
}

draw.bracket <- function(x0,y0,xl,xh) {
    llen <- 0.15
    lines(c(x0,x0), c(y0,y0+llen))
    lines(c(xl,xh), y0*c(1,1))
    lines(c(xl,xl), c(y0,y0-llen))
    lines(c(xh,xh), c(y0,y0-llen))
}

layout.matrix <- function(row.h, col.w) {
    layout.linear <- c()
    panel.number <- c(1:(length(row.h)*length(col.w)))
    panel.matrix <- matrix(panel.number, nrow=length(row.h), ncol=length(col.w))
    for (x in seq_along(col.w)) {
        column.linear <- c()
        for (y in seq_along(row.h)) {
            column.linear <- c(column.linear, rep(panel.matrix[y,x], row.h[y])) 
        }
        layout.linear <- c(layout.linear, rep(column.linear, col.w[x])) 
    }
    return(matrix(layout.linear, nrow=sum(row.h), ncol=sum(col.w)))
}

# To make log10 stacked barplot
logx <- function(m) {
    # Assuming that each row is a stack to be plotted in the barplot
    rnum <- nrow(m)

    # Find cumulative sum of rows
    mc <- m
    for (i in c(2:rnum)) mc[i,] <- mc[i,]+mc[i-1,]

    # Compute log10
    mcl <- log10(mc)
    # NOTE: log of 0 is -Inf; suppressing it internally
    mcl[mc == 0] <- 0

    # Find cumulative difference
    mcld <- mcl
    for (i in c(rnum:2)) mcld[i,] <- mcld[i,]-mcld[i-1,]

    # Return the result
    return(mcld)

}

# Select baseline date that is the lastest date before treatment
# Assuming that if the blood draw was on the same day as the treatment
# that it was collected before the treatment
get_baseline <- function(d, dtx) {

    del_d <- as.numeric(difftime(d, dtx, units = "days"))

    if(all(del_d>0)) {
        return(NA)
    } else {
        return(d[which.max(del_d<=0)])
    }

}

# Select the P2 / second visit timepoint
get_p2 <- function(d, v) {

    # If there is no second visit
    if(sum(!is.na(v)&(v==2)) == 0) {
        return(NA)
    # If there is a visit
    } else {
        return(d[!is.na(v)&(v==2)])
    }

}

bor_labels <- function(x) {

    case_when(x == 1 ~ "CR",
              x == 2 ~ "PR",
              x == 3 ~ "SD",
              x == 4 ~ "PD",
              x == 5 ~ "NE",
              .default = NA)

}

excel_date <- function(x) {

    as.Date(as.numeric(x), origin = "1899-12-30")

}
