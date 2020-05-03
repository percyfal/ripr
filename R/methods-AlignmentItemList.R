##' AlignmentItemList
##'
##' @export
##' @rdname AlignmentItemList-class
##'
##' @importFrom methods new
##'
setMethod("AlignmentItemList", "list",
          function(obj) new("AlignmentItemList", listData = obj))


##' as.data.frame
##'
##' @description Convert AlignmentItemList to data.frame.
##'
##' @param x AlignmentItemList object
##' @param ... additional parameters to lapply
##' @param .id name of id column added by dplyr::bind_rows
##'
##' @return data.frame
##'
##' @export
##'
##' @importFrom dplyr bind_rows
##'
setMethod("as.data.frame", signature = "AlignmentItemList",
          function(x, ..., .id="id") {
    data <- dplyr::bind_rows(lapply(x, as.data.frame, ...), .id=.id)
    ## Ensure factor order is preserved
    data[[.id]] <- factor(data[[.id]], levels=unique(data[[.id]]))
    data
})



##' autoplot.AlignmentItemList
##'
##' @param object AlignmentItemList
##' @param aes aes mapping
##' @param vars variable mapping to facet plots
##' @param ... additional parameters to ggplot function
##' @param which which plot to make. 'grid' option makes a scatter
##'     plot with marginal densities
##' @param xlim set x limits
##' @param ylim set y limits
##'
##'
##' @importFrom ggplot2 autoplot ggplot geom_point geom_boxplot
##'     geom_violin geom_density scale_color_viridis_d
##'     scale_fill_viridis_d quo_name facet_wrap element_blank
##'     coord_flip theme scale_y_reverse
##' @importFrom cowplot plot_grid
##'
##' @export
##'
autoplot.AlignmentItemList <- function(object, aes, vars, ..., which="point", xlim = NULL, ylim = NULL) {
    data <- as.data.frame(object)
    clabs <- NULL
    flabs <- NULL
    if ("colour" %in% names(aes))
        clabs <- unique(as.character(data[, quo_name(aes$colour)]))
    if ("fill" %in% names(aes))
        flabs <- unique(as.character(data[, quo_name(aes$fill)]))
    if (which != "grid") {
        p <- ggplot(data, {{ aes }})
        which <- match.arg(which, c("density", "point", "boxplot", "violin", "grid"))
        if (which == "point")
            p <- p + geom_point(...)
        else if (which == "boxplot")
            p <- p + geom_boxplot(...)
        else if (which == "violin")
            p <- p + geom_violin(...)
        else if (which == "density")
            p <- p + geom_density(...)
        if (!is.null(xlim))
            p <- p + xlim(xlim)
        if (!is.null(ylim))
            p <- p + ylim(ylim)
        p <- p + scale_color_viridis_d(labels=clabs) +
            scale_fill_viridis_d(labels=flabs)
    } else if (which == "grid") {
        p.scatter <- ggplot(data, {{ aes }}) + geom_point(...) +
            theme(axis.title.y = element_blank(),
                  panel.border = element_blank(),
                  axis.ticks.x = element_blank()) +
            scale_color_viridis_d(labels=clabs) +
            scale_fill_viridis_d(labels=flabs) + ylim(ylim)

        aes.y <- aes
        aes.y$x <- aes.y$y
        aes.y$y <- NULL
        p.y <- ggplot(data, {{ aes.y }} ) + geom_density(...) +
            coord_flip(xlim=xlim) + scale_y_reverse() +
            scale_color_viridis_d(labels=clabs) +
            scale_fill_viridis_d(labels=flabs) +
            theme(legend.position = "none", panel.border = element_blank(),
                  axis.text.x = element_blank(), axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())

        p <- cowplot::plot_grid(p.y + theme(panel.border = element_blank()),
                                p.scatter, ncol=2, align="hv", rel_widths=c(1, 4))
    }
    if (!missing(vars))
        p <- p + facet_wrap( {{ vars }} )
    p
}


##' plot
##'
##' @description plot an AlignmentItemList
##'
##' @param x object to plot
##' @param ... additional parameters for autoplot
##'
##' @export
##' @importFrom graphics plot
plot.AlignmentItemList <- function(x, ...) {
    print(autoplot(x, ...))
}
