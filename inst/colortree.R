#!/usr/bin/env Rscript

require(colortree)

#' Plot the output of \code{read_color_nexus}
#'
#' @param tre treedata object with data fields as specified in \code{read_color_nexus}
#' @return ggplot object
plot_color_nexus <- function(tre){
  # branch_color
  # label_color
  ggtree::ggtree(tre) + ggtree::geom_tree(color=tre@data$branch_color) 
}

#' Export a plot to a PDF
#'
#' @param g ggplot object
#' @param path filename for the output file
#' @param width width of the figure
#' @param height height of the figure
#' @param units units for width and height
#' @param device format of output (default="pdf")
#' @return nothing
export_color_nexus <- function(g, path, width=5, height=5, units="in", device="pdf"){
  ggplot2::ggsave(filename=path, plot=g, width=width, height=height, units="in", device=device)
}

#' Convert a nexus file to a pdf
#'
#' @param trefile the path to a nexus file
nexus2pdf <- function(trefile){
  pdffile <- file.path(sub("\\.[^.]*$", ".pdf", basename(trefile)))
  trefile %>% colortree::read_color_nexus(.) %>%
    plot_color_nexus(.) %>%
    export_color_nexus(., path=pdffile)
}

args = commandArgs(trailingOnly=TRUE)

nexus2pdf(args[1])
