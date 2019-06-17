#!/usr/bin/env Rscript

#' Read a nexus file with color annotations into a treedata object with a color column
#'
#' @param trefile path to a nexus file
#' @export
#' @examples
#' tre <- system.file("examples", "H1_new.tre", package="colortree")
read_color_nexus <- function(trefile){
  # read in the tree file as a character vector
  tre_str <- scan(file=trefile, what="", sep="\n", quiet=TRUE)

  # Extract leaf node colors from the newick tree. These represent the color of
  # the branches leading to the leaf.
  # I currently handle only nexus files with a single tree, with the pattern:
  #     tree tree_1 = [&R] <tree>
  newick_tree <- tre_str[grep("tree tree_1", tre_str, perl=TRUE)]
  newick_tree <- sub("^[^(]*", "", x=newick_tree, perl=TRUE)
  leaf_branch_colors <- newick_tree %>%
    # find all named nodes/tips in the newick tree that have meta data
    stringr::str_extract_all("'[^']+'\\[[^\\]]+\\]") %>% '[['(1) %>%
    # extract the node/tip name and the color field
    stringr::str_replace("'([^']+)'\\[.*color=(#[0-9a-f]+).*\\]", "\\1\t\\2") %>%
    # split the resulting list into a matrix where columns contain name and color
    stringr::str_split_fixed("\t", n=2) %>%
    tibble::as_tibble() %>%
    magrittr::set_colnames(c("name", "branch_color"))

  # extract the taxlabels column, this is not a very reliable method, I should
  # extract the lines bound by the expressions '^\ttaxlabels' and '^;'.
  taxlabels <- tre_str[stringr::str_detect(tre_str, "^\t'[^']+'")]

  # Extract label colors from the taxlabels block.
  leaf_label_colors <- taxlabels %>%
    stringr::str_replace("^\t*'([^']+)'\\[.*color=(#[0-9a-f]+).*\\]$", "\\1\t\\2") %>%
    # split the resulting list into a matrix where columns contain name and color
    stringr::str_split_fixed("\t", n=2) %>%
    tibble::as_tibble() %>%
    magrittr::set_colnames(c("name", "label_color"))

  leaf_colors <- merge(leaf_branch_colors, leaf_label_colors, by="name", all=TRUE)

  tre_str <- gsub("'\\[[^\\]]+\\]", "'", tre_str, perl=TRUE)

  # quote metadata
  # tre_str <- gsub("(\\[[\\]]+\\])", "'\\1'", tre_str, perl=TRUE)
  tre_str <- gsub("(\\[[^\\]]*\\])", "'\\1'", tre_str, perl=TRUE)

  # Remove any repeated quotes that are produced.
  # Single quotes occur if the node is named, for example:
  #   'foo'[&label=0.5,!color=#ffffff] --> 'foo''[&label=0.5,!color=#ffffff]'
  # Removing double quotes gives the desired form:
  #   'foo[&label=0.5,!color=#ffffff]'
  tre_str <- gsub("''", "", tre_str)

  # Nexus format does not require that names be quoted.
  # Unquoted names will be processed as follows:
  #   foo[&label=0.5,!color=#ffffff] --> foo'[&label=0.5,!color=#ffffff]'
  # This is not correct, we need to move the quote to include the name.  The
  # expression below SHOULD work, though I am not sure what characters are
  # legal in an unquoted name.
  tre_str <- gsub("([a-zA-Z0-9_-]+)'\\[&", "'\\1[&", tre_str)

  # The lines with the trees start with a little label, for example
  #    tree tree_1 = [&R] ((('KX0884
  # The expressions above quote the '[&R]', which is not what we want. I can
  # fix this immediate problem by simply replacing it with the unquoted string:
  tre_str <- gsub("'[&R]'", "[&R]", tre_str, fixed=TRUE)
  # Though ultimately a more general solution will probably be needed.

  # NOTE: currently I only support color, adding other node metadata is
  # possible, but don't use '=', since that raises an error.
  tre_str <- gsub("\\)'\\[&[^\\]]*color=#([0-9a-f]+)[^\\]]*\\]'", ")'#\\1'", x=tre_str, perl=TRUE)

  new_trefile <- tempfile(pattern=basename(trefile), fileext=".tre")
  write(tre_str, file=new_trefile)
  cat("written to: ", new_trefile, "\n")

  branch_colormap <- leaf_colors$branch_color
  names(branch_colormap) <- leaf_colors$name
  label_colormap <- leaf_colors$label_color
  names(label_colormap) <- leaf_colors$name

  tre2 <- treeio::as.treedata(ape::read.nexus(new_trefile))
  tre2@phylo$tip.label <- gsub("'", "", tre2@phylo$tip.label)
  tre2@phylo$node.label <- gsub("'", "", tre2@phylo$node.label)

  N <- tre2@phylo$Nnode
  d <- tibble::tibble(node = 1:(length(tre2@phylo$tip.label)+N))
  d$name <- c(tre2@phylo$tip.label, rep("", N))
  d$branch_color <- c(branch_colormap[tre2@phylo$tip.label], tre2@phylo$node.label)
  d$branch_color <- ifelse(d$branch_color == "", "black", d$branch_color) 
  d$label_color <- c(label_colormap[tre2@phylo$tip.label], rep("black", N))

  tre2@data <- d

  file.remove(new_trefile)

  tre2@phylo <- tre2@phylo %>% reorder(.) %>% midpoint(.)
  
  tre2
}

#' Plot the output of \code{read_color_nexus}
#'
#' @param tre treedata object with data fields as specified in \code{read_color_nexus}
#' @return ggplot object
#' @export
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
#' @export
export_color_nexus <- function(g, path, width=5, height=5, units="in", device="pdf"){
  ggplot2::ggsave(filename=path, plot=g, width=width, height=height, units="in", device=device)
}

#' Convert a nexus file to a pdf
#'
#' @param trefile the path to a nexus file
#' @export
nexus2pdf <- function(trefile){
  pdffile <- file.path(dirname(trefile), sub("\\.[^.]*$", ".pdf", basename(trefile)))
  trefile %>% read_color_nexus(.) %>% plot_color_nexus(.) %>% export_color_nexus(., path=pdffile)
  #export_color_nexus(plot_color_nexus(read_color_nexus(trefile)), path=pdffile)
}

args = commandArgs(trailingOnly=TRUE)

nexus2pdf(args[1])
