#' @importFrom magrittr "%>%"
utils::globalVariables(c("%>%", "%$%", "."))
NULL

#' colortree: Read colored nexus files and convert them to PDFs
#'
#' @docType package
#' @name colortree
NULL

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

  # Remove the metadata from the taxlabel entries (it is stored in leaf_colors)
  tre_str <- gsub("'\\[[^\\]]+\\]$", "'", tre_str, perl=TRUE)

  # Remove the leaf metadata from the newick entry
  tre_str <- gsub("([^)]'[^']+')\\[[^\\]]+\\]", "\\1", tre_str, perl=TRUE)
  #                 ^--^
  #                 | internal labels always start with ')', leaves never do

  # Quote the node metadata 
  tre_str <- gsub("\\[&([^\\]]+)\\]([^ ])", "'\\1'\\2", tre_str, perl=TRUE)
  # Replace metadata equal signs with ':'
  # This is really brittle
  tre_str <- gsub("([^ ])=([^ ])", "\\1oxxo\\2", tre_str, perl=TRUE)
  tre_str <- gsub("(oxxo[^,']+),", "\\1oolo", tre_str, perl=TRUE)

  # Remove any repeated quotes that are produced.
  # Single quotes occur if the node is named, for example:
  #   'foo'[&label=0.5,!color=#ffffff] --> 'foo''[&label=0.5,!color=#ffffff]'
  # Removing double quotes gives the desired form:
  #   'foo[&label=0.5,!color=#ffffff]'
  tre_str <- gsub("''", "", tre_str)

  # Create a temporary file to store the modified tree file
  # This is necessary since ape will read only from a file
  # The temporary file will be deleted after being used 
  new_trefile <- tempfile(pattern=basename(trefile), fileext=".tre")
  write(tre_str, file=new_trefile)
  cat("written to: ", new_trefile, "\n")

  # Read the temporary tree file and get the midpoint
  tre2 <- treeio::as.treedata(ape::read.nexus(new_trefile))
  tre2@phylo <- tre2@phylo %>%
    ape::reorder.phylo(.) %>%
    phangorn::midpoint(., node.labels="label")
  # Kill the pesky temporary file
  file.remove(new_trefile)

  # Get the number of nodes and leafs in the tree
  N <- tre2@phylo$Nnode
  L <- length(tre2@phylo$tip.label)

  # Remove quotes from names
  tre2@phylo$tip.label <- gsub("'", "", tre2@phylo$tip.label)
  tre2@phylo$node.label <- gsub("'", "", tre2@phylo$node.label)

  # If there are no node labels, initialize empty ones
  if(is.null(tre2@phylo$node.label)){
    tre2@phylo$node.label <- rep("", N)
  }

  # Parse the horribly mangled node names
  node_metadata <- stringr::str_replace(tre2@phylo$node.label, "!", "") %>%
    # Replace the previously generated delimiter with the forbidden '=' sign
    stringr::str_replace_all("oxxo", "=") %>%
    # Split metadata lists into key/value pairs
    stringr::str_split("oolo", simplify=TRUE) %>%
    # Turn this into a data frame. The columns do not necessarily correspond to
    # the a particular label.
    as.data.frame(stringsAsFactors=FALSE)

  # Createa a long table, for example:
  #   node   color label
  # 1 1151 #80cdc1 0.889
  # 2 1152 #80cdc1 0.819
  # 3 1153 #80cdc1 0.902
  # 4 1154 #80cdc1  0.93
  node_metadata$node <- 1:N+L
  node_metadata <- node_metadata %>%
    reshape2::melt(id.vars="node") %>%
    tidyr::separate(col="value", into=c("k", "v"), sep="=") %>%
    dplyr::select(-variable) %>%
    dplyr::filter(!is.na(v)) %>%
    reshape2::dcast(node ~ k)

  # Create a node branch color column. If there are no node colors, default to black.
  if("color" %in% names(node_metadata)){
    node_metadata <- dplyr::rename(node_metadata, branch_color = color)
  } else {
    node_metadata$branch_color <- "#000000"
  }
  node_branch_color_map <- node_metadata$branch_color
  names(node_branch_color_map) <- node_metadata$node
  node_metadata$branch_color <- NULL

  # Currently I do not preserve node labels. Or do I?
  tre2@phylo$node.label <- NULL

  # Create leaf branch colors. I assume these exist for now.
  branch_colormap <- leaf_colors$branch_color
  names(branch_colormap) <- leaf_colors$name
  label_colormap <- leaf_colors$label_color
  names(label_colormap) <- leaf_colors$name

  d <- tibble::tibble(node = 1:(L+N))
  d$name <- c(tre2@phylo$tip.label, rep("", N))
  d$branch_color <- c(branch_colormap[tre2@phylo$tip.label], node_branch_color_map[as.character((L+1):(L+N))])
  d$branch_color <- ifelse(d$branch_color == "" | is.na(d$branch_color), "#000000", d$branch_color)
  d$label_color <- c(label_colormap[tre2@phylo$tip.label], rep("#000000", N))

  # Merge in the non-color node metadata fields
  d <- dplyr::left_join(d, node_metadata, by="node")

  # Write the metadata field into the treedata object. These fields can now be
  # accessed as ggplot2/ggtree 'aes' fields.
  tre2@data <- d

  tre2
}

#' Plot the output of \code{read_color_nexus}
#'
#' @param tre treedata object with data fields as specified in \code{read_color_nexus}
#' @return ggplot object
#' @export
default_plot_color_nexus <- function(tre){
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
default_export_color_nexus <- function(g, path, width=5, height=5, units="in", device="pdf"){
  ggplot2::ggsave(filename=path, plot=g, width=width, height=height, units="in", device=device)
}

#' Convert a nexus file to a pdf
#'
#' @param trefile the path to a nexus file
#' @export
default_nexus2pdf <- function(trefile){
  pdffile <- file.path(dirname(trefile), sub("\\.[^.]*$", ".pdf", basename(trefile)))
  trefile %>%
    read_color_nexus %>%
    default_plot_color_nexus %>%
    default_export_color_nexus(path=pdffile)
}
