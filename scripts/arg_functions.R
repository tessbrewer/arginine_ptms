library(ggplot2)
library(stringr)
library(stringi)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ape)
library(ggtree)
library(ggnewscale)
library(caper)
library(ggsignif)
library(ggrepel)
library(Biostrings)
library(ggmsa)
library(phytools)
library(ggtext)

#Arginine project functions =================================

#return top n largest colony sizes
subsample <- function(x, to_split, n){
  collect <- list()
  unique_to_split = unique(x[[to_split]])
  count = 0
  for (i in unique_to_split){
    count = count + 1
    subbie <-  x[x[[to_split]] == i, ]
    if (nrow(subbie) > n){
      subbie <- subbie[1:n, ]
    }
    collect[[count]] <- subbie
  }
  finished <- do.call(rbind, collect)
  return(finished)
}


#only assign a name to n many taxonomic groups, otherwise label them other
#for assigning colors to n many groups
assign_color <- function(mapping, level, n){
  countz <- table(mapping[[level]])
  mapping$color <- ifelse(mapping[[level]] %in% names(tail(sort(countz), (n - 1))), mapping[[level]], 'Other')
}


#get node based on vector of tips 
#DomBennett MoreTreeTools
getParent <- function (tree, node=NULL, tips=NULL, edges=NULL) {
  if (!is.null (node) & length (node) == 1) {
    if (!is.numeric (node)) {
      stop ('Node must be numeric')
    }
    if (node > getSize (tree) + tree$Nnode) {
      stop ('Node not in tree')
    }
    if ((node == getSize (tree) + 1) & is.rooted (tree)) {
      # if node is root, return it
      return (node)
    }
    return (tree$edge[tree$edge[ ,2] == node, 1])
  } else if (!is.null (tips)) {
    if (is.character (tips)) {
      # if tips are labels
      edges <- match (match (tips, tree$tip.label), tree$edge[,2])
    } else {
      # ... else they're numbers
      edges <- match (tips, tree$edge[,2])
    }
  } else if (!is.null (node)) {
    edges <- which (tree$edge[ ,2] %in% node)
  } else if (!is.null (edges)) {
    if (is.character (edges) & !is.null (tree$edge.label)) {
      # assume they are labels
      edges <- match (edges, tree$edge.label)
    }
  } else {
    stop ('Must provide either edges, tips or nodes argument')
  }
  end.nodes <- tree$edge[edges, 1]
  term.node <- length (tree$tip.label) + 1
  while (TRUE){
    if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
      break
    }
    end.nodes <- sort (end.nodes, TRUE)
    start.node <- end.nodes[1]
    edge <- match (start.node, tree$edge[,2])
    end.node <- tree$edge[edge,1]
    edges <- c(edges, edge)
    end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
  }
  return (end.nodes[1])
}


#Arginine project common files =================================

#key species
species <- data.frame(code = c('heau', 'mepr', 'nico', 'dera', 'gefeR', 'deac', 'shon', 'eco', 'ppu'), 
                      species = c('Herpetosiphon aurantiacus', 'Mesotoga prima', 'Nitrosomonas communis',
                                  'Deinococcus radiodurans', 'Geoalkalibacter ferrihydriticus', 'Denitrovibrio acetiphilus',
                                  'Shewanella oneidensis', 'Escherichia coli', 'Psuedomonas putida'))
species$abbrev <- paste0(sapply(strsplit(as.character(species$species), ''), '[', 1), '. ', 
                         sapply(strsplit(as.character(species$species), ' '), '[', 2))




