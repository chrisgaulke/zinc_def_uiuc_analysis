options("stringsAsFactors" = F)

# FUNCTIONS ---------------------------------------------------------------

###
#        Function phylotype_analysis             #
###

phylotype_analysis <- function(obj, tax) {
  #obj: microbiome object with at least 1 slot (data)
  #tax: a tax object (named list taxa as names values in the list are seq ids)
  obj.out <- NULL
  for (h in 1:length(tax)) {
    df <- NULL
    #print(h)#debugging
    for (i in 1:length(tax[[h]])) {
      print(i)#debugging
      v1       <- obj$data[, unlist(tax[[h]][[i]])]
      v2       <- names(tax[[h]])[i]
      if (is.null(dim(v1))) {
        df[[v2]] <- v1
      } else{
        df[[v2]] <- rowSums(v1)
      }
    }
    obj.out[[names(tax)[h]]] <- as.data.frame(df)
  }
  return(obj.out)
}

make_taxa_df <- function(tax){

  kingdom.df <- replicate(length(unique(tax[, 2])), c())
  names(kingdom.df) <- unique(tax[, 2])
  phylum.df  <- replicate(length(unique(tax[, 3])), c())
  names(phylum.df) <- unique(tax[, 3])
  class.df   <- replicate(length(unique(tax[, 4])), c())
  names(class.df) <- unique(tax[, 4])
  order.df   <- replicate(length(unique(tax[, 5])), c())
  names(order.df) <- unique(tax[, 5])
  family.df  <- replicate(length(unique(tax[, 6])), c())
  names(family.df) <- unique(tax[, 6])
  genus.df   <- replicate(length(unique(tax[, 7])), c())
  names(genus.df) <- unique(tax[, 7])

  for (i in 1:nrow(tax)) {
    kingdom.df[[tax[i, 2]]] <-
      c(kingdom.df[[tax[i, 2]]], tax[i, 1])
    phylum.df[[tax[i, 3]]]  <-
      c(phylum.df[[tax[i, 3]]], tax[i, 1])
    class.df[[tax[i, 4]]]   <-
      c(class.df[[tax[i, 4]]], tax[i, 1])
    order.df[[tax[i, 5]]]   <-
      c(order.df[[tax[i, 5]]], tax[i, 1])
    family.df[[tax[i, 6]]]  <-
      c(family.df[[tax[i, 6]]], tax[i, 1])
    genus.df[[tax[i, 7]]]   <-
      c(genus.df[[tax[i, 7]]], tax[i, 1])
}

  tax.obj <- NULL
  tax.obj$kingdom <- kingdom.df
  tax.obj$phylum  <- phylum.df
  tax.obj$class   <- class.df
  tax.obj$order   <- order.df
  tax.obj$family  <- family.df
  tax.obj$genus   <- genus.df

return(tax.obj)
}

# DATA: AGGREGATE PHYLOTYPES -----------------------------------------------------

#Added
rownames(taxa.print) <- paste0("ASV", 1:738)

#we need to add an ASV identifier for the tax fucntion to work properly
tax <- taxa.print
tax <- data.frame(asv = rownames(tax), tax)

#make a data frame of taxonomy
zinc.tax <- make_taxa_df(tax = tax)

# aggregate phylotype counts for this we will need to provide an pseudoobject with a
# slot called data
zinc.obj <- NULL
zinc.obj$data <- seqtab_nochim.relabd

#classifying taxonomy as also called phyotyping
zinc_phylotype <- phylotype_analysis(zinc.obj, tax = zinc.tax)

# This is a little different than we have seen before because this is an list of
# data frames which means we can use the '$' operator to select different
# data frames of interest . For example, zinc_phylotype$phylum will allow us
# to access a data frame with phylum level abundance counts

View(zinc_phylotype$genus)
