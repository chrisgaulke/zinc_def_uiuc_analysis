# Title and author information --------------------------------------------
#!/usr/bin/R

#################################################
#                                               #
#          2021_08_26_zinc_def_uiuc_dada2.R     #
#                                               #
#################################################

#Title: Effects of moderate and severe zinc deficiency on gut microbiome
#
#Copyright (C) 2021-2022  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

#Purpose: The purpose of this script is to generate ASV profiles for
# zinc deficiency data


# SET ENVIRONMENT ---------------------------------------------------------

## Install Bioconductor

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")

## Install dada2
#BiocManager::install("dada2")
library(dada2)
library(ggplot2)
library(vegan)

options(stringsAsFactors = F)

#clever girl

#decompress the and untar the archive
#tar -xvf Zinc_16S_821.2021812.tar.bz2

#used this to move all empty files to a dump dir
#find . -maxdepth 1 -type f -size -100 -exec mv {} dumps/ \;


# Functions ---------------------------------------------------------------

###
#              Function filter_df                #
###

filter <- function(df) {
  df <- df[which(rowSums(df) > 0), which(colSums(df) > 0)]
  return(df)
}

###
#             Function normalize                 #
###

#normalize counts either by relative abundance of rarefying
#depends on vegan

normalize <- function(df, method="rel", depth=depth){
  #default method = relative abundance
  if(method == "rare"){
    if( is.null(depth)){
      depth <- min(rowSums(df))
    }
    ndf <- rrarefy(df, depth)
    ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }else{
    ndf <- sweep(df,
                 1,
                 rowSums(df),
                 `/`)
    #  ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }
  return(ndf)
}


# IMPORT DATA -------------------------------------------------------------

path <- "~/unsynced_projects/raw_data/2021_08_12_gaulke_zinc_mouse_uiuc_16S_raw/"
filt.path <- "/Users/cgaulke/Documents/research/zinc_def_uiuc/zinc_def_uiuc_analysis/data/filtered_data/" #filtered file directory make sure to update
#git ignore to ignore this
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


# ANALYSIS: QC ------------------------------------------------------------

#make sure the lengths are the same
length(fnFs) == length(fnRs)

#get sample names, in our case we have to leave some extra on the end
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

#preview
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

#aggregate all data together now
fnFs_qual.plot <- plotQualityProfile(fnFs,aggregate = T)
fnRs_qual.plot <- plotQualityProfile(fnRs,aggregate = T)

#set up for filtering
filtFs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_R_filt.fastq.gz"))

#make these named
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim
#by looking at the data we see a rapid drop in quality around 200bp. Since the
#average drops below ~25 around 230 we will truncate at 220 for the reverse. The
#forward looks better (this is usual) so we will truncate around 260

#note the original tutorial uses generic variable names

#stopped here

filter.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#take a look
View(filter.out)

colMeans(filter.out)
mean(1-(filter.out[,2]/filter.out[,1]))

#lets look at these numbers in a little more detail
fivenum(1-(filter.out[,2]/filter.out[,1]))
hist(1-(filter.out[,2]/filter.out[,1]))

# since some libraries look like there is a higher level of filtration
# than others lets take a closer look at this

sort(1-(filter.out[,2]/filter.out[,1]))

#let's keep this in mind moving forward

# ANALYSIS: ERROR ---------------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# ANALYSIS: MERGE AND FILTER -----------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#how much data was wrapped up in chimeras
sum(seqtab.nochim)/sum(seqtab)

write.table(seqtab.nochim,
            file = "data/dada2/seqtab_nochim.txt",
            quote = FALSE,
            sep = "\t",
            )

# ANALYSIS: TRACK READS ---------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(filter.out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# ANALYSIS: ADD TAX -------------------------------------------------------

#***Here is where you will need to go and download the silva databases.
#***Be sure to get the right ones (the names are the same as the ones below)
#***These files can be downloaded here:https://zenodo.org/record/4587955#.YSlzKC1h1hA


taxa <- assignTaxonomy(seqtab.nochim,
          "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_nr99_v138.1_train_set.fa",
          multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_species_assignment_v138.1.fa")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# IMPORT: METADATA --------------------------------------------------------

mouse_metadata.df <- read.table("/Users/cgaulke/Documents/research/zinc_def_uiuc/data/combined_mouse_metadata.txt",
                                sep="\t",
                                header = T)

rownames(mouse_metadata.df) <- mouse_metadata.df$id


mouse_ids.names <- sapply(rownames(seqtab.nochim),FUN = function(x){strsplit(x,"_")[[1]][1]} )
names(mouse_ids.names) <- NULL

master_metadata.df <- data.frame(matrix(nrow = 155,ncol = 8, "NA"))
colnames(master_metadata.df) <- c("id", "weight", "sex", "age", "group", "gf", "day", "tissue")

rownames(master_metadata.df) <- rownames(seqtab.nochim)
master_metadata.df$id <- mouse_ids.names

#must use the <<- instead of the <- operator
#note this code is not very readable

lapply(rownames(master_metadata.df),
       FUN= function(x){
         master_metadata.df[x,c("weight", "sex", "age", "group")] <<-
           mouse_metadata.df[master_metadata.df[x,"id"],
                             c("weight", "sex", "age", "group")] })

# now one by one we can add additional data with grep

master_metadata.df$gf <-
  as.numeric(grepl(pattern = "gZ",x = master_metadata.df$id))


master_metadata.df$tissue[grep(pattern = "cecum",
                               x = rownames(master_metadata.df),
                               ignore.case = T)] <- "cecum"

master_metadata.df$tissue[grep(pattern = "Duodenum",
                               x = rownames(master_metadata.df),
                               ignore.case = T)] <- "Duodenum"

master_metadata.df$tissue[grep(pattern = "Jejunum",
                               x = rownames(master_metadata.df),
                               ignore.case = T)] <- "Jejunum"

master_metadata.df$tissue[grep(pattern = "feces",
                               x = rownames(master_metadata.df),
                               ignore.case = T)] <- "feces"

master_metadata.df$day[grep(pattern = "Day_0_",
                               x = rownames(master_metadata.df),
                               ignore.case = T)] <- "0"

master_metadata.df$day[grep(pattern = "Day_21_",
                            x = rownames(master_metadata.df),
                            ignore.case = T)] <- "21"

master_metadata.df$day[grep(pattern = "Day_42_",
                            x = rownames(master_metadata.df),
                            ignore.case = T)] <- "42"

master_metadata.df$day[grep(pattern = "Day",
                            x = rownames(master_metadata.df),
                            ignore.case = T, invert = T)] <- "43"



# ANALYSIS: NORMALIZE  ----------------------------------------------------

all(rownames(seqtab.nochim) == rownames(master_metadata.df))

raresamp <- min(rowSums(seqtab.nochim))

#now we can plot the curves
rarecurve(seqtab.nochim, step = 1000,
          sample = 0,
          label = F,
          col = "black")

abline(v = 2000, col = "red")
abline(v = 5000, col = "yellow")
abline(v = 7500, col = "yellow")
abline(v = 10000, col = "black")


seqtab_nochim.filter <- seqtab.nochim[which(rowSums(seqtab.nochim) > 4999),]

seqtab_nochim.filter <- normalize(seqtab_nochim.filter,method = "rare", depth = 5000)
seqtab_nochim.filter <- filter(seqtab_nochim.filter)

filter_metadata.df <- master_metadata.df[which(rownames(master_metadata.df) %in% rownames(seqtab_nochim.filter)),]

#make a lookup so we can change header
seqtab_nochim.lookup <- data.frame(seqs = colnames(seqtab.nochim),
                                   asv = paste0("ASV", seq(from = 1, to =ncol(seqtab.nochim)))
)

colnames(seqtab.nochim) <- paste0("ASV", seq(from = 1, to =ncol(seqtab.nochim)))


seqtab_nochim.relabd <- sweep(seqtab.nochim,
                              MARGIN = 1,
                              rowSums(seqtab.nochim),
                              FUN = '/')

# ANALYSIS: PCA all -----------------------------------------------------------

#Not run
seqtab_nochim.pca <- prcomp(seqtab_nochim.filter,scale. = F, center = F)

all(rownames(seqtab_nochim.pca$x) == rownames(filter_metadata.df))


seqtab_nochim_pca.df <- as.data.frame(seqtab_nochim.pca$x[,1:5])

seqtab_nochim_pca.df$sex <- filter_metadata.df$sex
seqtab_nochim_pca.df$group <- filter_metadata.df$group
seqtab_nochim_pca.df$gf <- filter_metadata.df$gf
seqtab_nochim_pca.df$day <- filter_metadata.df$day
seqtab_nochim_pca.df$tissue <- filter_metadata.df$tissue

seqtab_nochim_pca.plot <- ggplot(data = seqtab_nochim_pca.df,
                                 aes(x = PC2,
                                     y = PC3,
                                     color = tissue,
                                     shape = group))

seqtab_nochim_pca.plot +
  geom_point()


# ANALYSIS: Adonis all --------------------------------------------------------

set.seed(731)
seqtab_nochim.adonis <- adonis(seqtab_nochim.filter~tissue+group+gf+id,
                               data = filter_metadata.df, permutations = 5000)

seqtab_nochim.adonis

# ANALYSIS: PCA Feces -----------------------------------------------------------

seqtab_nochim_feces.pca <- prcomp(subset(seqtab_nochim.filter,
                             subset = rownames(seqtab_nochim.filter) %in%
                               rownames(filter_metadata.df[which(filter_metadata.df$tissue == "feces"),]
                                        )
                             ),
                             scale. = F,
                             center = F)

all(rownames(seqtab_nochim_feces.pca$x) == rownames(filter_metadata.df[which(filter_metadata.df$tissue == "feces"),]))


seqtab_nochim_feces_pca.df <- as.data.frame(seqtab_nochim_feces.pca$x[,1:5])

seqtab_nochim_feces_pca.df$sex    <- filter_metadata.df[which(filter_metadata.df$tissue == "feces"),"sex"]
seqtab_nochim_feces_pca.df$group  <- filter_metadata.df[which(filter_metadata.df$tissue == "feces"),"group"]
seqtab_nochim_feces_pca.df$gf     <- filter_metadata.df[which(filter_metadata.df$tissue == "feces"),"gf"]
seqtab_nochim_feces_pca.df$day    <- filter_metadata.df[which(filter_metadata.df$tissue == "feces"),"day"]
seqtab_nochim_feces_pca.df$tissue <- filter_metadata.df[which(filter_metadata.df$tissue == "feces"),"tissue"]

seqtab_nochim_feces_pca.plot <- ggplot(data = seqtab_nochim_feces_pca.df,
                                 aes(x = PC1,
                                     y = PC2,
                                     color = group,
                                     shape = day))

seqtab_nochim_feces_pca.plot +
  geom_point()

# ANALYSIS: adonis Feces -----------------------------------------------------------

set.seed(731)
seqtab_nochim_feces.adonis <- adonis(subset(seqtab_nochim.filter,
                                            subset = rownames(seqtab_nochim.filter) %in%
                                              rownames(filter_metadata.df[which(filter_metadata.df$tissue == "feces"),]
                                              ))~group+day+sex+id,
                               data = filter_metadata.df[which(filter_metadata.df$tissue == "feces"),], permutations = 5000)

seqtab_nochim_feces.adonis

# ANALYSIS: BRAY ----------------------------------------------------------
seqtab_nochim.bray <- vegdist(seqtab_nochim.filter, method = "bray")
seqtab_nochim.bray <- as.matrix(seqtab_nochim.bray)

fecal.rnames <- rownames(filter_metadata.df[which(filter_metadata.df$tissue == "feces"),])
fecal.ids <- unique(filter_metadata.df[fecal.rnames, "id"])

seqtab_nochim_bray.feces <- seqtab_nochim.bray[which(rownames(seqtab_nochim.bray) %in% fecal.rnames),which(colnames(seqtab_nochim.bray) %in% fecal.rnames) ]



#get intra individual diversity across time points
x<- sapply(fecal.ids,
       FUN = function(x){
         y <- rownames(filter_metadata.df[which(filter_metadata.df$tissue == "feces" & filter_metadata.df$id == x ),]);
         seqtab_nochim_bray.feces[y,y]
        }
       )

#get rid of dups
y <- sapply(x, FUN = function(x){ unique(as.vector(x))})

y <- lapply(y, function(x){x[x>0]})


#boxplot
boxplot(unlist(y[1:8]), unlist(y[9:16]), unlist(y[17:24]),
        names = c("low", "moderate", "control"))

#is this sig
kruskal.test(c(unlist(y[1:8]), unlist(y[9:16]), unlist(y[17:24])),
             g = rep(c("low", "moderate", "high"), times= c(23,22,22)) )

#post hoc
pairwise.wilcox.test(c(unlist(y[1:8]), unlist(y[9:16]), unlist(y[17:24])),
                     g = rep(c("low", "moderate", "high"), times= c(23,22,22)) )

#both the kruskal and posthoc is significant. This means that temporal variability
#differs between groups. Increased temporal variability could indicate increased sensitivity
#to stimuli?

# ANALYSIS: PCA Jejunum -----------------------------------------------------------

seqtab_nochim_Jejunum.pca <- prcomp(subset(seqtab_nochim.relabd,
                                           subset = rownames(seqtab_nochim.relabd) %in%
                                             rownames(master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),]
                                             )
),
scale. = F,
center = F)

all(rownames(seqtab_nochim_Jejunum.pca$x) == rownames(master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),]))


seqtab_nochim_Jejunum_pca.df <- as.data.frame(seqtab_nochim_Jejunum.pca$x[,1:5])

seqtab_nochim_Jejunum_pca.df$sex    <- master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),"sex"]
seqtab_nochim_Jejunum_pca.df$group  <- master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),"group"]
seqtab_nochim_Jejunum_pca.df$gf     <- master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),"gf"]
seqtab_nochim_Jejunum_pca.df$day    <- master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),"day"]
seqtab_nochim_Jejunum_pca.df$tissue <- master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),"tissue"]

seqtab_nochim_Jejunum_pca.plot <- ggplot(data = seqtab_nochim_Jejunum_pca.df,
                                         aes(x = PC1,
                                             y = PC2,
                                             color = group,
                                             shape = day))

seqtab_nochim_Jejunum_pca.plot +
  geom_point()

# ANALYSIS: adonis Jejunum -----------------------------------------------------------

set.seed(731)
seqtab_nochim_Jejunum.adonis <- adonis(subset(seqtab_nochim.relabd,
                                              subset = rownames(seqtab_nochim.relabd) %in%
                                                rownames(master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),]
                                                ))~group* sex,
                                       data = master_metadata.df[which(master_metadata.df$tissue == "Jejunum"),], permutations = 5000)

seqtab_nochim_Jejunum.adonis

# ANALYSIS: PCA Duodenum -----------------------------------------------------------

seqtab_nochim_Duodenum.pca <- prcomp(subset(seqtab_nochim.relabd,
                                            subset = rownames(seqtab_nochim.relabd) %in%
                                              rownames(master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),]
                                              )
),
scale. = F,
center = F)

all(rownames(seqtab_nochim_Duodenum.pca$x) == rownames(master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),]))


seqtab_nochim_Duodenum_pca.df <- as.data.frame(seqtab_nochim_Duodenum.pca$x[,1:5])

seqtab_nochim_Duodenum_pca.df$sex    <- master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),"sex"]
seqtab_nochim_Duodenum_pca.df$group  <- master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),"group"]
seqtab_nochim_Duodenum_pca.df$gf     <- master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),"gf"]
seqtab_nochim_Duodenum_pca.df$day    <- master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),"day"]
seqtab_nochim_Duodenum_pca.df$tissue <- master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),"tissue"]

seqtab_nochim_Duodenum_pca.plot <- ggplot(data = seqtab_nochim_Duodenum_pca.df,
                                          aes(x = PC1,
                                              y = PC2,
                                              color = group,
                                              shape = day))

seqtab_nochim_Duodenum_pca.plot +
  geom_point()

# ANALYSIS: adonis Duodenum -----------------------------------------------------------

set.seed(731)
seqtab_nochim_Duodenum.adonis <- adonis(subset(seqtab_nochim.relabd,
                                               subset = rownames(seqtab_nochim.relabd) %in%
                                                 rownames(master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),]
                                                 ))~group * sex,
                                        data = master_metadata.df[which(master_metadata.df$tissue == "Duodenum"),], permutations = 5000)

seqtab_nochim_Duodenum.adonis



# ANALYSIS: PCA cecum -----------------------------------------------------------

seqtab_nochim_cecum.pca <- prcomp(subset(seqtab_nochim.relabd,
                                         subset = rownames(seqtab_nochim.relabd) %in%
                                           rownames(master_metadata.df[which(master_metadata.df$tissue == "cecum"),]
                                           )
),
scale. = F,
center = F)

all(rownames(seqtab_nochim_cecum.pca$x) == rownames(master_metadata.df[which(master_metadata.df$tissue == "cecum"),]))


seqtab_nochim_cecum_pca.df <- as.data.frame(seqtab_nochim_cecum.pca$x[,1:5])

seqtab_nochim_cecum_pca.df$sex    <- master_metadata.df[which(master_metadata.df$tissue == "cecum"),"sex"]
seqtab_nochim_cecum_pca.df$group  <- master_metadata.df[which(master_metadata.df$tissue == "cecum"),"group"]
seqtab_nochim_cecum_pca.df$gf     <- master_metadata.df[which(master_metadata.df$tissue == "cecum"),"gf"]
seqtab_nochim_cecum_pca.df$day    <- master_metadata.df[which(master_metadata.df$tissue == "cecum"),"day"]
seqtab_nochim_cecum_pca.df$tissue <- master_metadata.df[which(master_metadata.df$tissue == "cecum"),"tissue"]

seqtab_nochim_cecum_pca.plot <- ggplot(data = seqtab_nochim_cecum_pca.df,
                                       aes(x = PC1,
                                           y = PC2,
                                           color = group,
                                           shape = factor(gf)))

seqtab_nochim_cecum_pca.plot +
  geom_point()

# ANALYSIS: adonis cecum -----------------------------------------------------------

set.seed(731)
seqtab_nochim_cecum.adonis <- adonis(subset(seqtab_nochim.relabd,
                                            subset = rownames(seqtab_nochim.relabd) %in%
                                              rownames(master_metadata.df[which(master_metadata.df$tissue == "cecum"),]
                                              ))~group*gf ,
                                     data = master_metadata.df[which(master_metadata.df$tissue == "cecum"),], permutations = 5000)

seqtab_nochim_cecum.adonis

