devtools::load_all(".")
rm(list = ls())

library(dplyr)

# load in the raw data

otu.csv <- t(read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/Study_850_output/OTUs149.txt"))

bacteria.id <- otu.csv[1,]
otu <- otu.csv[-1,]

bact.nameid.csv <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/bacteria_id.csv", stringsAsFactors = FALSE) %>% filter(bact.id %in% bacteria.id)

bact.names.full <- sapply(bacteria.id, function(id){ # this should order the names with respect to the otu table
    name.id <- which(bact.nameid.csv$bact.id == id)
    bact.nameid.csv$tax[name.id]
})

genus <- strsplit(bact.names.full, ";") %>%
    sapply(function(s) substr(s[[6]], start=4, stop=nchar(s[[6]])))

colnames(otu) <- genus

clean.names <- sapply(rownames(otu), function(s) substr(s, start = 1, stop = nchar(s) - 7))
names(clean.names) <- NULL

otu.df <- as.data.frame(otu)
otu.df$id <- clean.names


ec.csv <- t(read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/kegg_ec_cleaned.csv"))

ec.df <- as.data.frame(ec.csv)
ec.df$id <- sapply(rownames(ec.csv), function(s) substr(s, start = 1, stop = nchar(s) - 7))

age.csv <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/age_id.csv", header=FALSE, stringsAsFactors = FALSE) %>% filter(V1 %in% clean.names) %>%
    transmute(id = V1, age = V2)


unifrac.csv <-read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/Study_850_output/DUniFrac_149.txt")[,-1] %>% apply(2, function(x) as.numeric(x))

unifrac.df <- as.data.frame(unifrac.csv)
colnames(unifrac.csv) <- sapply(colnames(unifrac.csv), function(s) substr(s, 1, nchar(s) - 7))
rownames(unifrac.csv) <- colnames(unifrac.csv)
unifrac.df$id <- colnames(unifrac.csv)


# create a data.frame to align the IDs correctly

temp.df <- merge(otu.df, age.csv, by="id") %>%
    merge(unifrac.df, by="id") %>%
    merge(ec.df, by="id")


# now extract the rearranged data

age <- as.vector(temp.df$age)
names(age) <- temp.df$id
Xraw <- select(temp.df, c(-id, -age)) %>% as.matrix %>% (function(x) x[,1:149])
rownames(Xraw) <- temp.df$id

ec.raw <- as.matrix(temp.df[,-(1:251)])
rownames(ec.raw) <- temp.df$id

ec.cols.to.ignore <- which(colMeans(ec.raw == 0) > 0.9) # check which pathways are not present in at lteast 90% of samples

ec.filtered <- ec.raw[,-ec.cols.to.ignore] # filter these pathways out

# we need to rearrange the unifrac distance matrix based on the `id` column in temp.df

unifrac <- sapply(1:100, function(i){

    cur.col.name <- temp.df$id[i]
    cur.col.id <- which(colnames(unifrac.csv) == cur.col.name)

    sapply(1:100, function(j){
        cur.row.name  <- temp.df$id[j]
        cur.row.id <- which(rownames(unifrac.csv) == cur.row.name)
        cur.dat <- unifrac.csv[cur.row.id ,cur.col.id]
        cur.dat
    })

})


# load in the Q matrix (we are assuming that this is ordered correctly)

Q <- as.matrix(read.table("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/Q.txt"))


geography.csv <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/geography.csv")

patristic.raw <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/Study_850_output/patristic149.txt", stringsAsFactors = FALSE)

patristic <- as.matrix(patristic.raw[,-1])
rownames(patristic) <- patristic.raw$X

brownian <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/Study_850_output/BrownianMotion149.txt")[,-1] %>% as.matrix


Q_ <- generateSimilarityKernel(patristic, squareValues = FALSE)

geography <- sapply(rownames(Xraw), function(s) filter(geography.csv, id == s)$geo)

E <- model.matrix(~ geography) %>% (function(mat) mat[,-1])

# now let's clean up the environment

rm(temp.df, unifrac.csv, unifrac.df,
   age.csv, otu.df, otu, clean.names, otu.csv,
   ec.csv, ec.df, ec.clr, ec.cols.to.ignore,
   geography.csv)


yatsunenko <- list(raw.counts = Xraw,
                   age = age,
                   patristic = patristic,
                   unifrac = unifrac,
                   ec = ec.filtered,
                   geography = geography)
