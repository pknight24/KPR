library(KPR)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(viridis)
rm(list = ls())

### Load data
data("yatsunenko")

counts <- yatsunenko$raw.counts
age <- yatsunenko$age
patristic <- yatsunenko$patristic
unifrac <- yatsunenko$unifrac
ec <- yatsunenko$ec

rm(yatsunenko)

### Prepare data

H <- generateSimilarityKernel(unifrac)
Q <- generateSimilarityKernel(patristic, squareValues = FALSE)

counts.clr <- log(counts + 1) - apply(log(counts + 1), 1, mean)
counts.final <- apply(counts.clr, 2, function(x) x - mean(x))

ec.clr <- log(ec + 1) - apply(log(ec + 1), 1, mean)
ec.final <- apply(ec.clr, 2, function(x) x - mean(x))

Y <- age - mean(age)

Y.perm <- Y[permute::shuffle(length(Y))]


bact.nameid.csv <- read.csv("C:/Users/pknight/Dropbox/Biplots_generalized/Data/Yatsunenko/used_in_KPR_paper/bacteria_id.csv", stringsAsFactors = FALSE)$tax

taxonomy.mat <- t(sapply(strsplit(bact.nameid.csv, ";"), function(x)
{
    sapply(1:6, function(i)
    {
        substr(x[[i]], start=4, stop = nchar(x[[i]]))

    })

}))

colnames(taxonomy.mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax <- as.data.frame(taxonomy.mat, stringsAsFactors = FALSE) %>%
    filter(Genus %in% colnames(counts)) %>% distinct(Genus, .keep_all = TRUE) %>%
    (function(x) x[match(colnames(counts), x$Genus),])

random <- matrix(rnorm(100 * 149), 100, 149)
H.rand <- random %*% t(random)

### Model fitting

kpr.out <- KPR(designMatrix = counts.final, Y = Y, Q = Q, H = H, scale=TRUE)

infer.out <- kpr.out$p.values

effects <- kpr.out$beta.hat

results.df <- data.frame(tax, effects, infer.out) %>% filter(effects < 0.8) %>%
  group_by(Class) %>% mutate(numInClass = n()) %>% filter(numInClass > 2) %>%
  select(-numInClass) %>% ungroup
idx <- 1:nrow(results.df)


### Plotting pvals and coefs

g <- ggplot(data = results.df) + theme_minimal()

effect.plot <- g + geom_point(aes(x = idx, y = effects, col=Class ), size=1.5) +
    xlab("Index") + ylab("Effect size") +
    geom_hline(yintercept = 0) +
    theme(legend.position = "none")

pval.plot <- g + geom_point(aes(x = idx, y = infer.out, col=Class ), size=1.5) +
    xlab("") + ylab("P-value") +
    theme(legend.position = "none")

cowplot::plot_grid(pval.plot, effect.plot, nrow=2)

### Tree plot

load("poster/subset110_genus.RData")

pruned <- prune_taxa(results.df$Genus, subset110_genus)

plot_tree(pruned, color="Class", nodelabf=nodeplotblank)

### Kernel heatmaps

K <- counts.final %*% t(counts.final)
K.sorted <- K[order(Y), order(Y)]
# heatmap(K.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))

K.Q <- counts.final %*% Q %*% t(counts.final)
K.Q.sorted <- K.Q[order(Y), order(Y)]
#heatmap(K.Q.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))

K.Z <- ec.final %*% t(ec.final)
K.Z.sorted <- K.Z[order(Y), order(Y)]
# heatmap(K.Z.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))
