library(KPR)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(viridis)
library(igraph)
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

# Phylogenetics

kpr.out <- KPR(designMatrix = counts.final, Y = Y, scale=FALSE)
kpr.out.Q <- KPR(designMatrix = counts.final, Y = Y, Q = Q, scale=FALSE)
infer.out <- kpr.out$p.values
infer.out.Q <- kpr.out.Q$p.values

effects <- kpr.out$beta.hat
effects.Q <- kpr.out.Q$beta.hat

results.df <- data.frame(tax, effects, infer.out) %>% filter(effects < 0.8) %>%
  group_by(Class) %>% mutate(numInClass = n()) %>% filter(numInClass > 2) %>%
  select(-numInClass) %>% ungroup
results.df.Q <- data.frame(tax, effects.Q, infer.out.Q) %>% filter(effects.Q < 0.8) %>%
  group_by(Class) %>% mutate(numInClass = n()) %>% filter(numInClass > 2) %>%
  select(-numInClass) %>% ungroup

idx <- 1:nrow(results.df)


### Plotting pvals and coefs

effect.plot <- ggplot(data = results.df) +
  geom_point(aes(x = idx, y = effects, col=Class, shape = infer.out < 0.05), size=5.5) +
    xlab("") + ylab("") +
    geom_hline(yintercept = 0) +
    ggtitle("without Q") +
    ylim(-.6, .6) +
    scale_shape_manual(values = c(1, 18)) +
    theme_minimal() +
    theme(legend.position = "none", axis.text=element_text(size=16),
          axis.title = element_text(size=24), plot.title = element_text(size=30, hjust=0.5))

effect.plot.Q <- ggplot(data = results.df.Q) +
  geom_point(aes(x = idx, y = effects.Q, col=Class, shape = infer.out.Q < 0.05), size=5.5) +
    xlab("") + ylab("Effect Size") +
    geom_hline(yintercept = 0) +
    ggtitle("with Q") +
    ylim(-.6, .6) +
    scale_shape_manual(values = c(1, 18)) +
    theme_minimal() +
    theme(legend.position = "none", axis.text=element_text(size=16),
          axis.title = element_text(size=24), plot.title = element_text(size=30, hjust=0.5))


##### Save plots
tiff("~/KPR/poster/figures/coef_withoutQ.tiff", width=1600, height=1600, res=300, compression = "none")
print(effect.plot)
dev.off()

tiff("~/KPR/poster/figures/coef_withQ.tiff", width=1600, height=1600, res=300, compression = "none")
print(effect.plot.Q)
dev.off()

tiff("~/KPR/poster/figures/biplot_withoutQ.tiff", width=1600, height=1600, res=300, compression = "none")
biplot(kpr.out)
dev.off()

tiff("~/KPR/poster/figures/biplot_withQ.tiff", width=1600, height=1600, res=300, compression = "none")
biplot(kpr.out.Q)
dev.off()
### Tree plot

load("poster/subset110_genus.RData")

pruned <- prune_taxa(results.df$Genus, subset110_genus)

tiff("~/KPR/poster/figures/phylo_tree.tiff", width=1600, height=1600, res=300, compression = "none")
print(plot_tree(pruned, color="Class", nodelabf=nodeplotblank))
dev.off()

### Graph generation
set.seed(42)
g <- erdos.renyi.game(75, 0.02, type = "gnp") %>%
  (function(x) delete.vertices(x, degree(x) == 0))
com <- walktrap.community(g)

tiff("~/KPR/poster/figures/graph.tiff", width=1600, height=1600, res=300, compression = "none")
plot(g, edge.arrow.size = 0, vertex.label = NA,
     vertex.color=membership(com), vertex.size = 10)
dev.off()

### Multi-omic integration

order <- findInterval(Y, sort(Y))
sample.col <- plasma(length(order))[order]

K <- counts.final %*% t(counts.final)
K.sorted <- K[order(Y), order(Y)]
# heatmap(K.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))
eigen.K <- eigen(K)
# plot(x = eigen.K$vectors[,1] * sqrt(eigen.K$values[1]),
#      y = eigen.K$vectors[,2] * sqrt(eigen.K$values[2]),
#      col = sample.col, pch=19)

K.Q <- counts.final %*% Q %*% t(counts.final)
K.Q.sorted <- K.Q[order(Y), order(Y)]
# heatmap(K.Q.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))
eigen.K.Q <- eigen(K.Q)
# plot(x = eigen.K.Q$vectors[,1] * sqrt(eigen.K.Q$values[1]),
#     y = eigen.K.Q$vectors[,2] * sqrt(eigen.K.Q$values[2]),
#     col = sample.col, pch=19)

K.Z <- ec.final %*% t(ec.final)
K.Z.sorted <- K.Z[order(Y), order(Y)]
# heatmap(K.Z.sorted, Colv="Rowv", Rowv=NA, symm=TRUE, col=magma(256))
eigen.K.Z <- eigen(K.Z)
# plot(x = eigen.K.Z$vectors[,1] * sqrt(eigen.K.Z$values[1]),
#     y = eigen.K.Z$vectors[,2] * sqrt(eigen.K.Z$values[2]),
#     col = sample.col, pch=19)
