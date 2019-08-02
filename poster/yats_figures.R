library(KPR)
library(dplyr)
library(ggplot2)
library(phyloseq)
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

kpr.out <- KPR(designMatrix = counts.final, Y = Y, Q = Q)

infer.out <- inference(kpr.out, method = "GMD")

effects <- kpr.out$beta.hat

results.df <- data.frame(tax, effects, infer.out) %>% filter(effects < 0.8)
idx <- 1:nrow(results.df)


# par(mfrow=c(2,1))
# plot(infer.out, col=as.factor(tax$Class))
# plot(effects, type="l")
# abline(a = 0, b = 0, col="red")
# par(mfrow=c(1,1))

g <- ggplot(data = results.df) + theme_minimal()

effect.plot <- g + geom_point(aes(x = idx, y = effects, col=Class ), size=3) +
    xlab("Index") + ylab("Effect size") +
    geom_hline(yintercept = 0)

pval.plot <- g + geom_point(aes(x = idx, y = infer.out, col=Class ), size=3) +
    xlab("Index") + ylab("P-value") +
    theme(legend.position = "none")

cowplot::plot_grid(pval.plot, effect.plot, nrow=2)
