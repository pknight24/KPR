rm(list = ls())
devtools::load_all(".")
library(ggplot2)
library(dplyr)
data(flax)
attach(flax)


dat <- data.frame(
  linear.fit.id = KPR(designMatrix = X, Y = Y, inference=FALSE)$beta.hat,
  linear.fit.Q = KPR(designMatrix = X, Y = Y, Q = Q, inference=FALSE)$beta.hat,
  linear.fit.H = KPR(designMatrix = X, Y = Y, H = H, inference = FALSE)$beta.hat,
  linear.fit.H.Q = KPR(designMatrix = X, Y = Y, Q = Q, H = H, inference = FALSE)$beta.hat,
  normal.fit.id = KPR(designMatrix = X, Y = Y, inference=FALSE, linear_solve = FALSE)$beta.hat,
  normal.fit.Q = KPR(designMatrix = X, Y = Y, Q = Q, inference=FALSE, linear_solve=FALSE)$beta.hat,
  normal.fit.H = KPR(designMatrix = X, Y = Y, H = H, inference = FALSE, linear_solve = FALSE)$beta.hat,
  normal.fit.H.Q = KPR(designMatrix = X, Y = Y, Q = Q, H = H, inference = FALSE, linear_solve = FALSE)$beta.hat
)
dat %>% mutate(index = 1:n()) %>%
  tidyr::pivot_longer(cols = -index,
                            names_to = c("method", "penalty"), values_to="betahat",
                            names_pattern="(.+).fit.(.+)") %>%
  ggplot(mapping = aes(x = index, y = betahat)) +
    geom_point(aes(col=method)) + facet_wrap(~penalty)
