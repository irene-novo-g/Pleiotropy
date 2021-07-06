
require(ggplot2)
require(dplyr)

load("final.rda")

gene <- dplyr::filter(final, pleidegree == 10)

#NCP = r2 * n * h2 = r2 * n * 2pq * beta2
#r2 = 1

gene$NCP <- gene$N * 2*gene$q*(1-gene$q) * (gene$effect^2)
gene$NCPlowef <- gene$N * 2*gene$q*(1-gene$q) * ((gene$effect/2)^2)

gene$nforNCP80 <- 80 / (2*gene$q*(1-gene$q) * (gene$effect^2))
gene$nforNCP80lowef <- 80 / (2*gene$q*(1-gene$q) * ((gene$effect/2)^2))

gene$ind <- 1:nrow(gene)
