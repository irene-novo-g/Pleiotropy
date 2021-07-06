
set.seed(1234)

require(plyr)
require(dplyr)
load("final.rda")


final$pleidegree <- factor(ifelse(as.numeric(final$pleidegree) >= 6, 6, final$pleidegree))
final <- dplyr::filter(final, !((final$chr == 6) & (as.numeric(as.character(final$pos37)) >= 25e6) & (as.numeric(as.character(final$pos37)) <= 34e6)))


mean.effect <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  mean.effect[[i]] <- mean(data$effect)
  print(i/6)
}
sd.effect <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.effect[[i]] <- sd(data$effect)
  print(i/6)
}

pleidegree <- 1:6
mean.effect <- data.frame(pleidegree, mean.effect, sd.effect)

gen <- as.character(unique(final$gene))
mean.effect2 <- vector()
for (i in 1:length(gen)) {
  data <- dplyr::filter(final, gene == gen[[i]])
  mean.effect2[[i]] <- mean(data$effect)
  print(i/length(gen))
}

mean.effect2 <- data.frame(gen, mean.effect2)
colnames(mean.effect2) <- c("gene", "mean.effect")

final <- merge(final, mean.effect2, by = "gene")

mean.effect3 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  mean.effect3[[i]] <- mean(data$mean.effect)
}

sd.effect3 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.effect3[[i]] <- sd(data$mean.effect)
}

mean.effect3 <- data.frame(pleidegree, mean.effect3, sd.effect3)

gen <- as.character(unique(final$gene))
sd.effect1 <- vector()
for (i in 1:length(gen)) {
  data <- dplyr::filter(final, gene == gen[[i]])
  sd.effect1[[i]] <- sd(data$effect)
  print(i/length(gen))
}

sd.effect1 <- data.frame(gen, sd.effect1)
colnames(sd.effect1) <- c("gene", "sd.effect")

final <- merge(final, sd.effect1, by = "gene")

sd.effect2 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.effect2[[i]] <- mean(data$sd.effect)
}

mean.effect3 <- data.frame(mean.effect3, sd.effect2)

ggplot(data = mean.effect3, aes(x = pleidegree, y = sd.effect2)) + 
  geom_point(colour = "black") +
  geom_smooth(method = "lm", se=FALSE, formula = y ~ x) +
  labs(x = "Pleiotropy degree", y = "\u03C3 (effect)") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 1:6, labels = c("1", "2", "3", "4", "5", "\u2265 6"))

x = mean.effect3$pleidegree
y = mean.effect3$sd.effect2

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)



lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*"", 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ggplot(data = mean.effect3, aes(x = pleidegree, y = sd.effect3)) + 
  geom_point(colour = "black") +
  geom_smooth(method = "lm", se=FALSE, formula = y ~ x) +
  labs(x = "Pleiotropy degree", y = "\u03C3 (effect)") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 1:6, labels = c("1", "2", "3", "4", "5", "\u2265 6")) 

x = mean.effect3$pleidegree
y = mean.effect3$sd.effect3

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)


data1 <- subset(final, select = c(pleidegree, mean.effect))
data1 <- dplyr::mutate(data1, duplicated(data1))
data1 <- subset(data1, data1$`duplicated(data1)` == TRUE)
fig1a <- ggplot() + 
  geom_boxplot(data = data1, aes(x = pleidegree, y = log10(mean.effect))) +
  geom_smooth(data = data1, aes(x = as.numeric(pleidegree), y = log10(mean.effect)), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "log10(effect)", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = c(1:6), labels = c("1", "2", "3", "4", "5", "\u2265 6"))

x <- as.numeric(data1$pleidegree)
y <- log10(data1$mean.effect)

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)


mean.h <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  mean.h[[i]] <- mean(data$Va)
  print(i/6)
}
sd.h <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.h[[i]] <- sd(data$Va)
  print(i/6)
}

pleidegree <- 1:6
mean.h <- data.frame(pleidegree, mean.h, sd.h)

gen <- as.character(unique(final$gene))
mean.h2 <- vector()
for (i in 1:length(gen)) {
  data <- dplyr::filter(final, gene == gen[[i]])
  mean.h2[[i]] <- mean(data$Va)
  print(i/length(gen))
}

mean.h2 <- data.frame(gen, mean.h2)
colnames(mean.h2) <- c("gene", "mean.h")

final <- merge(final, mean.h2, by = "gene")

mean.h3 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  mean.h3[[i]] <- mean(data$mean.h)
  print(i/6)
}

sd.h3 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.h3[[i]] <- sd(data$mean.h)
  print(i/6)
}

pleidegree <- 1:6
mean.h3 <- data.frame(pleidegree, mean.h3, sd.h3)

gen <- as.character(unique(final$gene))
sd.h1 <- vector()
for (i in 1:length(gen)) {
  data <- dplyr::filter(final, gene == gen[[i]])
  sd.h1[[i]] <- sd(data$Va)
  print(i/length(gen))
}

sd.h1 <- data.frame(gen, sd.h1)
colnames(sd.h1) <- c("gene", "sd.h")

final <- merge(final, sd.h1, by = "gene")

sd.h2 <- vector()
for (i in 1:6) {
  data <- dplyr::filter(final, pleidegree == i)
  sd.h2[[i]] <- mean(data$sd.h)
  print(i/6)
}

mean.h3 <- data.frame(mean.h3, sd.h2)

data2 <- subset(final, select = c(pleidegree, mean.h))
data2 <- dplyr::mutate(data2, duplicated(data2))
data2 <- subset(data2, data2$`duplicated(data2)` == TRUE)

ggplot() + 
  geom_boxplot(data = data2, aes(x = pleidegree, y = mean.h)) +
  geom_smooth(data = data2, aes(x = as.numeric(pleidegree), y = mean.h), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "h\u00B2", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = c(1:6), labels = c("1", "2", "3", "4", "5", "\u2265 6")) +
  scale_y_continuous(limits = c(0, 0.005))

x = as.numeric(data2$pleidegree)
y = data2$mean.h

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)


uniquegene <- dplyr::select(final, gene, pleidegree)
uniquegene <- dplyr::distinct(uniquegene)
nrow(uniquegene) == length(unique(uniquegene$gene))

freq <- vector()
for (i in 1:6) {
  data <- dplyr::filter(uniquegene, pleidegree == i)
  freq[[i]] <- nrow(data)/nrow(uniquegene)
  print(i/6)
}

freq <- data.frame(pleidegree, freq)

ggplot(data = freq, aes(x = pleidegree, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 1:6, labels = c("1", "2", "3", "4", "5", "\u2265 6"))


data <- dplyr::mutate(final, as.numeric(pleidegree))
cor <- vector()
for (i in 1:5) {
  data2 <- dplyr::filter(data, data$`as.numeric(pleidegree)` >= i)
  cor[[i]] <- cor(data2$effect, data2$`as.numeric(pleidegree)`)
  print(i/5)
}

cor.mean <- vector()
for (i in 1:5) {
  data2 <- subset(data, data$`as.numeric(pleidegree)` >= i, select = c(mean.effect, gene, `as.numeric(pleidegree)`))
  data2 <- dplyr::mutate(data2, duplicated(data2$gene))
  data2 <- subset(data2, data2$`duplicated(data2$gene)` == FALSE)
  cor.mean[[i]] <- cor(data2$mean.effect, data2$`as.numeric(pleidegree)`)
  print(i/5)
}

cor.df <- data.frame(1:5, cor, cor.mean)
colnames(cor.df) <- c("pleidegree", "cor", "cor.mean")

ggplot(data = cor.df, aes(x = pleidegree, y = cor)) + 
  geom_point(colour = "black") + 
  geom_smooth(method = "lm", se=FALSE, formula = y ~ x) +
  labs(x = "Pleiotropy degree", y = "r (pleiotropy degree vs. effect) ") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 1:5, labels = c("\u2265 1", "\u2265 2", "\u2265 3", "\u2265 4", "\u2265 5"))

x <- cor.df$pleidegree
y <- cor.df$cor

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)


MAF <- vector()
for(i in 1:nrow(final)){
  if(final$q[i] > 0.5){
    MAF[i] <- 1-final$q[i]
  } else {
    MAF[i] <- final$q[i]
  }
}

final$MAF <- MAF

data1 <- subset(final, select = c(pleidegree, MAF))
ggplot() + 
  geom_boxplot(data = data1, aes(x = pleidegree, y = MAF)) +
  geom_smooth(data = data1, aes(x = as.numeric(pleidegree), y = MAF), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "MAF", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = c(1:6), labels = c("1", "2", "3", "4", "5", "\u2265 6"))

x <- as.numeric(data1$pleidegree)
y <- data1$MAF

cor <- cor(x, y, method = "pearson", use = "na.or.complete")
cors <- vector()
for(i in 1:10000){
  reorder <- sample(y, replace = F)
  cors[i] <- cor(x, reorder, method = "pearson", use = "na.or.complete")
  cors <- sort(cors)
}
pos <- which.min(abs(cors-cor))
P <- pos/10000
if(P > 0.5){
  P <- 1 - P
}
cor
P

linearMod <- lm(y ~ x)
summary(linearMod)


trait.type <- ifelse(as.character(final$trait.type) == "Neurological", "Neurological/Psychiatric", ifelse(as.character(final$trait.type) == "Psychiatric", "Neurological/Psychiatric", as.character(final$trait.type)))
final$trait.type <- trait.type

plei <- dplyr::filter(final, pleidegree != 1)
final$pleiotropy <- ifelse(final$gene %in% plei$gene, TRUE, FALSE)

traits.genes <- ddply(final, .(cluster, pleiotropy, trait.type), nrow)
noplei.trait <- subset(traits.genes, pleiotropy == F, select = c(cluster, trait.type, V1))
colnames(noplei.trait) <- c("cluster", "fd", "noplei.genes")
plei.trait <- subset(traits.genes, pleiotropy == T, select = c(cluster, V1))
colnames(plei.trait) <- c("cluster", "plei.genes")
traits.genes <- merge(plei.trait, noplei.trait, by = "cluster")

traits.genes <- dplyr::mutate(traits.genes, 100*(plei.genes/(plei.genes+noplei.genes)))
colnames(traits.genes)[5] <- "plei.perc"


colores <- c("Dermatological" = "plum1", 
             "Cardiovascular" = "darkgoldenrod1", 
             "Metabolic" = "khaki1", 
             "Neurological/Psychiatric" = "lightskyblue3", 
             "Skeletal" = "grey90", 
             "Cancer" = "darkseagreen1", 
             "Gastrointestinal" = "darkgoldenrod", 
             "Hematological" = "brown1", 
             "Immunological" = "burlywood1", 
             "Endocrine" = "gold")

traits.genes$cluster <- factor(traits.genes$cluster, levels = c("Basal cell carcinoma","Chronic lymphocytic leukemia","Lung cancer","Prostate cancer","Prostate-specific antigen levels","Testicular germ cell tumor",
                                                                   "Atrial fibrillation","Coronary artery disease","Coronary heart disease","Myocardial infarction",
                                                                   "Atopic dermatitis","Psoriasis","Vitiligo",
                                                                   "Body mass index","Menarche (age at onset)","Obesity","Type 2 diabetes",
                                                                   "Digestive disease","Ulcerative colitis",
                                                                   "Glycated hemoglobin levels","Mean platelet volume","Monocyte count","Neutrophil traits","Red blood cell traits",
                                                                   "Primary biliary cholangitis","Rheumatoid arthritis","Systemic lupus erythematosus","Type 1 diabetes",
                                                                   "Cholesterol","HDL","Triglycerides","Urate levels",
                                                                   "Migraine","Parkinson's disease",
                                                                   "Schizophrenia",
                                                                   "Bone mineral density","Height","Waist circumference","Waist-related traits","Waist-hip ratio","Waist-to-hip-related traits"))




ggplot(data = traits.genes, aes(x = cluster, weight = plei.perc, fill = fd)) + 
  geom_bar(colour = "black") + 
  labs(fill = "Functional domain", x = "Trait", y = "Pleiotropic genes (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 18), plot.title = element_text(size = 25, hjust = 0.5)) +
  scale_fill_manual(values = colores) 


