
require(ggplot2)
require(dplyr)


load("dfgenes_37.rda")
dfgenes <- dplyr::filter(dfgenes, !is.na(B))

HLAgenes <- vector()
for(i in 1:nrow(dfgenes)){
  if(is.na(dfgenes$posf[i])){
    if(as.numeric(dfgenes$chr[i]) == 6){
      if((as.numeric(dfgenes$posi[i]) >= 25e6) & (as.numeric(dfgenes$posi[i]) <= 34e6)){
        HLAgenes <- c(HLAgenes, i)
      }
    }
  } else {
    if(as.numeric(dfgenes$chr[i]) == 6){
      if((as.numeric(dfgenes$posi[i]) <= 25e6) & (as.numeric(dfgenes$posf[i]) >= 34e6)){
        HLAgenes <- c(HLAgenes, i)
      } else if (((as.numeric(dfgenes$posi[i]) >= 25e6) & (as.numeric(dfgenes$posi[i]) <= 34e6)) | ((as.numeric(dfgenes$posf[i]) >= 25e6) & (as.numeric(dfgenes$posf[i]) <= 34e6))){
        HLAgenes <- c(HLAgenes, i)
      }
    }
  }
}

dfgenes <- dfgenes[-HLAgenes,]

plei <- vector()
freq <- vector()
for(i in 1:max(dfgenes$plei)){
  a <- dplyr::filter(dfgenes, plei == i)
  freq[i] <- nrow(a)/nrow(dfgenes)
  plei[i] <- i
}

freq <- data.frame(plei = plei, freq = freq)


ggplot(data = freq, aes(x = plei, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", family = "serif", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = c(1:10))

dfgenes$plei2 <- ifelse(dfgenes$plei >= 6, 6, dfgenes$plei)

ggplot() + 
  geom_boxplot(data = dfgenes, aes(x = as.factor(plei2), y = B)) +
  geom_smooth(data = dfgenes, aes(x = as.numeric(as.factor(plei2)), y = B), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "B", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = 1:6, labels = c("1", "2", "3", "4", "5", "\u2265 6"))

x <- as.numeric(dfgenes$plei)
y <- dfgenes$B
z <- log10(dfgenes$RR)
a <- data.frame(x,y,z)
a <- a[is.finite(rowSums(a)),]

linearMod <- lm(y ~ x)
summary(linearMod)

linearMod <- lm(a$y ~ a$x + a$z)
summary(linearMod)




ggplot() + 
  geom_boxplot(data = dfgenes, aes(x = as.factor(plei2), y = log10(RR))) +
  geom_smooth(data = dfgenes, aes(x = as.numeric(as.factor(plei2)), y = log10(RR)), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "log10(Rec. Rate (cM/Mb))", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = 1:6, labels = c("1", "2", "3", "4", "5", "\u2265 6"))


x <- as.numeric(dfgenes$plei)
y <- log10(dfgenes$RR)
z <- dfgenes$B
a <- data.frame(x,y,z)
a <- a[is.finite(rowSums(a)),]

linearMod <- lm(a$y ~ a$x)
summary(linearMod)

linearMod <- lm(a$y ~ a$x + a$z)
summary(linearMod)


ggplot() + 
  geom_point(data = dfgenes, aes(x = B, y = log10(RR))) +
  geom_smooth(data = dfgenes, aes(x = B, y = log10(RR)), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "log10(Recombination Rate (cM/Mb))", x = "B") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", family = "serif", size = 15), plot.title = element_text(size = 16, hjust = 0.5))


x <- dfgenes$B
y <- log10(dfgenes$RR)
z <- as.numeric(dfgenes$plei)
a <- data.frame(x,y,z)
a <- a[is.finite(rowSums(a)),]


linearMod <- lm(y ~ x)
summary(linearMod)

#REGRESSION
linearMod <- lm(a$y ~ a$x + a$z)
summary(linearMod)
