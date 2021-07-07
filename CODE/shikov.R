
library(readxl)
library(dplyr)
library(ggplot2)

data <- read_excel("Shikov et al. 2020 Table S1.xlsx")

load("bsc.1.rda")
bsc1 <- bsc
load("bsc.2.rda")
bsc2 <- bsc
load("bsc.3.rda")
bsc3 <- bsc
load("bsc.4.rda")
bsc4 <- bsc
load("bsc.5.rda")
bsc5 <- bsc
load("bsc.6.rda")
bsc6 <- bsc
bsc6 <- dplyr::filter(bsc6, chr == 6)
load("bsc.7.rda")
bsc7 <- bsc
load("bsc.8.rda")
bsc8 <- bsc
load("bsc.9.rda")
bsc9 <- bsc
load("bsc.10.rda")
bsc10 <- bsc
load("bsc.11.rda")
bsc11 <- bsc
load("bsc.12.rda")
bsc12 <- bsc
load("bsc.13.rda")
bsc13 <- bsc
load("bsc.14.rda")
bsc14 <- bsc
load("bsc.15.rda")
bsc15 <- bsc
load("bsc.16.rda")
bsc16 <- bsc
load("bsc.17.rda")
bsc17 <- bsc
load("bsc.18.rda")
bsc18 <- bsc
load("bsc.19.rda")
bsc19 <- bsc
load("bsc.20.rda")
bsc20 <- bsc
load("bsc.21.rda")
bsc21 <- bsc
load("bsc.22.rda")
bsc22 <- bsc

BSC <- rbind(bsc1,bsc2,bsc3,bsc4,bsc5,bsc6,bsc7,bsc8,bsc9,bsc10,bsc11,bsc12,bsc13,bsc14,bsc15,bsc16,bsc17,bsc18,bsc19,bsc20,bsc21,bsc22)

meanB <- vector()
for(i in 1:nrow(data)){
  chrdata <- data$chr[i]
  posidata <- data$start[i]
  posfdata <- data$stop[i]
  bdata <- dplyr::filter(BSC, (chr == chrdata) & (posi >= posidata) & (posf <= posfdata))
  meanB[i] <- mean(bdata$B)
  print(i/nrow(data))
}

data$B <- meanB


data <- data[-(507:515),]

######################
## PLEIDEGREE = MED ##
######################

meddata <- data.frame(med = as.factor(ifelse(data$clust_med == 2.5, 3, ifelse(data$clust_med == 3.5, 4, ifelse(data$clust_med >= 4.5, 5, data$clust_med)))), B = data$B)

freq <- vector()
for (i in 2:5) {
  dat <- dplyr::filter(meddata, med == i)
  freq <- c(freq, nrow(dat)/nrow(meddata))
}

freq <- data.frame(med = 2:5, freq)


ggplot(data = freq, aes(x = med, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 2:5, labels = c("2", "3", "4", "\u2265 5"))


ggplot() + 
  geom_boxplot(data = meddata, aes(x = med, y = B)) +
  geom_smooth(data = meddata, aes(x = as.numeric(med), y = B), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "B", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = 2:5, labels = c("2", "3", "4", "\u2265 5"),
                   limits = c("2", "3", "4", "5"))

x <- as.numeric(meddata$med)
y <- meddata$B

linearMod <- lm(y ~ x)
summary(linearMod)


######################
## PLEIDEGREE = MAX ##
######################

maxdata <- data.frame(max = as.factor(ifelse(data$clust_max >= 10, 10, data$clust_max)), B = data$B)

freq <- vector()
for (i in 2:10) {
  dat <- dplyr::filter(maxdata, max == i)
  freq <- c(freq, nrow(dat)/nrow(maxdata))
}

freq <- data.frame(max = 2:10, freq)


ggplot(data = freq, aes(x = max, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = 2:10, labels = c("2", "3", "4", "5", "6", "7", "8", "9", "\u2265 10"))


ggplot() + 
  geom_boxplot(data = maxdata, aes(x = max, y = B)) +
  geom_smooth(data = maxdata, aes(x = as.numeric(max), y = B), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "B", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = 2:10, labels = c("2", "3", "4", "5", "6", "7", "8", "9", "\u2265 10"),
                   limits = c("2", "3", "4", "5", "6", "7", "8", "9", "10"))

x <- as.numeric(maxdata$max)
y <- maxdata$B

linearMod <- lm(y ~ x)
summary(linearMod)

