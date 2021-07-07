
library(readxl)
library(dplyr)
library(ggplot2)
library(readr)

data <- read_excel("Watanabe.xlsx", skip = 1, sheet = "ST 12")

SNPs_plei <- data.frame(rs = data$rsID, traits = data$`#traits`, domains = data$`#domains`)
SNPs_pos <- read_table2("SNPs_watanabe_positions_GRCh37")

plei <- vector()
plei2 <- vector()
for(i in 1:nrow(SNPs_pos)){
  ind <- which(as.character(SNPs_plei$rs) == as.character(SNPs_pos$ID[i]))
  plei[i] <- SNPs_plei$traits[ind]
  plei2[i] <- SNPs_plei$domains[ind]
  print(i/nrow(SNPs_pos))
}

colnames(SNPs_pos) <- c("chr", "pos", "rs")
SNPs_pos$traits <- plei
SNPs_pos$domains <- plei2

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

B <- vector()
for(i in 1:nrow(SNPs_pos)){
  chrdata <- SNPs_pos$chr[i]
  posdata <- SNPs_pos$pos[i]
  bdata <- dplyr::filter(BSC, chr == chrdata)
  bdata <- bdata[which.min(abs(bdata$posi - posdata)),]
  if(nrow(bdata) == 0){
    B[i] <- NA
  } else {
    B[i] <- bdata$B
  }
  print(i/nrow(SNPs_pos))
}

SNPs_pos$B <- B

meanB <- vector()
for(i in 1:nrow(SNPs_pos)){
  chrdata <- SNPs_pos$chr[i]
  posidata <- SNPs_pos$pos[i] - 250000
  posfdata <- SNPs_pos$pos[i] + 250000
  bdata <- dplyr::filter(BSC, (chr == chrdata) & (posi >= posidata) & (posf <= posfdata))
  meanB[i] <- mean(bdata$B)
  print(i/nrow(SNPs_pos))
}

SNPs_pos$meanB <- meanB


chr1 <- read.csv(paste(PATH, "genetic_map_chr1.txt", sep = ""), sep="")
chr2 <- read.csv(paste(PATH, "genetic_map_chr2.txt", sep=""), sep="")
chr3 <- read.csv(paste(PATH, "genetic_map_chr3.txt", sep=""), sep="")
chr4 <- read.csv(paste(PATH, "genetic_map_chr4.txt", sep=""), sep="")
chr5 <- read.csv(paste(PATH, "genetic_map_chr5.txt", sep=""), sep="")
chr6 <- read.csv(paste(PATH, "genetic_map_chr6.txt", sep=""), sep="")
chr7 <- read.csv(paste(PATH, "genetic_map_chr7.txt", sep=""), sep="")
chr8 <- read.csv(paste(PATH, "genetic_map_chr8.txt", sep=""), sep="")
chr9 <- read.csv(paste(PATH, "genetic_map_chr9.txt", sep=""), sep="")
chr10 <- read.csv(paste(PATH, "genetic_map_chr10.txt", sep=""), sep="")
chr11 <- read.csv(paste(PATH, "genetic_map_chr11.txt", sep=""), sep="")
chr12 <- read.csv(paste(PATH, "genetic_map_chr12.txt", sep=""), sep="")
chr13 <- read.csv(paste(PATH, "genetic_map_chr13.txt", sep=""), sep="")
chr14 <- read.csv(paste(PATH, "genetic_map_chr14.txt", sep=""), sep="")
chr15 <- read.csv(paste(PATH, "genetic_map_chr15.txt", sep=""), sep="")
chr16 <- read.csv(paste(PATH, "genetic_map_chr16.txt", sep=""), sep="")
chr17 <- read.csv(paste(PATH, "genetic_map_chr17.txt", sep=""), sep="")
chr18 <- read.csv(paste(PATH, "genetic_map_chr18.txt", sep=""), sep="")
chr19 <- read.csv(paste(PATH, "genetic_map_chr19.txt", sep=""), sep="")
chr20 <- read.csv(paste(PATH, "genetic_map_chr20.txt", sep=""), sep="")
chr21 <- read.csv(paste(PATH, "genetic_map_chr21.txt", sep=""), sep="")
chr22 <- read.csv(paste(PATH, "genetic_map_chr22.txt", sep=""), sep="")
chr1$chr = 1
chr2$chr = 2
chr3$chr = 3
chr4$chr = 4
chr5$chr = 5
chr6$chr = 6
chr7$chr = 7
chr8$chr = 8
chr9$chr = 9
chr10$chr = 10
chr11$chr = 11
chr12$chr = 12
chr13$chr = 13
chr14$chr = 14
chr15$chr = 15
chr16$chr = 16
chr17$chr = 17
chr18$chr = 18
chr19$chr = 19
chr20$chr = 20
chr21$chr = 21
chr22$chr = 22
map <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)

RR <- vector()
for(i in 1:nrow(SNPs_pos)){
  chrdata <- as.numeric(SNPs_pos$chr[i])
  posdata <- as.numeric(SNPs_pos$pos[i])
  RRdata <- dplyr::filter(map, (chr == chrdata) & (position >= posdata))
  if(nrow(RRdata) == 0){
    RR[i] <- NA
    print("RRdata has 0 lines")
  } else {
    RRdata <- RRdata[1,]
    RR[i] <- RRdata$COMBINED_rate.cM.Mb.
  }     
  print(i/nrow(SNPs_pos))
}

SNPs_pos$RR <- RR


SNPs_pos <- dplyr::filter(SNPs_pos, !((chr == 6) & (pos >= 25e6) & (pos <= 36e6)))

plei <- ifelse(SNPs_pos$domains >= 4, 4, SNPs_pos$domains)
SNPs_pos$domains2 <- plei

plei <- vector()
freq <- vector()
for(i in 1:max(SNPs_pos$domains)){
  a <- dplyr::filter(SNPs_pos, domains == i)
  freq[i] <- nrow(a)/nrow(SNPs_pos)
  plei[i] <- i
}

freq <- data.frame(plei = plei, freq = freq)


ggplot(data = freq, aes(x = plei, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous(breaks = c(1:12))


ggplot() + 
  geom_boxplot(data = SNPs_pos, aes(x = as.factor(domains2), y = B)) +
  geom_smooth(data = SNPs_pos, aes(x = as.numeric(as.factor(domains2)), y = B), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "B", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = 1:4, labels = c("1", "2", "3", "\u2265 4"))

x <- as.numeric(SNPs_pos$domains)
y <- log10(SNPs_pos$RR)
z <- SNPs_pos$B
a <- data.frame(x,y,z)
a <- a[is.finite(rowSums(a)),]

linearMod <- lm(a$y ~ a$x)
summary(linearMod)

linearMod <- lm(a$y ~ a$x + a$z)
summary(linearMod)



plei <- ifelse(SNPs_pos$traits >= 20, 20, SNPs_pos$traits)
SNPs_pos$traits2 <- plei

plei <- vector()
freq <- vector()
for(i in 2:max(SNPs_pos$traits)){
  a <- dplyr::filter(SNPs_pos, traits == i)
  freq[i] <- nrow(a)/nrow(SNPs_pos)
  plei[i] <- i
}

freq <- data.frame(plei = plei, freq = freq)


ggplot(data = freq, aes(x = plei, weight = freq)) + 
  geom_bar(colour = "black", fill = "grey69") + 
  labs(x = "Pleiotropy degree", y = "Frequency") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_continuous()

ggplot() + 
  geom_boxplot(data = SNPs_pos, aes(x = as.factor(traits2), y = B)) +
  geom_smooth(data = SNPs_pos, aes(x = as.numeric(as.factor(traits2)), y = B), method = "lm", se=FALSE, formula = y ~ x) +
  labs(y = "Regional B", x = "Pleiotropy degree") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "grey60"), text = element_text(face = "bold", size = 15), plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_x_discrete(breaks = seq(from = 2, to = 20, by = 2), labels = c("2", "4", "6", "8", "10", "12", "14", "16", "18", "\u2265 20"))

x <- as.numeric(SNPs_pos$traits)
y <- SNPs_pos$B
z <- log10(SNPs_pos$RR)
a <- data.frame(x,y,z)
a <- a[is.finite(rowSums(a)),]


linearMod <- lm(y ~ x)
summary(linearMod)

linearMod <- lm(a$y ~ a$x + a$z)
summary(linearMod)




