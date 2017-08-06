# This is code to replicate the analyses in Hartmann (subm.)

#################
# PRELIMINARIES #
#################

# load package
if(!is.element("dplyr", installed.packages())) { install.packages("dplyr") }
library(dplyr)

# set options
options(stringsAsFactors = F)
Sys.setlocale("LC_ALL", "de_DE")


############################
# read and manipulate data #
############################

# read data
ungbaby_annotated <- read.csv("ungbaby5.csv", 
                              sep="\t", head=T, quote="'", encoding = "UTF-8")

ungbaby_annotated <- filter(ungbaby_annotated, Compound!="y" & Syntagma!="ja")

# make list of all tokens
alltokens_freq <- c()
for(i in 1:6) {
  period_current <- read.csv(paste("pos_period", i, ".txt", sep="", collapse=""), sep="\t", encoding = "UTF-8", head=F, quote="")
  period_current <- period_current[-grep("\\$", period_current$V2),] #remove punctuation tokens
  alltokens_freq[i] <- sum(as.numeric(period_current$V1))
  rm(period_current)
}

# 50-year periods
periods50 <- seq(1600, 1900, 50)

# 50-year-period-vector
periods50_print <- sapply(1:6, function(i) paste(periods50[i], "-", periods50[i+1], sep=""))

###############################################
# constructions in which ung-nominals appear: #
# P V-ung, plural and determiner construction #
###############################################

# rel. freq. of PREP V-ung construction
prep.ung <- c()

for(i in 1:length(levels(factor(ungbaby_annotated$Period)))) {
  prep.ung[i] <- length(which(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])$Prep_V_ung!="")) /
    nrow(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i]))
}


# plot PREP V-ung
# png("prep_v_ung.png", width=5, height=5, un="in", res=300)
plot(1:6, prep.ung, ylim=c(0,0.2), type="b", xaxt="n", pch=18,
     main="[P V-ung] construction", ylab="Relative Frequency", xlab="Period")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.7)
grid(ny=10, nx=0, col="darkgrey")
# dev.off()

# png("prep_v_ung_DE.png", width=5, height=5, un="in", res=300)
plot(1:6, prep.ung, ylim=c(0,0.2), type="b", xaxt="n", pch=18,
     main="[P V-ung]-Konstruktion", ylab="Relative Frequenz", xlab="Periode", lwd=2,
     col="darkblue")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.7)
grid(ny=10, nx=0, col="darkgrey")
# dev.off()

cor.test(1:6, prep.ung, method="kendall")


# determiners
det.ung <- c()

for(i in 1:length(levels(factor(ungbaby_annotated$Period)))) {
  det.ung[i] <- length(which(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])$Determiner!="")) /
    nrow(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i]))
}

# plot rel. freq. of determiners
# png("determiners_ungbaby.png", width=5, height=5, un="in", res=300)
plot(1:6, det.ung*100, ylim=c(0,100), type="b", xaxt="n", pch=18,
     main="ung-nominals with determiner", ylab="Relative Frequency", xlab="Period", yaxt="n")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.5)
axis(2, at=c(0:10)*10, labels=paste((0:10)*10, " %", sep=""), las=2, cex.axis=0.7)
grid(ny=10, nx=0, col="darkgrey")
# dev.off()

cor.test(1:6, det.ung, method="kendall")


# pluralization
ungbaby_annotated$Number <- NA
for(i in 1:nrow(ungbaby_annotated)) {
  if(ungbaby_annotated$Key[i] %in% grep(".*ungen$", ungbaby_annotated$Key, value=T)) {
    ungbaby_annotated$Number[i] <- "plural"
  } else {
    ungbaby_annotated$Number[i] <- "singular"
  }
  print(paste(i, " of ", nrow(ungbaby_annotated), sep="", collapse=""))
}


pl.ung <- sapply(1:6, function(i) nrow(filter(ungbaby_annotated, Number=="plural" & Period==levels(factor(ungbaby_annotated$Period))[i])) /
  nrow(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])))

# plot frequency of determiner & pluralization

# png("det_plur.png", width=10, height=5, un="in", res=300)
par(mfrow=c(1,2))
plot(1:6, det.ung*100, ylim=c(0,100), type="b", xaxt="n", pch=18,
     main="ung-nominals with determiner", ylab="Relative Frequency", xlab="Period", yaxt="n")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.5)
axis(2, at=c(0:10)*10, labels=paste((0:10)*10, " %", sep=""), las=2, cex.axis=0.7)
grid(ny=10, nx=0, col="darkgrey")

plot(1:6, pl.ung*100, xaxt="n", main="Pluralization", pch=20, type="b", ylim=c(0, 25),
     xlab="Period", ylab="Relative Frequency", yaxt="n")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.5)
axis(2, at=c(0:20)*5, labels=paste((0:20)*5, " %", sep=""), las=2, cex.axis=0.7)
grid(ny=10, nx=0, col="darkgrey")
# dev.off()
par(mfrow=c(1,1))

cor.test(1:6, pl.ung, method="kendall")


#####################################################
# type and token frequency & potential productivity #
#####################################################

# token and type frequency
tok.freq.ungbaby <- type.freq.ungbaby <- c()

for(i in 1:length(levels(factor(ungbaby_annotated$Period)))) {
  tok.freq.ungbaby[i] <-  length(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])$Lemma) / 
    alltokens_freq[i]
  
  type.freq.ungbaby[i] <- length(levels(factor(filter(
      ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])$Lemma))) / 
    alltokens_freq[i]
}

# png("rel_type_tok_freq_ungbaby.png", width=13, height=6.5, un="in", res=300)
par(mfrow=c(1,2))
plot(1:6, tok.freq.ungbaby, ylim=c(0,0.03), xaxt="n", 
     main=expression(paste("Token Frequency, DTA", italic("baby"), sep="", collapse="")),
     xlab="Period", ylab="Relative Frequency", type="b", lwd=2, pch=18)
axis(1, at=c(1:6), labels=periods50_print, cex.axis=0.5)
grid(ny=10, nx=0, col="darkgrey")

plot(1:6, type.freq.ungbaby, xaxt="n", 
     main=expression(paste("Type Frequency, DTA", italic("baby"), sep="", collapse="")),
     xlab="Period", ylab="Relative Frequency", type="b", lwd=2, pch=18, ylim=c(0,0.005))
axis(1, at=c(1:6), labels=periods50_print, cex.axis=0.5)
grid(ny=10, nx=0, col="darkgrey")
# dev.off()
 par(mfrow=c(1,1))

cor.test(1:6, tok.freq.ungbaby, method="kendall")
cor.test(1:6, type.freq.ungbaby, method="kendall")

# potential productivity
pot.prod.ungbaby.annotated <- c()
hapaxes_ungbaby <- as.character(filter(as.data.frame(table(ungbaby_annotated$Lemma)), 
                                       Freq==1)$Var1)

for(i in 1:length(levels(factor(ungbaby_annotated$Period)))) {
  pot.prod.ungbaby.annotated[i] <- length(filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i] &
                  Lemma %in% hapaxes_ungbaby)$Lemma) / length(
                    filter(ungbaby_annotated, Period==levels(factor(ungbaby_annotated$Period))[i])$Lemma)
}


plot(1:6, pot.prod.ungbaby.annotated, xaxt="n", main="Potential Productivity",
     pch=18, type="b", ylim=c(0,0.1), ylab="Potential Productivity", xlab="Period")
axis(1, at=c(1:6), labels = periods50_print, cex.axis=0.5)
grid(ny=10, nx=0, col="darkgrey")

