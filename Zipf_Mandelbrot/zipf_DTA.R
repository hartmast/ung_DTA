
# This is code for replicating the results reported on in Hartmann (forthc.):
# Derivational morphology in flux. A case study on word-formation change in German,
# to appear in: Cognitive Linguistics.

# This script dates back to a period when I knew next to nothing about big data
# handling in R, which is why it takes considerable time to run.
# As it works with randomized datasets, there will be differences between your
# results and the results reported in the paper.


#################
# Preliminaries #
#################

# install and load packages
sapply(c("dplyr", "reshape2", "lme4", "zipfR", "data.table", "compiler", "ggplot2", "gridExtra", "psych"), 
       function(x) if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))
lapply(list("dplyr", "reshape2", "lme4", "zipfR", "data.table", "compiler", "ggplot2", "gridExtra", "psych"), 
       require, character.only=T)


Sys.setlocale("LC_ALL", "de_DE")


###############################################
# Zipf-Mandelbrot model without bootstrapping #
###############################################

# read in files with all ung-nouns in the DTA (semi-automatically corrected)

centurydta1 <- read.table("centurydta1.csv", sep = "\t", head = T, quote = "",
                          fill = T, encoding = "UTF-8")
centurydta2 <- read.table("centurydta2.csv", sep = "\t", head = T, quote = "",
                          fill = T, encoding = "UTF-8")
centurydta3 <- rbind(read.table("centurydta3a.csv", sep = "\t", head = T, quote = "",
                          fill = T, encoding = "UTF-8"),
                     read.table("centurydta3b.csv", sep = "\t", head = T, quote = "",
                                fill = T, encoding = "UTF-8"))


# randomize datasets
centurydta1_random <- centurydta1[sample(1:nrow(centurydta1)),]
centurydta2_random <- centurydta2[sample(1:nrow(centurydta2)),]
centurydta3_random <- centurydta3[sample(1:nrow(centurydta3)),]



# preparing emp.vgc for non-randomized version

for(j in 1:3) {
  current_sequence <- seq(1000, nrow(get(paste("centurydta", j, sep="", collapse=""))), 1000)
  vgc_df <- as.data.frame(matrix(nrow=length(current_sequence), ncol=3))
  colnames(vgc_df) <- c("N", "V", "V1")
  
  current_vector <- tolower(get(paste("centurydta", j, sep="", collapse=""))$V2)
  
  for(i in 1:length(current_sequence)) {
    vgc_df$N[i] <- current_sequence[i]
    vgc_df$V[i] <- length(levels(factor(current_vector[1:current_sequence[i]])))
    vgc_df$V1[i] <- length(as.data.frame(table(current_vector[1:current_sequence[i]]))[which(as.data.frame(table(current_vector[1:current_sequence[i]]))$Freq==1),1])
    
  }
  
  write.table(vgc_df, file=paste("dtacentury", j, ".emp.vgc", sep="", collapse=""),
              row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
  
  
}


# preparing emp.vgc for randomized version

for(j in 1:3) {
  current_sequence <- seq(1000, nrow(get(paste("centurydta", j, "_random", 
                                               sep="", collapse=""))), 1000)
  vgc_df <- as.data.frame(matrix(nrow=length(current_sequence), ncol=3))
  colnames(vgc_df) <- c("N", "V", "V1")
  
  current_vector <- tolower(get(paste("centurydta", j, "_random", sep="", collapse=""))$V2)
  
  for(i in 1:length(current_sequence)) {
    vgc_df$N[i] <- current_sequence[i]
    vgc_df$V[i] <- length(levels(factor(current_vector[1:current_sequence[i]])))
    vgc_df$V1[i] <- length(as.data.frame(table(current_vector[1:current_sequence[i]]))[which(as.data.frame(table(current_vector[1:current_sequence[i]]))$Freq==1),1])
    
  }
  
  write.table(vgc_df, file=paste("dtacentury_random", j, ".emp.vgc", sep="", collapse=""),
              row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
  
  
}


# create spc files for non-randomized version
for(j in 1:3) {
  spec_df <- as.data.frame(table(tolower(get(paste("centurydta", j, 
                                                   sep="", collapse=""))$V2)))
  spec_df <- as.data.frame(table(spec_df$Freq))
  colnames(spec_df) <- c("m", "Vm")
  
  write.table(spec_df, file=paste("dtacentury", j, ".spc", sep="", collapse=""),
              row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
}

# create spc files for randomized version
for(j in 1:3) {
  spec_df <- as.data.frame(table(tolower(get(paste("centurydta", j, "_random", sep="", collapse=""))$V2)))
  spec_df <- as.data.frame(table(spec_df$Freq))
  colnames(spec_df) <- c("m", "Vm")
  
  write.table(spec_df, file=paste("dtacentury_random", j, ".spc", sep="", collapse=""),
              row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
}


# read in files
century1.emp.vgc <- read.vgc("dtacentury1.emp.vgc")
century2.emp.vgc <- read.vgc("dtacentury2.emp.vgc")
century3.emp.vgc <- read.vgc("dtacentury3.emp.vgc")

century1.spc <- read.spc("dtacentury_random1.spc")
century2.spc <- read.spc("dtacentury_random2.spc")
century3.spc <- read.spc("dtacentury_random3.spc")



# Zipf-Mandelbrot model
fzm_all1 <- lnre("fzm", century1.spc, m.max=1)
fzm_all2 <- lnre("fzm", century2.spc, m.max=1)
fzm_all3 <- lnre("fzm", century3.spc, m.max=1)



# extrapolated vgcs
century_all1.vgc <- lnre.vgc(fzm_all1, 1:5000000, m.max = 1)
century_all2.vgc <- lnre.vgc(fzm_all2, 1:5000000, m.max = 1)
century_all3.vgc <- lnre.vgc(fzm_all3, 1:5000000, m.max = 1)


# hapax growth curve
plot(century_all1.vgc$N, century_all1.vgc$V1, type="l", col="blue", lwd=2, main="Hapax Growth Curve",
     ylab="Extrapolated number of hapaxes", xlab="N")
points(century1.emp.vgc$N, century1.emp.vgc$V1, pch=10, col="blue", cex=0.5)
lines(century_all2.vgc$N, century_all2.vgc$V1, type="l", col="red", lwd=2)
points(century2.emp.vgc$N, century2.emp.vgc$V1, pch=12, col="red", cex=0.5)
lines(century_all3.vgc$N, century_all3.vgc$V1, type="l", col=rgb(red=0.1, green=1.0, blue=0.1, alpha=0.2), lwd=2)
points(century3.emp.vgc$N, century3.emp.vgc$V1, pch=15, col=rgb(red=0.1, green=1.0, blue=0.1, alpha=0.2), cex=0.5)
legend("topleft", inset=c(0.01,0.01), pch=c(10,12,15),
       col=c("blue", "red", "green"), lty=1, lwd=2,
       legend=c("17th century", "18th century", "19th century"),
       cex=0.6)



# potential productivity
plot(century_all1.vgc$N, (century_all1.vgc$V1/century_all1.vgc$N), type="l", col="grey75", cex=0.5, xlim=c(0,100000))
lines(century_all2.vgc$N, (century_all2.vgc$V1/century_all2.vgc$N), type="l", col="grey55", cex=0.5)
lines(century_all3.vgc$N, century_all3.vgc$V1/century_all3.vgc$N, type="l", col="grey35", cex=0.5)

# extrapolated productivity at N=5,000,000
png("extrapolated_prod_col.png", width=5, height=5, un="in", res=300)
plot(1:3, c(century_all1.vgc$V1[5000000] / century_all1.vgc$N[5000000],
       century_all2.vgc$V1[5000000] / century_all2.vgc$N[5000000],
       century_all3.vgc$V1[5000000] / century_all3.vgc$N[5000000]),
     type="b", ylab="V1/N (N=5,000,000)", xlab="Century", xaxt="n",
     main="Extrapolated Potential Productivity", pch=18, lwd=2, col="darkblue")
axis(1, at=c(1:3), labels=c("17th century", "18th century", "19th century"), cex.axis=.7)
dev.off()


############################################
# Zipf-Mandelbrot model with bootstrapping #
############################################

# check extrapolation quality via bootstrapping -------------------------------

fzm_boot_list01 <- fzm_boot_list02 <- fzm_boot_list03 <- list() # lists for storing the results 
                                                                # (one for each century)


for(current_boot in 1:100) {
  # preparing emp.vgc for bootsrapped version
  
  centurydta1_boot <- centurydta1[sample(1:nrow(centurydta1), 100000),]
  centurydta2_boot <- centurydta2[sample(1:nrow(centurydta2), 100000),]
  centurydta3_boot <- centurydta3[sample(1:nrow(centurydta3), 100000),]
  
  for(j in 1:3) {
    current_sequence <- seq(1000, nrow(get(paste("centurydta", j, "_boot", sep="", collapse=""))), 1000)
    vgc_df <- as.data.frame(matrix(nrow=length(current_sequence), ncol=3))
    colnames(vgc_df) <- c("N", "V", "V1")
    
    current_vector <- tolower(get(paste("centurydta", j, sep="", collapse=""))$V2)
    
    for(i in 1:length(current_sequence)) {
      vgc_df$N[i] <- current_sequence[i]
      vgc_df$V[i] <- length(levels(factor(current_vector[1:current_sequence[i]])))
      vgc_df$V1[i] <- length(as.data.frame(table(current_vector[1:current_sequence[i]]))[which(as.data.frame(table(current_vector[1:current_sequence[i]]))$Freq==1),1])
      
    }
    
    write.table(vgc_df, file=paste("dtacentury", j, "_boot.emp.vgc", sep="", collapse=""),
                row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
    
    
  }
  
  
  # create spc files for bootstrapped version
  for(j in 1:3) {
    spec_df <- as.data.frame(table(tolower(get(paste("centurydta", j, "_boot", sep="", collapse=""))$V2)))
    spec_df <- as.data.frame(table(spec_df$Freq))
    colnames(spec_df) <- c("m", "Vm")
    
    write.table(spec_df, file=paste("dtacentury_random", j, "_boot.spc", sep="", collapse=""),
                row.names=F, quote=F, sep="\t", fileEncoding = "UTF-8")
  }
  
  
  ###read in files
  century1.boot.emp.vgc <- read.vgc("dtacentury1_boot.emp.vgc")
  century2.boot.emp.vgc <- read.vgc("dtacentury2_boot.emp.vgc")
  century3.boot.emp.vgc <- read.vgc("dtacentury3_boot.emp.vgc")
  
  century1.boot.spc <- read.spc("dtacentury_random1_boot.spc")
  century2.boot.spc <- read.spc("dtacentury_random2_boot.spc")
  century3.boot.spc <- read.spc("dtacentury_random3_boot.spc")
  
  ###Zipf-Mandelbrot model for bootstrapping
  fzm_boot1 <- lnre("fzm", century1.boot.spc, m.max=1)
  fzm_boot2 <- lnre("fzm", century2.boot.spc, m.max=1)
  fzm_boot3 <- lnre("fzm", century3.boot.spc, m.max=1)
  
  
  ###extrapolated vgcs
  fzm_boot_list01[[current_boot]] <- lnre.vgc(fzm_boot1, 1:500000, m.max = 1)
  fzm_boot_list02[[current_boot]] <- lnre.vgc(fzm_boot2, 1:500000, m.max = 1)
  fzm_boot_list03[[current_boot]] <- lnre.vgc(fzm_boot3, 1:500000, m.max = 1)
 
  print(current_boot)
   
}



png("extrapolated_hgc_boot_and_prod_bw.png", height=5, width=10, un="in", res=300)
par(mfrow=c(1,2))
plot(fzm_boot_list01[[1]]$N, fzm_boot_list01[[1]]$V1, type="l", col=rgb(red = 0, green = 0.3, blue = 0.3, alpha = 0.3), ylab="V1", xlab="N", main="Hapax Growth Curves\n(100 bootstrapping iterations)", cex.main=0.7)
sapply(2:100, function(i) lines(fzm_boot_list01[[i]]$N, fzm_boot_list01[[i]]$V1, type="l", col=rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 0.5)))
lines(seq(1, 500000, 10000), sapply(seq(1, 500000, 10000), function(k) mean(sapply(1:100, function(i) fzm_boot_list01[[i]]$V1[k]))), col="black", 
      type="l", lwd=2)

lines(fzm_boot_list02[[1]]$N, fzm_boot_list02[[1]]$V1, type="l", col=rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5))
sapply(2:100, function(i) lines(fzm_boot_list02[[i]]$N, fzm_boot_list02[[i]]$V1, type="l", col=rgb(red=0.5, green=0.5, blue=0.5, alpha=0.5)))
lines(seq(1, 500000, 10000), sapply(seq(1, 500000, 10000), function(k) mean(sapply(1:100, function(i) fzm_boot_list02[[i]]$V1[k]))), col="black", 
      type="l", lwd=3, lty=2)

lines(fzm_boot_list03[[1]]$N, fzm_boot_list03[[1]]$V1, type="l", col=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.5))
sapply(2:100, function(i) lines(fzm_boot_list03[[i]]$N, fzm_boot_list03[[i]]$V1, type="l", col=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.5)))
lines(seq(1, 500000, 10000), sapply(seq(1, 500000, 10000), function(k) mean(sapply(1:100, function(i) fzm_boot_list03[[i]]$V1[k]))), col="black", 
      type="l", lwd=3, lty=3)

legend("topleft", inset=c(0.01,0.01), lty=c(1,2,3), col=c(rgb(red=0.3, green=0.3, blue=0.3), 
                                                          rgb(red=0.5, green=0.5, blue=0.5), 
                                                          rgb(red=0.5, green=0.5, blue=0.5)), 
       legend=c("17th century", "18th century", "19th century"), cex=0.6)

boxplot(matrix(c(sapply(1:100, function(i) fzm_boot_list01[[i]]$V1[500000] / fzm_boot_list01[[i]]$N[500000]),
                    sapply(1:100, function(i) fzm_boot_list02[[i]]$V1[500000] / fzm_boot_list02[[i]]$N[500000]),
                    sapply(1:100, function(i) fzm_boot_list03[[i]]$V1[500000] / fzm_boot_list03[[i]]$N[500000])), 
                  ncol=3), notch=T,
        xaxt="n", main="Extrapolated Potential Productivity", ylab="Potential productivity (N=500,000)", xlab="Century")
axis(1, at=c(1:3), labels=c("17th century", "18th century", "19th century"), cex.axis=0.7)
dev.off()
par(mfrow=c(1,1))


