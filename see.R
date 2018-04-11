args <- commandArgs(TRUE)
srcFile <- args[1]

d <- read.table(srcFile,header=T)

c=d$par1
I=d$t
plot(c~I,ann=FALSE,type="n",xlim=c(min(I),max(I)),ylim=c(min(c),max(c)))
lines(c~I,lwd=2)
title("par1",xlab="X",ylab="Y")


