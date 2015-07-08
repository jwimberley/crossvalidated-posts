# Settings
numGroups <- 20
filename <- "toystudy.csv"
#filename <- "realstudy.csv"

# Load tagdata
tagdata <- read.csv(filename)
if (dim(tagdata)[[2]] == 2)
    tagdata$W <- rep(1,nrow(tagdata))
colnames(tagdata) <- c("Tag","DLL","W")
print(summary(tagdata))

# Do fit
library(rms)
fit <- lrm(Tag ~ DLL, weights = W, data = tagdata, x=TRUE, y=TRUE)
print(fit)
print(latex(residuals(fit,"gof")))

# Bin data, recognizing weights
library(Hmisc)
tagdata$group <- as.numeric(cut2(tagdata$DLL,g=numGroups))
Tag <- vector("numeric",numGroups)
TagError <- vector("numeric",numGroups)
DLL <- vector("numeric",numGroups)
DLLError <- vector("numeric",numGroups)
for (i in 1:numGroups) {
    Tag[i] <- weighted.mean(tagdata[tagdata$group==i,]$Tag,tagdata[tagdata$group==i,]$W)
    DLL[i] <- weighted.mean(tagdata[tagdata$group==i,]$DLL,tagdata[tagdata$group==i,]$W)
    N <- length(tagdata[tagdata$group==i,]$W)
    TagError[i] <- sqrt(weighted.mean((tagdata[tagdata$group==i,]$Tag-Tag[i])^2)/(N-1))
    DLLError[i] <- sqrt(weighted.mean((tagdata[tagdata$group==i,]$DLL-DLL[i])^2)/(N-1))
}
tagdata.binned <- data.frame(Tag,TagError,DLL,DLLError)
print(summary(tagdata.binned))

# Make plots
png(filename="plotWithWeights.png",width=800,height=400)
with(tagdata.binned,plot(DLL,Tag,ylim=range(c(Tag-TagError,Tag+TagError)),xlim=range(c(DLL-DLLError,DLL+DLLError)),pch=19))
with(tagdata.binned,arrows(DLL,Tag-TagError,DLL,Tag+TagError,length=0.05,angle=90,code=3))
with(tagdata.binned,arrows(DLL-DLLError,Tag,DLL+DLLError,Tag,length=0.05,angle=180,code=3))
lines(tagdata.binned$DLL,predict(fit,newdata=tagdata.binned,type="fitted"))
dev.off()

# Bin data ignoring weights
library(Hmisc)
tagdata$group <- as.numeric(cut2(tagdata$DLL,g=numGroups))
Tag <- vector("numeric",numGroups)
TagError <- vector("numeric",numGroups)
DLL <- vector("numeric",numGroups)
DLLError <- vector("numeric",numGroups)
for (i in 1:numGroups) {
    Tag[i] <- mean(tagdata[tagdata$group==i,]$Tag)
    DLL[i] <- mean(tagdata[tagdata$group==i,]$DLL)
    N <- length(tagdata[tagdata$group==i,]$W)
    TagError[i] <- sqrt(var(tagdata[tagdata$group==i,]$Tag)/(N-1))
    DLLError[i] <- sqrt(var(tagdata[tagdata$group==i,]$DLL)/(N-1))
}
tagdata.binned.unw <- data.frame(Tag,TagError,DLL,DLLError)
print(summary(tagdata.binned.unw))

# Make plots
png(filename="plotWithoutWeights.png",width=800,height=400)
with(tagdata.binned.unw,plot(DLL,Tag,ylim=range(c(Tag-TagError,Tag+TagError)),xlim=range(c(DLL-DLLError,DLL+DLLError)),pch=19))
with(tagdata.binned.unw,arrows(DLL,Tag-TagError,DLL,Tag+TagError,length=0.05,angle=90,code=3))
with(tagdata.binned.unw,arrows(DLL-DLLError,Tag,DLL+DLLError,Tag,length=0.05,angle=180,code=3))
lines(tagdata.binned.unw$DLL,predict(fit,newdata=tagdata.binned.unw,type="fitted"))
dev.off()

