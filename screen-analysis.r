bcscreen <- function(y,dmsocols=grep("DMSO",names(y))) {
	outlist <- list("reads" = y)
	wnorm <- apply(y[,-1],2,function(x) (x/sum(x))*10000000)
	outlist[["norm"]] = wnorm
	wlog <- log2(wnorm)
	outlist[["log"]] = wlog
	wmedc <- apply(wlog,2,function(x) x-median(x))
	outlist[["medc"]] = wmedc
	wrescaled <- data.frame(apply(wmedc,2,function(x) x/mad(x)))
	wrescaled$dmsomed <- apply(wrescaled[,dmsocols],1,median)
	outlist[["rescaled"]] = wrescaled
	wde <- apply(wrescaled,2,function(x,y) x-y,y=wrescaled$dmsomed)
	outlist[["de"]] = wde
	outlist
}

degraphs <- function(x,y="degraphs.pdf") {
pdf(y)
wde <- x[["de"]]
for(i in 1:ncol(wde)) { plot(wde[,i][order(wde[,i])],pch=20,main=dimnames(wde)[[2]][i],ylab="Drug Effect")} 
dev.off()
}
