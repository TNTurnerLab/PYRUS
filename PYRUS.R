#!/usr/bin/Rscript 
library(optparse)
library(dplyr)
library(seqminer)

option_list<-list(make_option(c("-f", "--file"), default=NULL, type="character",help='Input an Inital Bed File'),
                  make_option(c("-d", "--dir"), default=NULL, type="character",help='Input Directory of Bed Files'), 
                  make_option(c("-p", "--dirpattern"), default=".bed.gz", type="character",help='Pattern of files in Directory, defualt .bed'),
                  make_option(c("-a", "--annotation"), default=NULL, type="character",help='Add an Annotation File '),
                  make_option(c("-c", "--cordfile"), default=NULL, type="character",help='Input a Chr Coordinates Bed File'), 
                  make_option(c("-l", "--lineColor"), default="blue,red", type="character",help='Input a Color(s) ex: blue,red'),
                  make_option(c("-t", "--plotTogether"), default=NULL, type="character",help='Graph Together. ex: BARX1,BARX1-DT'),
                  make_option(c("-r", "--HighRes"), default=NULL, type="character",help='Save Plot(s) as PNG Image, Defualt is off'),
                  make_option(c("-b", "--box"), default=NULL, type="character",help='Plot Box-plot Window'),
                  make_option(c("-s", "--singleplots"), default=NULL, type="character",help='plots each file in dirpath individually'),
                  make_option(c("-y", "--rename"), default=NULL, type="character",help='Rename file -f on plot'),
                  make_option(c("-n", "--annInput"), default=NULL, type="character",help='Input For Annotation (fill,boarder,name), ex : "red,blue,Exons'),
                  make_option(c("-o", "--outfile"), default=NULL, type="character",help='Outputs The Chromosome Plots Into One PDF, Defualt Is Plots Each Cordinate To Its Own PDF File'))
opt<-parse_args(OptionParser(option_list=option_list))
options(scipen = 999)
outfile=opt$outfile
colorofann<-"red"
bordcolofann<-"black"
nameofann<-"Exons"

annlen<-(-0.12)
colofline<-strsplit(toString(opt$lineColor), ',', fixed=TRUE)[[1]][1]
getcol<-strsplit(toString(opt$lineColor), ',', fixed=TRUE)[[1]][2]

if(!is.null(opt$cordfile)){
  coord <- read.delim(gzfile(opt$cordfile), header=F)
  if(is.null(opt$plotTogether)){
    COORD_INPUT<- list(paste0(coord[,4],",",coord[,1], ":", coord[,2], "-", coord[,3], sep=""))
  }
  if(!is.null(opt$plotTogether)){
    getnameInitial<-strsplit(toupper(opt$plotTogether), split=',', fixed=TRUE)[[1]][1]
    getnewlistname<-strsplit(toupper(opt$plotTogether), split=',', fixed=TRUE)[[1]][-1]
    coord1<- coord %>% filter(coord[,4] %in% getnameInitial)
    COORD_INPUT<- list(paste0(coord1[,4],",",coord1[,1], ":", coord1[,2], "-", coord1[,3], sep=""))
    NEWLISTPlot <- tryCatch({ 
      x<- lapply(getnewlistname,function(x){
        coord<- coord %>% filter(coord[,4] %in% x)
        NEWLISTPlot<- list(paste0(coord[,4],",",coord[,1], ":", coord[,2], "-", coord[,3], sep=""))
        NEWLISTPlot<- unlist(NEWLISTPlot)
        return(NEWLISTPlot)
      })},
      error = function(e){
        y<- lapply(getnewlistname,function(x){
          stop("Missing Input: Chromosome Position Bed File")}) })}
}

filename<-toString(rbind(strtrim(basename(opt$file),15)))
if(!is.null(opt$rename)){
  filename=toString(rbind(strtrim(opt$rename,15)))}

if(is.null(opt$highRes)){
    outfilename=opt$outfile
    pdf(outfilename)}

old.par <- par(no.readonly = TRUE)

if(!is.null(opt$annInput)){
  annsplit<-strsplit(opt$annInput, ',', fixed=TRUE)
  colorofann<-annsplit[[1]][1]
  bordcolofann<-annsplit[[1]][2]
  nameofann<-toString(rbind(strtrim(annsplit[[1]][3],15)))
  annlen<-(-0.2)
}
if(!is.null(opt$dir)){
  pattern<-toString(opt$dirpattern)
  pattern1<-paste0(pattern,"$")
  directory<-list.files(path=toString(opt$dir) ,pattern=pattern1, full.names=TRUE)
  if(length(directory) >= 10){
    opt$box="yes"}
}

func<-function(t){
  lapply(COORD_INPUT[[1]], function(x){
    par(mar=c(5.1, 4.1, 4.1, 6.1))
    rangeofit<-toString(strsplit(x, ',', fixed=TRUE)[[1]][2])
    chr<-strsplit(rangeofit, ':', fixed=TRUE)[[1]][1]
    nameofit<-strsplit(x, ',', fixed=TRUE)[[1]][1]
    initialfile<-file.path(opt$file)
    fr<-tabix.read.table(initialfile,rangeofit,col.names = TRUE,stringsAsFactors = FALSE)
    wait<-1
    if(length(fr)>=1){
      par(old.par)
      if(is.null(opt$outfile)){
        if(!is.null(opt$HighRes) && is.null(opt$singleplots) &&is.null(opt$plotTogether)){
          png(paste(nameofit, ".png", sep=""), width = 600, height = 600, units = "px", res=100) }
        if(is.null(opt$HighRes) && is.null(opt$singleplots) && is.null(opt$plotTogether) ){
          pdf(paste0(nameofit,".pdf")) } 
        if(!is.null(opt$HighRes) && !is.null(opt$singleplots) ){
            png(paste(nameofit, "_",basename(t),".png", sep=""), width = 600, height = 600, units = "px", res=100) }
        if(is.null(opt$HighRes) && !is.null(opt$singleplots) ){
            pdf(paste0(nameofit,"-",basename(t),".pdf")) }
        }
    initialmean<-mean(fr$V4)
    par(mar=c(5.1, 4.1, 4.1, 6.1))
    if(is.null(opt$plotTogether)){
      xaxisname<-paste(chr, "position")
      plot(as.numeric(fr$V2), y=as.numeric(fr$V4), col=("transparent"), xlab=xaxisname ,
           ylab="Estimated Copy Number", main=(bquote(paste(italic(.(nameofit))))), pch=19,
           ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1,las=1)
      abline(h=2, col="grey", lty=2)

      legend("topright", inset=c(-0.2,0.5), legend=filename,pch=23,col=colofline,lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
      abline(h=2.7, col="grey", lty=2)
      abline(h=2, col="grey", lty=1)
      abline(h=1.3, col="grey", lty=2)
      lines(as.numeric(fr$V2), y=as.numeric(fr$V4), type="l", col=colofline)
      old.par <- par(mar = c(0, 0, 0, 0))
    }
      if(!is.null(opt$plotTogether)){
          secondmean<-lapply(NEWLISTPlot, function(p){
          plottogethernameofit<-strsplit(p, ',', fixed=TRUE)[[1]][1]
          plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
          splitfromchrplot<-strsplit(plottogetherrange, ':', fixed=TRUE)
          chrplot<-strsplit(rangeofit, ':', fixed=TRUE)[[1]][1]
          plotlinetogether1<-tabix.read.table(initialfile,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE)
          
          if(length(plotlinetogether1)>=1 ){
            if(is.null(opt$outfile) && is.null(opt$singleplots)){
            pdf(paste0(nameofit,"_and_",plottogethernameofit,".pdf"))}
            xaxisname<-paste("positions",chr, "and" ,chrplot)
            secmean<-mean(plotlinetogether1$V4)
            plotlinetogether<-rbind(fr,plotlinetogether1)
            par(mar=c(5.1, 4.1, 4.1, 6.1))
            
             plot(as.numeric(plotlinetogether$V2), y=as.numeric(plotlinetogether$V4),  xlab=xaxisname,col="transparent",
                 ylab="Estimated Copy Number", main=(bquote(paste(italic(.(nameofit)),paste(" and "),italic(.(plottogethernameofit))))), pch=19,
                 ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1,las=1)
            legend("topright", inset=c(-0.2,0.5), legend=strtrim(nameofit,15),col=colofline,lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5) 
            legend("topright", inset=c(-0.2,0.45), legend=strtrim(plottogethernameofit,15),col=toString(getcol),lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
           return(secmean) }
          else{ plotlinetogether<-NULL }
        })
      }
    printline<-function(){
      abline(h=2.7, col="grey", lty=2)
      abline(h=2, col="grey", lty=1)
      abline(h=1.3, col="grey", lty=2)
      lines(as.numeric(fr$V2), y=as.numeric(fr$V4), type="l", col=colofline)
    }
        if(!is.null(opt$annotation)){
          par(mar=c(5.1, 4.1, 4.1, 6.1))
          annot<-file.path(opt$annotation)
          annota<-tabix.read.table(annot,rangeofit,col.names = TRUE,stringsAsFactors = FALSE)
          if(length(annota)>=1){
            rect(annota$V2, 0, annota$V3, 0.35, col= colorofann, border =bordcolofann)
          }
          if(!is.null(opt$plotTogether)){
            lapply(NEWLISTPlot, function(p){
              plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
              annopt<-tabix.read.table(annot,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
             
              if(length(annopt)>=1){
                rect(annopt$V2, 0, annopt$V3, 0.35, col= colorofann, border =bordcolofann) }
              else{ annopt<-NULL }
            })
          }
          legend("topright", inset=c(annlen,0.6), fill=colorofann,border=bordcolofann,legend=nameofann,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
        }
    if(!is.null(opt$dir)){
        if(!is.null(opt$singleplots)){ 
          par(mar=c(5.1, 4.1, 4.1, 6.1))
          dirpass<-tabix.read.table(t,rangeofit,col.names = TRUE,stringsAsFactors = FALSE) 
            if(length(dirpass)>=1){
              getmeans<-mean(dirpass$V4)
              lines(as.numeric(dirpass$V2), y=as.numeric(dirpass$V4), type="l")
              dirnames<-toString(rbind(strtrim(basename(t),15)))
              legend("topright", inset=c(as.numeric(-0.2),0.55), legend=dirnames,col="grey24",lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
              if(!is.null(opt$plotTogether)){
                  meansget<-lapply(NEWLISTPlot, function(p){
                    plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
                    listdir1<-tabix.read.table(t,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
                    plotlinetogether1<-tabix.read.table(initialfile,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE)
                    if(length(plotlinetogether1)>=1){
                      lines(as.numeric(plotlinetogether1$V2), y=as.numeric(plotlinetogether1$V4), type="l", col=getcol)
                    }
                    if(length(listdir1)>=1){
                      par(new=TRUE)
                      lines(as.numeric(listdir1$V2), y=as.numeric(listdir1$V4), type="l",col="grey24")
                      getmeans1<-cbind(listdir1$V4)
                      return(getmeans1) }
                    else{listdir1<-NULL }
                })
                  printlines<-printline()
              }
              if(!is.null(opt$box)){
                if(!is.null(opt$plotTogether)){
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(fig=c(0.47,0.9,0.6,1), new=TRUE)
                  
                  boxplot(getmeans, col="white",ylim=c(0,6),axes=FALSE,xpd = TRUE)
                  points(initialmean, col=colofline, pch=23)
                  par(fig=c(0.5,1,0.6,1), new=TRUE)
                  axis(side = 4,las=1,cex.axis=0.5)
                  abline(h=2, col="grey", lty=1)
                  abline(h=2.7, col="grey", lty=2)
                  abline(h=1.3, col="grey", lty=2)
                  
                  par(fig=c(0.57,1,0.6,1), new=TRUE)
                  boxplot(meansget, col="white",ylim=c(0,6),axes=FALSE,xpd = TRUE)
                
                if(!is.null(opt$plotTogether)){
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(fig=c(0.57,1,0.6,1), new=TRUE)
                  lapply(secondmean,function(m){
                    points(m, col=getcol, pch=23)})  
                  }
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(mfrow=c(2,1))
                  par(old.par) }
              
              if(is.null(opt$plotTogether)){
                par(fig=c(0.5,1,0.6,1), new=TRUE)
                boxplot(getmeans, col="white",ylim=c(0,6),axes=FALSE,xpd = TRUE)
                axis(side = 4,las=1,cex.axis=0.5)
                points(initialmean, col=colofline, pch=23)
                abline(h=2, col="grey", lty=1)
                abline(h=2.7, col="grey", lty=2)
                abline(h=1.3, col="grey", lty=2)
                par(old.par) }
              }
            }
        }
        if(is.null(opt$singleplots)){
        x<-lapply(directory,function(y){
          if(y != initialfile){
            listdir<-tabix.read.table(y,rangeofit,col.names = TRUE,stringsAsFactors = FALSE) 
            if(length(listdir)>=1){
              getmeans<-mean(listdir$V4)
              lines(as.numeric(listdir$V2), y=as.numeric(listdir$V4), type="l")
              printlines<-printline()
              dirnames<- "Files From Directory"
              par(new=TRUE)
              legend("topright", inset=c(as.numeric(-0.225),0.55), legend=dirnames,col="grey24",lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
              if(!is.null(opt$plotTogether)){
                lapply(NEWLISTPlot, function(p){
                  plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
                  listdir1<-tabix.read.table(y,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
                  plotlinetogether1<-tabix.read.table(initialfile,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE)
                  if(length(plotlinetogether1)>=1){
                    lines(as.numeric(plotlinetogether1$V2), y=as.numeric(plotlinetogether1$V4), type="l", col=getcol,lwd=1.0)
                  }
                  if(length(listdir1)>=1){
                    getmeans<-cbind(listdir1$V4)
                    par(new=TRUE)
                    lines(as.numeric(listdir1$V2), y=as.numeric(listdir1$V4), type="l",col="grey24")
                    return(getmeans) }
                  else{listdir1<-NULL}
                  })
              }
              return(getmeans) }
            else{getmeans<-NULL }
          }
        })
      if(!is.null(opt$box)){
          if(is.null(opt$plotTogether)){
          x<-do.call(rbind,x) 
          }
          par(fig=c(0.5,1,0.6,1), new=TRUE)
          boxplot(x, col="white",ylim=c(0,6),axes=FALSE,xpd = TRUE)
          axis(side = 4,las=1,cex.axis=0.5)
          points(initialmean, col=colofline, pch=23)
        if(is.null(opt$singleplots) && !is.null(opt$plotTogether)){
          lapply(secondmean,function(m){
          points(3,m, col=getcol, pch=23)}) }
          abline(h=2, col="grey", lty=1)
          abline(h=2.7, col="grey", lty=2)
          abline(h=1.3, col="grey", lty=2)
          par(old.par) }
        }
      }
    }
  else{fir<-NULL}
  })
}
  if(!is.null(opt$singleplots)){ 
  lapply(directory,func) }
  if(is.null(opt$singleplots)){ 
    lapply(1,func)}
dev.off()