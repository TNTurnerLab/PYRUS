#!/usr/bin/Rscript 
library(optparse)
library(dplyr)
library(seqminer)
library(data.table)

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
                  make_option(c("-g", "--regions"), default="1.3,2.7", type="character",help='Plot only regions above 2.7 CNV or below 1.3 CNV'),
                  make_option(c("-v", "--minmax"), default=NULL, type="character",help='print between regions'),
                  make_option(c("-j", "--Trio"), default=NULL, type="character",help='Plot trio with files from parent:A(Male) and parent:B(Female) alongside the initial file -f'),
                  make_option(c("-k", "--subtractbynum"), default=1, type="character",help='Subtract values given in -g flag by input for only chrY and chrX if karyotype is M, Default is 1 '),
                  make_option(c("-u", "--topylim"), default=6, type="integer",help='change top -y lim value, defualt is 6'),
                  make_option(c("-x", "--sex"), default=F, type="character",help='M or F sample, if male sample, -g flag max and mins will be subtracted by 1'),
                  make_option(c("-n", "--annInput"), default=NULL, type="character",help='Input For Annotation (fill,boarder,name), ex : "red,blue,Exons'),
                  make_option(c("-o", "--outfile"), default=NULL, type="character",help='Outputs The Chromosome Plots Into One PDF, Defualt Is Plots Each Cordinate To Its Own PDF File'))
opt<-parse_args(OptionParser(option_list=option_list))
options(scipen = 999)
outfile=opt$outfile
addname<-''
colorofann<-"red"
bordcolofann<-"black"
nameofann<-"Exons"
max<-(as.numeric(strsplit(toString(opt$regions), ',', fixed=TRUE)[[1]][2]))
min<-(as.numeric(strsplit(toString(opt$regions), ',', fixed=TRUE)[[1]][1]))
annlen<-(-0.12)
colofline<-strsplit(toString(opt$lineColor), ',', fixed=TRUE)[[1]][1]
getcol<-strsplit(toString(opt$lineColor), ',', fixed=TRUE)[[1]][2]

sampleFile = opt$file
if(is.null(opt$Trio)){
  inputdata<-sampleFile
 
}
if(!is.null(opt$Trio)){
  inputdata<-append(toString(opt$file),unlist(strsplit(toString(opt$Trio), ',', fixed=TRUE)))
  if(is.null(opt$outfile)){
    warning("PLEASE ADD (-o) OUTPUTFILE")
    stop()}
  if(!is.null(opt$plotTogether)){
    warning("REMOVE -t FLAG")
    stop()}
}
suppressWarnings({
if(inputdata==sampleFile){
infile <- as.data.frame(fread(opt$file), header=F)
infile <- infile[which(infile$V1 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")),]
chromCount <- sum(round(tapply(infile$V4, infile$V1, mean)))
chrX <- paste0(rep("X", round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrX"]), collapse="")
chrY <- paste0(rep("Y", round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrY"]), collapse="")
if((as.integer(round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrY"])) >=1 ){
  opt$sex<-'M'}
}})
if(!is.null(opt$cordfile)){
  coord <- read.delim(gzfile(opt$cordfile), header=F)
  if(is.null(opt$plotTogether) ){
    COORD_INPUT<- list(paste0(coord[,4],",",coord[,1], ":", coord[,2], "-", coord[,3], sep=""))
  }
  if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
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

if(is.null(opt$highRes)){
    outfilename=opt$outfile
    pdf(outfilename)}
if(!is.null(opt$annInput)){
  annsplit<-strsplit(opt$annInput, ',', fixed=TRUE)
  colorofann<-annsplit[[1]][1]
  bordcolofann<-annsplit[[1]][2]
  nameofann<-toString(rbind(strtrim(annsplit[[1]][3],15)))
  annlen<-(-0.2)}

if(!is.null(opt$dir)){
  pattern<-toString(opt$dirpattern)
  pattern1<-paste0(pattern,"$")
  directory<-list.files(path=toString(opt$dir) ,pattern=pattern1, full.names=TRUE)
  if(length(directory) >= 10){
    opt$box="yes"}
}
if(is.null(opt$Trio)){
  old.par <- par(no.readonly = TRUE)}
if(!is.null(opt$Trio)){
layout(matrix(c(2,3,1,1), 2, 2, byrow = TRUE))
  old.par<-NULL}

func<-function(t){
  holding<-lapply(COORD_INPUT[[1]], function(x){

    par(mar=c(5.1, 4.1, 4.1, 6.1))
    rangeofit<-toString(strsplit(x, ',', fixed=TRUE)[[1]][2])
    chr<-strsplit(rangeofit, ':', fixed=TRUE)[[1]][1]
    stsp<-strsplit(rangeofit, ':', fixed=TRUE)[[1]][2]
    nameofit<-strsplit(x, ',', fixed=TRUE)[[1]][1]
    startcord<-strsplit(stsp, '-', fixed=TRUE)[[1]][1]
    stopcord<-strsplit(stsp, '-', fixed=TRUE)[[1]][2]
    filehold<-lapply(inputdata, function(k){
      inputfile<-k
      filename<-toString(rbind(strtrim(strsplit(toString(basename(inputfile)), '.bed', fixed=TRUE)[[1]][1],15)))
     
    if(!is.null(opt$Trio)){
      if(inputfile==opt$file){
        triomember<- "Child"
      ledgendvalue<-(-0.15)}
      if(inputfile!=opt$file){
        triomember<-"Parent"
        ledgendvalue<-(-0.5)
        annlen<-(-0.5)}
      infile <- as.data.frame(fread(inputfile), header=F)
      infile <- infile[which(infile$V1 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")),]
      
      chromCount1 <- sum(round(tapply(infile$V4, infile$V1, mean)))
      chrX1 <- paste0(rep("X", round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrX"]), collapse="")
      chrY1 <- paste0(rep("Y", round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrY"]), collapse="")
      triomember<-paste0(triomember,"(",chromCount1,chrX1,chrY1,")")
      addname<- paste0(triomember,":")
      if((as.integer(round(tapply(infile$V4, infile$V1, mean))[names(round(tapply(infile$V4, infile$V1, mean))) == "chrY"])) >=1 ){
        opt$sex<-'M'
      }
    }
      
    initialfile<-file.path(inputfile)
    fr<-tabix.read.table(initialfile,rangeofit,col.names = TRUE,stringsAsFactors = FALSE)
    
    if(length(fr)>=1){
      par(old.par)
      initialmean<-mean(fr$V4)
      
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
    if(chr == "chrY"){
      min<-(min-as.numeric(opt$subtractbynum))
      max<-(max-as.numeric(opt$subtractbynum))}
    if(opt$sex == "M" && chr == "chrX"){
      max<-(max-as.numeric(opt$subtractbynum))
      min<-(min-as.numeric(opt$subtractbynum))}
      
    if(is.null(opt$minmax) || (!is.null(opt$minmax) && (as.numeric(initialmean)<=min || as.numeric(initialmean)>=max))){
    par(mar=c(5.1, 4.1, 4.1, 6.1))
    bgcolor<-function(){
      par(fig=c(0.5,1,0.6,1),new=TRUE)
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white", border="transparent")}
    
    if(is.null(opt$plotTogether) ){
      xaxisname<-paste(chr, "position")
      plot(as.numeric(fr$V2), y=as.numeric(fr$V4), col=("transparent"), xlab=xaxisname ,
           ylab="Estimated Copy Number", main=(bquote(paste(paste(.(addname)),italic(.(nameofit))))), pch=23,
           ylim= c(0,opt$topylim), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1,las=1)
      abline(h=2, col="grey", lty=2)
      
      if(is.null(opt$Trio)){
      legend("topright", inset=c(-0.2,0.5), legend=filename,pch=23,col=colofline,lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)}
      
      abline(h=2.7, col="grey", lty=2)
      abline(h=2, col="grey", lty=1)
      abline(h=1.3, col="grey", lty=2)
      lines(as.numeric(fr$V2), y=as.numeric(fr$V4), type="l", col=colofline)}
    
      if(!is.null(opt$plotTogether) && is.null(opt$Trio) ){
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
                 ylim= c(0,opt$topylim), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1,las=1)
            
            if(is.null(opt$Trio)){
            legend("topright", inset=c(-0.2,0.5), legend=strtrim(nameofit,15),col=colofline,lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5) 
            legend("topright", inset=c(-0.2,0.45), legend=strtrim(plottogethernameofit,15),col=toString(getcol),lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)}
            lines(as.numeric(plotlinetogether1$V2), y=as.numeric(plotlinetogether1$V4), type="l", col=getcol)
            lines(as.numeric(fr$V2), y=as.numeric(fr$V4), type="l", col=colofline)
           return(secmean) }
          else{ plotlinetogether<-NULL }
        })
      }
    printline<-function(){
      abline(h=2.7, col="grey", lty=2)
      abline(h=2, col="grey", lty=1)
      abline(h=1.3, col="grey", lty=2)
      lines(as.numeric(fr$V2), y=as.numeric(fr$V4), type="l", col=colofline)}
    
        if(!is.null(opt$annotation)){
          par(mar=c(5.1, 4.1, 4.1, 6.1))
          annot<-file.path(opt$annotation)
          annota<-tabix.read.table(annot,rangeofit,col.names = TRUE,stringsAsFactors = FALSE)
          
          if(length(annota)>=1){
            rect(annota$V2, 0, annota$V3, 0.35, col= colorofann, border =bordcolofann)}
          
          if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
            lapply(NEWLISTPlot, function(p){
              plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
              annopt<-tabix.read.table(annot,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
             
              if(length(annopt)>=1){
                rect(annopt$V2, 0, annopt$V3, 0.35, col= colorofann, border = bordcolofann) }
              else{ annopt<-NULL }
            })}
          
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
              legend("topright", inset=c(ledgendvalue,0.55), legend=dirnames,col="grey24",lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
              if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
                  meansget<-lapply(NEWLISTPlot, function(p){
                    plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
                    listdir1<-tabix.read.table(t,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
                    plotlinetogether1<-tabix.read.table(initialfile,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE)
                    
                    if(length(plotlinetogether1)>=1){
                      lines(as.numeric(plotlinetogether1$V2), y=as.numeric(plotlinetogether1$V4), type="l", col=getcol)}
                    if(length(listdir1)>=1){
                      par(new=TRUE)
                      lines(as.numeric(listdir1$V2), y=as.numeric(listdir1$V4), type="l",col="grey24")
                      getmeans1<-cbind(listdir1$V4)
                      return(getmeans1) }
                    else{listdir1<-NULL }
                })
                  printlines<-printline()}
              
              if(!is.null(opt$box) && is.null(opt$Trio)){
                colobgit<-bgcolor()
                if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(fig=c(0.47,0.9,0.6,1), new=TRUE)
                  boxplot(getmeans, col="white",ylim=c(0,opt$topylim),axes=FALSE,xpd = TRUE)
                  points(initialmean, col=colofline, pch=23)
                  par(fig=c(0.5,1,0.6,1), new=TRUE)
                  par(bg = "blue")
                  axis(side = 4,las=1,cex.axis=0.5)
                  abline(h=2, col="grey", lty=1)
                  abline(h=2.7, col="grey", lty=2)
                  abline(h=1.3, col="grey", lty=2)
                  
                  par(fig=c(0.57,1,0.6,1), new=TRUE)
                  boxplot(meansget, col="white",ylim=c(0,opt$topylim),axes=FALSE,xpd = TRUE)
                
                if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(fig=c(0.57,1,0.6,1), new=TRUE)
                  lapply(secondmean,function(m){
                  points(m, col=getcol, pch=23)})
                  }
                  par(mar=c(5.1, 4.1, 4.1, 6.1))
                  par(old.par) }
              
              if(is.null(opt$plotTogether)){
                par(fig=c(0.5,1,0.6,1), new=TRUE)
                boxplot(getmeans, col="white",ylim=c(0,opt$topylim),axes=FALSE,xpd = TRUE)
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
              getmeans<-listdir$V4
              lines(as.numeric(listdir$V2), y=as.numeric(listdir$V4), type="l")
              printlines<-printline()
              dirnames<- "Files From Directory"
              par(new=TRUE)
              if(is.null(opt$Trio)){legend("topright", inset=c(as.numeric(-0.225),0.55), legend=dirnames,col="grey24",lty=1,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)}
              if(!is.null(opt$plotTogether) && is.null(opt$Trio)){
                y<-lapply(NEWLISTPlot, function(p){
                  plottogetherrange<-strsplit(p, ',', fixed=TRUE)[[1]][2]
                  listdir1<-tabix.read.table(y,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE) 
                  plotlinetogether1<-tabix.read.table(initialfile,plottogetherrange,col.names = TRUE,stringsAsFactors = FALSE)
                  if(length(plotlinetogether1)>=1){
                    lines(as.numeric(plotlinetogether1$V2), y=as.numeric(plotlinetogether1$V4), type="l", col=getcol,lwd=1.0)
                  }
                  if(length(listdir1)>=1){
                    getmeans1<-listdir1$V4
                    n<-max(length(getmeans),length(getmeans1))
                    getmeans<-data.table(getmeans=getmeans[1:n],getmeans1=getmeans1[1:n])
                  
                    par(new=TRUE)
                    lines(as.numeric(listdir1$V2), y=as.numeric(listdir1$V4), type="l",col="grey24")
                    return(getmeans) }
                  else{listdir1<-NULL}
                  })
                getmeans<-do.call(rbind,y)
              }
              return(getmeans) }
            else{getmeans<-NULL }
          }
        })
        x<-do.call(rbind,x)
      if((!is.null(opt$box) && is.null(opt$Trio)) ){
          
          if(is.null(opt$plotTogether)){
          x<-mean(x)}
          else{ 
            x<-cbind(mean(na.omit(x$getmeans)),mean(na.omit(x$getmeans1)))}
          print(x)
          colobgit<-bgcolor()
          par(fig=c(0.5,1,0.6,1),new=TRUE)
          boxplot(x, col="white",ylim=c(0,opt$topylim),axes=FALSE,xpd = TRUE)
          axis(side = 4,las=1,cex.axis=0.5)
          points(initialmean, col=colofline, pch=23)
          abline(h=2, col="grey", lty=1)
          abline(h=2.7, col="grey", lty=2)
          abline(h=1.3, col="grey", lty=2)
        if(is.null(opt$singleplots) && !is.null(opt$plotTogether) && is.null(opt$Trio)){
          lapply(secondmean,function(m){
          points(2,m, col=getcol, pch=23)}) 
          }
      }
        if(!is.null(opt$Trio)){
          dirmean<-strtrim(toString(mean(na.omit(x))),5)
          printmeans<-paste("D:",dirmean, '\n', triomember,":",strtrim(toString(initialmean),5))
          legend("topright", inset=c(as.numeric(ledgendvalue),0.3), legend=printmeans,col="grey24",pch=23,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)
          legend("topright", inset=c(as.numeric(ledgendvalue),0.5), legend=filename,pch=23,col=colofline,xpd = TRUE,horiz = TRUE,bty = "n",cex=0.5)}
        }
    }
    else{fir<-NULL}
    }
    if((as.numeric(initialmean)<=min || as.numeric(initialmean)>=max) ){
      
      if(is.null(opt$Trio)){
      holdthething<-paste0(chr,"\t",startcord,"\t",stopcord,"\t",nameofit,"\t",initialmean)
      return(holdthething)}
      if(!is.null(opt$Trio)){
      holdthething<-paste0(chr,"\t",startcord,"\t",stopcord,"\t",nameofit,"\t",initialmean, "\t", triomember)
      return(holdthething)}
      }
    }
  else{fir<-NULL}
    })
   fileholdagain<-filehold[!sapply(filehold, is.null)]
    return(fileholdagain)
    })
  if(is.null(opt$plotTogether)){
    filenameisthis<-paste0(chromCount, chrX, chrY,"_",toString(basename(opt$file)),"_chromosomal_positions_outside_of_",min,"_and_",max,".txt")
    if(is.null(opt$Trio)){
      header<-paste0("chrom","\t","start","\t","stop","\t","name","\t","WTC_mean","\t")
      karyotype<-paste0("The karyotype in ", inputdata, " is: ", paste0(chromCount, chrX, chrY), "\n")
      header<-append(karyotype,header)}
    if(!is.null(opt$Trio)){
      header<-paste0("chrom","\t","start","\t","stop","\t","name","\t","WTC_mean","\t","Parent/Child")
      filenameisthis<-paste0('TRIO:',filenameisthis)}
    a.list<-append(header,unlist(holding[!sapply(holding, is.null)]))
    lapply(a.list, write, filenameisthis, append=TRUE, ncolumns=5)}
}
  if(!is.null(opt$singleplots)){ 
  lapply(directory,func) }
  if(is.null(opt$singleplots)){ 
    lapply(1,func)}
dev.off()
