#!/usr/bin/Rscript 
library(optparse)
library(base)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(stringi)
library(randomcoloR)
library(pdftools)
library(R.utils)
library(patchwork)
library(png)
library(grid)

#library(seqminer)
#library(future)
#plan(multisession)
#library("parallel")
#no_of_cores = detectCores()

start1 <- Sys.time() 
option_list<-list(make_option(c("-f", "--file"), default=NULL, type="character",help='Input an Inital Quick-mer2 Bed File'),
                    make_option(c("-d", "--dir"), default=NULL, type="character",help='Input Directory of Bed Files'), 
                    make_option(c("-p", "--dirpattern"), default=NULL, type="character",help='Pattern of files in Directory, defualt .bed'),
                    make_option(c("-a", "--annotation"), default=NULL, type="character",help='Input an Annotation'),
                    make_option(c("-c", "--cordfile"), default=NULL, type="character",help='Input a Chr Coordinates Bed File'), 
                    make_option(c("-l", "--lineColor"), default=NULL, type="character",help='Input a Color(s) ex: blue,red'),
                    make_option(c("-i", "--lineMult"), default=NULL, type="character",help='Input a Color for -f line when using -d flag ex: blue'),
                    make_option(c("-t", "--plotTogether"), default=NULL, type="character",help='Plot Chromosomes Together. ex: NGF,POGZ'),
                    make_option(c("-r", "--png"), default=NULL, type="character",help='Save Plot(s) as PNG Image'),
                    make_option(c("-b", "--box"), default=NULL, type="character",help='Plot Box-plot Window'),
                    make_option(c("-o", "--multifile"), default=NULL, type="character",help='Outputs The Chromosome Plots Into One PDF, Defualt Is Plots Each Cordinate To Its Own PDF File'))
opt<-parse_args(OptionParser(option_list=option_list))

if(!is.null(opt$annotation)){
    
    typeline<-function(msg="Enter text: "){
    
        if(!interactive()){
            cat(msg);
            exonisbed<-readLines("stdin",n=1);
        } 
        
        else{
            exonstartcol<-readLines("stdin",n=1);
            exonendcol<-readLines("stdin",n=1);
            exonalpha<-readLines("stdin",n=1);
            exoncolor<-readLines("stdin",n=1);
            exonborder<-readLines("stdin",n=1);
            exonledgend<-readLines("stdin",n=1);
        
        }

        return(exonisbed)
        return(exonstartcol)
        return(exonsendcol)
        return(exonalpha)
        return(exoncolor)
        return(exonborder)
        return(exonledgend)
    }
    print(" Building custom annotation track...")
    print("See Documentation for more info")
 
    exonisbed=tolower(toString(typeline("If input file a bed file (y/n):")))
  
    if(exonisbed %like% "y" || exonisbed %like% "yes"){

        exonchrcol=toString(typeline("What column number is your chr in (value): "))
        exonstartcol=toString(typeline("What column is your start position in (value): "))
        exonendcol=toString(typeline("What column is your end position in (value): "))

        exonalpha=toString(typeline("Alpha (value): "))
        exoncolor=toString(typeline("Color (char): "))
        exonborder=toString(typeline("Outline Color (char): "))
        exonledgend=tolower(toString(typeline("Include a legend (y/n):")))
    }

    else{
        warning("FILE MUST BE BEDFILE")
    }

    if((exonledgend %like% "y") || (exonledgend %like% "yes")){
        
        type<-function(msg="Enter text: "){
            if(!interactive()){
                cat(msg);
                nameexon<-readLines("stdin",n=1);
            }
            return(nameexon)
        }
        nameexon=toString(type("Name Of Legend: "))
    }
}

pf<-NULL
merge.png.pdf <- function(pdfFile, pngFiles, deletePngFiles=FALSE) {
  
  #### Package Install ####
  pngPackageExists <- require ("png")
  
  if ( !pngPackageExists ) {
    install.packages ("png")
    library ("png")
    
  }
  #########################
  
  pdf(pdfFile)
  
  n <- length(pngFiles)
  
  for( i in 1:n) {
    
    pngFile <- pngFiles[i]
    
    pngRaster <- readPNG(pngFile)
    
    grid.raster(pngRaster, width=unit(0.8, "npc"), height= unit(0.8, "npc"))
    
    if (i < n) plot.new()
    
  }
  
  dev.off()
  
  if (deletePngFiles) {
    
    unlink(pngFiles)
    
  }
  
}
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################

getChrom<-function(){ 

    if(!is.null(opt$cordfile)){
        lm<-as.data.frame(fread(opt$cordfile))
        colnames(lm)<-c('Chr','start','stop', 'GeneName')
        rownames(lm)<-lm$GeneName
        #print(lm)
        
        lm.list<-split(lm, seq(nrow(lm)))
        
    }

    else{
        lm.list<-NULL
        stop('Missing Input: Chromosome Position Bed File')
    }
    return(lm.list)
    
}  
chrom<-getChrom()


print(chrom)
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
system.time(
getInput<-function(){ 
    df=NULL
    

    if(!is.null(opt$file)){
      ucdf<-NULL
      
      l <- 1
      while(l <= length(chrom)){
        if (l <= length(chrom[[l]]) ){
        chrv1<-(chrom[[l]]$Chr)
        ucdf<- rbind(ucdf, as.data.frame(chrv1))
        }
        l <- l+1
      }
      
      listchr<-c(unique(ucdf))
      

        if(!is.null(opt$dir)){
            
            if(is.null(opt$lineMult)){
                colofline=toString("blue")
            }

            if(!is.null(opt$lineMult)){
                getcol<-toString(opt$lineMult)
                colofline<-toString(getcol) 
            }
            
            fr<-as.data.frame(fread(opt$file))
            
            
            colnames(fr)<-c('Chr','start','stop', 'Ecopynum')
            fr<- fr %>% filter(fr$Chr %in% listchr)
            namef<-rep('Inital', times= length(fr$stop))
            colit<-rep(colofline, times= length(fr$stop))
            fr$FILENAME<-namef
            fr$COLOR<-colit
            

            dir<-toString(opt$dir)

            if(!is.null(opt$dirpattern)){
                pattern<-toString(opt$dirpattern)
                pattern1<-paste0(pattern,"$")
                c<-list.files(path=dir , pattern=pattern1, full.names=TRUE)
            }
            
            if(is.null(opt$dirpattern)){
                pattern1<-paste0(".bed","$")
                c<-list.files(path=dir , pattern=pattern1, full.names=TRUE)
            }
            

            
            readdata <- function(x)
            {
              
              mkrd <- fread(x)
              mkrd<-as.data.frame(mkrd)
              
              colnames(mkrd)<-c('Chr','start','stop', 'Ecopynum')
              mkrd<- mkrd %>% filter(mkrd$Chr %in% listchr)
              
              df=NULL
              dt=NULL
              exonfilt=NULL
              
              
              colit<-rep('black', times= length(mkrd$stop))
              mkrd$COLOR<-colit
              return(mkrd)
            }
            
            f<-lapply(c,readdata)
            f<-mapply(cbind, f, "FILENAME"=1:length(f), SIMPLIFY=F)
            
            d<- do.call("rbind", f)
            df<-rbind(fr, data.frame(d))
        }
            
        
  
        if(is.null(opt$dir)){
            df<-as.data.frame(fread(opt$file))
            colnames(df)<-c('Chr','start','stop', 'Ecopynum')
            df<- df %>% filter(df$Chr %in% listchr)
        }
    }

    else{
        df<-NULL
        stop('Missing Input: Quickmer-2 Bed File')
    }
    
    return(df)
    
   
}
)

dataf<-getInput()
data<-as.data.frame(dataf)


######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
lista<-NULL
holda<-NULL
if(is.null(opt$plotTogether) && is.null(opt$dir) ){
  # Opt -l flag
  
  if(is.null(opt$lineColor)){
    lim=toString("black")
  }
  
  if(!is.null(opt$lineColor)){
    getcol<-toString(opt$lineColor)
    getcol<-strsplit(getcol, ',', fixed=TRUE)
    
    if(as.numeric(length(getcol[[1]])) <2){
      lim<-toString(opt$lineColor)
    }
    
    else{
      warning("TOO MANY ARGUMENTS ALLOWED FOR -l FLAG, NO -t FLAG CALLED")
    }
  }
  i<-1
  
  while(i <= length(chrom)){
   
    if(i >= 0){
      # Filter Based on Each Input
      
      test1<-data %>% filter(data$Chr %in% chrom[[i]]$Chr)
      #a %<-% {
      #print(test2$Chr)
      test1<-test1 %>% filter(as.numeric(test1$start) <= as.numeric(chrom[[i]]$stop))
      #}
      #b %<-% {
      test1<-test1 %>% filter(as.numeric(test1$stop) >= as.numeric(chrom[[i]]$start))
      #}
      #system.time(test2<-cbind(a,b))
      
      if(nrow(test1)>=1){
      test2<-as.data.frame(test1)
      
      }
      
      # Create String For Title With Gene Name
      getnamegene<-toString(rownames(chrom[[i]]))
      titlestring<-toString("Gene")
      titlename<-paste(getnamegene,titlestring)
      if(is.null(opt$png)){
       
      #Rename File 
      endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
      nameoffile<-paste(getnamegene,endfilename, sep="")
      }
      
      if(is.null(opt$multifile) && !is.null(opt$png)){
        endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
        nameoffile<-paste(getnamegene,endfilename, sep="")
      }
      if(!is.null(opt$multifile) && !is.null(opt$png)){
        endfilename<-toString(".png")
        nameoffile<-paste(i,endfilename, sep="")
      }
      # Create String For X Axis Name
      xstring<-chrom[[i]]$Chr
      str<-toString("Chromosome Position (")
      str1<-paste(str,xstring)
      str2<-toString(")")
      xaxisname<-paste(str1,str2)
    
      if((is.null(opt$multifile)) && (is.null(opt$png))){
        
          pdf("Rplots.pdf")
        }
        if((is.null(opt$multifile) )&& (!is.null(opt$png))){
          png(nameoffile)
          
        }
        if((!is.null(opt$multifile)) && (!is.null(opt$png))){
        png(nameoffile)
       }
      if(is.null(opt$multifile)){
        
        plot(as.numeric(test2$start), y=as.numeric(test2$Ecopynum), col=lim, xlab=xaxisname ,
             ylab="Estimated Copy Number", main=titlename , pch=19,
             ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1)
        abline(h=2, col="grey", lty=2)
        
        #######WORKING########
        ###ANNOTATION TRACK###
        if(!is.null(opt$annotation)){
          
          exondf<-as.data.frame(fread(opt$annotation))
          exonstartcol<-as.numeric(exonstartcol)
          exonendcol<-as.numeric(exonendcol)
          exonchrcol<-as.numeric(exonchrcol)
          exoncolor<-toString(exoncolor)
          exonalpha<-as.numeric(exonalpha)
          exonborder<-tolower(toString(exonborder))
          
          names(exondf)[exonstartcol]<-"start"
          names(exondf)[exonendcol]<-"end"
          names(exondf)[exonchrcol]<-"Chr"
          exondf<-exondf %>% filter(exondf$Chr %in% chrom[[i]]$Chr)
          
          exondf<-exondf %>% filter(as.numeric(exondf$start) <= as.numeric(chrom[[i]]$stop))
          exondf<-exondf %>% filter(as.numeric(exondf$end) >= as.numeric(chrom[[i]]$start))
          
          if(1 <= length(exondf$start)){
            j<-1
            
            while(j <= length(exondf$start)){
              
              for(j in 1:nrow(exondf)){
                
                strte<-exondf[j, "start"]
                stope<-exondf[j, "end"]
                rect(strte, 0, stope, 0.2, col=exoncolor, border=exonborder)
              }
              j<-j+1
            }
          }
          
          if(exonledgend %like% "y" || exonledgend %like% "yes"){
            legend("topright", fill=exoncolor, legend=nameexon, bty="n", border=exonborder)
          }
        }
        ######################
        #########END##########
        #########END##########
        if(!is.null(opt$box)){
          if(nrow(test2)>=2){
          box<- test2 %>% select('Ecopynum')
          print(box)
          plotdim <- par("plt")
          xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * 0.25
          xright   = plotdim[2]  #
          ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.25  #
          ytop     = plotdim[4]
          par(
            fig = c(xleft, xright, ybottom, ytop)
            , mar=c(0,0,0,0)
            , new=TRUE
          )
          boxplot(box$Ecopynum, col="white")
          points(mean(box$Ecopynum), col="Blue", pch=23)
          }
        }
        
      }
      
      if(!is.null(opt$multifile)){
        
        
        
        plot(test2$start, y=test2$Ecopynum, col=lim, xlab=xaxisname ,
             ylab="Estimated Copy Number", main=titlename , pch=19,
             ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1)
        abline(h=2, col="grey", lty=2)
        
        #######WORKING########
        ###ANNOTATION TRACK###
        
        if(!is.null(opt$annotation)){
          
          exondf<-as.data.frame(fread(opt$annotation))
          exonstartcol<-as.numeric(exonstartcol)
          exonendcol<-as.numeric(exonendcol)
          exonchrcol<-as.numeric(exonchrcol)
          exoncolor<-toString(exoncolor)
          exonalpha<-as.numeric(exonalpha)
          exonborder<-tolower(toString(exonborder))
          
          names(exondf)[exonstartcol]<-"start"
          names(exondf)[exonendcol]<-"end"
          names(exondf)[exonchrcol]<-"Chr"
          exondf<-exondf %>% filter(exondf$Chr %in% chrom[[i]]$Chr)
          
          exondf<-exondf %>% filter(as.numeric(exondf$start) <= as.numeric(chrom[[i]]$stop))
          exondf<-exondf %>% filter(as.numeric(exondf$end) >= as.numeric(chrom[[i]]$start))
          
          if(1 <= length(exondf$start)){
            j<-1
            
            while(j <= length(exondf$start)){
              
              for(j in 1:nrow(exondf)){
                
                strte<-exondf[j, "start"]
                stope<-exondf[j, "end"]
                rect(strte, 0, stope, 0.2, col=exoncolor, border=exonborder)
              }
              j<-j+1
            }
          }
          
          if(exonledgend %like% "y" || exonledgend %like% "yes"){
            legend("topright", fill=exoncolor, legend=nameexon, bty="n", border=exonborder)
          }
        }  
        ######################
        #########END##########
        if(!is.null(opt$box)){
          if(nrow(test2)>=2){
          box<- test2 %>% select('Ecopynum')
          
          plotdim <- par("plt")
          xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * 0.25
          xright   = plotdim[2]  #
          ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.25  #
          ytop     = plotdim[4]
          par(
            fig = c(xleft, xright, ybottom, ytop)
            , mar=c(0,0,0,0)
            , new=TRUE
          )
          boxplot(box$Ecopynum, col="white")
          points(mean(box$Ecopynum), col="Blue", pch=23)
          }
        }
      
      
      
      
    
    if( (!is.null(opt$multifile) ) && !(is.null(opt$png)) ){
      list1<-c(nameoffile)
      lista<-rbind(lista,(list1))
      working<-getwd()
      
      
      
    }
      }
    }
    i<-i + 1 
  }
  
  dev.off()
  graphics.off()

  if(!is.null(opt$multifile) && (is.null(opt$png)) ){
    file.rename("Rplots.pdf",shortname)
  }

  if(!is.null(opt$multifile) && !(is.null(opt$png)) ){
    
    list1<-lista
    print(length(list1))
    i<-1
    while(i <= length(list1)){
      if(i <= length(list1)){
        
        print(i)
        print(paste0(working,'/',lista[i,]))
        file.exists(paste0(working,'/',lista[i,]))
        merge.png.pdf (pdfFile = paste0(i,".pdf"), pngFiles = paste0(working,'/',lista[i,]), deletePngFiles = F)
        hold1<-c(paste0(i,".pdf"))
        holda<- rbind(holda,hold1)
        if(!is.null(opt$multifile)){
          filename<-toString(opt$multifile)
          pdf_combine(holda,output=filename)
          
        }
        i<-i+1
      }
      
    }
    
    
  }
 
lapply(holda,function(x){file.remove(x)})

}
lapply(lista,function(x){unlink(paste0(working,"/",x))})
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################

if(!is.null(opt$cordfile) && !is.null(opt$plotTogether) &&  is.null(opt$dir)){
  
  colored<-randomColor()
  
  getname<-strsplit(toupper(opt$plotTogether), split=',', fixed=TRUE)
  getname<-as.data.frame(getname)
  colnames(getname)<-'GeneName'
  getname
  
  newdata<-as.data.frame(fread(opt$cordfile))
  colnames(newdata)<-c('Chr','start','stop','GeneName')
  rownames(newdata)<-toupper(newdata$GeneName)
  
  for(test in 1:nrow(getname)){
    testingit<-toupper(getname[test,"GeneName"])
    
    if(testingit %in% row.names(newdata)){
      newdataset<-newdata %>% filter(row.names(newdata) %in% getname$GeneName)
    }
    
    else{
      warning("ONE OR MORE CHROMOSOMES ARE MISSING FROM -c CHR.BED FILE")
    }
  }
  
  test<-data %>% filter(data$Chr %in% newdataset$Chr)
  
  gg<-as.character(newdataset$Chr)
  names(gg)<-as.character(newdataset$GeneName)
  
  gh<-as.character(newdataset$GeneName)
  names(gh)<-as.character(newdataset$Chr)
  hh<-toString(paste(gg,":",gh))
  tname<-paste("Genes:", toString(gh))
  chromname<-paste("Chromosome Position(s) (" ,hh ,")")
  chromname<-str_replace(chromname,fixed(",")," , ")
  
  df=NULL
  dt=NULL
  exonfilt=NULL
  
  for(row in 1:nrow(newdataset)){
    
    colored<-randomColor()
    genenames<-newdataset[row,"GeneName"]
    startit<-newdataset[row,"start"]
    stopit<-newdataset[row,"stop"]
    
    testx<-test %>% filter(as.numeric(test$start) <= as.numeric(stopit))
    testz<-testx %>% filter(as.numeric(testx$stop) >= as.numeric(startit))
    
    geneinputname<-rep(genenames, times= length(testz$stop))
    
    newdataframe<-testz %>% select("start","Ecopynum")
    newdataframe$GeneNames<-geneinputname
    
    testp<-testz %>% select("Chr","start","stop")
    newdataframe$Color<-colored  
    
    df<-rbind(df, data.frame(newdataframe))
    exonfilt<-rbind(exonfilt,data.frame(testp))
      
  }
  if(is.null(opt$png)){
  shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.pdf")
  shortname<-str_replace(shortname,fixed(", "),"_")
  shortname<-str_replace(shortname,fixed(","),"_")
  shortname<-str_replace(shortname,fixed(" "),"")
  }
  if(!is.null(opt$png)){
    shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.png")
    shortname<-str_replace(shortname,fixed(", "),"_")
    shortname<-str_replace(shortname,fixed(","),"_")
    shortname<-str_replace(shortname,fixed(" "),"")
  }

  if(!is.null(opt$plotTogether) && (is.null(opt$lineColor))){
    
    col<-as.character(df$Color)
    names(col)<-as.character(df$GeneNames)
    
    plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=GeneNames)) + geom_line(aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1))  + ggtitle(tname) +
      xlab(chromname) + ylab("Estimated Copy Number") + ggtitle(tname) + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
    plot<-plot + geom_hline(yintercept=2, color="light grey")
    
    #######WORKING########
    ###ANNOTATION TRACK###
    
    if(!is.null(opt$annotation)){
      
      exondf<-as.data.frame(fread(opt$annotation))
      exonstartcol<-as.numeric(exonstartcol)
      exonendcol<-as.numeric(exonendcol)
      exonchrcol<-as.numeric(exonchrcol)
      exoncolor<-toString(exoncolor)
      exonalpha<-as.numeric(exonalpha)
      exonborder<-tolower(toString(exonborder))
      
      names(exondf)[exonstartcol]<-"start"
      names(exondf)[exonendcol]<-"end"
      names(exondf)[exonchrcol]<-"Chr"
      
      exondf<-exondf %>% filter(exondf$Chr %in% exonfilt$Chr)
      
      max<-max(exonfilt$stop)
      min<-min(exonfilt$start)
      
      exondf<-exondf %>% filter(as.numeric(exondf$start) <= max)
      exondf<-exondf %>% filter(as.numeric(exondf$end) >= min)
      
      if(1 <= length(exondf$start)){
        j<-1
        
        while(j <= length(exondf$start)){
          
          for(j in 1:nrow(exondf)){
            
            strte<-exondf[j, "start"]
            stope<-exondf[j, "end"]
            
            if(exonledgend %like% "y" || exonledgend %like% "yes"){
              plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon))+ guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(tname) + xlab(chromname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
              plot<-plot1 + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
              
            }
            
            else{
              plot<-plot + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
              
            }
          }
          j<-j+1
        }
      }
    }
    #####################
    #########END#########
    if(!is.null(opt$box)){
      box<- df %>% select("GeneNames",'Ecopynum','Color')
      
      gd<- box %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
      
      p1<-box %>% ggplot(aes(x=GeneNames,y=Ecopynum)) + geom_boxplot() + stat_boxplot(geom = "errorbar", width = 0.5) + geom_point(data=gd,color="blue", pch=23) + geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
      
      plot<- plot + inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
    }
    ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
  }
  
  if(!is.null(opt$plotTogether) && (!is.null(opt$lineColor))){
    
    getlength<-strsplit(opt$plotTogether, split=',', fixed=TRUE)
    getcol<-toString(opt$lineColor)
    getcol<-strsplit(getcol, ',', fixed=TRUE)
    
    if((as.numeric(length(getcol[[1]])) >= 2) && (as.numeric(length(getcol[[1]])) == as.numeric(length(getlength[[1]])))){
      colorit<-c(getcol[[1]])
      
      plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=GeneNames)) + geom_line(aes(color=GeneNames), alpha=.5) +scale_alpha_manual(values=c(1,0.4,1)) +scale_color_manual(values=colorit) + ggtitle(tname) +
        xlab(chromname) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
      plot<-plot + geom_hline(yintercept=2, color="light grey")
      
      #######WORKING########
      ###ANNOTATION TRACK###
      
      if(!is.null(opt$annotation)){
        exondf<-as.data.frame(fread(opt$annotation))
        exonstartcol<-as.numeric(exonstartcol)
        exonendcol<-as.numeric(exonendcol)
        exonchrcol<-as.numeric(exonchrcol)
        exoncolor<-toString(exoncolor)
        exonalpha<-as.numeric(exonalpha)
        exonborder<-tolower(toString(exonborder))
        
        names(exondf)[exonstartcol]<-"start"
        names(exondf)[exonendcol]<-"end"
        names(exondf)[exonchrcol]<-"Chr"
        
        exondf<-exondf %>% filter(exondf$Chr %in% exonfilt$Chr)
        
        max<-max(exonfilt$stop)
        min<-min(exonfilt$start)
        
        exondf<-exondf %>% filter(as.numeric(exondf$start) <= max)
        exondf<-exondf %>% filter(as.numeric(exondf$end) >= min)
        
        if(1 <= length(exondf$start)){
          j<-1
          
          while(j <= length(exondf$start)){
            
            for(j in 1:nrow(exondf)){
              
              strte<-exondf[j, "start"]
              stope<-exondf[j, "end"]
              
              if(exonledgend %like% "y" || exonledgend %like% "yes"){
                plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(tname) + xlab(chromname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
                plot<-plot1 + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
              
              }
              
              else{
                plot<-plot + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
                 
              }    
            }
            j<-j+1
          }
        }  
      }
      #######################
      ##########END##########
      if(!is.null(opt$box)){
        box<- df %>% select("GeneNames",'Ecopynum','Color')
        
        gd<- box %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
        
        p1<-box %>% ggplot(aes(x=GeneNames,y=Ecopynum)) + geom_boxplot() + stat_boxplot(geom = "errorbar", width = 0.5) + geom_point(data=gd,color="blue", pch=23) + geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
        
        plot<- plot + inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
      }
      ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
    }
    
    if(as.numeric(length(getcol[[1]])) <2){
      warning("TOO FEW ARGUMENTS ALLOWED FOR -l FLAG")
    }
    
    if(as.numeric(length(getcol[[1]])) > as.numeric(length(getlength[[1]]))){
      warning("TOO MANY ARGUMENTS ALLOWED FOR -l FLAG")
    }
    
    if(as.numeric(length(getcol[[1]])) < as.numeric(length(getlength[[1]]))){
      warning("TOO MANY ARGUMENTS ALLOWED FOR -t FLAG")
    }
  }
}

#--------------------#
######################
######################
######################
######################
######################
######################
######################
######################
######################
######################
#--------------------#

if(!is.null(opt$dir) && !is.null(opt$cordfile)){
  
  if(is.null(opt$plotTogether)){
    
    i<-1
    list<-NULL
    hold<-NULL
    while(i <= length(chrom)){
      
      if(i >= 0){
        
        # Filter Based on Each Input
        test2<-data %>% filter(data$Chr %in% chrom[[i]]$Chr)
        test2<-test2 %>% filter(as.numeric(test2$start) <= as.numeric(chrom[[i]]$stop))
        test2<-test2 %>% filter(as.numeric(test2$stop) >= as.numeric(chrom[[i]]$start))
        # for (r in 1:nrow(test2)){
        #   check<- test2[r,"Ecopynum"]
        #   name<-test2[r,"FILENAME"]
        #   x<-as.numeric(check)
        #   if((x > 3)==TRUE){
        #     print(name,x)
        #   }
        # }
        
        if(is.null(opt$lineColor)){
          #Create String For Title With Gene Name
          getnamegene<-toString(rownames(chrom[[i]]))
          titlestring<-toString("Gene")
          titlename<-paste(getnamegene,titlestring)
          
          #Rename File 
          if(is.null(opt$png)){
          endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
          nameoffile<-paste(getnamegene,endfilename, sep="") 
          }
          if(!is.null(opt$png)){
            endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
            nameoffile<-paste(getnamegene,endfilename, sep="") 
          }
          xstring<-chrom[[i]]$Chr
          str<-toString("Chromosome Position (")
          str1<-paste(str,xstring)
          str2<-toString(")")
          xaxisname<-paste(str1,str2)
          
          col<-test2 %>% select("COLOR","FILENAME")
          col<-unique(col)
          
          col1<-as.character(col$COLOR)
          names(col1)<-as.character(col$FILENAME)
          h<-toString(opt$file)
          a<-toString(rbind(str_trunc(h,20,"right")))
          
          test2<-test2 %>% group_by(FILENAME)
          if(is.null(opt$lineMult)){
            colofline=toString("blue")
          }
          
          if(!is.null(opt$lineMult)){
            getcol<-toString(opt$lineMult)
            colofline<-toString(getcol)  
          }
          
          
          plot<-test2 %>% ggplot(aes(x=start, y=Ecopynum, group=FILENAME,color=COLOR)) +scale_color_manual(values=c("black",colofline), labels=c("Files From Directory",a), name='Data')+ geom_line() + ggtitle(titlename) + xlab(xaxisname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
          plot<-plot + geom_hline(yintercept=2, color="light grey")
        
          
          
          #######WORKING########
          ###ANNOTATION TRACK###
          
          if(!is.null(opt$annotation)){
            
            exondf<-as.data.frame(fread(opt$annotation))
            exonstartcol<-as.numeric(exonstartcol)
            exonendcol<-as.numeric(exonendcol)
            exonchrcol<-as.numeric(exonchrcol)
            exoncolor<-toString(exoncolor)
            exonalpha<-as.numeric(exonalpha)
            exonborder<-tolower(toString(exonborder))
            
            names(exondf)[exonstartcol]<-"start"
            names(exondf)[exonendcol]<-"end"
            names(exondf)[exonchrcol]<-"Chr"
            
            exondf<-exondf %>% filter(exondf$Chr %in% chrom[[i]]$Chr)
            exondf<-exondf %>% filter(as.numeric(exondf$start) <= as.numeric(chrom[[i]]$stop))
            exondf<-exondf %>% filter(as.numeric(exondf$end) >= as.numeric(chrom[[i]]$start))
            
            if(1 <= length(exondf$start)){
              j<-1
              
              while(j <= length(exondf$start)){
                
                for(j in 1:nrow(exondf)){
                  
                  strte<-exondf[j, "start"]
                  stope<-exondf[j, "end"]
                  
                  if(exonledgend %like% "y" || exonledgend %like% "yes"){
                    plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder) + scale_fill_discrete(labels=c(nameexon))+ scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon))
                    plot<-plot1 + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder)
                    #plot<-plot1 
                  }
                  
                  else{
                    plot<-plot + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
                    #plot<-plot1 
                  }
                }
                j<-j+1
              } 
            }
          }
          #########################
          ###########END###########
          if(!is.null(opt$box)){
            if(nrow(test2)>=2){
            a<-as.data.frame(sapply(split(test2$Ecopynum,test2$FILENAME),mean))
            colnames(a)<-"means"
            
            
            fun_mean <- function(x){
              return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}
            
            
            p1<-a %>% ggplot(aes(x='',y=means)) + geom_boxplot() + stat_summary(fun= mean, geom="point",colour="black", size=1) + geom_point(aes(x='', y=a["Inital",1]), color="blue",pch=23)
            
            #stat_summary(fun.data = fun_mean, geom="text",colour="black",size=1, vjust=-0.7)
            p1<- p1 + scale_color_manual(values=c("black",colofline),)+ geom_hline(yintercept=2, color="light grey")+ geom_hline(yintercept=2.7, color="light grey") + geom_hline(yintercept=1.3, color="light grey")  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
            plot<-plot+  inset_element(p1, left=0.6, bottom=0.6, right=1.05, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE) 
            }
            
          }
          ggsave(plot=plot, width=7, height=6, dpi=300, filename=nameoffile)
          
          if(!is.null(opt$png) && !is.null(opt$multifile)){
            list1<-c(nameoffile)
            
            pngFiles<-rbind(list, (list1))
            merge.png.pdf (pdfFile = paste0(i,".pdf"), pngFiles = pngFiles, deletePngFiles = T)
            
            
            hold1<-c(paste0(i,".pdf"))
            hold<- rbind(hold,hold1)
            if(!is.null(opt$multifile)){
            filename<-toString(opt$multifile)
            pdf_combine(hold,output=filename)
            }
            
          }
          if(is.null(opt$png) && !is.null(opt$multifile)){
            list1<-c(nameoffile)
            list<-rbind(list, (list1))
            if(!is.null(opt$multifile)){
            filename<-toString(opt$multifile)
            pdf_combine(list, output=filename)
            }
            
          }
        }
        i<-i + 1
        
        if(!is.null(opt$lineColor)){
          warning("-l FLAG MUST BE NOT BE USED")
        }
      }
      
    }
    
    if(!is.null(opt$multifile)){
      lapply(list,function(x){file.remove(x)})
    }
    
    if(!is.null(opt$png)){
      lapply(hold,function(x){file.remove(x)})
    } 
  }
  
  #--------------------#
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  ######################
  #--------------------#
  
  if(!is.null(opt$plotTogether)){
    
    colored<-randomColor()
    
    getname<-strsplit(toupper(opt$plotTogether), split=',', fixed=TRUE)
    getname<-as.data.frame(getname)
    colnames(getname)<-'GeneName'
    getname
    
    newdata<-as.data.frame(fread(opt$cordfile))
    colnames(newdata)<-c('Chr','start','stop','GeneName')
    rownames(newdata)<-toupper(newdata$GeneName)
    
    for(test in 1:nrow(getname)){
      
      testingit<-toupper(getname[test,"GeneName"])
      
      if(testingit %in% row.names(newdata)){
        newdataset<-newdata %>% filter(row.names(newdata) %in% getname$GeneName)   
      }
      
      else{
        warning("ONE OR MORE CHROMOSOMES ARE MISSING FROM -c CHR.BED FILE")
      }
    }
    
    test<-data %>% filter(data$Chr %in% newdataset$Chr)
    
    gg<-as.character(newdataset$Chr)
    names(gg)<-as.character(newdataset$GeneName)
    
    gh<-as.character(newdataset$GeneName)
    names(gh)<-as.character(newdataset$Chr)
    hh<-toString(paste(gg,":",gh))
    tname<-paste("Genes:", toString(gh))
    chromname<-paste("Chromosome Position(s) (" ,hh ,")")
    chromname<-str_replace(chromname,fixed(",")," , ")
    
    df=NULL
    dt=NULL
    exonfilt=NULL
    
    for(row in 1:nrow(newdataset)){
      
      colored<-randomColor()
      genenames<-newdataset[row,"GeneName"]
      startit<-newdataset[row,"start"]
      stopit<-newdataset[row,"stop"]
      
      testx<-test %>% filter(as.numeric(test$start) <= as.numeric(stopit))
      testz<-testx %>% filter(as.numeric(testx$stop) >= as.numeric(startit))
      
      geneinputname<-rep(genenames, times=length(testz$stop))
      
      newdataframe<-testz %>% select("start","Ecopynum")
      newdataframe$GeneNames<-geneinputname
      
      
      testp<-testz %>% select("Chr","start","stop")
      newdataframe$Color<-colored
      newdataframe$dog<-paste(newdataframe$GeneNames, " ", testz$FILENAME)
      df<-rbind(df, data.frame(newdataframe))
      exonfilt<-rbind(exonfilt,data.frame(testp))
      
    }
    if(!is.null(opt$png)){
      shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.png")
      shortname<-str_replace(shortname,fixed(", "),"_")
      shortname<-str_replace(shortname,fixed(","),"_")
      shortname<-str_replace(shortname,fixed(" "),"")
    }
    if(!is.null(opt$png)){
    shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.pdf")
    shortname<-str_replace(shortname,fixed(", "),"_")
    shortname<-str_replace(shortname,fixed(","),"_")
    shortname<-str_replace(shortname,fixed(" "),"")
    }
    
    if(!is.null(opt$plotTogether) && (is.null(opt$lineColor))){
      
      if(is.null(opt$lineMult)){
        colofline=toString("blue")
      }
      
      if(!is.null(opt$lineMult)){
        getcol<-toString(opt$lineMult) 
        colofline<-toString(getcol) 
      }
      
      col<-as.character(df$Color)
      names(col)<-as.character(df$GeneNames)
      `%notlike%`<-Negate(`%like%`)
      
      plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog)) + geom_line(data=subset(df,dog %notlike% "Inital"),aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1))  + ggtitle(tname) +
        xlab(chromname) + ylab("Estimated Copy Number") + ggtitle(tname) + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
      
      plot<-plot + geom_hline(yintercept=2, color="light grey")
      plot<-plot
      plot<-plot + geom_line(data=subset(df,dog %like% "Inital"),color=colofline)
      
      #######WORKING########
      ###ANNOTATION TRACK###
      
      if(!is.null(opt$annotation)){
        
        exondf<-as.data.frame(fread(opt$annotation))
        exonstartcol<-as.numeric(exonstartcol)
        exonendcol<-as.numeric(exonendcol)
        exonchrcol<-as.numeric(exonchrcol)
        exoncolor<-toString(exoncolor)
        exonalpha<-as.numeric(exonalpha)
        exonborder<-tolower(toString(exonborder))
        
        names(exondf)[exonstartcol]<-"start"
        names(exondf)[exonendcol]<-"end"
        names(exondf)[exonchrcol]<-"Chr"
        
        exondf<-exondf %>% filter(exondf$Chr %in% exonfilt$Chr)
        
        max<-max(exonfilt$stop)
        min<-min(exonfilt$start)
        
        exondf<-exondf %>% filter(as.numeric(exondf$start) <= max)  
        exondf<-exondf %>% filter(as.numeric(exondf$end) >= min)
        
        if(1 <= length(exondf$start)){
          j<-1
          
          while(j <= length(exondf$start)){
            
            for(j in 1:nrow(exondf)){
              
              strte<-exondf[j, "start"]
              stope<-exondf[j, "end"]
              
              if(exonledgend %like% "y" || exonledgend %like% "yes"){
                plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon))+ guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(tname) + xlab(chromname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
                plot<-plot1 + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
                 
              }
              
              else{
                plot<-plot + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
                
              }        
            }
            j<-j+1
          }  
        }
      }
      #########################
      ###########END###########
      
      
      #stat_summary(fun.data = fun_mean, geom="text",colour="black",size=1, vjust=-0.7)
      if(!is.null(opt$box)){
        box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
        box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
        gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
        
        p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
        
        plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
      }
     
      ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
    }
    
    if(!is.null(opt$plotTogether) && (!is.null(opt$lineColor))){
      
      if(is.null(opt$lineMult)){
        colofline=toString("blue")
      }
      
      if(!is.null(opt$lineMult)){
        getcol<-toString(opt$lineMult)
        colofline<-toString(getcol) 
      }
      
      getlength<-strsplit(opt$plotTogether, split=',', fixed=TRUE)
      getcol<-toString(opt$lineColor)
      getcol<-strsplit(getcol, ',', fixed=TRUE)
      
      if((as.numeric(length(getcol[[1]])) >= 2) && (as.numeric(length(getcol[[1]])) == as.numeric(length(getlength[[1]])))){
        
        colorit<-c(getcol[[1]])
        
        plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog)) + geom_line(aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1)) +scale_color_manual(values=colorit) + ggtitle(tname) +
          xlab(chromname) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) 
        plot<-plot + geom_line(data=subset(df,dog %like% "Inital"),color=colofline)  
        
        #######WORKING########
        ###ANNOTATION TRACK###
        
        if(!is.null(opt$annotation)){
          
          exondf<-as.data.frame(fread(opt$annotation))
          exonstartcol<-as.numeric(exonstartcol)
          exonendcol<-as.numeric(exonendcol)
          exonchrcol<-as.numeric(exonchrcol)
          exoncolor<-toString(exoncolor)
          exonalpha<-as.numeric(exonalpha)
          exonborder<-tolower(toString(exonborder))
          
          names(exondf)[exonstartcol]<-"start"
          names(exondf)[exonendcol]<-"end"
          names(exondf)[exonchrcol]<-"Chr"
          
          exondf<-exondf %>% filter(exondf$Chr %in% exonfilt$Chr)
          
          max<-max(exonfilt$stop)
          min<-min(exonfilt$start)
          
          exondf<-exondf %>% filter(as.numeric(exondf$start) <= max)
          exondf<-exondf %>% filter(as.numeric(exondf$end) >= min)
          
          if(1 <= length(exondf$start)){
            j<-1
            
            while(j <= length(exondf$start)){
              
              for(j in 1:nrow(exondf)){
                
                strte<-exondf[j, "start"]
                stope<-exondf[j, "end"]
                
                if(exonledgend %like% "y" || exonledgend %like% "yes"){
                  plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(tname) + xlab(chromname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
                  plot1<-plot1 + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder)  
                  plot<-plot1 
                }
                
                else{
                  plot1<-plot + annotate("rect", xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor, alpha=exonalpha, color=exonborder) 
                  plot<-plot1 
                } 
              }
              j<-j+1
            }
          }    
        }
        #######################
        ##########END##########
        if(!is.null(opt$box)){
          box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
          box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
          gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
          
          p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
          
          plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
        }
        ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
      }
      
      if(as.numeric(length(getcol[[1]])) <2){
        warning("TOO FEW ARGUMENTS ALLOWED FOR -l FLAG")
      }
      
      if(as.numeric(length(getcol[[1]])) > as.numeric(length(getlength[[1]]))){
        warning("TOO MANY ARGUMENTS ALLOWED FOR -l FLAG")
      }
      
      if(as.numeric(length(getcol[[1]])) < as.numeric(length(getlength[[1]]))){
        warning("TOO MANY ARGUMENTS ALLOWED FOR -t FLAG")
      }
    }
  } 
}

end1 <- Sys.time()  
time1 <- end1 - start1 
print(time1)

