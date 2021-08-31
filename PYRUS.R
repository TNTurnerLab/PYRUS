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
                    make_option(c("-s", "--singleplots"), default=NULL, type="character",help='plots each file in dirpath individually'),
                    make_option(c("-y", "--rename"), default=NULL, type="character",help='Rename file -f on plot'),
                    make_option(c("-n", "--non"), default=NULL, type="character",help='Non-interactive Annotation'),
                    make_option(c("-o", "--multifile"), default=NULL, type="character",help='Outputs The Chromosome Plots Into One PDF, Defualt Is Plots Each Cordinate To Its Own PDF File'))
opt<-parse_args(OptionParser(option_list=option_list))



#########Non-interactive annotation track###########
if(!is.null(opt$non) & !is.null(opt$annotation)){ 
  annotate<-strsplit(opt$non, split=',', fixed=TRUE)
  
  if(((as.numeric(length(annotate[[1]])) >=8)) & (as.numeric(length(annotate[[1]])) <=9)){
  
    isbed<-annotate[[1]][1]
    isbed<-as.character(isbed)
  
    if(isbed =="y" |isbed == "yes"){
    
      exonchrcol<-as.numeric(annotate[[1]][2])
      exonstartcol<- as.numeric(annotate[[1]][3])
      exonendcol<- as.numeric(annotate[[1]][4])
      exonalpha<- as.numeric(annotate[[1]][5])
  
      exoncolor<- as.character(annotate[[1]][6])
      exonborder<- as.character(annotate[[1]][7])
      exonledgend<- as.character(annotate[[1]][8])
    
      if((as.numeric(length(annotate[[1]])) ==9)){
        nameexon<- as.character(annotate[[1]][9])
      }
      if(exonledgend=="y"| exonledgend=="yes"){
        if((as.numeric(length(annotate[[1]])) !=9)){
          warning("NO EXON TITLE NAME GIVEN")
          stop
        }
      }
  
    }
    else{
      warning("FILE MUST BE BEDFILE")
    }
  }
  else{
    if((as.numeric(length(annotate[[1]])) <=8)){
      warning("TOO FEW INPUTS")
      stop()
    }
    if((as.numeric(length(annotate[[1]])) >=9)){
      warning("TOO MANY INPUTS")
      stop()
    }
  }
}
################################

###########Interactive annotation track###########
if(is.null(opt$non)){
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
        stop()
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

}

#####################################

##############Print PNG#################
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
    
    grid.raster(pngRaster, width=unit(18, "cm"), height= unit(15, "cm"))
    
    if (i < n) plot.new()
    
  }
  
  dev.off()
  
  if (deletePngFiles) {
    
    unlink(pngFiles)
    
  }
  
}

#####################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#

###########Read and get data from Chr file##################
getChrom<-function(){ 

    if(!is.null(opt$cordfile)){
        lm<-as.data.frame(fread(opt$cordfile))
        colnames(lm)<-c('Chr','start','stop', 'GeneName')
        rownames(lm)<-lm$GeneName

        lm.list<-split(lm, seq(nrow(lm)))
        
    }

    else{
        lm.list<-NULL
        stop('Missing Input: Chromosome Position Bed File')
    }
    return(lm.list)
    
}  
chrom<-getChrom()

#####################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#

##########Read and get data from -f flag file##############
system.time(
getInput<-function(){ 
    df=NULL
    fr<-NULL

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
      
      listchr<-unique(ucdf)
      lschr<-as.data.frame(listchr)
      colnames(lschr)<- "name"
     
      

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
           
            fr<- fr %>% filter(fr$Chr %in% lschr$name)
            
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
            
            
            
            if(is.null(opt$singleplots)){
            readdata <- function(x)
            {

              mkrd <- fread(x)
              mkrd<-as.data.frame(mkrd)

              colnames(mkrd)<-c('Chr','start','stop', 'Ecopynum')
              mkrd<- mkrd %>% filter(mkrd$Chr %in% lschr$name)

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
            
            if(!is.null(opt$singleplots)){

              readdatasingle <- function(x)
              {

                mkrd <- fread(x)
                mkrd<-as.data.frame(mkrd)

                colnames(mkrd)<-c('Chr','start','stop', 'Ecopynum')
                mkrd<- mkrd %>% filter(mkrd$Chr %in% lschr$name)
                #mkrd<- mkrd[mkrd$Chr %in% listchr,]
               

                name<-toString(basename(x))
                FILENAME<- rep(name, times= length(mkrd$stop))
                mkrd$FILENAME<-FILENAME

                df=NULL
                dt=NULL
                exonfilt=NULL


                colit<-rep('black', times= length(mkrd$stop))
                mkrd$COLOR<-colit

               return(mkrd)
              }
              f<-lapply(c,readdatasingle)

              d<- do.call("rbind", f)
              df<-rbind(fr,data.frame(d))

            }
            
            }
      
      if(is.null(opt$dir)){
        df<-as.data.frame(fread(opt$file))
        colnames(df)<-c('Chr','start','stop', 'Ecopynum')
        
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




#####################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#

#########Basic plotter without -t or -d flag################


lista<-NULL
holda<-NULL

if(is.null(opt$plotTogether) && is.null(opt$dir) ){

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
     
      test2<-data %>% filter(data$Chr %in% chrom[[i]]$Chr)
      
      test2<-test2 %>% filter(as.numeric(test2$start) <= as.numeric(chrom[[i]]$stop))
     
      test2<-test2 %>% filter(as.numeric(test2$stop) >= as.numeric(chrom[[i]]$start))
     
      
      test2<-as.data.frame(test2)
      

      getnamegene<-toString(unique(chrom[[i]]$GeneName))
      titlestring<-toString("Gene")
  
     
      if(is.null(opt$png)){

      endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
      nameoffile<-paste(getnamegene,endfilename, sep="")
      }
      
      if(is.null(opt$multifile) && !is.null(opt$png)){
        endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
        nameoffile<-paste(getnamegene,endfilename, sep="")
      }
      
      
      if(is.null(opt$multifile) && !is.null(opt$png)){
        endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
        nameoffile<-paste(getnamegene,endfilename, sep="")
      }
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
             ylab="Estimated Copy Number", main=(bquote(paste(italic(.(getnamegene)),paste(" "),paste(.(titlestring))))), pch=19,
             ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1)
        abline(h=2, col="grey", lty=2)
       

        ###ANNOTATION TRACK###
        if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
          
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
          
          if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
            legend("topright", fill=exoncolor, legend=nameexon, bty="n", border=exonborder)
          }
        }
        ######################

        ###########Box###########
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
        ##########Plot###########
      }
      
      if(!is.null(opt$multifile)){
        
        plot(test2$start, y=test2$Ecopynum, col=lim, xlab=xaxisname ,
             ylab="Estimated Copy Number", main=(bquote(paste(italic(.(getnamegene)),paste(" "),paste(.(titlestring))))), pch=19,
             ylim= c(0,6), type= 'l' , cex.main=1, cex.lab=1, cex.axis=1)
        abline(h=2, col="grey", lty=2)

        ###ANNOTATION TRACK###
        
        if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
          
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
          
          if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
            legend("topright", fill=exoncolor, legend=nameexon, bty="n", border=exonborder)
          }
        }  
        ######################
        
        ###########Box###########
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
        ##########Plot###########
      
        if( (!is.null(opt$multifile) ) && !(is.null(opt$png)) ){
          list1<-c(nameoffile)
          lista<-rbind(lista,(list1))
          working<-getwd()
      
        }
      }
      i<-i + 1
      if(is.null(opt$multifile) && (is.null(opt$png)) ){
        shortname<-nameoffile
        file.rename("Rplots.pdf",nameoffile)
      }
     
  }
  }
  dev.off()
  graphics.off()

  if(!is.null(opt$multifile) && (is.null(opt$png)) ){
    shortname<- toString(opt$multifile)
    file.rename("Rplots.pdf",shortname)
  }
  
  if(!is.null(opt$multifile) && !(is.null(opt$png)) ){
    
    list1<-lista
    
    i<-1
    
    while(i <= length(list1)){
      
      if(i <= length(list1)){
        
        
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


#####################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#

#########Plotter with -t flag and without -d ############

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
      stop()
    }
  }
  
  test<-data %>% filter(data$Chr %in% newdataset$Chr)
  
  gg<-as.character(newdataset$Chr)
  names(gg)<-as.character(newdataset$GeneName)
  
  gh<-as.character(newdataset$GeneName)
  names(gh)<-as.character(newdataset$Chr)
  hh<-toString(paste(gg,":",gh))
  tname<-paste(toString(gh))
  chr1<-paste(hh)
  
  
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
    
    plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=GeneNames)) + geom_line(aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) + xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
    plot<-plot + geom_hline(yintercept=2, color="light grey")

    ###ANNOTATION TRACK###
    
    if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
      
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
            
            if((exonledgend %like% "y" || exonledgend %like% "yes")| (as.numeric(length(annotate[[1]])) ==9)){
              plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon))+ guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +           xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
    
    ###########Box###########
    if(!is.null(opt$box)){
      box<- df %>% select("GeneNames",'Ecopynum','Color')
      
      gd<- box %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
      
      p1<-box %>% ggplot(aes(x=GeneNames,y=Ecopynum)) + geom_boxplot() + stat_boxplot(geom = "errorbar", width = 0.5) + geom_point(data=gd,color="blue", pch=23) + geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
      
      plot<- plot + inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
    }
    ##########Plot###########
    
    ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
  }
  
  if(!is.null(opt$plotTogether) && (!is.null(opt$lineColor))){
    
    getlength<-strsplit(opt$plotTogether, split=',', fixed=TRUE)
    getcol<-toString(opt$lineColor)
    getcol<-strsplit(getcol, ',', fixed=TRUE)
    
    if((as.numeric(length(getcol[[1]])) >= 2) && (as.numeric(length(getcol[[1]])) == as.numeric(length(getlength[[1]])))){
      colorit<-c(getcol[[1]])
      
      plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=GeneNames)) + geom_line(aes(color=GeneNames), alpha=.5) +scale_alpha_manual(values=c(1,0.4,1)) +scale_color_manual(values=colorit) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) + xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
      plot<-plot + geom_hline(yintercept=2, color="light grey")

      
      
      
      ###ANNOTATION TRACK###
      
      if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
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
              
              if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
                plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +           xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
      
      ###########Box###########
      if(!is.null(opt$box)){
        box<- df %>% select("GeneNames",'Ecopynum','Color')
        
        gd<- box %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
        
        p1<-box %>% ggplot(aes(x=GeneNames,y=Ecopynum)) + geom_boxplot() + stat_boxplot(geom = "errorbar", width = 0.5) + geom_point(data=gd,color="blue", pch=23) + geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
        
        plot<- plot + inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
      }
      
      ##########Plot###########
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

##############










#####################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#
if(is.null(opt$singleplots)){
#########Plotting with -d flag##############
if(!is.null(opt$dir) && !is.null(opt$cordfile)){
  
  #########plotting without -t flag################
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
          h<-toString(basename(opt$file))
          a<-toString(rbind(str_trunc(h,20,"right")))
          
          
          
          if(is.null(opt$lineMult)){
            colofline=toString("blue")
          }
          
          if(!is.null(opt$lineMult)){
            getcol<-toString(opt$lineMult)
            colofline<-toString(getcol)  
          }
          if(!is.null(opt$rename)){
            a<-toString(opt$rename)
          }
          if(is.null(opt$rename)){
            a<-toString(rbind(str_trunc(h,20,"right")))
          }
          
          
          plot<-test2 %>% ggplot(aes(x=start, y=Ecopynum, group=FILENAME,color=COLOR)) +scale_color_manual(values=c("black",colofline), labels=c("Files From Directory",a), name='Data')+ geom_line() + ggtitle(bquote(paste(italic(.(getnamegene)),paste(" "),paste(.(titlestring))))) + xlab(xaxisname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
          plot<-plot + geom_hline(yintercept=2, color="light grey")
          
        
          ###ANNOTATION TRACK###
          
          if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
            
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
                  
                  if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
                    plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder) + scale_fill_discrete(labels=c(nameexon))+ scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon))
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
          
          ###########Box###########
          if(!is.null(opt$box)){
            if(nrow(test2)>=2){
              a<-as.data.frame(sapply(split(test2$Ecopynum,test2$FILENAME),mean))
              colnames(a)<-"means"
            
            
              fun_mean <- function(x){
                return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

              p1<-a %>% ggplot(aes(x='',y=means)) + geom_boxplot() + stat_summary(fun= mean, geom="point",colour="black", size=1) + geom_point(aes(x='', y=a["Inital",1]), color="blue",pch=23)
            
              p1<- p1 + scale_color_manual(values=c("black",colofline),)+ geom_hline(yintercept=2, color="light grey")+ geom_hline(yintercept=2.7, color="light grey") + geom_hline(yintercept=1.3, color="light grey")  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
              plot<-plot+  inset_element(p1, left=0.6, bottom=0.6, right=1.05, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE) 
            }
          }
          ##########Plot###########
          
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
  

  
  ##########plotting -d and -t flags#################
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
    #tname<-paste("Genes:", toString(gh))
    chr1<-paste(hh)
    #chromname<-str_replace(chromname,fixed(",")," , ")
    tname<-paste(toString(gh))
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
      
      plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog)) + geom_line(data=subset(df,dog %notlike% "Inital"),aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1))  + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +
        xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
      
      plot<-plot + geom_hline(yintercept=2, color="light grey")
      plot<-plot
      plot<-plot + geom_line(data=subset(df,dog %like% "Inital"),color=colofline)

      ###ANNOTATION TRACK###
      
      if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
        
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
              
              if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
                plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon))+ guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +           xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
     
      ###########Box###########
      if(!is.null(opt$box)){
        box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
        box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
        gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
        
        p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
        
        plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
      }
      ##########Plot###########
      
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
        
        plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog)) + geom_line(aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1)) +scale_color_manual(values=colorit) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) + xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) 
        plot<-plot + geom_line(data=subset(df,dog %like% "Inital"),color=colofline)  
        plot<-plot + geom_hline(yintercept=2, color="light grey")
        
  #       
  #       bquote(paste(paste("Genes: "),italic(.(tname))))) + xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")"))))
        ###ANNOTATION TRACK###
        
        if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
          
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
                
                if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
                  plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) + xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
        
        
        ###########Box#########
        if(!is.null(opt$box)){
          
          box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
          box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
          gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
          
          p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
          
          plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
        }
        ##########Plot#########
        
        ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
      }
      
      if(as.numeric(length(getcol[[1]])) <2){
        warning("TOO FEW ARGUMENTS ALLOWED FOR -l FLAG")
        stop()
      }
      
      if(as.numeric(length(getcol[[1]])) > as.numeric(length(getlength[[1]]))){
        warning("TOO MANY ARGUMENTS ALLOWED FOR -l FLAG")
        stop()
      }
      
      if(as.numeric(length(getcol[[1]])) < as.numeric(length(getlength[[1]]))){
        warning("TOO MANY ARGUMENTS ALLOWED FOR -t FLAG")
        stop()
      }
    }
  } 
}

}

############################################

#-----------------------------------#
#####################################
#####################################
#-----------------------------------#

###########Plot Directory Full##############

if(!is.null(opt$singleplots)){
#########Plotting with -d flag##############
if(!is.null(opt$dir) && !is.null(opt$cordfile)){
  
  #########plotting without -t flag################
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
        test2<-test2 %>% group_by(FILENAME)
        
        if(nrow(test2)>=1){
        `%notlike%`<-Negate(`%like%`)
        tes1<-unique(test2$FILENAME)
        tes1<-as.data.frame(tes1)
        colnames(tes1)<-"names"
        tes1<- tes1 %>% filter(tes1$names %notlike% "Inital")
        
        
        for (r in 1:nrow(tes1)){
          name<-tes1[r, "names"]
          
        if(is.null(opt$lineColor)){
          #Create String For Title With Gene Name
          getnamegene<-toString(rownames(chrom[[i]]))
          #getnamegene<-italic(getnamegene)
          titlestring<-toString("Gene")
          titlename<-paste(getnamegene,titlestring)
          
          #Rename File 
          if(is.null(opt$png)){
            endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
            nameoffile<-paste0(name,"_",getnamegene,endfilename, sep="") 
          }
          
          if(!is.null(opt$png)){
            endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
            nameoffile<-paste0(name,"_",getnamegene,endfilename, sep="") 
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
          h<-toString(basename(opt$file))
          if(!is.null(opt$rename)){
            a<-toString(opt$rename)
          }
          if(is.null(opt$rename)){
          a<-toString(rbind(str_trunc(h,20,"right")))
          }
          name1<-toString(rbind(str_trunc(name,15,"right")))
          
          if(is.null(opt$lineMult)){
            colofline=toString("blue")
          }
          
          if(!is.null(opt$lineMult)){
            getcol<-toString(opt$lineMult)
            colofline<-toString(getcol)  
          }
         
              plot<-test2 %>% ggplot(aes(x=start, y=Ecopynum), group=FILENAME,color=COLOR)+ geom_line(data=subset(test2,FILENAME %like% name),aes(color="black"))+ geom_line(data=subset(test2,FILENAME %like% "Inital"),aes(color=colofline))  + scale_color_manual(values=c("black",colofline), labels=c(name1,a), name='Data') + ggtitle(bquote(paste(italic(.(getnamegene)),paste(" "),paste(.(titlestring))))) + xlab(xaxisname) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
              
              plot<-plot + geom_hline(yintercept=2, color="light grey")
              
              nameoffile<-paste0(toString(name),nameoffile)

        ###ANNOTATION TRACK###
        
        if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
          
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
                
                if((exonledgend %like% "y" || exonledgend %like% "yes") || (as.numeric(length(annotate[[1]])) ==9)){
                  plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder) + scale_fill_discrete(labels=c(nameexon))+ scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon))
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
        
        ###########Box###########
        if(!is.null(opt$box)){
          
            a<-as.data.frame(sapply(split(test2$Ecopynum,test2$FILENAME),mean))
            colnames(a)<-"means"
            a<-subset(a,rownames(a) %in% name | rownames(a) %in% "Inital")

            fun_mean <- function(x){
              return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}
            
            p1<-a %>% ggplot(aes(x='',y=means)) + geom_boxplot() + stat_summary(fun= mean, geom="point",colour="black", size=1) + geom_point(aes(x='', y=a["Inital",1]), color="blue",pch=23)
            
            p1<- p1 + scale_color_manual(values=c("black",colofline),)+ geom_hline(yintercept=2, color="light grey")+ geom_hline(yintercept=2.7, color="light grey") + geom_hline(yintercept=1.3, color="light grey")  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
            plot<-plot+  inset_element(p1, left=0.6, bottom=0.6, right=1.05, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE) 
            
        }
        ##########Plot###########
        
        ggsave(plot=plot, width=7, height=6, dpi=300, filename=nameoffile)
        
        if(!is.null(opt$png) && !is.null(opt$multifile)){
          list1<-c(nameoffile)
          
          pngFiles<-rbind(list, (list1))
          merge.png.pdf (pdfFile = paste0(nameoffile,".pdf"), pngFiles = pngFiles, deletePngFiles = T)
          
          hold1<-c(paste0(nameoffile,".pdf"))
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

##########plotting -d and -t flags#################
if(!is.null(opt$plotTogether)){
  i<-1
  list<-NULL
  hold<-NULL
  
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
  
  tname<-paste(toString(gh))
  chr2<-paste("Chromosome Position(s) (" )
  chr1<-paste(hh)
  #chromname<-str_replace(chromname,fixed(",")," , ")
 
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
    newdataframe$FILENAME<-testz$FILENAME
    df<-rbind(df, data.frame(newdataframe))
    
    exonfilt<-rbind(exonfilt,data.frame(testp))
    
  }
  
  if(!is.null(opt$png)){
    shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.png")
    shortname<-str_replace(shortname,fixed(", "),"_")
    shortname<-str_replace(shortname,fixed(","),"_")
    shortname<-str_replace(shortname,fixed(" "),"")
    shortname<-paste0(shortname)
  }
  
  if(!is.null(opt$png)){
    shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.pdf")
    shortname<-str_replace(shortname,fixed(", "),"_")
    shortname<-str_replace(shortname,fixed(","),"_")
    shortname<-str_replace(shortname,fixed(" "),"")
    shortname<-paste0(shortname)
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
    tes1<-unique(df$FILENAME)

    tes1<-as.data.frame(tes1)
    colnames(tes1)<-"names"
    tes1<- tes1 %>% filter(tes1$names %notlike% "Inital")
    getnamegene<-toString(rownames(chrom[[i]]))
    
    for (r in 1:nrow(tes1)){

      name<-tes1[r, "names"]
      
      
        if(is.null(opt$png)){
          endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
          shortname<-paste0(name,"_",getnamegene,endfilename, sep="") 
        }
        
        if(!is.null(opt$png)){
          endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
          shortname<-paste0(name,"_",getnamegene,endfilename, sep="") 
        }
      
      name1<-toString(rbind(str_trunc(name,10,"right")))
      
      plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog)) + geom_line(data=subset(df,dog %like% name),color="black")  + geom_line(data=subset(df,dog %like% "Inital"),aes(group=dog,color=GeneNames)) +scale_alpha_manual(values=c(1,0.4,1))  + scale_color_discrete(name=name1)+ ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +
      xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number")  + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
    
      plot<-plot + geom_hline(yintercept=2, color="light grey")
      plot<-plot 
   
 
    
    ###ANNOTATION TRACK###
    
    if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
      
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
            
            if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
              plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon))+ guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +           xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
    
    ###########Box###########
    if(!is.null(opt$box)){
      box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
      box<-subset(box,dog %like% name | dog %like% "Inital")
      
      box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
      gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
      
      p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
      
      plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
    }
    ##########Plot###########
    
    ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
    
    
    if(!is.null(opt$png) && !is.null(opt$multifile)){
      list1<-c(shortname)
      
      pngFiles<-rbind(list, (list1))
      merge.png.pdf (pdfFile = paste0(shortname,".pdf"), pngFiles = pngFiles, deletePngFiles = T)
      
      hold1<-c(paste0(shortname,".pdf"))
     
      hold<- rbind(hold,hold1)
      
      if(!is.null(opt$multifile)){
        filename<-toString(opt$multifile)
        pdf_combine(hold,output=filename)
        
      }
      
    }
    if(is.null(opt$png) && !is.null(opt$multifile)){
      
      list1<-c(shortname)
      list<-rbind(list, (list1))
      
      if(!is.null(opt$multifile)){
        filename<-toString(opt$multifile)
        pdf_combine(list, output=filename)
      }
    }
    
    }
  }
  if(!is.null(opt$multifile)){
    lapply(list,function(x){file.remove(x)})
  }
  
  if(!is.null(opt$png)){
    lapply(hold,function(x){file.remove(x)})
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
      
      `%notlike%`<-Negate(`%like%`)
      tes1<-unique(df$FILENAME)
      
      tes1<-as.data.frame(tes1)
      
      colnames(tes1)<-"names"
      
      tes1<- tes1 %>% filter(tes1$names %notlike% "Inital")
      
      getnamegene<-toString(rownames(chrom[[i]]))
      for (r in 1:nrow(tes1)){
        
        name<-tes1[r, "names"]
        
        
        if(is.null(opt$png)){
          endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
          shortname<-paste0(name,"_",getnamegene,endfilename, sep="") 
        }
        
        if(!is.null(opt$png)){
          endfilename<-toString("_CHR_Position_vs_Copy_Num.png")
          shortname<-paste0(name,"_",getnamegene,endfilename, sep="") 
        }
        
        name1<-toString(rbind(str_trunc(name,10,"right")))
        colorit<-c(getcol[[1]])
        
        plot<-df %>% ggplot(aes(x=start, y=Ecopynum, group=dog),color=GeneNames)+ geom_line(data=subset(df,dog %like% name),color="Black") + geom_line(data=subset(df,dog %like% "Inital"),aes(color=GeneNames)) + scale_alpha_manual(values=c(1,0.4,1))  +scale_color_manual(values=colorit, name=name1) +  ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +
          xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw() + theme(plot.title=element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))
        
        plot<-plot + geom_hline(yintercept=2, color="light grey")
        plot<-plot 
        
      
      ###ANNOTATION TRACK###
      
      if((!is.null(opt$annotation) | !is.null(opt$non) ) | (!is.null(opt$non) & !is.null(opt$annotation))){
        
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
              
              if((exonledgend %like% "y" || exonledgend %like% "yes") | (as.numeric(length(annotate[[1]])) ==9)){
                plot1<-plot  + geom_rect(aes(xmin=as.numeric(strte), xmax=as.numeric(stope), ymin=0, ymax=0.2, fill=exoncolor), alpha=exonalpha, color=exonborder)  + scale_fill_discrete(labels=c(nameexon)) + scale_fill_manual(values=exoncolor, labels=c(nameexon)) + guides(fill=guide_legend(title=nameexon)) + theme(legend.justification=c(1, 1), legend.position=c(1, 1),  title=element_blank(),legend.key=element_rect(colour=NA, fill=NA)) + ggtitle(bquote(paste(paste("Genes: "),italic(.(tname))))) +           xlab(bquote(paste(paste("Chromosome Position(s) ("),italic(.(chr1)),paste(")")))) + ylab("Estimated Copy Number") + scale_x_continuous(n.breaks=6) + scale_y_continuous(limits=c(0, 6), n.breaks=6)  + theme_bw()
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
      
      
      ###########Box#########
      if(!is.null(opt$box)){
        
        box<- df %>% select("GeneNames",'Ecopynum','Color','dog')
        box<-subset(box,dog %like% name | dog %like% "Inital")
        box2<- subset(box,dog %like% "Inital", select=c(GeneNames,Ecopynum))
        gd<- box2 %>% group_by(GeneNames) %>% summarise(Ecopynum=mean(Ecopynum))
        
        p1<-box %>% ggplot(aes( x=GeneNames,y=Ecopynum)) + geom_boxplot(notch = TRUE) + stat_boxplot(geom = "errorbar", width = 0.5) + scale_color_manual(values=c("black",colofline),)+ geom_point(data=gd,color="blue", pch=23) +geom_hline(yintercept=2, color="light grey") + geom_hline(yintercept=2.7, color="light grey",lty=3) + geom_hline(yintercept=1.3, color="light grey",lty=3)  + scale_y_continuous(limits=c(0, 6), n.breaks=6,position="right")+ theme_bw() + theme(legend.position = "none",rect = element_rect(fill = "transparent"),panel.border = element_blank(),axis.line.y = element_line(colour = "dark grey"),plot.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(color="black", size=5), axis.ticks.x=element_blank())
        
        plot<-plot+ inset_element(p1, left=0.6, bottom=0.6, right=1.0425, top=0.999, align_to = "plot", on_top=TRUE, clip=FALSE)
      }
      ##########Plot#########
      
      ggsave(plot=plot, width=7, height=6, dpi=300, filename=shortname)
        if(!is.null(opt$png) && !is.null(opt$multifile)){
          list1<-c(shortname)
          
          pngFiles<-rbind(list, (list1))
          merge.png.pdf (pdfFile = paste0(shortname,".pdf"), pngFiles = pngFiles, deletePngFiles = T)
          
          hold1<-c(paste0(shortname,".pdf"))
          
          hold<- rbind(hold,hold1)
          
          if(!is.null(opt$multifile)){
            filename<-toString(opt$multifile)
            pdf_combine(hold,output=filename)
            
          }
          
        }
        if(is.null(opt$png) && !is.null(opt$multifile)){
         
          list1<-c(shortname)
          list<-rbind(list, (list1))
          
          if(!is.null(opt$multifile)){
            filename<-toString(opt$multifile)
            pdf_combine(list, output=filename)
          }
        }
        
      }
    
    if(!is.null(opt$multifile)){
      lapply(list,function(x){file.remove(x)})
    }
    
    if(!is.null(opt$png)){
      lapply(hold,function(x){file.remove(x)})
    }
        
    
  
    if(as.numeric(length(getcol[[1]])) <2){
      warning("TOO FEW ARGUMENTS ALLOWED FOR -l FLAG")
      stop()
    }
    
    if(as.numeric(length(getcol[[1]])) > as.numeric(length(getlength[[1]]))){
      warning("TOO MANY ARGUMENTS ALLOWED FOR -l FLAG")
      stop()
    }
    
    if(as.numeric(length(getcol[[1]])) < as.numeric(length(getlength[[1]]))){
      warning("TOO MANY ARGUMENTS ALLOWED FOR -t FLAG")
      stop()
    }
  }
  }
}
}

}

#####################################
end1 <- Sys.time()  
time1 <- end1 - start1 
print(time1)