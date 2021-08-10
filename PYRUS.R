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

option_list<-list(make_option(c("-f", "--file"), default=NULL, type="character",help='Input an Inital Quick-mer2 Bed File'),
                    make_option(c("-d", "--dir"), default=NULL, type="character",help='Input Directory of Bed Files'), 
                    make_option(c("-p", "--dirpattern"), default=NULL, type="character",help='Pattern of files in Directory, defualt .bed'),
                    make_option(c("-a", "--annotation"), default=NULL, type="character",help='Input an Annotation'),
                    make_option(c("-c", "--cordfile"), default=NULL, type="character",help='Input a Chr Coordinates Bed File'), 
                    make_option(c("-l", "--lineColor"), default=NULL, type="character",help='Input a Color(s) ex: blue,red'),
                    make_option(c("-i", "--lineMult"), default=NULL, type="character",help='Input a Color for -f line when using -d flag ex: blue'),
                    make_option(c("-t", "--plotTogether"), default=NULL, type="character",help='Plot Chromosomes Together. ex: NGF,POGZ'),
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
        lm.list<-split(lm, seq(nrow(lm)))
    }

    else{
        lm.list<-NULL
        stop('Missing Input: Chromosome Position Bed File')
    }
    return(lm.list)
}  
chrom<-getChrom()

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

getInput<-function(){ 
    df=NULL

    if(!is.null(opt$file)){
    
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
            namef<-rep('Inital', times= length(fr$stop))
            colit<-rep(colofline, times= length(fr$stop))
            fr$FILENAME<-namef
            fr$COLOR<-colit
        
            dir<-toString(opt$dir)

            if(!is.null(opt$dirpattern)){
                pattern<-toString(opt$dirpattern)
                pattern1<-paste0(pattern,"$")
            }
            
            if(is.null(opt$dirpattern)){
                pattern1<-paste0(".bed","$")
            }

            c<-list.files(path=dir , pattern=pattern1)
            c<-as.data.frame(c)
            dd<-NULL
    
            for(row in 1:nrow(c)){

                fileplot<-toString(c[row,'c'])
                fileplot1<-paste0(dir,'/',fileplot)
                
                d<-as.data.frame(fread(fileplot1))
                colnames(d)<-c('Chr','start','stop', 'Ecopynum')
                name<-rep(fileplot, times= length(d$stop))
                colit<-rep('black', times= length(d$stop))
                d$COLOR<-colit
                d$FILENAME<-name
                
                dd<-rbind(dd, data.frame(d))
            }
            df<-rbind(fr, data.frame(dd))
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

if(is.null(opt$plotTogether) && is.null(opt$dir)){
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
            test2<-data %>% filter(data$Chr %in% chrom[[i]]$Chr)
            test2<-test2 %>% filter(as.numeric(test2$start) <= as.numeric(chrom[[i]]$stop))
            test2<-test2 %>% filter(as.numeric(test2$stop) >= as.numeric(chrom[[i]]$start))

            test2<-as.data.frame(test2)

            # Create String For Title With Gene Name
            getnamegene<-toString(rownames(chrom[[i]]))
            titlestring<-toString("Gene")
            titlename<-paste(getnamegene,titlestring)
            
            #Rename File 
            endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
            nameoffile<-paste(getnamegene,endfilename, sep="")
            
            # Create String For X Axis Name
            xstring<-chrom[[i]]$Chr
            str<-toString("Chromosome Position (")
            str1<-paste(str,xstring)
            str2<-toString(")")
            xaxisname<-paste(str1,str2)
            
            if(is.null(opt$multifile)){
                pdf("Rplots.pdf")
            
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
                shortname<-nameoffile
                file.rename("Rplots.pdf",nameoffile)
            }
            
            if(!is.null(opt$multifile)){

                shortname<-toString(opt$multifile)
                
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
            }
        }
        i<-i + 1  
    }
    dev.off()

    if(!is.null(opt$multifile)){
        file.rename("Rplots.pdf",shortname)
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
        datatableout<-testz %>% select("start","stop","Ecopynum","Chr")
        datatableout$GeneNames<-geneinputname
        dt<-rbind(dt, data.frame(datatableout))    
    }

    shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.pdf")
    shortname<-str_replace(shortname,fixed(", "),"_")
    shortname<-str_replace(shortname,fixed(","),"_")
    shortname<-str_replace(shortname,fixed(" "),"")
    
    tablename<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.txt")
    tablename<-str_replace(tablename,fixed(", "),"_")
    tablename<-str_replace(tablename,fixed(" "),"")
    
    tablename<-str_replace(toString(tablename),fixed(","),"_")
    dt<-dt %>% select("GeneNames","Chr","start","stop","Ecopynum")
    write.table(dt, file=tablename, row.names=FALSE)
    
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
        #####################
        #########END#########
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
                    endfilename<-toString("_CHR_Position_vs_Copy_Num.pdf")
                    nameoffile<-paste(getnamegene,endfilename, sep="") 

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
                    #########################
                    ###########END###########
                    ggsave(plot=plot, width=7, height=6, dpi=300, filename=nameoffile)
                    list1<-c(nameoffile)
                    list<-rbind(list, (list1))
                    pdf_combine(list, output="joined.pdf")
                }
                i<-i + 1
                    
                if(!is.null(opt$lineColor)){
                    warning("-l FLAG MUST BE NOT BE USED")
                }
            }
        }
        
        if(!is.null(opt$multifile)){
            
            filename<-toString(opt$multifile)
            tab<-str_replace(filename,fixed(".pdf"),"_table.txt")
            file.rename("joined.pdf",filename)

            for(l in list){
                file.remove(l)
            }
        }
    
        if(is.null(opt$multifile)){
            file.remove("joined.pdf")
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
            datatableout<-testz %>% select("start","stop","Ecopynum","Chr")
            datatableout$GeneNames<-geneinputname
            dt<-rbind(dt, data.frame(datatableout))
        
        }
        
        shortname<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.pdf")
        shortname<-str_replace(shortname,fixed(", "),"_")
        shortname<-str_replace(shortname,fixed(","),"_")
        shortname<-str_replace(shortname,fixed(" "),"")
        
        tablename<-paste0("Genes_",toString(gh),"_CHR_Position_vs_Copy_Num.txt")
        tablename<-str_replace(tablename,fixed(", "),"_")
        tablename<-str_replace(tablename,fixed(" "),"")
        
        tablename<-str_replace(toString(tablename),fixed(","),"_")
        dt<-dt %>% select("GeneNames","Chr","start","stop","Ecopynum")
        write.table(dt, file=tablename, row.names=FALSE)
        
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
            #########################
            ###########END###########
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