# Quick2Plot

This script offered a quick, and customization way to create Ling Graphs based on data produced by Quick-mer2.

##Requirements

###R studio with the following packages:

optparse
base
data.table
ggplot2
dplyr
tibble
stringr
stringi
randomcoloR
pdftools

Quick-mer2 Bed file

Chromosome file

###Flag Options:

"-f": -REQUIRED

  default=NULL
  
    Input the initial/main Quick-mer2 bed file
    
"-d":

  default=NULL
  
    Input the directory pathway to bed files, if multiple bed files are wanted to be plotted onto one plot. 
    
"-a":

  default=NULL
  
    Flag used if user would like to add an annotation track using a bedfile. This option will is interactive and will require answers to the following questions:
    
    If input file a bed file (y/n):
    What column number is your chr in (value): 
    What column is your start position in (value): 
    What column is your end position in (value): 
    Alpha (value): 
    Color (char): 
    Outline Color (char): 
    Include a legend (y/n):
    Name Of Legend: 
    
"-c": -REQUIRED

  default=NULL
  
    Input a Chromosome Coordinates Bed File
    
"-l":

  default=NULL
  
    Input a Color For Lines 
    
    ex: -l blue,red
    
"-t":

  default=NULL
  
    Flag to Plot Certain Chromosomes, From Chromosome Bed File, Together.
    
      ex: -t NGF,POGZ
    
"-o"

  default=NULL
  
    Allows User to Name, and Place All Chromosome Plots Outputs Into One PDF, Defualt Is Plots Each Cordinate To Its Own PDF File
    
    
Example Runs:

  Rscript Quick2Plot.R -f myInitialQuick-merBedFile.bed -c myChrFile.bed -a MyAnnotationFile.bed  -d /Users/elvisa/Desktop/Rtesting/fake -o MyMultiplePDFOutputFile.pdf
