<img align="right" src="https://user-images.githubusercontent.com/77067493/129043760-26d57868-5f15-46bb-a00e-e82f8c50fa9c.png" height="20%" width="20%" >

<h1>PYRUS</h1>

  PYRUS is an R script that utilizes flags to create a quick, and customizable way to create line graphs based on bed file data. Utilation of PYRUS includes having a quick plotter for copy number variation within a given gene domain range of a referenced bed file. Along with the estimated number of copies portrayed on the y axis, the x axsis will show the gene chromosome and location. If given a directory of files, and a main reference file, PYRUS has the ability to plot a gene's CNVs from all bed files inside a directory. This program will create a single plot for every gene inside a chromosome bed file, but can also plot multiple genes CNVs onto a single plot given the correct flag. It is recommened to use this feature for genes that are at a reasonable distance from one another to generate the best visual plot. In addition, there is an option to add an annotation track to the plot. The annotation track is in the lower quadrent of the graph and can be used for plotting exons in relation to the gene's CNVs. The script also gives the user to option to plot a box-plot window of the average copy number values of all files fed into the script. The user will see a blue, rhombus denoting the average copy number of the file represented by the -f flag. If the user wishes, they can also print out a png image in place of a pdf file. A combination of these options can be used to generate a CNVs plot tailored to the needs of the user.

## Requirements

### R studio with the following packages:

* optparse
* base
* data.table
* ggplot2
* dplyr
* tibble
* stringr
* stringi
* randomcoloR
* pdftools
* grid
* png
* patchwork

### File Requirments

Bed files are the only acceptable file types for this script

### Flag Options:

"-f": REQUIRED

  Input the initial/main Quick-mer2 bed file
  
    default=NULL
    
      ex: -f Quick-mer2BedFile.bed

    
"-d": 

  Input the directory pathway to bed files, if multiple bed files are wanted to be plotted onto one plot.
    
    default=NULL
    
      ex: -d /path/to/a/directory/with/bed/files
"-p": 

  Input a ending pattern for files found inside the directory in "-d" flag. Must be the endding pattern of the file.
    
    default=.bed
    
      ex: -p .t1.cn.bed.gz  
"-i":

  Input a Color for initial line when using -d line, default color is blue.
  
    default=NULL
    
      ex: -i green
    
"-a": 

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
    
    default=NULL
      
      ex: -a exonbedfile.bed
    
"-c": REQUIRED

  Input a Chromosome Coordinates Bed File
    
    default=NULL
    
      ex: myChrbedfile.bed
  
    
"-l":

  Input a Color For Lines 
    
    default=NULL
    
      ex: -l blue,red

"-t":

  Flag to Plot Certain Chromosomes, From Chromosome Bed File, Together.
    
    default=NULL
    
      ex: -t NGF,POGZ

"-o":

  
  Allows User to Name, and Place All Chromosome Plots Outputs Into One PDF, Default Is Plots Each Coordinate To Its Own PDF File
  
    default=NULL
     
     ex: -o multiplePagePDFFileofMyPlots.pdf
     
"-r":

  Flag to change the output from a pdf file to a png file. When creating multiple plots at a single time, if the -o flag is used in combination with the -r flag, output will be png files combined into one pdf file.
    
    defualt=NULL

      ex: -r y or -o output.pdf -r y
      
"-b":

   This flag will create a box-plot window in the top right corner of the existing line plot. The box plot will display the individual, marked with the -f flag, with a blue dot.

    defualt=NULL
      
      ex: -b y 
    
Example Runs:

  Rscript Quick2Plot.R -f myInitialQuick-merBedFile.bed -c myChrFile.bed -a MyAnnotationFile.bed  -d /Users/elvisa/Desktop/Rtesting/fake -o MyMultiplePDFOutputFile.pdf
  
  
 ## Plotting combinations
 
  * -f -c

    * Plots all chromosomes in -c flag's file with a black line and in seprate PDF files
      
      * Add -l {COLORNAME} to change line color. 
       
        * ex : -f qm.bed -c chr.bed -l blue
      
      * Add -o {OUTPUTFILENAME.PDF} to make a multiple page PDF. 
       
        * ex : -f qm.bed -c chr.bed -o myoutput.pdf 
      
      * Add -d {/path/to/directory} to plot multiple lines from other bedfiles on top of the initial file given with -f flag.  
        
        * ex : -f qm.bed -c chr.bed -d /path/to/directory
      
      * Add -d {/path/to/directory} and -i {COLORNAME} to change the initial -f file line color when plotting against all bed files in directory. 
        
          * ex :  -f qm.bed -c chr.bed -d /path/to/directory -i green
       * Add -d {/path/to/directory} and -p {PATTERN} to grab only files with the inpputed ending inside directory to be used in plot. 
        
          * ex :  -f qm.bed -c chr.bed -d /path/to/directory -i green -p .bed.gz

       * Add -a {AnnotationBedFile.bed} to add an annotation track to all plots. 
      
          * ex : -f qm.bed -c chr.bed -d /path/to/directory -a myexonstoplot.bed
          
          * ex : -f qm.bed -c chr.bed -a myexonstoplot.bed 
         
      * Add -b y to create a box-plot window on the line graph displaying the average copy number estimation and in a blue rhombus, the average copy number of the file presented in -f flag will be displayed on top of tye box plot. 
	        
          * ex: -f qm.bed -c chr.bed -b y
      
      * Add -r y to create a png file output. 
	        
          * ex: -f qm.bed -c chr.bed -r y -b y
       * Add -r y in combination with the -o flag to make a pdf of png images.
          
          * ex :  -f qm.bed -c chr.bed -o output.pdf -r y -b y


  * -f -c -t 
   
   * Plots only named chromosomes from -c flag's file. Can not use -o flag when -t flag is present.

      * Add -t {CHR1,CHR2} to plot certain chromosomes together on a single plot. -f qm.bed -c chr.bed -t XYZ,ABC
      
        * Add -t {CHR1,CHR2} -l {COLORNAME1,COLORNAME2} to change the line color of plotted chromosomes. 
        
          *  ex: -f qm.bed -c chr.bed -t XYZ,ABC -l blue,red

        * Add -a {AnnotationBedFile.bed} to add an annotation track to all plot. 
      
          * ex : -f qm.bed -c chr.bed -t XYZ,ABC -l blue,red -a myexonstoplot.bed 
      
        * Add -d {/path/to/directory} to plot multiple lines from other bedfiles on top of the initial file given with -f flag.  
        
          * ex : -f qm.bed -c chr.bed -d /path/to/directory -t XYZ,ABC -l blue,red -a myexonstoplot.bed 
          
          * ex : -f qm.bed -c chr.bed -d /path/to/directory -i green -p .bed.gz -t XYZ,ABC -l blue,red -a myexonstoplot.bed 
   
      
        * Add -r y to create a png file output. 
	        
          * ex: -f qm.bed -c chr.bed -r y
      
        * Add -b y to create a box-plot window with box-plots labled to each gene name inputted with the -t flag. Each box plot will contain the average copy number and a blue, rhombus showing the average copy number for each respective gene.   
	        
          * ex: -f qm.bed -c chr.bed -d /path/to/directory -t XYZ,ABC -l blue,red -a myexonstoplot.bed -b y 
