<img align="right" src="https://user-images.githubusercontent.com/77067493/129043760-26d57868-5f15-46bb-a00e-e82f8c50fa9c.png" height="20%" width="20%" >

<h1>PYRUS</h1>

PYRUS is a plotting tool that uses tabix files to create line graphs from bed file data. Utilization of PYRUS includes having a quick plotter for copy number variation within a given chromosomal range of a referenced bed file. 


## Running Pyrus

### Docker:
Clone the PYRUS repository:

	git clone git@github.com:TNTurnerLab/PYRUS.git

Pull the following Docker image:

	docker pull tnturnerlab/pyrus
	
Run Docker image :

	docker run -v /path/to/working/directory:/path/to/working/directory/ -w /path/to/cloned/Repo/ -it tnturnerlab/pyrus Rscript /path/to/cloned/Repo/PYRUS/PYRUS.R -f /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}

To run on LSF with Docker:
	
	bsub -G name -q name -R 'rusage[mem=20GB]' -a 'docker(tnturnerlab/pyrus)' LC_ALL=C.UTF-8  Rscript /path/to/cloned/Repo/PYRUS/PYRUS.R /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}

### Runing on local command line with R:

	Rscript /path/to/cloned/Repo/PYRUS.R -f /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}
 
 ## Requirements
Users have the option to run PYRUS with the given Docker image, or via a command line if there is a configured R application.
 
#### R studio must contain the following packages:

* optparse
* dplyr
* seqminer

#### File Requirements

* Bed files are the only acceptable file types for this script







## Required Flags :

-f	:

Input the initial/main bed file
  
    	default=NULL	Ex: -f BedFile.bed.gz

-c	:

Input a chromosome coordinates bed file.
    
    	default=NULL	Ex: myChrbedfile.bed.gz



## Optional Flags:

-o	:
  
 Allows user to name, and place all outputputted plots into one PDF, default is plots each coordinate to its own PDF file
  
   * default=NULL	
   * ex: -o multiplePagePDFFileofMyPlots.pdf

-r	:

Flag to change the outputs from a pdf file to a low resolution png file.
    
   * default=NULL	
   * Ex: -r y 
    
-l	:

  Input a color for the line(s) created from the file given to the -f flag. 
    
   * default=blue,red	
   * Ex: -l blue,red

-y	: 

To change the legend name of the file -f when plotting. The default name will always be the first 15 characters of the file's name.
    
  * default=NULL	
  * Ex: -y MyNewName


-u:

Change the ylim max values plot window. 
    
   * default=6	
   * Ex: -u 7
    
-d	: 
Input the directory pathway to bed files, if multiple bed files are wanted to be plotted onto one plot. The program will automatically search for .tbi files and require them to run. 
    
   * default=NULL	
   * Ex: -d /path/to/a/directory/with/bed/file
   	
	-s: 
   		Flag used if user would like to plot every file in the given directory 
		against the file denoted by the -f flag. Must be used along with the -d flag.
 			default=NULL	Ex: -d /path/to/dir/ -s y

	-p: 

  		Input an ending pattern for files found inside the directory in "-d" flag. Must be the 
		endding pattern of the file. Can also be used to specify a specific file from a directory.
    
    			default=.bed.gz	Ex: -d /path/to/dir/ -p .t1.cn.bed.gz	Ex: -d /path/to/dir/ -p exact.t1.cn.bed.gz 

-a	: 

Flag used if user would like to add an annotation track using a bedfile. The default color settings or the annotation track is a red fill, black boarder, and Exon as the name,

* default=NULL	
* ex: -a exonbedfile.bed
 	
		-n: 

  		Flag must be used along with "-a" to change the fill,boarder, and name of annotation track. 
  
    
    			default=NULL	Ex: -a exonbedfile.bed.gz -n "color,color,name"
    
-b:

This flag will create a box-plot window in the top right corner of the existing line plot. The box plot will display the individual, marked with the -f flag, with a blue dot. This will automattically turn on if there are more than 10 files when multiplotting.

* default=NULL	
* Ex: -b y 

-t:

This flag takes two chromosomes, found in the chromosome bed file, and plots them side by side on a single plot. Data from the fourth column of the chromosome bed file must be used, typically it is a name/identifier. Default color of these lines are blue and red, but can be modified with the -l flag.
    
* default=NULL	
* Ex: -t BARX1,BARX1-DT

-g: 

Threshold values for printing out names that are with below the 1st value, or above the second. Note that X and Y chromosomes will be 1 less than the values given if the sex is known. 

* default="1.3,2.7"	
* Ex: -g 2,3

	
-x: 

If sex is known of the file called in the -f flag. Default is F. This flag only has to be called if user wants to change sex to M.
    
* default=F	
* Ex: -x M

-v: 

 Plot only those whose threshold values are below the 1st value of -g flag, or above the second of the -g flag.Note that X and Y chromosomes will be 1 less than the values given if the sex is known. 
    
 * default=NULL	
 * Ex: -v y

  	
  
