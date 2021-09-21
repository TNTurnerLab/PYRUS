<img align="right" src="https://user-images.githubusercontent.com/77067493/129043760-26d57868-5f15-46bb-a00e-e82f8c50fa9c.png" height="20%" width="20%" >


Maintainer: Elvisa Mehinovic

Laboratory of Dr. Tychele N. Turner, Ph.D.

Washington University in St. Louis

<h1>PYRUS</h1>
PYRUS is a tool for plotting copy number estimate data, from an individual, for user-specified regions of the genome. It has several options including plotting other individuals in the same region, plotting an annotation track, and writing out specific regions where the individuals have a copy number below or above given values. The input to the tool are bgzipped and tabix indexed bed files, which enables rapid plotting of the data.

### ADD IMAGE HERE ###

### Get Started
Clone the PYRUS repository:

	git clone git@github.com:TNTurnerLab/PYRUS.git

 ## Running Pyrus Local Requirements
Users have the option to run PYRUS using Docker, or via the command line using Rscript.
 
#### The following R packages are required to run PYRUS:

* optparse
* dplyr
* seqminer
* data.table

Running on local command line with R:

	Rscript /path/to/cloned/Repo/PYRUS.R -f /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}
	
## Running Pyrus with Docker:

Pull the following Docker image:

	docker pull tnturnerlab/pyrus
	
Run Docker image :

	docker run -v /path/to/working/directory:/path/to/working/directory/ -w /path/to/cloned/Repo/ -it tnturnerlab/pyrus Rscript /path/to/cloned/Repo/PYRUS/PYRUS.R -f /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}

To run on LSF with Docker:
	
	bsub -G name -q name -R 'rusage[mem=20GB]' -a 'docker(tnturnerlab/pyrus)' LC_ALL=C.UTF-8  Rscript /path/to/cloned/Repo/PYRUS/PYRUS.R /path/to/input/file/My_Bed_File.bed.gz -c /path/to/input/file/My_Chr_file.bed.gz {INSERT ANY MORE FLAGS IF NECESSARY}



#### File Requirements

* bgzipped and tabix-indexed bed files are the only acceptable file types for this script. For more information about how to prepare these files see [bgzip](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html)

An example of the input bed file prior to bgzipping and tabix indexins is shown below:

```
chr1	0	54484	1.394791
chr1	54484	60739	1.258411
chr1	60739	68902	1.155175
chr1	68902	82642	2.633766
```

## Required Flags :

-f	:

Input the initial/main bgzipped and tabix indexed copy number estimate bed file for the individual
  
   * default=NULL	
   * Ex: -f bedFile.bed.gz

-c	:

Input a bgzipped and tabix indexed chromosome coordinate file containing the regions for plotting with PYRUS.
    
   * default=NULL	
   * Ex: myChrbedfile.bed.gz

## Optional Flags:

-o	:
  
 Enables the user to combine all plots into one PDF with user-specified name, the default setting is to plot each coordinate to its own PDF file
  
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

To change the sample name in the legend for the file specified with the `-f` argument. The default name will always be the first 15 characters of the filename.
    
  * default=NULL	
  * Ex: -y MyNewName

-u:

Change the ylim maximum value for the plot window. The default is a copy number of 6. 
    
   * default=6	
   * Ex: -u 7
    
-d	: 

Input the directory path to bgzipped and tabix indexed bed files from additional samples. This option is used when the user wants to plot multiple samples. Note: the program will automatically search for .tbi files and requires them to run. 
    
   * default=NULL	
   * Ex: -d /path/to/a/directory/with/bed/file
   	
	-s: 
   		Flag used if user would like to plot every file in the given directory 
		against the file denoted by the -f flag. Must be used along with the -d flag.
 			default=NULL	Ex: -d /path/to/dir/ -s y

	-p: 

  		Input an ending pattern for files found inside the directory in "-d" flag. Must be the 
		ending pattern of the file. Can also be used to specify a specific file from a directory.
    
    			default=.bed.gz	Ex: -d /path/to/dir/ -p .t1.cn.bed.gz	Ex: -d /path/to/dir/ -p exact.t1.cn.bed.gz 

-a	: 

Flag used if user would like to add an annotation track (e.g., exons) using a bgzipped and tabix indexed bed file. The default color settings for the annotation track included a red fill, black border, and Exon as the name,

* default=NULL	
* ex: -a exonbedfile.bed
	
	  -n: 

  		Flag must be used along with "-a" to change the fill,boarder, and name of annotation track. 
  
    			default=NULL	Ex: -a exonbedfile.bed.gz -n "color,color,name"
    
-b:

This flag will create a box-plot window in the top right corner of the existing line plot. The box plot will display the individual, marked with the -f flag, with a blue dot. This will automatically turn on if there are more than 10 files when multiplotting.

* default=NULL	
* Ex: -b y 

-t:

This flag takes two names, found in the chromosome coordinate file, and plots them side by side on a single plot. The two names must be from regions on the same chromosome. Data from the fourth column of the chromosome bed file must be used, typically it is a name/identifier. The default colors of these lines are blue and red, but can be modified with the -l flag.
    
* default=NULL	
* Ex: -t BARX1,BARX1-DT

-g: 

Threshold values for printing out names that are below the first value or above the second value. Note that if the user specificies the individual is male (-x M), the threshold for the X and Y chromosomes will be 1 less than the values indicated with this option. 

* default="1.3,2.7"	
* Ex: -g 2,3
	
-x: 

If the sex is known for the individual in the file specified with the -f flag enter the value with this option. Default is F. This flag only has to be called if user wants to change sex to M.
    
* default=F	
* Ex: -x M

-v: 

Plot only the regions whose threshold values are below the first value of -g flag, or above the second of the -g flag. Note that the X and Y chromosomes will have a threshold one less than the values given if the sex is male. 
    
 * default=NULL	
 * Ex: -v y
