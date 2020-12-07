Preprocessing raw NGS data I
-----------------------------------------------------------------------

### Description

Preprocessing raw NGS data I is built to process raw reads with one adaptor sequence. Thereby it is focused on mapping miRNA and possible isomiR's by constructing an miRNA isomiR DB. For comparison samples can be grouped (e.g. healthy and diseased) and results will be pictured through an extensive number of plots.  

The rules of this snakemake wrapper repository do following steps:

***1. Adaptor Trimming***

Sequence of the adaptor is trimmed by btrim32-static (software is in tools folder). Therefore the adaptor sequence is stored in a text file (refs/linker3.txt). 
If another adaptor sequence is used than here provided textfile has to be adjusted. 


***2. Quality Control of the reads***

Quality control of the sequenced reads is performed with fastqc. Report is summarized as .html an .zip file. FastQC is embedded in the environment of the snakmake wrapper repository.


***3. Mapping reads to annotated references of rRNA, tRNA, snoRNA and snRNA***

Mapping is carried out with bowtie. Therefore indexfiles for mature rRNA, tRNA, snoRNA and snRNA are already included in the refs/RNAcentral folder. 
It is possible to use any other reference files for mapping by creating equivalent index files with bowtie-built. The corresponding snakemake rule has to be adjusted by path, therefore.  

***4. Construction of an miRNA isomir DB and mapping reads to this one***

 Mapping to miRNA is extended by creating first an isomir DB from mature miRNA's. Therefore the isomiRROR pipeline is implemented.

***5. Vizualisation of reads content and mapping results***

Vizualisation of the results is performed by a R script. An extensive number of plots gives an overview about mapping distributions. Thereby samples can be grouped together by creating an TAB seperated Groupnames.txt file with partial match to sample names (Groupnames.txt has to be located in the /Seq_data_and_Output folder/).





Installing and running "Preprocessing raw NGS data I"
---------------------------------------------

Requirements to run this code is conda and bowtie, everything else will be delivered by the repository.

### Step 1: Installing bowtie and conda
First install bowtie globally on your system. The following link will give you access for necessary files:
	
	http://bowtie-bio.sourceforge.net/tutorial.shtml

if you work on a linux system following command will install bowtie:

	sudo apt-get install bowtie

Next you need to download and install miniconda3:
for linux
    
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

for mac os

    curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

### Step 2: Clone the workflow
Clone the worflow to the folder where you run your pipelines

	git clone https://gitlab.lrz.de/ImmunoPhysio_private/Preprocessing_raw_NGS_data_I.git

### Step 3: Create an environment based on environment.yaml

    cd Preprocessing_raw_NGS_data_I
    conda update conda
    conda env create --file environment.yaml

### Step 4: Supply reference and input files
Preprocessing raw NGS data I needs reference files for mapping reads and for construction of the miRNA isomir DB. Index files for the mapping of rRNA, tRNA, snRNA and snoRNA can be built with the bowtie-built tool and have to be located in: /tools/bowtie_references. Please make sure using reference files for the desired species. Index files for the construction of miRNA isomirs are built by running code. Therefor a fasta file of mature miRNA is provided (refs/mature.fa). This mature miRNA file can be customized (has to be located in the folder /refs/). 
Furthermore the header line of each mature miRNA in that file should only contain the name of the miRNA. Everything after the first space has to be removed. Provide hairpin sequences of miRNAs as well (hairpin.fa file must be located in :/refs/hairpin.fa folder, header line of each sequence must be the same as for mature.fa file). Make sure mature.fa and hairpin.fa files contain only sequences of the species of interest.
Sample fastq files have to be located in the /Seq_data_and_Output/ folder togethere with the Groupnames.txt file. Groupnames must be tab seperated.     

### Step 5: Complete config.yaml, samples.csv and plotting.R with the missing information
In order to run the pipeline you will need to complete the config.yaml file and the samples.csv file. 

1. config.yaml – reference files and isomiR length limits
The config.yaml stores the information of where to find the fasta files for the isomiR reference generation as well the minimum and maximum length of isomiRs to consider for alignment.

References
* reference_folder: refs or PATH/TO/LOCAL/FOLDER
* mature_db_file: mature
* hairpin_ref_file: hairpin

By default isomiRROR will look  for the reference sequences in the /refs subfolder but users can also specify a different folder on their workstation. In both cases the name of the mature miRNA fasta file (mature_db_file) and their corresponding precursor hairpin fasta file (hairpin_ref_file) need to be given without the .fa ending.

Input variables
* min_length: 16
* max_length: 30

Since sequences below a certain length are very prone to align to multiple locations and represent most likely degraded RNA from other RNA subspecies, it is recommended to limit the range of isomiRs before alignment. Likewise the number of isomiRs above a certain lenths is rapidly declining and computational resources and time can be saved. Out of personal experience I recommend using a minimum length (min_length) of 16 or higher  and a maximum length (max_length) of around 30, but users are encouraged to experiment and find a fitting length range for their own experiments. Since isomiR lengths are stored in the ID found in the final readcount file as well, filtering of certain lengths post isomiRROR is possible as well. 

2. samples.csv -

*  The samples.csv simply holds a list of the sample names that should be alignend to the isomiR reference with each sample on a new line. Sample names should be given without the .fastq ending. First line should always have the string "samples"

    samples  
    samplename_one  
    samplename_two  
    ...

3. plotting.R 

*  The plotting.R file is to visualize the results through plots and is located in the /src/ folder. The working directory has to be adjusted before the script can run. Set as working directory the /Seq_Data_and_Output/ folder. Here also the Groupnames.txt file should be located. This file is necessary for the R-script and should include the grounames which are TAB seperated.

### Step 6: Run "Preprocessing raw NGS data I"!
Simply run :

    cd Preprocessing_raw_NGS_data_I
    snakemake --use-conda

together with your choice of parametery (i.e --cores) in the root folder of your isomiRROR folder and it will be running from beginning to end.

### Output files

Preprocessing raw NGS data I will generate different output folders located in /Seq_Data_and_Output/.
The FastQC folder includes the qualtiy report of the sequenced reads. The Plots folder includes all plots generated by the R-script (plotting.R) and text files which summarise the necessary information for the plots. The Readcounts folder includes all Readcount files generated by Preprocessing raw NGS data I.
The Trimming folder includes files with information about the adaptor trimming part. The Unmapped folder includes files which summarise the unmapped reads.
Beside this the miRNA isomiR part will generate a new isomiR reference fasta file, an alignment index corresponding to the mapping tool of your choice splitted by lengths, a text file containing the possible multireads (isomiRs with identical sequences but originating from different mature miRNA sequences) and of course a readcount file containing a table with isomiR ID in the first column and corresponding readcounts for each sample in the following.