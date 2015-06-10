#SnIPRE PREP SCRIPTS

The objective of the two scripts (snipre_prep_bash.sh and snipre_prep.R) is to create an input file for SnIPRE program (http://med.stanford.edu/bustamantelab/software.html).

##snipre_prep_bash.sh

This script takes in as input all the .bam files used to create the .vcf file, a fasta file and a gff file. The output is a .coverage file with the
following format... 

Group:Pos	Total_Depth	Average_Depth	Depth_one	Depth_two	Depth_three	Depth_four	Depth_five

Group1.1:2	2	0.20    0      1	0       0	1

The bash script creates a .coverage file using GATK command, and then splits up the .coverage file based on unique Group name. 

Note : This bash script must be run prior to snipre_prep.R as its output files are used by the R script.

The .jar file needed for this bash script (GenomeAnalysisTK.jar) can be downloaded from https://www.broadinstitute.org/gatk/.


##snipre_prep.R

This R script takes in eight command-line arguments to create an output file that can be used by the SnIPRE program. The eight command-line arguments
are (in order)...

	1) .vcf file
	2) starting bee column number (positive integer) of the vcf file (1-based indexing)
	3) ending bee column number (positive integer) of the vcf file (1-based indexing)
	4) SnpEff file (for Format Specifications, see Part 2)	
	5) gff file (for Format Specifications, see Part 5)
	6) Folder where output files from snipre_prep_bash.sh are located (path name)
	7) nout : outgroup x 2 (positive integer)
	8) npop : population size x 2 (positive integer)
	
The output file created from this R script has the following format...
	
Group					PR		FR		PS		FS		tsil	trepl		nout	npop		
XP_003398050.1			10		3		3 		10 		445		483			20		16	
XP_003398051.00			3		10		10		3 		544		384			20		16
	
### Usage of the R Script:
		Rscript		snipre_prep.R		sample.vcf	10	17	syn_nonsyn_sample.eff	sample.gff	/coverage/Groups		16	20

				
##### Sani Patel, Winter 2015
