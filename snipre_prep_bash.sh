# In this script out goal is to take the bam files used to create the vcf file and run Depth of Coverage (feature of GATK) command on it. The second part
# of the script creates an output file for each unique gene due to the huge size of the output file created from GATK command.

# The output file has the following format (for Group1.1 let's say)...
	# Group:Pos			Total_Depth	Average_Depth_Sample	Depth_for_one	Depth_for_two	Depth_for_three	Depth_for_four	Depth_for_five	Depth_for_six	Depth_for_seven
	# Group1.1:2    	2       	0.20    				0       		1       		0       		0       		1       		0       		0
	# Group1.1:4    	4       	0.40    				1       		1       		1       		1       		1       		1       		1
	
# For the below bash command, you need GenomeAnalysisTK.jar file. 
# sample.fa = Fasta file for the genome
# output file is called sample.coverage in the example  
java -Xmx2g -jar GenomeAnalysisTK.jar -R sample.fa -T DepthOfCoverage  -I one.bam -I two.bam -I three.bam -I four.bam -I five.bam -I six.bam -I seven.bam --out sample.coverage --minBaseQuality 20

#To run it without java command (may get memory usage issue):
	# gatk -R sample.fa -T DepthOfCoverage  -I one.bam -I two.bam -I three.bam -I four.bam -I five.bam -I six.bam -I seven.bam --out sample.coverage --minBaseQuality 20

# Now break down the above output file based on unique Group name.  
awk '{ print $1 }' sample.gff | uniq > uniq_groups.txt
while read line; do grep $line: sample.coverage >> $line.txt ; done < uniq_groups.txt