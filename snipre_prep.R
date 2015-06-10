# There are several parts to this script (eight to be exact). This R script uses output files from snipre_prep_bash.sh script, 
# so must be ran after it.  

# The following are the command-line arguments that must be provided with the script...
#	1) .vcf file
# 	2) starting bee column number (positive integer) of the vcf file (1-based indexing)
#	3) ending bee column number (positive integer) of the vcf file (1-based indexing)
#	4) SnpEff file (for Format Specifications, see Part 2)	
#	5) gff file (for Format Specifications, see Part 5)
#	6) Folder where output files from snipre_prep_bash.sh are located (path name)
#	7) nout : outgroup x 2
#	8) npop : population size x 2

# Usage : Rscript	snipre_prep.R		sample.vcf	10	17	syn_nonsyn_sample.eff	sample.gff	/media/data1/bombus/coverage/Groups		16	20

# The final output that it creates is the snipre input file required to run the SnIPRE program in R mode. Download : http://med.stanford.edu/bustamantelab/software.html

# The final output file has the following format.. 
	# Group					PR		FR		PS		FS		tsil	trepl		nout	npop		
	# XP_123456789.1		10		3		3 		10 		445		483			20		16	
	# XP_123456789.0		3		10		10		3 		544		384			20		16
	
	
# Sani Patel, Winter 2015
	

								#################################
#################################			PART 1 				######################################################################
								#################################

# This part of the script takes in a VCF file and user-specified Bee columns (hard-coded on lines 91 and 92). The format of the VCF file is as follows...
	#	V1    	 V2 	V3 V4 V5    V6   	V7 		V8  	V9								V10        	               V11						 	V12                    	    V13 ...
	#	Group1.1 8873  	.  A  G 	5360.35 PASS  	. 		GT:AD:DP:GQ:MLPSAC:MLPSAF:PL	1:0,15:15:99:1:1.00:623,0  1:0,13:13:99:1:1.00:549,0	1:0,11:11:99:1:1.00:470,0   1:0,9:9:99:1:1.00:371,0 ... 
		
#The script looks for fixed as well as polymorphic differences between bee columns (columns V10 to end), and creates a column in the output file
#denoting either an "F" or a "P" corresponding to "fixed" and "polymorphic" differences, respectively. The script also calculates the number of
#missing data columns, and outputs that number in the last column of output file. 

#The format of the output file is as follows...
	#CHROM 		POS REF ALT P/F 0
	#Group1.1 	123 A   C   F 	0								
					
print("Starting Part 1... ")
					
args <- commandArgs(TRUE)
						
source("VCFFunctions.r")

i= args[1]
		
vcf = Read.VCF(i)

#The function that calculates P/F column of the output file
myfunc_FP <- function(x, num1, num2) {
	
	#The following two lines check if the number is positive, else halts!!!
	stopifnot(num1 >= 0)
	stopifnot(num2 >= 0)
	
	len_one <- length(grep("^1[:]", x[num1:num2]))
	len_zero <- length(grep("^0[:]", x[num1:num2]))
	
	right = (num2 - num1) + 1
	what = ""
	
	if (len_one == 0 & len_zero == right) {
		what <- "F"
	}
	else if (len_zero == 0 & len_one == right) {
		what <- "F"
	}
	else if ((len_zero != 0 & len_one != right) |  (len_zero != right & len_one != 0)) {
		what <- "P"
	}
	else {
		what <= "?" # shouldn't exist. If you see a column with "?" in P/F column, check input file.
	}
	return(what) 
}

#The function that calculates missing data column in output file
myfunc_missing_data <- function(x, num1, num2) {
	
	#The following two lines check if the number is positive, else halts!!!
	stopifnot(num1 >= 0)
	stopifnot(num2 >= 0)
	
	#dot represents missing information
	len_dot <- length(grep("^[.:]", x[num1:num2]))
	return (len_dot)
}

# To myfunc_FP and myfunc_missing_data functions, pass in the column ranges of the bees you're interested in
# (bee columns in vcf file, let's say columns 10 to 19...)

# The column parameters must be positive.
startpos = as.numeric(args[2])
endpos = as.numeric(args[3])

F_or_P = sapply(vcf, function(x) myfunc_FP(x, startpos, endpos))
Missing_Data = sapply(vcf, function(x) myfunc_missing_data(x, startpos, endpos))

#The following two lines check if the number of rows are the same 
stopifnot(length(vcf) == length(Missing_Data))
stopifnot(length(vcf) == length(F_or_P))

#x[1] = CHROM
#x[2] = POS
#x[3] = REF
#x[4] = ALT
one = sapply(vcf, function(x) x[1])
two = sapply(vcf, function(x) x[2])
four = sapply(vcf, function(x) x[4])
five = sapply(vcf, function(x) x[5])

#append the columns vertically
aa <- cbind(one, two)
bb <- cbind(aa, four)
cc <- cbind(bb, five)
dd <- cbind(cc, F_or_P)
ee <- cbind(dd, Missing_Data)

write.table(ee, file="output1.txt", col.names=F, row.names=F, quote=F)								
							
print("... Ending Part 1")							

								#################################
################################# 			PART 2 				############################################################
								#################################
								
# This part of the script takes in the output file from last part (PART 1) as well as a SnpEFF file to output a column for either synonymous
# or non-synonymous mutation. The format of SnpEff file is as follows...  
	# Group1.1        46     C       T       SNP             3442.32 0               Gene_rna0       XP_012345678.1          rna0    Exon_Group1.1_337_468   5       NON_SYNONYMOUS_CODING   G/S     Ggc/Agc 180     1
	
#	** Note snpeff file must have column number 15 allocated for "SYNONYMOUS_CODING" or "NON_SYNONYMOUS_CODING" column. Also, no row should be missing "SYNONYMOUS_CODING" 
# 	or "NON_SYNONYMOUS_CODING" column. For this you can use Bash Command : 
# 	grep "SYNONYMOUS_CODING" SnpEff_File.eff > sample.eff 

#The format of the output file is as follows...
	# XP_003398050.1		123 A   G	SYNONYMOUS_CODING	P 	1
	# XP_003398050.1		124 A   G	SYNONYMOUS_CODING	P 	1	

print("Starting Part 2... ")	
	
#read in the SnpEff file ("\t" delimited)
input1 = args[4]
snpeff = read.table(input1, header=FALSE, sep="\t")

#read in output file from last part 
fixpoly <- read.table("output1.txt", header=FALSE) 

#To see which positions within a group are common between the two files, we combine Group and Pos column to make easier the use of "%in%" function. 
scaff = paste(fixpoly$V1, fixpoly$V2, sep="_")
scaff_snp = paste(snpeff$V1, snpeff$V2, sep="_")

matches = fixpoly[which(scaff %in% scaff_snp),]
new = snpeff[which(scaff_snp %in% scaff),]

new$scaff= paste(new$V1, new$V2, sep="_")
matches$scaff= paste(matches$V1, matches$V2, sep="_")

# Merge the two data frames based on commonality in Group+Pos column. 
test = new
test = test[which(test$scaff %in% matches$scaff),]
test = merge(new, matches, by="scaff")

# Output just the needed columns: Column 12 = Gene ID, Column 17 : SYN/NONSYN... 
final = test[,c(12,3,4,5,17,30,31)]

write.table(final, file="output2.txt", col.names=F, row.names=F, quote=F)

print("... Ending Part 2")


								#################################
################################# 			PART 3 				######################################################################
								#################################						

# This part of the script takes in the output file from last step and creates an output file with a column for type mutation. 
# The type mutation column (called "Mix") contains following values 1) PS, 2) PR, 3) FS, and 4) FR depending on:
	#PS = polymorphic + synonymous
	#FS = fixed + synonymous
	#PR = polymorphic + non-synonymous
	#FR = fixed + non-synonymous
	
# The output file format ...
	# 	Group 			Pos 		Ref Alt Syn 					P/F Missing Mix
	#	XP_012345678.1 	10022757 	G 	A 	SYNONYMOUS_CODING 		F 	0 		FS
	# 	XP_012345678.2 	10044759 	C 	A 	SYNONYMOUS_CODING 		P 	0 		PS
	#	XP_012345678.3 	10055780	G	T	NON_SYNONYMOUS_CODING	F	0		FR 

print("Starting Part 3... ")	
	
z <- read.table("output2.txt", head=FALSE)

# This function aims to combine the P/F and Synonymous/ Non-synonymous column of the input file and 
# create a column with four possible values: 1) FS, 2) FR, 3) PS, and 4) PR. 
myfunc <- function(x, num1, num2) {

	value = ""
	syn = length(grep("^SYNONYMOUS_CODING", x[num1]))
	non = length(grep("^NON_SYNONYMOUS_CODING", x[num1]))
	
	polymo = length(grep("P", x[num2]))
	fixed = length(grep("F", x[num2]))
	
	if(syn > 0) {
		if (polymo > 0) {
			value = "PS"
		} else if (fixed > 0) {
			value = "FS"
		}	
	} else if (non > 0) {
		if (polymo > 0) {
			value = "PR"
		} else if (fixed > 0) {
			value = "FR"
		}		
	}
	return (value)	
}

# The Parameter "5" is column for "SYN..." or "NONSYN..." and parameter "6" is column for P/F (supplied to function `myfunc')
y <- apply(z, 1, function(x) myfunc(x, 5, 6))
z$Mix <- y

write.table(z, file="output3.txt", col.names=F, row.names=F, quote=F)

print("... Ending Part 3")


								#################################
################################# 			PART 4 				######################################################################
								#################################


# This part of the script takes in output file from last part (PART 4) and aggregates for each unique gene the number of 
# PS, FS, FR, and PR sites using "length" function. 

# The output file created has the following format...
	# 	Group 			PR 	FR 	PS 	FS
	#	XP_012345678.1 	1 	0 	3 	1
	#	XP_012345678.2 	1 	0 	3 	1

print("Starting Part 4... ")	
	
#read file obtained from last Part	
x <- read.table("output3.txt", head=TRUE)

colnames(x) = c("Group", "Pos", "Ref", "Alt", "Syn", "P.F", "Missing", "Mix")

#obtain unique groups and make it a matrix
q = unique(x$Group)
q1 = as.matrix(q)

# This for loop creates a data frame subset for each unique group name in input file, and 
# then it calculates the number of occurrence for each unique value in "Mix" column 
# (how many "PS", "PR", "FS", and "FR").
for (i in 1:length(q)) {

	a <- x[x$Group == q[i],]
	
	pr = length(which(a$Mix == "PR"))
	ps = length(which(a$Mix == "PS"))
	fs = length(which(a$Mix == "FS"))
	fr = length(which(a$Mix == "FR"))
		
	#horizontal binding	
	aa = cbind(q1[i], pr)	
	bb = cbind(aa, fr)
	cc = cbind(bb, ps)
	dd = cbind(cc, fs)
			
	write.table(dd, file = "output4.txt", quote = F, append = TRUE, row.names=F, col.names=F)
}

print("... Ending Part 4")


								#################################			
################################# 			PART 5 				######################################################################
								#################################

# This part of the script takes in a gff file of the following format... 
	# Group1.1        RefSeq  CDS     1       75      .       -       0       "ID=cds0;Name=XP_012345678.1;Parent=rna0;Dbxref=GeneID:123456,Genbank:XP_012345678.1;gbkey=CDS;product=random odorant receptor 63a-like;protein_id=XP_012345678.1"
	# Group1.2        RefSeq  CDS     114     694     .       -       1       "ID=cds0;Name=XP_012345678.2;Parent=rna0;Dbxref=GeneID:123456,Genbank:XP_012345678.2;gbkey=CDS;product=amazing odorant receptor 61a-like;protein_id=XP_012345678.2"

# ... and calculate the length of each XP (subtract columns 4 from 5 and then aggregate for each unique gene)

# The output file has the following format (no header)...
	# XP_012345678.1 972
	# XP_012345678.2 1716
	# XP_012345678.3 1044

print("Starting Part 5... ")	
	
# Read the gff file 	
input2 = args[5]
bter = read.table(input2, head=FALSE)

# Get right format for XP
z = strsplit(as.character(bter$V9), ";", fixed=TRUE)
t = sapply(z, "[[", 2)

#calculate length (column5 - column4 + 1)
bter$xp = substr(t, 6, 19)
bter$sum = bter$V5 - bter$V4 + 1
bter = bter[,c(10, 11)]
 
n = aggregate(bter$sum, list(bter$xp), sum)
 
write.table(n, file="output5.txt", col.names=F, row.names=F, quote=F)	
 
print("... Ending Part 5")

 
								################################# 
################################# 			PART 6 				######################################################################
								################################# 
 
# This part of the script takes in the output file from previous part (PART 5) as input file, as well as the output file from Part 4
# to create an output file of the following format... 
	# Group 			PR	FR	PS	FS	sum
	# XP_012345678.1 	1	0 	0 	2	483
	# XP_012345678.2 	3 	4 	0 	0 	484
	# XP_012345678.3 	1 	0 	7 	0 	485
	
# Below shown snippet of code "merges" the two files based on Group column ("sum" column refers to length of each gene). 
	
print("Starting Part 6... ")	
	
x <- read.table("output4.txt")
n <- read.table("output5.txt", head=FALSE)

colnames(n) = c("Group", "sum")
colnames(x) = c("Group", "PR", "FR", "PS", "FS")

z = merge(x, n, by="Group")

write.table(z, file="output6.txt", col.names=T, row.names=F, quote=F)

print("... Ending Part 6")


								#################################
################################# 			PART 7				######################################################################
								#################################

# This part of the script aims to use the output files created from Bash script snipre_prep_bash.sh and use sample.gff file to 
# output a "tsil" column for SnIPRE. Here, the cutoff for "tsil" is average_depth >= 8.0 (can be modified on line 370)

# Format : sample.gff ...
	# Group1.1        RefSeq  CDS     7       45      .       -       0       "ID=cds0;Name=XP_012345678.1;Parent=rna0;Dbxref=GeneID:12345,Genbank:XP_012345678.1;gbkey=CDS;product=amazing odorant thing 6-like;protein_id=XP_012345678.1"

# Format : Output File from snipre_prep_bash.sh ...	
	# Group1.1:129    80      8.00    2       11      4       12      2       12      11
	# Group1.1:130    80      8.00    2       11      4       12      2       12      11
	
	
# The format of the output file is as follows (second column is "tsil")...
	# XP_012345678.1 534
	# XP_012345678.2 438
	# XP_012345678.3 4755
	
# This part of the script is the most time-consuming since we have an for loop within a for loop. Polynomial running time. 

print("Starting Part 7... ")

# Read the gff file
input2 = args[5]
bter = read.table(input2, head=FALSE)

z = strsplit(as.character(bter$V9), ";", fixed=TRUE)
t = sapply(z, "[[", 2)

bter$xp = substr(t, 6, 19)

x <- bter[,c(1, 4, 5, 10)]
colnames(x) = c("Group", "start", "end", "xp")

# Select each file made from snipre_prep_bash.txt and use "outer" function to count how many positions have Average_Depth >= 8.

# path to files created in snipre_prep_bash.sh
# make sure to not include "/" at last index
path = args[6] 

file.names = dir(path, pattern=".txt")

# For each GroupX.txt file, the following code creates data frame with Average_Depth >= 8.0. Then for each unique Gene ID in the 
# sample.gff file, determine the start and end of the gene. Within this gene region, let's say from 122 to 254 in Group1.1, 
# count how many positions from the "Group1.1.txt" file (created from snipre_prep_bash.sh) fall in it. Lastly, Use aggregate 
# function to sum up the "tsil" count. 

for (i in 1:length(file.names)) {

	fullname = paste(path, file.names[i], sep="")
		
	y = read.table(fullname)

	y = y[,c(1,3)]
	
	z = strsplit(as.character(y$V1), ":", fixed=TRUE)

	z1 = lapply(z, function(x) x[2])
	z2 = lapply(z, function(x) x[1])

	y$Group = z2
	y$pos = z1

	y = y[,c(3,4,2)]
	y1 = subset(y, y[3] >= 8.0)

	grp = unlist(unique(y1$Group))
	print(grp)	
			
	tp1 = x[which(x$Group==grp),]; yp1=y1[which(y1$Group==grp),]
	xps = unique(tp1$xp)
	
	for (xp in xps) {
	
			print(xp)
	
			tp2 = tp1[which(tp1$xp == xp),]
	
			#blah = positions in Group file greater than start position in gff file 
			blah=outer(as.numeric(unlist(tp2$start)), as.numeric(unlist(yp1$pos)), "<=") 
			
			#blah1 = positions in Group file less than the end position in gff file 
			blah1=outer(as.numeric(unlist(tp2$end)), as.numeric(unlist(yp1$pos)), ">=") 
			
			blah=(which(blah1=="TRUE" & blah=="TRUE", arr.ind=T))
			blah=data.frame(blah)
			
			#aggregate based on Gene ID (nrow(blah) > 0 used to make sure use of aggregate function is valid)
			if (nrow(blah) > 0) {
				test=with(blah, aggregate(col, by=list(row), length))
				tp2$Group.1=seq(1,nrow(tp2),1)
				test=merge(test,tp2,by="Group.1",all.y=T)
				
				final = aggregate(x ~ xp, test, FUN = sum)	
				write.table(final, file="output7.txt", row.names=F, col.names=F, quote=F, append=T)
			}
	}
}

print("... Ending Part 7")

								#################################
################################# 			PART 8				######################################################################
								#################################	

# This part of the script "merges" files from Part 6 and Part 7 to obtain "tsil" column for each gene along with the rest of the file.
# The missing value for "nout" (outgroup) and "npop" (population size) column is calculated based on how many bam files were used in vcf file. 

# This final output file has the following format:
	# Group				PR		FR		PS		FS		tsil	trepl		nout	npop		
	# XP_012345678.1	10		3		3 		10 		445		483			20		16
	
print("Starting Part 8... ")	

x = read.table("output7.txt")
y = read.table("output6.txt", head=TRUE)

colnames(x) = c("Group", "tsil")
z = merge(x, y, by="Group")

# Obtain Gene ID, PR, FR, PS, FS, Tsil, and TrepL columns 
z1 = z[,c(1,3,4,5,6,2,7)]

#nout : outgroup for A.cerana 1 bee x 2 
#npop : population size for A.scutellata 11 bees x 2 
z1$nout = args[7]
z1$npop = args[8]

colnames(z1) = c("geneID", "PR", "FR", "PS", "FS", "Tsil", "Trepl", "nout", "npop")
write.table(z1, file="snipre_input.txt", row.names=F, col.names=T, quote=F)

print("... The End")