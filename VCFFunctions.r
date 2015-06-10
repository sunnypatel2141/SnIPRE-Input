Read.VCF=function(file=i){
	#To load a VCF file set i="FILENAME"
	#then type "vcf=Read.VCF()
	x=readLines(i);has=grep(pattern="##", x)
	x=x[-has]
	x=strsplit(x,split="\t")
	return(x)	
	}