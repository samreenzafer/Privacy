########################################################################################################
########################################################################################################
#
#  Encryption tool to run association testing on collaborated data
#
########################################################################################################
########################################################################################################

We have developed a tool which performs association testing on collaborative data from multiple investigators.
Since patient privacy and confidentiality are a concern, we use public key security paradigm to encrypt data from investigators before collaborating them; following this association testing can be done and results can be decrypted by investigators.

The application is currently under development and not complete and tested.

Given multiple collaborators ("users"), each user is required to register with their username and password. Following registration by all users, a public key is generated and shared with all users. Users must submit this public key along with their data, to enable encryption and data processing, to make their data ready (in appropriate format ) for association testing. Following this, processed-encrypted data is sent to a trusted third party, for accumulating data, and running the required analysis. The results may be shared with all collaborators thereafter. These steps and format of input data are detailed below.




########################################################
##                                                    ##
## 1. Input data format : VCF files    	              ##
##                                                    ##
########################################################

1. It is necessary that the data (case or control) submitted by any user must be in Variant Call Format (VCF) and have case and control samples in different VCF files. The VCF files must be valid VCF files.

2. Format for a valid VCF file can be seen on http://www.1000genomes.org/node/101

3. Validation of VCF files can be done using tools like vcftools_0.1.6 ( http://vcftools.sourceforge.net/perl_module.html ).

########################################################
##   TODO                                             ##
## 2. Output format		                      ##
##                                                    ##
########################################################



########################################################
##                                                    ##
## 3. Public Key Generation			      ##
#	(I) User Registeration 			      ##
##                                                    ##
########################################################

Main code file 
	**  engine.pl
Supporting code file 
	**  decidekey1.pl

1. Each user must register with their username and password. Run the engine.pl script with flag "--register" and required fields USERNAME and PASSWORD;
	** Run the engine as ./engine.PL --register --user=USERNAME --pwd=PASSWORD

2. OUTPUT: A temporary key is generated for each user and stored in a file named USERNAME_keypart1 under directory "users".



########################################################
##                                                    ##
## 3. Public Key Generation			      ##
#	(II) create public key 			      ##
##                                                    ##
########################################################

Main code file
	**  engine.pl
Supporting code file 
	**  decidekey2.pl

1. Any user may fire the key creation step. The public key is stored in a file "key" under subdirectory "users" under directory containing code files. If the "key" file already exists then firing the key creation step will fail. This prevents the key creating step to be performed multiple times. Run the engine.pl with flag "--genkey" ;
	** Run the engine as  ./engine.PL --genkey --user=USERNAME

2. OUTPUT: A public key stored in file  "key" under subdirectory "users". i.e. "code_directory/users/key"


########################################################
##                                                    ##
## 4. Data submission by user			      ##
#	(I) Annotation of input VCF		      ##
##                                                    ##
########################################################
http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/inputs.shtml

1. Input VCF files should be annotated by vcfCodingSNPS v1.5 ( http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/index.shtml ).
2. Download and Install the vcfCodingSNPS v1.5 tool using steps indicated http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/Installation.shtml
3. The tool requires a reference genome and a list of genes downloaded from the UCSC genome browser (http://genome.ucsc.edu/) for doing the annotation. Some reference genomes and gene lists for humans are available in the tool itself and maybe used drectly. The user may download other reference genomes and genelists as per their need. 
	* Reference genome must be in fasta format
	* Gene list must be in GenePred format (http://genome.ucsc.edu/FAQ/FAQformat#format9) described in brief;
		The first11 fields of gene file are required tab delimited fields and must be put in the order as following:
			1st	string	name	"Name of gene"
			2nd	string	chrom	"Chromosome name"
			3rd	char[1]	strand	"+ or - for strand"
			4th	uint	txStart	"Transcription start position"
			5th	uint	txEnd	 "Transcription end position"
			6th	uint	cdsStart	"Coding region start"
			7th	uint	cdsEnd	"Coding region end"
			8th	uint	exonCount	"Number of exons"
			9th	uint[exonCount]	exonStarts	"Exon start positions"
			10th	uint[exonCount]	exonEnds	"Exon end positions"
			11th	string	gene symbol	"Standard gene symbol"
		Note: the 11th field is a mandatory field for running vcfCodingSnps. 
		
		Here is an example of input gene file headlines:

			##ucscname  chrom  strand  txStart  txEnd  cdsStart  cdsEnd  exonCount  exonStarts  exonEnds  genename
			uc001aaa.2  chr1  +  1115  4121  1115  1115  3  1115,2475,3083, 2090,2584,4121,  BC032353
			uc009vip.1  chr1  +  1115  4272  1115  1115  2  1115,2475,  2090,4272,  AX748260

	* For detailed description of input files (vcf, referencegenome and gene list) to vcfCodingSNPS tool, refer to http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/inputs.shtml

4. VCF annotation can be fired using the "annotate_script.sh" as. 
	* qsub -o out -e err annotate_script.sh PATH_INPUT_FILE  OUTPUT_DIRECTORY RUN_NAME_PREFIX
	* Default reference genome we use is : NCBI released B36 reference genome file	(this can be changed in the script)
	* Default gene list we use is : UCSC known genes for B36 (this can be changed in the script)

NOTE: TODO: Add the annotation code inside the main engine.pl .

########################################################
##                                                    ##
## 4. Data submission by user			      ##
#	(I) Encryption of Annotated VCFs	      ##
##                                                    ##
########################################################

Main code file
	**  engine.pl
Supporting code file 
	**  encrypt.pl

1. A user submits annotated vcf files for encryption. Required inputs to the engine.pl are described below:
	(1) --user = Username 
	(2) --keyfile = Path of a file containing the pulic key. 
	(3) --genelist = a genelist file, exactly the same that was used during annotation.
	(4) --weight = a tab delimited file, with atleast one line, where each line has two columns , described below:
			1st column = name of functional category (any categories that are handled by vcfCodingSNPs tools. A list of these can be seen at http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/Tutorial.shtml )
			2nd column = a floating point weight between 0 and 1 that is to be assigned to the particular category.

		An example weight file looks like "weights.txt" and is shown below. This contains the exhaustive list of categories=>
			5'UTR	0.35
			3'UTR   0.35
			INTRONIC	0.1
			SYNONYMOUS_CODING	0.15
			NON_SYNONYMOUS_CODING	0.2
			SPLICE_SITE	0.15
			STOP_GAINED	0.02
			STOP_LOST	0.03
			UPSTREAM	0.2
			DOWNSTREAM	0.2

	(5) --case = one or more input case files.
	(6) --control = one or more input control files.
	(7) --out = output directory. (the output directory must be empty)
	(8) --non_syn , a flag to indicate processing only NON_SYNONYMOUS_CODING SNVs. This means that all other categories mentioned in the weight file will be ignored and only NON_SYNONYMOUS_CODING	will be processed for. Also if the --non_syn flag is provided then care must taken that the weight file being inputted must atleast contain the functional category NON_SYNONYMOUS_CODING .

2. Encryption can be fired as 
		** ./engine.pl --encrypt  --user=S --keyfile=S --case=S{,} --control=S{,} --genelist=S --weight=S --out=S [--non_syn]"
	OR directly through the encrypt.pl as,
		** ./encrypt.pl  --user=S --keyfile=S --case=S{,} --control=S{,} --genelist=S --weight=S --out=S [--non_syn]"

3. OUTPUT: The output directory will contain 
	(1) A directory GENES -> this contains vcf files per gene encountered in all input nvcf files cumulatively. The gene files are named GENEXXX, where XXX is a number assigned to the gene based on counter (incremented whenever a new gene is encountered). The gene files are split into subdirectories EVEN and ODD to avoid running into system dependent file handlig issues.
	(2) GeneInfo.txt -> this tab-delimited file contains the name (encrypted) of all genes encountered along with the path to their vcf files (into the GENES directory)
		Example:
			@#%#&&!      output/GENES/EVEN/DIR_4/GENE4.vcf
			@$^@fdf      output/GENES/ODD/DIR_3/GENE3.vcf
			@$^@3df      output/GENES/EVEN/DIR_2/GENE2.vcf
			@$^@8df      output/GENES/EVEN/DIR_0/GENE0.vcf
			@$^dsgf      output/GENES/ODD/DIR_1/GENE1.vcf
			
	(3) Pheno.txt -> This tab-delimited file contains all the samples IDs (encrypted) encoutered in all the input vcf files along with their phenotype ( 1 or -1) 
		Example:

			@%^%^#$^@	1
			!$#^#^#^#	1
			^$#!#!$#%	-1
			&@!!@#$#g	-1

##########################################################
##
##  Association testing
##
##########################################################



##########################################################
##
##  Decryption of results
##
##########################################################



##########################################################
##
##  Pre-requisites for running the application
##
##########################################################

Pre-requisites:

1. Download and Install vcfCodingSNPs as mentioned in http://www.sph.umich.edu/csg/liyanmin/vcfCodingSnps/Installation.shtml
2. Download and Install the perl module Crypt.
	* Module can be downloaded from http://search.cpan.org/~dparis/Crypt-DES-2.05/
	* If crypt is installed in a directory say "/root/temp_lib" , then add this to the perl5lib environment variable as,
		export PERL5LIB=/root/temp_lib/lib64/perl5/site_perl/5.10.0/x86_64-linux-thread-multi/Crypt
3. If you put all source code in a directory say "working_dir" , then make note of the following:
	* Create a subdirectory "working_dir/users", this must be empty before the 1st user registers. All temporary keys from all users will go inside "users". Also the final public key will be created inside "users/key". 

4. Outputs from encryption.pl will be directed to the output directory specified by each user.






##########################################################
##
##  TODO
#
##########################################################


* Synchronize output vcf filenames for genes between multiple users, currently this is synchronized for multiple input vcf files from one user only.
* 
*

