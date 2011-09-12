#!/usr/bin/perl -w

## Requires installing the modules Crypt and its supporting modules and adding to the path.

use Crypt::CBC;
use Getopt::Long;
use Switch;
use POSIX qw/floor/;


## Other assumptions ->
# Input vcf files are annotated using vcfCodingSNPSv1.5 .
# All input vcf files are have been validated, for e.g. using vcf-tools-validate
# Control files given Phenotype = -1 & Case files give Phenotype = 1
# REquires the genelist used for annotating the vcf files as input
# Requires a weights files, containing function category for annotation with their weights.
#
# VCF input : assumption that a line in VCF file has  9 fixed  fields i.e. chro, pos, ....  info, format.
#

{ package DATA;
	my @input_file=();  # array of hashes {NAME, HANDLE, PHENOTYPE, Flag_read, @line, #IDstart, #IDend}
	my $output_dir="";		# Input parameter
	my $master_gene;		# Reference pointing to structure 'user' containing user's information, key and cipher - Input parameter
	my $user;			# Reference pointing to structure 'gene' containing list of genes, functional weights - Input parameter

	my %curr_gene_pool=();		# a pool which maintains list of genes currently having open file handles. Once a gene is completely analysed, its file handle is closed and removed from this hash
	my $geneCount=0;		# a counter of total number of genes encountered, aids in naming new vcf files as GeneXXX.vcf
	my @IDs=();			# Array of Individuals IDs filled while reading headers of all input vcf files
	my %IDhash=();			# Hash with key=> Individual ID, Value=> Phenotype 
	my @header=();			# Common header to be printed to new vcf files.


sub new {       	#subroutine to create new object. Paramters set are 'user' reference and 'gene' reference
	my $self = shift;
	$user=shift;
	$master_gene=shift;
	bless {} }#end-sub-new

sub set_out_dir {	#subroutine to set output directory
	my $self = shift;
	$output_dir=shift;
	$output_dir=~s/\/$//;
	}#end-sub-set_out_dir

sub open_input_files{	#subroutine to open all input files and store file handles and other info in an array of hashes indexed by input-file name '@input-file'.
        my $self = shift;
	my(@case_files) = @{$_[0]};
        my(@control_files) = @{$_[1]};
        my $count=0;  my $file;
        if($#case_files > -1) { #signals a case file
        foreach $file (@case_files[0..$#case_files]) {
                if($file =~ m/vcf/){
                open("FH$count", "<", "$file") or die "Can't open case file $file\n";
                push @input_file, {NAME=> $file, HANDLE=> "FH$count", PHENO=> 1};
 #push(@case, "FH$count");
                $count++;
                }
        }}
        if($#control_files > -1) { #signals a control file
        foreach $file (@control_files[0..$#control_files]) {
                if($file =~ m/vcf/){
                open("FH$count", "<", "$file") or  die "Can't open control file $file\n";
                push @input_file, {NAME=> $file, HANDLE=> "FH$count", PHENO=> -1};
 #push(@control, "FH$count");
                $count++;
                }
        }}
        return ($count);
}#end-sub-open_input_files

sub extractIndividualsIDs{		#subroutine to extract all Individual IDs, modify them, and store in corresponding data structures with relevant information like Phenotype.
        my $self = shift;
	my $param1=shift;
	my @parts = split(/\t/, $param1);
	my $aFile= shift;
	my $VCFfile = $aFile->{NAME}; #my $VCFfile = $_[1];
	my $phenotype = $aFile->{PHENO}; #my $phenotype = $_[2];
	
	$aFile->{ID_START}= scalar( @IDs);
	my $mark;
	my $partnum = scalar(@parts);
  	if($partnum>10){ #if there is more than one individual in the file
		if($phenotype==-1){  #if it is a control file
                	for($i=9;$i<$partnum;$i++){
                        	if($parts[$i]!~/\W/ && !exists($IDhash{$parts[$i]})){ #if the id hasn't been seen before and isn't blank
                                	$IDhash{$user->{USERNAME}.'control'.$parts[$i]} = -1;
                                        push(@IDs,$user->{USERNAME}.'control'.$parts[$i]); #add the id to the array
                                        $mark =1;
                                }
                                elsif($parts[$i]!~/\W/  && $mark!=2){ #if it is blank, use the filehandle as the ID unless the file handle has already been used
                                	push(@IDs,$user->{USERNAME}.'control'.$VCFfile);
                                        $mark=2;
                                }
                        } #end for
                        if($mark==1) {chomp($IDs[(scalar(@IDs)-1)])};   #removes the end line char from the last ID if taken from file input.
                } #end if pheno is control file
                else{   #this does the same as above, but for case files
                	for($i=9;$i<$partnum;$i++){
                        	if($parts[$i]!~/\W/ && !exists($IDhash{$parts[$i]})){
                                	push(@IDs,$user->{USERNAME}.'case'.$parts[$i]);
                                        $IDhash{ $user->{USERNAME}.'case'.$parts[$i]} = 1;
                                        $mark =1;
                                }
                                elsif($parts[$i]!~/\W/ && $mark!=2){
                                        push(@IDs,$user->{USERNAME}.'case'.$VCFfile);
                                	$mark=2;
                                }
                        } #end for
                        if($mark==1) {chomp($IDs[(scalar(@IDs)-1)])};   #removes the end line char from the last ID
                } #end else
        }
	else{  #if no id is specified in the file, use the filehandle
        	if($phenotype == -1 ){ 
		  	push(@IDs,$user->{USERNAME}.'control'.$VCFfile);   
			$IDhash{ $user->{USERNAME}.'control'.$VCFfile } = -1;
		}
                elsif($phenotype==1){
			push(@IDs,$user->{USERNAME}.'case'.$VCFfile);      
			$IDhash{ $user->{USERNAME}.'case'.$VCFfile } = 1;
		}
        }
        $aFile->{ID_END}= scalar(@IDs);
}#end-sub-extractIndividualsIDs



sub read_all_headers{		#subroutine to read the headers of all files, calls subroutine extractIndividualsIDs() for extracting IDs of individuals from the headers.
        my $self = shift;
        my $line = readline $input_file[0]->{HANDLE};
        push (@header, $line);
	$line= "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	push (@header, $line );
	$line= "##INFO=<ID=.,Number=.,Type=Float,Description=\"Weight of functional category of annotation\">\n";
	push (@header, $line );

	foreach $aFile(@input_file){
		while ($line = readline $aFile->{HANDLE}) {
#			print "$line \n";		#test printing
	        	if ($line =~ m/^##/)    # if line begins with a hash it is header
		                { next;	}
		        elsif  ($line =~ m/^#CHROM/)       { 

				chomp $line;
				$self->extractIndividualsIDs($line,$aFile);
        			$aFile->{FLAG_READ}= 1;			#0=don't read, 1=allow read, 2=EOD reached, don't read
				if ($header[(scalar(@header))-1]  =~ m/^#CHROM/){}
				else{
					my @items= split(/\t/,$line);
					push (@header, join("\t", @items[0..8]));
				}
				last;
			}
		}
	}
	$header[(scalar(@header))-1].="\t";
	$header[(scalar(@header))-1].= (join("\t", @IDs)."\n");
}#end-sub-read_all_headers



sub LHSinfo($) {  	# extract the LHS of a function-gene tuple from info field, i.e. - functional name
        my $self = shift;
	my @temp = split (/=/, $_[0]);
	return ($temp[0]);
}#end-LHSinfo

sub RHSinfo($) {   	# extract the RHS of a function-gene tuple from info field. i.e. - gene name
        my $self = shift;
        my @temp = split (/=/, $_[0]);
        if (scalar(@temp) < 2) {
		return ("");
	}
	else {
		$temp[1] =~ m/(\(.*\))/;
		my $geneName = substr($1, 1, (length($1)-2));
		return ($geneName);		
	}
}#end-RHSinfo


sub create_gene_file{ 		#subroutine to create a gene file called "geneXXX" where XXX is the current gene count. Its put in appropriate directory heirachy to incorporate tens of thousands of genes. Directory root is user-inputted $output_dir
        my $self = shift;
	my $name_key=shift;
	my $count = $geneCount;
	my $suffix_digit =0;
        my $path=$output_dir."/GENES"; if (! -d $path) { mkdir $path, 0777; chmod 0777, $path;}
	$path.= ($count%2==1) ? "/ODD" : "/EVEN"; 
	if (! -d $path) { mkdir $path, 0777; chmod 0777, $path;}
        $suffix_digit = $count%400;
	$path.="/DIR_$suffix_digit";
	if (! -d $path)  { mkdir $path, 0777;  chmod 0777, $path;}
	$path.="/GENE$count.vcf";
        open("FH_$name_key", ">>", $path) or die "Cannot create gene file for $name_key called $path";
        chmod 0777, $path;
	$current_gene_pool{$name_key} = ["FH_$name_key", $path];

	#TODO add in the @header a line for info
	foreach $line(@header){			# print common header to all output geneXXX.vcf files
		print {"FH_$name_key"} $line;
	}

	$geneCount++;
}#end-sub-create_gene_file

sub empty_gene_pool{ 		#subroutine to close a geneXXX.vcf file, print the gene file's path in a common GeneInfo.txt file, and finally empty the current gene pool
        my $self = shift;
	open FH_info, ">>".$output_dir."/GeneInfo.txt" or die "Cannot open ".$output_dir."/GeneInfo.txt file ";
	for $name_key ( keys %current_gene_pool ) {
#		print FH_info  $user->{CIPHER}->encrypt($name_key)."\t". $current_gene_pool{$name_key}[1]."\n";   #Update GeneInfo.txt with Gene name(encrypted) & its path on system
		 print FH_info  $name_key."\t". $current_gene_pool{$name_key}[1]."\n";   #Update GeneInfo.txt with Gene name (unencrypted) & its path on system
	 	close $current_gene_pool{$name_key}[0];  
	}
	close FH_info;
	for (keys %current_gene_pool){		  # completely empty %current_gene_pool
        	delete $current_gene_pool{$_};
	}	
}#end-sub-empty_gene_pool


sub get_genotype($){ 		#subroutine - extracts genotype from the Individual's data fields of a vcf file. Genotypes are assigned as follows; 0|0 - 0, 1|1 - 2, anything else - 1.
        my $self = shift;
	my $item = $_[0];
	if($item =~ /0\|0/ )		#Assumes GT format in vcf file is 0|0
		{ return '0';}
	elsif($item =~ /1\|1/ )		#Assumer GT format in vcf file is 1|1
		{ return '2';}
	else 				#Any other GT is genotyped with "1" for e.g. ./. or 0|1
		{return '1';}
}#end-get_genotype


sub  append_gene_file{  	#subroutine to print to an output geneXXX.vcf file. 
        my $self = shift;
	my $buffer = shift;

	for (keys %$buffer) {
		my $temp_handle= $current_gene_pool{$_}[0];
		my @arr = @{$buffer->{$_}};
#		my $string = $user->{CIPHER}->encrypt($arr[0]) ."\t";
#		$string.= $user->{CIPHER}->encrypt( $arr[1] - $master_gene->{GENE_HASH}->{$_} ) ."\t";
#		$string.= $user->{CIPHER}->encrypt( $arr[2] ) ."\t";
               my $string = $arr[0] ."\t";					#these 3 lines for printing unencrypted data - for testing.
               $string.= ($arr[1] - $master_gene->{GENE_HASH}->{$_} ) ."\t";
		$string.= $arr[2]  ."\t";
		$string .= join ("\t", @arr[3.. ((scalar @arr )-1) ]);
		$string.="\n";
		print  {$temp_handle} $string;
	}
}#end-sub-append_gene_file



sub extract_ID_data{ 	#subroutine to extract all the data of individuals from a line from an input VCF file. 
        my $self = shift;
	my @line = @_;

	my @IDdata=(("") x ((scalar @line)-9));
	for ($ind = 9; $ind < (scalar @line); $ind++) {
        	@IDdata [($ind -9)] = $self->get_genotype( $line[$ind] );
	}
	#print "sub-extract_ID_data IDdata  @IDdata \n";	#test printing
	return (@IDdata);
}#end-extract_ID_data

sub update_op_buffer{ 	#arguments - $object_ref (discarded), $op_buffer, $afile , $curr_chr, $curr_pos;
        my $self = shift;
	my $buffer_ref = shift;
	my $aFile=shift;

	if ($aFile->{FLAG_READ} == 0){
	        my $the_chr=shift;
	        my $the_pos=shift;
		if ( (@{$aFile->{LINE}}[0] == $the_chr ) &&  ( @{$aFile->{LINE}}[1] == $the_pos )){
			my (@IDdata) = $self->extract_ID_data(@{$aFile->{LINE}});
			for (keys %$buffer_ref){
				@{$buffer_ref->{$_}}[ (9+$aFile->{ID_START}) .. (9+ $aFile->{ID_END} -1) ] = @IDdata;
			}
			return 1;
		}
	}
	return 0;
}#end-update_OP_buffer


sub extractInfo {		#arguments - $object_ref (discarded), $file handle
	#subroutine to extract all information from a line from an input vcf file. Extracted things are chromose, position, functional annotaion + gene name and genotype data for all individuals.

 my $self =  shift;
 my $aFile=shift;
 my @line = @{$aFile->{LINE}}; 			
 my @info = split(/;/, $line[7]); 		
 my $output_buffer={};			

 my (@IDdata) = $self->extract_ID_data(@line);

#print "IDdata @IDdata \n";		# test printing following 5 lines
#print "line @line \n";
#print "info @info \n";
#print "IDend $aFile->{ID_END} \n";
#print "IDstart $aFile->{ID_START} \n";

 my @buff;

 foreach $info_item (@info){
	my $functionName = $self->LHSinfo($info_item);
	if (exists $master_gene->{WEIGHTS}->{$functionName}) {
	  my ($theGene) =  $self->RHSinfo($info_item);
	    if (($theGene ne "") && (exists $master_gene->{GENE_HASH}->{$theGene})){
		if (! exists $current_gene_pool{$theGene}){
#      print "inside create_gene_file theGene \n";		#test printing
			$self->create_gene_file($theGene);
		}
		if (! exists $output_buffer->{$theGene}){
#	print "inside   \@ {\$ output_buffer->{$ theGene}} = @ buff   \n";		#test printing
	                @buff = ((0) x (9+scalar(@IDs)));
			@buff[0..6] = @line[0..6];
			$buff[7] = $master_gene->{WEIGHTS}->{$functionName};
			$buff[8] = "GT";
			@buff[ (9 + $aFile->{ID_START}) .. (9 + $aFile->{ID_END} -1) ] = @IDdata;
			$output_buffer->{$theGene} = ();
			@{$output_buffer->{$theGene}} = @buff;
			undef @buff;
		}
		else {  
#	print "inside update functional weight \t Old weight= $output_buffer->{$theGene}[7] \t New weight= $master_gene->{WEIGHTS}->{$functionName}  \n"; #test printing
	                if($output_buffer->{$theGene}[7] < $master_gene->{WEIGHTS}->{$functionName} ) {
	                        $output_buffer->{$theGene}[7] = $master_gene->{WEIGHTS}->{$functionName} ;
	                }
 		}
            }
	}
 }
 return ($output_buffer);
}#end-extractInfo

sub find_next_SNP{  		#subroutine - from all the open files handles (input vcf files), this subroutine decides which file will have the next line processed. All input files are processed in the order of increasing SNP positions. The file which is decided to be processed next will have its FLAG_READ flag set to 1, to indicate that a file read is to be done afterwards.

        my $self = shift;
        my $min_chr="23";
        my $min_pos="18446744073709551615";
	my $minFile={};

        foreach $aFile(@input_file){
                if ($aFile->{FLAG_READ} != 2){
                        if ( @{$aFile->{LINE}}[0] < $min_chr ) {
                                $minFile=$aFile;
                                $min_chr= @{$aFile->{LINE}}[0];
                                $min_pos= @{$aFile->{LINE}}[1];
 			}
			elsif (( @{$aFile->{LINE}}[0] = $min_chr ) && ( @{$aFile->{LINE}}[1] < $min_pos)) {
                                $minFile=$aFile;
                                $min_chr= @{$aFile->{LINE}}[0];
                                $min_pos= @{$aFile->{LINE}}[1];
                        }
                }
        }
        if (defined $minFile) {
                foreach $aFile(@input_file){
                        if ( ($aFile->{FLAG_READ} != 2) && ($aFile ne $minFile) ){
                                $aFile->{FLAG_READ} = 0;
                        }
                }
	        $minFile->{FLAG_READ}=1;
        }
	return $minFile;
}#end-find_next_SNP


sub read_VCF_line{	#this sub performs a readline operation on those input vcf files which have their FLAG_READ set to 1. This is because only those files had been selected for processing. Once a readline operation cannot be done because of EOF the FLAG_READ is set to 2. If atleast one file has been read subroutine returns 1 (true) else 0 (false). 0 means no more files have any lines left to be read.
        my $self = shift;
	my $flag=0;
	foreach $aFile(@input_file){
		if ($aFile->{FLAG_READ} == 1){		#read a line if it was marked for reading i.e. FLAG_READ =1
			my $line = readline $aFile->{HANDLE};
			if (defined $line) {
				chomp $line;
			        @{$aFile->{LINE}} = split(/\t/, $line);
				$flag=1;
			}
			else {
				$aFile->{FLAG_READ} = 2;		#if EOF is reached
			}
		}
	}
	return $flag;
}#end-sub-read_VCF_line

sub shatterVCF{		#main function which reads all input vcf files, extracts data in increasing order of SNPs and outputs all encrypted information gene-wise, in files named as geneXXX.vcf .

        my $self = shift;
	my $curr_chr="0"; my $curr_pos="0";
	my $OPbuffer={};	

	$self->read_all_headers();

	my $isContinue = $self->read_VCF_line();	#reads a new line from any file whose FLAG_READ  is set to 1.

	while ($isContinue == 1){
		my $aFile = $self->find_next_SNP();	#gets the VCF file handle which has the next closest SNP.
		if (! defined $aFile) { last;}

		my $this_chr = @{$aFile->{LINE}}[0];
		my $this_pos = @{$aFile->{LINE}}[1];
		if ((scalar keys %current_gene_pool)==0) {  # if gene_pool is empty
			$curr_chr= $this_chr ; 
		}
		elsif ($curr_chr ne $this_chr){  # check for new chromosome, if so close all open geneXXX.vcf file handles in current gene pool 
			$self->empty_gene_pool();
			$curr_chr= $this_chr ;
		}

		($OPbuffer)= $self->extractInfo($aFile);
		
		foreach $otherFile(@input_file)  {
		   if ($otherFile->{FLAG_READ} == 0){
			my $isUpdate = $self->update_op_buffer($OPbuffer, $otherFile, $this_chr, $this_pos);	#from remaining files.. i.e. if any chr& pos overlap with this
			if ($isUpdate) {
				$otherFile->{FLAG_READ} = 1;
			}
		   }
		}
#		for (keys %$OPbuffer){	print " GENE= $_ \t LINE= @{$OPbuffer->{$_}} \n";}		#test printing
#		foreach $fil(@input_file){	print $fil->{NAME} . "\t" . $fil->{FLAG_READ}. "\n";	}	#test printing

		$self->append_gene_file($OPbuffer);

		for (keys %$OPbuffer){            # empty the buffer after printing to files.
                	delete $OPbuffer->{$_};
	        }

		$isContinue = $self->read_VCF_line();	#reads a new line from all those files whose flag is set.
	}
	$self->empty_gene_pool();
}#end-sub-shatterVCF

sub write_pheno_file{		#prints the phenotypes of all individuals from all input vcf files, into a common Pheno.txt file
        my $self = shift;
	open FH_id, ">>".$output_dir."/Pheno.txt" or die "Cannot create $output_dir/Pheno.txt ";
	for $id ( keys %IDhash ) {
#		print FH_id  $user->{CIPHER}->encrypt($id) . "\t"  . $IDhash{$id} . "\n";	# prints Individual ID (encrypted) & phenotype
		print FH_id  $id. "\t" . $IDhash{$id} . "\n";   		# prints Individual ID (not encrypted) & phenotype  for testing
        }
        close FH_id;
        for (keys %IDhash){      
                delete $IDhash{$_};
        }
}#end-write_pheno_file


} #Package close


package main;

my $options={INITIALIZE => sub {
	my $self= shift;
	my $programOptions = GetOptions (
                        "case=s{,}" => \@{$self->{CASE}},
                        "control=s{,}" => \@{$self->{CONTROL}},
                        "weight=s" => \$self->{WEIGHT_FILE},
                        "genelist=s" => \$self->{GENE_FILE},
                        "user=s" => \$self->{USERNAME},
                        "keyfile=s" => \$self->{KEY_FILE},
                        "out=s" => \$self->{OUTPUT_DIR},
                        "non_syn" => \$self->{FLAG_NON_SYN},     
                        "help|?" => \$self->{FLAG_HELP});
        }};

my $user={
	GENERATE_CIPHER => sub {
		my $self=shift; $self->{KEY_FILE}=shift;
	        my $key;

	        open keyFH, "<".$self->{KEY_FILE} or die "Can't open key file\n";
        	read(keyFH, $key, 8);
	        $self->{CIPHER} = Crypt::CBC->new(-key         => "$key",
                                'regenerate_key'  => 0,
                                -iv => "$key",
                                'prepend_iv' => 0 ); #generates the cipher
	        close keyFH;
	}};
	
my $gene={
	INITIALIZE => sub {
		my $self=shift;
		$self->{WEIGHT_FILE}= shift;
		$self->{GENE_FILE}=shift;
		$self->{FLAG_NON_SYN} = shift() ? 1: 0;
		$self->{WEIGHTS} = {};
		$self->{GENE_HASH} = {};
		},

	INIT_WEIGHTS => sub {
		my $self=shift;
        	open funcFH, "<".$self->{WEIGHT_FILE} or die "Can't open functional weights file\n";
	        while (<funcFH>){
	                chomp; my @temp = split(/\t/, $_);
			if ($self->{FLAG_NON_SYN}) {
				if ( $temp[0] eq "NON_SYNONYMOUS") {
		        	        $self->{WEIGHTS}->{$temp[0]} = $temp[1];
				}
			}
			else{
				 $self->{WEIGHTS}->{$temp[0]} = $temp[1];
			}
	        	}
        	close funcFH;
		#if ((scalar keys %$self->{WEIGHTS}) == 0 ) { die "Weights file is not valid\n";}
		},
	READ_GENE_FILE => sub {
                my $self=shift;
		open geneFH, "<". $self->{GENE_FILE} or die "Can't open gene_list file\n";

		my $temp = <geneFH>; chomp($temp); my @line = split(/\t/, $temp);
		my $indexGeneName=-1; my $indexTxstart=-1; my $i;
		for $i (0 .. $#line){
		        if ($line[$i] eq "\#name") {$indexGeneName=$i; }      
			if ($line[$i] eq "txStart") {$indexTxstart=$i; }
			}
		# Fill geneList with names of genes in order of appearance in the geneFile
		while (<geneFH>)  {
		        chomp;  @line = split(/\t/, $_);
		        $self->{GENE_HASH}->{$line[$indexGeneName]} = $line[$indexTxstart];
		        }
		close geneFH;
		}
};
		


if (@ARGV > 0 ) {
	&{$options->{INITIALIZE}}($options);

	 if ($options->{FLAG_HELP}) {print "\nHelp Options..\n";
                print "\n (1) Encrypt Data: \n\t ./encrypt.pl --user=S --keyfile=S --case=S{,} --control=S{,} -genelist=S -weight=S --out=S [--non_syn]";
                print "\n (2) Help: \n\t ./encrypt.pl --help\n\n";
                exit;
        }

	if($options->{USERNAME}){ 
                $user->{USERNAME} = $options->{USERNAME};}
		else { die "Username must be provided. \n"; }
	if(($options->{WEIGHT_FILE} ) && ($options->{GENE_FILE} )) {
		&{$gene->{INITIALIZE}}($gene, $options->{WEIGHT_FILE}, $options->{GENE_FILE}, $options->{FLAG_NON_SYN}); }
                else { die "Weights file and Gene list file must be provided. \n"; }
 
	if($options->{KEY_FILE} ne ""){ #signals the encryption key
		&{$user->{GENERATE_CIPHER}}($user, $options->{KEY_FILE}); }
                else { die "An encryption key must be provided in a file. \n"; }
 
	$VCFobj= DATA->new($user, $gene) ;

	$VCFobj->set_out_dir($options->{OUTPUT_DIR});
	
	&{$gene->{INIT_WEIGHTS}}($gene);
	&{$gene->{READ_GENE_FILE}}($gene);

	$VCFobj->open_input_files($options->{CASE},$options->{CONTROL});

	$VCFobj->shatterVCF();

	$VCFobj->write_pheno_file();

}

exit;




