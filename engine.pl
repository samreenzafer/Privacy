#!/usr/bin/perl -w
use strict;
# use warnings;
use Getopt::Long;
use Switch;
use File::Spec;
#  All things marked with TODO are to incomplete or need modifications.
#
### TODO : ADD NON_SYN flag..
#
my $CURR_DIR = dirname (File::Spec->rel2abs( __FILE__ ));

require "$CURR_DIR/decidekey1.pl";
my $DECIDEKEY2="$CURR_DIR/decidekey2.pl";
my $ENCRYPT="$CURR_DIR/encrypt.pl";
my $MERGE="$CURR_DIR/mergeVtFiles.pl";
my $ASSOC="$CURR_DIR/runAssociation.pl";
my $DECRYPT="$CURR_DIR/decryptScores.pl";

# require "decidekey2.pl";

my $programOptions; 
my $flag_register; 
my $flag_genkey; 
my $flag_encrypt; 
my $flag_assoc; 
my $flag_decrypt; 
my $flag_help;
my $username=""; 
my $password=""; 
my $case=""; 
my $control=""; 
my $output=""; 
my @src_dir; 
my $dest_dir; 
my $info=""; 
my $key="";
my $genelist="";
my $weightfile="";

if ( @ARGV > 0 ) {
	$programOptions = GetOptions (
			"register" => \$flag_register,
			"genkey" => \$flag_genkey,
			"encrypt" => \$flag_encrypt,
			"assoc" => \$flag_assoc,
			"decrypt" => \$flag_decrypt,
			"user:s" => \$username,
			"pwd:s" => \$password,
			"case:s" => \$case,
			"control:s" => \$control,
			"genelist:s" => \$genelist,
			"weight:s" => \$weightfile,
			"out:s" => \$output,
			"key:s" => \$key,
			"src:s{2}" => \@src_dir,
			"dest:s" => \$dest_dir,
			"score:s" => \$info,
			"help|?" => \$flag_help);
	
	if ($flag_help) {print "\nHelp Options..\n";
		print "\n (1) New user Registration: \n\t ./engine.PL --register --user=S --pwd=S";
		print "\n (2) Generate Key: \n\t ./engine.PL --genkey --user=S";
 		print "\n (3) Encrypt Data: \n\t ./engine.pl --encrypt --user=S --keyfile=S --case=S{,} --control=S{,} --genelist=S --weight=S --out=S [--non_syn]";
		print "\n (4) Run Association Tests: \n\t ./engine.PL --assoc --src=S --dest=S [--out=S]";
		print "\n (5) Decrypt Scores: \n\t ./engine.PL --decrypt --key=S --out=S --score=S";
		print "\n (6) Help: \n\t ./engine.PL --help\n\n";
	}
	elsif ($flag_register)
	{
		if (($username eq "") || ($password eq ""))
			{die "\nIncorrect options for -register\n";}
		else
			{   user_key ($username, $password); }  ##subroutine from decidekey1.pl called to make part 1 of public key per user.
	}
	elsif ($flag_genkey)
	{
		if ($username eq "")
			{die "\nIncorrect options for -genkey\n";}
		else
		{
			if (-e  "$CURR_DIR/users/key") {
				print "\nThe Primary key is stored in $CURR_DIR/users/key :\n";
			open FH, "<", "$CURR_DIR/users/key" or die "Can't open $CURR_DIR/users/key\n"; 
			print <FH>;
			close (FH);
			}
			else {
				system ($DECIDEKEY2);  #calls decidekey2.pl script
			}
		}
	}
	elsif ($flag_encrypt)
	{

		##### TODO ANNOTATE input files using VCFcodingSNPS
		#
		#
		#
		#

		##### Submit annotated files for encryption, public key, username, genelist (same as that used in annotation) and weights files MUST be provided.
		my @arguments = ($ENCRYPT); # calls encrypt.pl script
		if ($username ne "")
			{ @arguments = (@arguments, "-user", $username);}
		if ($key ne "")
			{ @arguments = (@arguments, "-keyfile", $key);}
		if ($case ne "")
			{ @arguments = (@arguments, "-case", $case);}
		if ($control ne "")
			{ @arguments = (@arguments, "-control", $control);}
		if ($output ne "")
			{ @arguments = (@arguments, "-out", $output);}
		if ($genelist ne "")
                        { @arguments = (@arguments, "-genelist", $output);}
 		if ($weightfile ne "")
                        { @arguments = (@arguments, "-weight", $output);}
 
		system(@arguments);
#			or die "\n\nEncryption Failed\n";
		print "\n\nEncryption Successful\n";
	}
	elsif ($flag_assoc)		#TODO Change below for new association tests
	{
		if ((scalar(@src_dir) < 2) || ($dest_dir eq ""))
#		|| ($output eq "") )
		# || ($info eq ""))
			{die "\nIncorrect options for -assoc\n";}
		else {
			$dest_dir=~s/\/$//;
			my @arguments = ($MERGE);  #"/ifs/scratch/c2b2/ip_lab/sz2317/privacy/workingdir/mergeVtFiles.pl",
					"-d1", $src_dir[0],
					"-d2", $src_dir[1],
					"-out", $dest_dir);
			system(@arguments);
#				or die "\n\nMerging Failed\n";
			print "\n\nMerging Successful\n";

			@arguments = ($ASSOC);  #"/ifs/scratch/c2b2/ip_lab/sz2317/privacy/workingdir/runAssociation.pl",
					"-out", "Scores.txt",
					"-dir", $dest_dir,
					"-info", $dest_dir."/info.txt");
			system(@arguments);
#				or die "\n\nAssociation Testing Failed\n";
			print "\n\nAssociation Testing Successful\n";
		}
	}
	elsif ($flag_decrypt)		#TODO change below for new decryption
	{
		if (($key eq "") || ($output eq "")|| ($info eq ""))
			{die "\nIncorrect options for -decrypt\n";}
		else {
			my @arguments = ($DECRYPT);  # "/ifs/scratch/c2b2/ip_lab/sz2317/privacy/workingdir/decryptScores.pl",
					"-key", $key,
					"-out", $output,
					"-scores", $info);
			system(@arguments);
#				or die "\n\nDecryption Failed\n";
			print "\n\nDecryption Successful\n";
		}
	}
	else
	{	die "\nRun engine with options or -help\n";}

}
else {
	die "\nRun engine with options or -help\n";
}


exit 0;




