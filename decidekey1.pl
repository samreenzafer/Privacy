#!/usr/bin/perl -w
#Rajan Banerjee
#reb2143
#This perl script reads in a username and password from a user and performs MD5 hexes to hash the password.
#That hash is then bitwise XOR'ed to the hashes of the current date, time, and a random number.
#The result is saved to a file "users/USERNAME_KeyPart1"  in directory "users"

use Digest::MD5 qw(md5_hex);
use File::Spec;
my $dirname = dirname( File::Spec->rel2abs( __FILE__ ) );

sub user_key  #parameters...
{
my @ARGS = @_;

	if(scalar (@ARGS)<2){ die "not enough arguments have been passed\n";} #checks that user gave an input file

my $localtime = localtime;		
my @localtime =  split(/\s/,$localtime);
my $date = "$localtime[1]" . ' ' . "$localtime[2]" . ' '."$localtime[4]"; #creates a string for the date as in mmm dd yyyy
my $user = $ARGS[0];  #takes the user's name, so that can add a check if needed
my $pw = $ARGS[1];	#takes pw from command line args to generate key



#generates md5_hex hashes of the password, current time, the date, and a random number

my $passHash = md5_hex($pw);  

my $dateHash = md5_hex($date);

my $timeHash = md5_hex($localtime[3]);

my $randHash = md5_hex(rand);


#performs a bitwise XOR on the hashes and prints them to a file.

my $keyPart1 = ((($passHass^$dateHash)^$timeHash)^$randHash);

if (-e "$dirname/users"){
	open OUT1, ">", "$dirname/users/$ARGS[0]_KeyPart1" or die "Can't open $ARGS[2]\n"; }
else { die "\n $dirname/users directory does not exist\n"; }

$keyPart1 =~ s#/#//#;

print OUT1 "$keyPart1";

close OUT1;

}
1;
