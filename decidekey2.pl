#!/usr/bin/perl -w
#Rajan Banerjee
#reb2143
#This perl script reads in part1 of key of all registerd users from the directoy "users" and generates a public key stored in file users/key
#The public key is a bitwise XOR operation performed on part1 from all users.

use Digest::MD5 qw(md5_hex);

opendir(DIR, "users");
@FILES= readdir(DIR); 
my $count = scalar(@FILES) -2;
my %handles=();
my @user_keys =();

# TODO directory "users" must atleast have keypart1 from 2 users.
# if ($count<2) {die "not enough users registered for key generation"; }

for ($index=0; $index < $count; $index++ )
{
        open($fixed_split{$index}, "<users/" . $FILES[$index+2] )|| die("Could not open users/$FILES[$index+2] file!");
	my $temp_handle = $fixed_split{$index};
	$user_keys[$index] = <$temp_handle>;
#	print "\nkeys.. $user_keys[$index] ";
}

closedir (DIR);
for ($index=0; $index < $count; $index++ )
{
        close ($fixed_split{$index});
}

open out, ">", "users/key" or die "Can't open key file!";

my $temp_key = $user_keys[0]^$user_keys[1];	# XOR keys from 1st two users and store in $temp_key
my $key="";
if($count>2) { 		#if there are more users, XOR temp_key with each one by one, finally store in $key
	for ($index=2; $index < $count; $index++ ) {
		$key = $temp_key^$user_keys[$index];
		$temp_key = $key;
	}
}
else
{	$key = $temp_key; }

print out "$key";
close out;
print "\nThe Primary key is stored in users/key :\n$key";
exit;
