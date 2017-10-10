#!/usr/bin/perl 
# Write Gaussian input files, read Gaussian output files 
package Chemistry::File::Gaussian;
use base qw(Chemistry::File);
use Chemistry::Mol;
use FileHandle;
use POSIX;
use strict;
use warnings;

Chemistry::Mol->register_format(log => __PACKAGE__); 
Chemistry::Mol->register_format(com => __PACKAGE__); 

# Read a Gaussian *.log file, returning a 
# Molecule object of the output with the Name set to total energy | filename 
#
sub read_mol{
	my ($self,$fh,%opts)=@_;
	return if $fh->eof;

	my $mol_class=$opts{mol_class}||"Chemistry::Mol";
	my $mol=$mol_class->new;
	$mol->name("10000.|uninitialized.log");
	my $finished=0; my $Natoms=0;

	my @lines=<$fh>;
	foreach my $l(@lines){
		$finished++ if $l=~m/Normal termination/;
		if($l=~m/^ NAtoms=/){my @f=split(/ +/,$l); $Natoms=$f[2];}
	}
	if($finished){
		# Read charge and final energy 
		foreach my $line(@lines){
			if($line=~m/^ Charge/){
				my @fields=split(/ +/,$line);
				$mol->charge($fields[2]);
			}
			if($line=~m/SCF Done/){
				my @fields=split(/ +/,$line);
				$mol->name($fields[5])."|".$fh;
			}
		}
		# Find the last "rientation" 
		my $last=0;
		foreach my $il(0..scalar(@lines)-1){
			$last=$il if $lines[$il]=~m/rientation/;
		}
		# Read the Cartesian coordinates 
		foreach my $iat(0..$Natoms-1){
			my @f=split(/ +/,$lines[$iat+5+$last]);
			$mol->new_atom(Z=>$f[2],coords=>[$f[4]*1.,$f[5]*1.,$f[6]*1.],name=>$iat);
		}
	}
	return($mol);
}

sub mult{
	# Given a Molecule, return its lowest possible multiplicity 
	my $mol= $_[0];
	my $nelec=-$mol->charge();
	foreach my $a($mol->atoms()){$nelec=$nelec+$a->Z;}
	my $m=$nelec%2+1;
	return($m);
}

sub write_string{
	my ($class, $mol, %opts)=@_;
	%opts=(head=>"# PM6 ",tail=>"\n",title=>"Dummy title",%opts);
	my $ret= $opts{head}."\n\n".$opts{title}."\n\n";

	my $m=mult($mol);
	if(defined($opts{multiplicity})){$m=$opts{multiplicity};}
	$ret=$ret.floor($mol->charge())." ".$m."\n";
	foreach my $a($mol->atoms()){
		#$ret=$ret.sprintf(" %5d %12.6f %12.6f %12.6f %20s  \n",$a->Z,$a->x3,$a->y3,$a->z3,$a->name);
		$ret=$ret.sprintf(" %5d %12.6f %12.6f %12.6f \n",$a->Z,$a->x3,$a->y3,$a->z3);
	}
	$ret=$ret."\n".$opts{tail};
	return($ret);
}
	
	
	
			
1;
