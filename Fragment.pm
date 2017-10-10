#!/usr/bin/perl 
package Fragment;
use strict;
use Exporter;
use Chemistry::Mol;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS); 

# Modfied 2016.02.05 
# Each Fragment may have several Isomers, each being a Molecule with a different geometry
# Each Fragment also has a list of parent Heavy atoms, and an array of pointers to sibling Fragments
# Intra-Fragment hydrogen transfers simply add isomeric Molecules 

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT_OK   = qw();
@EXPORT      = ();
%EXPORT_TAGS = (
all=>[@EXPORT, @EXPORT_OK],
); 

sub new{
	my $class=shift;
	my $self = {
		heavy => shift,
	};
	bless $self, $class;
	return($self);
}

sub addIsomer{
	# Push a pointer to an isomeric Molecule
	my ($self,$i_ref)=@_;
	push(@{$self->{isomers}},$i_ref);
}
sub addSibling{
	# Push a pointer to a sibling Fragment
	my ($self,$c_ref)=@_;
	push(@{$self->{siblings}},$c_ref);
}
sub numIsomers{
	my $self=shift;
	return scalar(@{$self->{isomers}}); 
}
sub numSiblings{
	my $self=shift;
	return  scalar(@{$self->{siblings}}); 
}
sub getIsomer{
	my ($self,$i)=@_;
	my @is=@{$self->{isomers}};
	return($is[$i]);
}
sub getFile{
	# Return the filename of the Isomers with minimum energy 
	my $self=shift;
	my $ret='';
	my $E0=100000000000.;
	foreach my $i(@{$self->{isomers}}){
		#print "TEH NAME IS TEH ".$i->name()."\n";
		#my @array = split(/,/,$i->name());
		#print "TEH VALUES ARE: ";
		#foreach my $s(@array){
		#	print "!$s!  ";
		#}
		#print "\n";
		my ($E,$file)=split(/,/,$i->name());
		#print "	TEH E IS TEH $E \n";
		#print "	TEH FILE IS TEH $file \n";
		if($E<$E0){$E0=$E;$ret=$file;}
	}
	#print "TEH RET IS TEH $ret \n";
	return $ret;
}
sub getE{
	# Return the minimum energy among all Isomers 
	my $self=shift;
	my $ret=100000000.;
	foreach my $i(@{$self->{isomers}}){
		my ($E,$file)=split(/,/,$i->name());
		if($E<$ret){$ret=$E;}
	}
	return 1.0*$ret;  
}
sub getAllFiles{
	# Return the filenames of my and my siblings' lowest energies
	my $self=shift;
	my $ret = $self->getFile();
	if(defined $self->{siblings}){
		foreach my $s(@{$self->{siblings}}){
			$ret = $ret .",".$s->getFile();
		}
	}
	return($ret);
}
sub getEtot{
	# Return the sum of my energy and my siblings' energies
	my $self=shift;
	my $ret = $self->getE();
	if(defined $self->{siblings}){
		foreach my $s(@{$self->{siblings}}){
			$ret = $ret + $s->getE();
		}
	}
	return($ret);
}
