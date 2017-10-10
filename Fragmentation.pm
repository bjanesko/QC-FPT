#!/usr/bin/perl 
# Ben Janesko 2015.07.29 
# Work routines for mass spec fragmentation 

package Fragmentation;
use strict;
use Exporter; 
use List::Util qw(min max);
use Chemistry::Mol;
use Chemistry::InternalCoords;
use Chemistry::File::SMILES;
use Chemistry::Bond::Find ':all';
use Fragment; 
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(unsaturatedNeighbors findTargetFragmentations Isomerize 
BoltzmannFactor mult moleculeName heavyAtoms);
%EXPORT_TAGS = (
All=>[qw(&unsaturatedNeighbors &findTargetFragmentations &Isomerize
&BoltzmannFactor &mult &moleculeName &heavyAtoms)]
); 

# 2016.04.28: Global value for keeping fragments 
my $MassThreshold = 0.9;

sub findTargetFragmentations{
	# Work routine to generate all fragmentations from breaking all bonds in a list of parent fragments 
	my ($parent,$q0,$Mtar,$frags_ref,$umols_ref,$unames_ref) = @_;
	my @bs=$parent->getIsomer(0)->bonds(); 
	my $Nbond=scalar(@bs);
	for(my $ibond=$Nbond-1;$ibond>=0;$ibond=$ibond-1){
		my $bi=$bs[$ibond];my @ais=$bi->atoms();
		if($ais[0]->Z>1 and $ais[1]->Z>1){
		for(my $jbond=$ibond;$jbond>=0;$jbond=$jbond-1){
			my $bj=$bs[$jbond];my @ajs=$bj->atoms();
			if($ajs[0]->Z>1 and $ajs[1]->Z>1){
			for(my $kbond=$jbond;$kbond>=0;$kbond=$kbond-1){
				my $bk=$bs[$kbond];my @aks=$bk->atoms();
				if($aks[0]->Z>1 and $aks[1]->Z>1){
				for(my $lbond=$kbond;$lbond>=0;$lbond=$lbond-1){
					my $bl=$bs[$lbond];my @als=$bl->atoms();
					  if($als[0]->Z>1 and $als[1]->Z>1 ){

						# The atom numbers of all of the atoms involved in the bond breaking 
						my $atomNumbers="";
						foreach my $i(@ais){$atomNumbers=$atomNumbers.sprintf("%3d,",$i->name+1);}
						foreach my $i(@ajs){$atomNumbers=$atomNumbers.sprintf("%3d,",$i->name+1);}
						foreach my $i(@aks){$atomNumbers=$atomNumbers.sprintf("%3d,",$i->name+1);}
						foreach my $i(@als){$atomNumbers=$atomNumbers.sprintf("%3d,",$i->name+1);}

						# The names of all of the atoms involved in the bond breaking
						my %brokeatoms;
						foreach my $i(@ais){$brokeatoms{$i}=1;}
						foreach my $i(@ajs){$brokeatoms{$i}=1;}
						foreach my $i(@aks){$brokeatoms{$i}=1;}
						foreach my $i(@als){$brokeatoms{$i}=1;}
						#push(@ais,push(@ajs,push(@aks,@als)));
						#my %brokeatoms; @brokeatoms{@als}=();
		
						# Delete bonds, generate new molecule(s) from separated parent 
						$parent->getIsomer(0)->delete_bond($bi);
						if($jbond ne $ibond){$parent->getIsomer(0)->delete_bond($bj);}
						if($kbond ne $ibond and $kbond ne $jbond){$parent->getIsomer(0)->delete_bond($bk);}
						if($lbond ne $ibond and $lbond ne $jbond and $lbond ne $kbond){$parent->getIsomer(0)->delete_bond($bl);}
			 			my @newmols=$parent->getIsomer(0)->separate; 
	
						# Rebuild parent 
						$parent->getIsomer(0)->add_bond($bi);
						if($jbond ne $ibond){$parent->getIsomer(0)->add_bond($bj);}
						if($kbond ne $ibond and $kbond ne $jbond){$parent->getIsomer(0)->add_bond($bk);}
						if($lbond ne $ibond and $lbond ne $jbond and $lbond ne $kbond){$parent->getIsomer(0)->add_bond($bl);}
	
						my $Nfrags=scalar(@newmols);
						if($Nfrags>1){

							# Find fragments that approach the target 
							my $keep=-1;my $keptMass=-1;
							foreach my $if(0..scalar(@newmols)-1){
								my $f=$newmols[$if];
								my $fm=$f->mass;
								if(abs($Mtar-$fm)<5){
									$keep=$if;
									$keptMass=$fm;


									print "Testing fragmentation bonds $ibond,$jbond,$kbond,$lbond atoms $atomNumbers fragment ".moleculeName($newmols[$keep])." of ".$Nfrags." mass ".$newmols[$keep]->mass."\n";
									# Get lists of each fragment's heavy atoms 
									my @heavies=();
									foreach my $m(@newmols){push(@heavies,heavyAtoms($m));}
		
									# Set the charges on the molecules in $newmols 
									foreach my $m(@newmols){$m->charge(0);}
									$newmols[$keep]->charge($q0);

									# Add this fragmentation pattern 
									if(abs($keptMass-$Mtar)<$MassThreshold){
										print "Keeping ".moleculeName($newmols[$keep])." mass ".$newmols[$keep]->mass."\n";
										updateSiblingFragments(\@newmols,$umols_ref,$unames_ref,\@heavies,$frags_ref);
									}
	
									# Do all possible single alpha hydrogen transfers 
									print "Starting H-atom transfer with kept fragment $keep mass ".$newmols[$keep]->mass."\n";
									my $newmols_ref = \@newmols;
									my $brokeatoms_ref= \%brokeatoms;
									transferAlphaHydrogens($newmols_ref,\@heavies,$frags_ref,$umols_ref,$unames_ref,$brokeatoms_ref,$keep,$Mtar);
								}
							}
						}
				}}
			}}
		}}
	}}
}
sub addSiblingFragments{
	# Add a set of sibling Molecules to the current Fragment list
	my($nm_ref,$heavy_ref,$frags_ref)=@_;
	my @nm=@$nm_ref; my @heavy=@$heavy_ref;
	my $keep=1;
	foreach my $m(@nm){if(scalar($m->atoms())<1){$keep=0;}}
	foreach my $h(@heavy){if(length($h)<2){$keep=0;}}
	if($keep){
		my @newfrags=();
		foreach my $im(0..scalar(@nm)-1){
			my $f=new Fragment($heavy[$im]);
			$f->addIsomer($nm[$im]);
			push(@newfrags,$f);
		}
		foreach my $if(0..scalar(@newfrags)-1){
			foreach my $jf(0..scalar(@newfrags)-1){
				if($if ne $jf){$newfrags[$if]->addSibling($newfrags[$jf]);}
			}
		}
		foreach my $f(@newfrags){push(@$frags_ref,$f);}
	}
}

sub transferAlphaHydrogens{
	# Do all possible single alpha hydrogen transfers between molecules in @newmols 
	# that give the correct kept mass 
	my ($newmols_ref,$heavies_ref,$frags_ref,$umols_ref,$unames_ref,$brokeatoms_ref,$keep,$Mtar)=@_;
	my @newmols=@$newmols_ref;
	my %brokeatoms=%$brokeatoms_ref;
	#print "Starting H-atom transfer with kept fragment $keep mass ".$newmols[$keep]->mass."\n";
	my $Nfrags=scalar(@newmols);
	foreach my $imol(0..$Nfrags-1){
		foreach my $jmol(0..$Nfrags-1){
			if($jmol ne $imol){
				#print "Testing transfer from molecules $jmol ".$newmols[$jmol]->sprintf("%s ")." to $imol ".$newmols[$imol]->sprintf("%s ")."\n";
				foreach my $iatom($newmols[$imol]->atoms()){
					if(exists $brokeatoms{$iatom}){
						#print "Fragment $imol atom $iatom is in broken \n";
						foreach my $jatom($newmols[$jmol]->atoms()){
							if(exists $brokeatoms{$jatom}){
								#print "Fragment $jmol atom $jatom is in broken \n";
								foreach my $katom($jatom->neighbors()){
									if($katom->Z eq 1){
										# Transfer atom katom from newmols(jmol) atom katom to newmols(imol) atom iatom 
										my $transName = $katom->name;
										my @temp= $katom->bonds;
										my $transBond = $temp[0];

										# Get the hydrogen added to iatom 
										my($capAtom,$capBond)=getAddedHydrogen($iatom,$transName);
					
										$newmols[$jmol]->delete_bond($transBond);
										$newmols[$jmol]->delete_atom($katom);
										$newmols[$imol]->add_atom($capAtom);
										$newmols[$imol]->add_bond($capBond);

										# See if the charged fragment hits our mass range
										if(abs($newmols[$keep]->mass -$Mtar)<$MassThreshold){
											print "Keeping H-transferred ".moleculeName($newmols[$keep])." mass ".$newmols[$keep]->mass."\n";
	
											# See if this molecule is unique , clone if needed 
											updateSiblingFragments(\@newmols,$umols_ref,$unames_ref,$heavies_ref,$frags_ref);
										}

										# Undo our changes to input 
										$newmols[$imol]->delete_bond($capBond);
										$newmols[$imol]->delete_atom($capAtom);
										$newmols[$jmol]->add_atom($katom);
										$newmols[$jmol]->add_bond($transBond);
									}
								}
							}
						}
					}
				}
			}
		}
	} 
}

sub getAddedHydrogen{
	my($iatom,$transName)=@_;
	my @nabor=$iatom->neighbors();
	my $cap=Chemistry::Atom->new(symbol=>'H',name=>$transName);
	my $ic=Chemistry::InternalCoords->new($cap,$iatom,1.1);
	if(scalar(@nabor)>0){$ic=Chemistry::InternalCoords->new($cap,$iatom,1.1,$nabor[0],120.0);}
	if(scalar(@nabor)>1){$ic=Chemistry::InternalCoords->new($cap,$iatom,1.1,$nabor[0],109.5,$nabor[1],120.0);}
	$cap->coords($ic->cartesians);
	#print "Made new atom name ".$cap->name." id ".$cap->id." \n";
	if(scalar(@nabor)>2){
		my $v=$cap->coords - $nabor[2]->coords;
		if($v->length<0.7){
			#print "Cap atom clash with ".$nabor[2]->coords->stringify." \n";
			#print "Moving cap atom from ".$cap->coords->stringify." \n";
			$ic=Chemistry::InternalCoords->new($cap,$iatom,1.1,$nabor[0],109.5,$nabor[1],-120.0);
			$cap->coords($ic->cartesians);
			#print "to ".$cap->coords->stringify." \n";
		}
	}
	my $capBond=Chemistry::Bond->new(atoms=>[$iatom,$cap]);
	#print "Made new bond id ".$capBond->id." \n";
	return($cap,$capBond);
}

sub Isomerize{
	# Isomerize all molecules in each fragment 
	my ($frags_ref,$umols_ref,$unames_ref) = @_;
	foreach my $frag(@$frags_ref){
		my $Nis=$frag->numIsomers();
		foreach my $imol(0..$Nis-1){
			#print "Isomerizing $imol of $Nis fragment ".moleculeName($frag->getIsomer(0))."\n";
			my $mol=$frag->getIsomer($imol);

			# Does this molecule have an atom with >1 implicit hydrogens? 
			foreach my $iatom($mol->atoms()){
				if($iatom->calc_implicit_hydrogens ge 2){
	
					# Does an adjacent atom have 0 implicit hydrogens? 			
					my @inabor=$iatom->neighbors();
					my $hatom = undef; 
					my $hbond = undef; 
					my $transName = ''; 
					my $un1 = undef; 
					foreach my $jatom(@inabor){
						if($jatom->Z>1 and $jatom->calc_implicit_hydrogens eq 0){
	
							# Does atom jatom have a hydrogen katom to transfer? 
							my @jnabor=$jatom->neighbors($iatom);
							foreach my $katom(@jnabor){
								if($katom->Z eq 1 ){ 
									# Transfer atom katom from jatom to iatom 
									#print "Molecule ".$mol->sprintf("%f")." unsat hopping\n";
									#print "Old name: ".moleculeName($mol)."\n";
									my @temp = $katom->bonds;
									my $oldbond = $temp[0];
									my($newatom,$newbond)=getAddedHydrogen($iatom,$katom->name);
							
									$mol->delete_bond($oldbond);
									$mol->delete_atom($katom);
									$mol->add_atom($newatom);
									$mol->add_bond($newbond);

									#print "New name: ".moleculeName($mol)."\n";
									my $mol2 = updateUniqueMolecules($mol,$umols_ref,$unames_ref,1);
									$frag->addIsomer($mol2);

									$mol->delete_bond($newbond);
									$mol->delete_atom($newatom);
									$mol->add_atom($katom);
									$mol->add_bond($oldbond);
								}
							}
						}
					}
				}
			}
			# Does this molecule have a saturated atom, containing a hydrogen, between two unsaturated atoms? 
			foreach my $iatom($mol->atoms()){
				if($iatom->Z>1 and $iatom->calc_implicit_hydrogens eq 0){
					my @nabor=$iatom->neighbors();
					my $hatom = undef; 
					my $hbond = undef; 
					my $transName = ''; 
					my $un1 = undef; 
					my $un2 = undef; 
					foreach my $jatom(@nabor){
						if($jatom->Z>1 and $jatom->calc_implicit_hydrogens >0){
							if(!(defined $un1)){$un1 = $jatom;}
							else{$un2 = $jatom;}
						}
						if($jatom->Z eq 1){
							$hatom = $jatom;
							$transName = $jatom->name;	
							my @temp= $jatom->bonds;
							$hbond= $temp[0];
						}
					}
					if(defined $hatom and defined $un1 and defined $un2){
						print "Molecule ".$mol->sprintf("%f")." hopping \n";
	
						# We have to be careful: the new bond is automatically added to the molecule! 
						 
						# Move to atom un1 
						my($cap1,$bond1)=getAddedHydrogen($un1,$transName);
						$mol->delete_bond($hbond);
						$mol->delete_atom($hatom);
						$mol->add_atom($cap1);
						$mol->add_bond($bond1);

						# Update 
						my $mol2 = updateUniqueMolecules($mol,$umols_ref,$unames_ref,1);
						$frag->addIsomer($mol2);

						# Undo changes 
						$mol->delete_bond($bond1);
						$mol->delete_atom($cap1);
						$mol->add_atom($hatom);
						$mol->add_bond($hbond);
	
						# Move to atom un2 
						my($cap2,$bond2)=getAddedHydrogen($un2,$transName);
						$mol->delete_bond($hbond);
						$mol->delete_atom($hatom);
						$mol->add_atom($cap2);
						$mol->add_bond($bond2);

						# Update 
						my $mol2 = updateUniqueMolecules($mol,$umols_ref,$unames_ref,1);
						$frag->addIsomer($mol2);

						# Undo changes 
						$mol->delete_bond($bond2);
						$mol->delete_atom($cap2);
						$mol->add_atom($hatom);
						$mol->add_bond($hbond);
	
					}
				}
			}
		}
	}
}
				

sub OldtransferAlphaHydrogens{
	# For each current fragment containing unsaturated heavy atoms, transfer a
	# hydrogen atom from each sibling fragment. 
	my ($pmass,$fragments_ref,$umols_ref,$unames_ref) = @_;
	my  @fragments=@$fragments_ref;
	my $Nstart = scalar(@fragments);
	print "Trasferring alpha hydrogens from the first $Nstart fragments \n";
	foreach my $if1(0..$Nstart-1){
		my $f1=$fragments[$if1];
		my @ns1=unsaturatedNeighbors($f1->getIsomer(0),0);
		if(scalar(@ns1)>0 and exists $f1->{siblings} and $pmass-$f1->getIsomer(0)->mass()>1.1){ 
			my @f2s = @{$f1->{siblings}};

			# Fragment f1 has at least one unsaturated atom. Let's find unsaturated atoms on its siblings 
			foreach my $if2(0..scalar(@f2s)-1){
				my $f2=$f2s[$if2];
				my @ns2=unsaturatedNeighbors($f2->getIsomer(0),1);
				if(scalar(@ns2)>1){

					# Fragment f2 has an alpha hydrogen
					# adjacent to an unsaturated atom.
					print "Hydrogen transfer to fragment $if1 from its $if2 th sibling\n";

					# Prepare an H atom to transfer from f2 to f1
					my $a=$ns1[0];
					my $transName = $ns2[1]->name;
					my @temp= $ns2[1]->bonds;
					my $transBond = $temp[0];
					my $cap=Chemistry::Atom->new(symbol=>'H',name=>$transName);
					my $ic=Chemistry::InternalCoords->new($cap,$a,1.0);
					if(scalar(@ns1)>1){$ic=Chemistry::InternalCoords->new($cap,$a,1.0,$ns1[1],120.0);}
					if(scalar(@ns1)>2){$ic=Chemistry::InternalCoords->new($cap,$a,1.0,$ns1[1],120.0,$ns1[2],109.5);}
					$cap->coords($ic->cartesians);
					my $capBond=Chemistry::Bond->new(atoms=>[$a,$cap]);

					# Transfer H atom from molecule m1 to molecule m2 
					$f2->getIsomer(0)->delete_bond($transBond);
					$f2->getIsomer(0)->delete_atom($ns2[1]);
					$f1->getIsomer(0)->add_atom($cap);
					$f1->getIsomer(0)->add_bond($capBond);

					# Check if these molecules are unique, deep copy them if so
					my $m1b = updateUniqueMolecules($f1->getIsomer(0),$umols_ref,$unames_ref,1);
					my $m2b = updateUniqueMolecules($f2->getIsomer(0),$umols_ref,$unames_ref,1);

					# Undo our changes to molecules m1 and m2: if necessary, those 
					# changes are preserved by the deep copy 
					$f1->getIsomer(0)->delete_bond($capBond);
					$f1->getIsomer(0)->delete_atom($cap);
					$f2->getIsomer(0)->add_atom($ns2[1]);
					$f2->getIsomer(0)->add_bond($transBond);

					# Build lists of molecules and heavyAtoms for this fragmentation pattern
					# This fragmentation pattern is identical to that producing f1, but with 
					# m1 and m2 updated by hydrogen transfer 
					my @newmols=($m1b,$m2b);
					my @newHeavies=(heavyAtoms($m1b),heavyAtoms($m2b));
					foreach my $jf2(0..scalar(@f2s)-1){
						if($jf2 ne $if2){
							push(@newmols,$f2s[$jf2]->getIsomer(0));
							push(@newHeavies,$f2s[$jf2]->{heavy});
						}
					}
					# Update the sibling fragments 
					addSiblingFragments(\@newmols,\@newHeavies,$fragments_ref);
				}
			}
		}
	}
}

sub updateSiblingFragments{
	my ($newmols_ref,$umols_ref,$unames_ref,$heavies_ref,$frags_ref)=@_;
	my @newmols=@$newmols_ref;
	my @newmols2=();
	foreach my $i(0..scalar(@newmols)-1){
		push(@newmols2,updateUniqueMolecules($newmols[$i],$umols_ref,$unames_ref,1));
	}
	addSiblingFragments(\@newmols2,$heavies_ref,$frags_ref);
}
	
					
					
			

sub unsaturatedNeighbors{
	# Return the first coordinatively unsaturated Atom in a Molecule and a list of neighboring Atoms
	# Type Behavior 
	# 0    All nearest neighbors 
	# 1    H atom nearest neighbors 
	# 2    H atom second nearest neighbors 
	my ($m,$type)=@_;
	if ($type<0 or $type>2){die "Bad type in unsaturatedNeighbors";}
	my @ret=();
	my $done=0;
	foreach my $a($m->atoms()){ 
		if($done<1 and $a->calc_implicit_hydrogens>0){
			$done=1;
			push(@ret,$a);
			foreach my $b($a->neighbors()){
				if($type eq 0 or ($type eq 1 and $b->Z eq 1)){
					push(@ret,$b);
				}
				if($type eq 2){
					foreach my $c($b->neighbors($a)){
						if($c->Z eq 1){
							push(@ret,$c);
						}
					}
				}
			}
		}
	}
	return(@ret);
}

sub updateUniqueMolecules{
	# Determine if a new molecule is in a set of unique molecules.
	# Add to the unique molecules if so, replace the "new" molecule if not 
	# Add option to add a clone 
	my($moli,$umols_ref,$unames_ref,$clone)=@_; 
	my $ret = $moli;
	if(scalar($moli->atoms())>0){
		my $found=0;
		my $namei=moleculeName($moli);
		foreach my $j(0..scalar(@$umols_ref)-1){
			if($namei eq ${$unames_ref}[$j]){
				$found=1;
				$ret = ${$umols_ref}[$j];
			}
		}
		if($found<1){
			if($clone>0){
				#print "Cloning molecule ".$moli->print("%f\n");
				$ret = $moli->safe_clone;
				#print "Result is ".$ret->print("%f\n");
				push(@$umols_ref,$ret);
			}
			else{
				push(@$umols_ref,$moli);
			}
			push(@$unames_ref,$namei);
		}
	}
	return($ret);
}

sub BoltzmannFactor{
	# Given an energy difference in kcal/mol , and either a temperature in
	# Kelvin or a collision energy in eV, return the corresponding
	# Boltzmann factor 
	my ($Ekcalmol,$TK,$EceV)=@_;
	if($TK<=0 and $EceV<=0){die "BoltzmannFactor requires either a temperature or a collision energy\n";}
	#my $Ekcalmol=$Eau*627.5095;
	my $EkJmol = $Ekcalmol*4.185; 
	my $EJmol = $EkJmol*1000.;
	my $RJmolK = 8.314 ; 
	my $exponent = 0.0; 
	if($TK>0){$exponent=($EJmol/($RJmolK*$TK))};
	if($EceV>0){
		my $EcHartree= $EceV/27.211;
		my $Eckcalmol= $EcHartree*627.5095;
		$exponent = $Ekcalmol/$Eckcalmol;
	}
	return(exp(-$exponent));
}

sub mult{
	# Given a Molecule, return its lowest possible multiplicity 
	my $mol= $_[0];
	my $nelec=-$mol->charge();
	foreach my $a($mol->atoms()){$nelec=$nelec+$a->Z;}
	my $m=$nelec%2+1;
	return($m);
}
sub moleculeName{
	# Return a unique name for each Molecule: 
	# Charge_Stoichiometry_SMILES 
	# 2015.10.06 SMILES strings are super expensive, use 
	# numbers of CC,CH,CO bonds  
	my $m=shift;
	my $ret=$m->sprintf("%q_%f");

	# Count each number of bond types 
	my $nCC=0;my $nCH=0;my $nCO=0; my $nOH=0; my $nCN=0; my $nNH=0;
#	foreach my $bond($m->bonds()){
#		my @as=$bond->atoms();
#		my @Zs=sort($as[0]->Z,$as[1]->Z);
#		if($Zs[0] eq 1 and $Zs[1] eq 6){ $nCH+=1;}
#		if($Zs[0] eq 1 and $Zs[1] eq 8){ $nOH+=1;}
#		if($Zs[0] eq 6 and $Zs[1] eq 6){ $nCC+=1;}
#		if($Zs[0] eq 6 and $Zs[1] eq 8){ $nCO+=1;}
#		if($Zs[0] eq 6 and $Zs[1] eq 7){ $nCN+=1;}
#		if($Zs[0] eq 1 and $Zs[1] eq 7){ $nNH+=1;}
#	}
#	# Count coordination numbers of each C atom 
#	my @nCs=(0.,0.,0.,0.,0.);
#	foreach my $a($m->atoms()){
#		if($a->Z eq 6){
#			my $v=$a->valence;
#			$nCs[$v]=$nCs[$v]+1;
#		}
#	}
#	my $name = $ret.sprintf("_%d_%d_%d_%d_%d_%d__%d_%d_%d",$nCC,$nCH,$nCO,$nOH,$nCN,$nNH,$nCs[4],$nCs[3],$nCs[2]);

	#my $s =$m->sprintf("%s");
	my $s =$m->print(format=>'smiles',unique=>1,auto_number=>1); 
	$s=~s/:[0-9]*\]/]/g;
	my $name = $ret."_".$s;
	
	return($name);
}
sub heavyAtoms{
	# Return a sorted list of the non-hydrogen atom Names in a Molecule
	# Used to determine degeneracy of Fragments 
	# For simplicity, use the Gaussian atom ordering starting from index 1, not 0 
	my $mol = $_[0];
	my @ret=(); 
	foreach my $atom($mol->atoms()){
		if($atom->Z>1){
			push(@ret,($atom->name)+1);
		}
	}
	my @ret2 = sort(@ret);
	my $ret3='';
	foreach my $a(@ret2){$ret3=$ret3."$a-";}
	return $ret3;
}

1;
