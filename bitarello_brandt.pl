#!/usr/bin/perl


#########################################################################################################################################
#  Program: final project (BIO5793)       Created in: 17.10.2012									#		
#																	#
#                             HLA consensus: a perl script which creates consensus sequences for HLA allelic lineages and other things. #
#                                                                                                                                       #
#  Author(s): Barbara Bitarello & Debora Brandt                                                                                         #
#                                                                                                                                       #
#  Description: A script to read in an alignment file of an HLA loci and a table with allele information for the bi-allelic             #
#				locus. HLA alleles have a very well established nomenclature method. The first two digits indicate the allelic          #
#				lineage to which an allele belongs. Differences in the next two digits inform that a certain group of alleles           #
#				are from the same lineages, but have non-synonymous substitutions among them. Finally, the next couple of digits        #
#				inform that a group of alleles only have differences that are synonymous.                                               #
#                                                                                                                                       #
#				Since common sequencing of exons 2 and 3 of HLA loci only allows the discrimination of these alleles up to              #
#				the second set of digits, our script reads in a table provided by the user, informing the allele names up to            #
#				the 2nd set of digitis. Our script then gathers all IMGT alleles compatible with each group contemplated in the         #
#				user table and produces a consensus sequence. Differences among alleles from the same group are then masked as          #
#				'?' in the consensus for each group. Finally, these consensus sequences are joined as an alignment in fasta             #
#				format. According to the cut-off value informed by the user, the script will eliminate the positions where              #
#				'?' is above that cut-off (considering all sequences). By varying the cutoff values and analysing the results           #
#				of donwstream population genetic analyses with these different cut-offs, the user can make a better informed            #
#				decision as to how to deal with the HLA allele ambiguity.                                                               #
#                                                                                                                                       #
#                                                                                                                                       #
#                                                                                                                                       #
#  Usage: bitarello_brandt.pl                                                                                                           #
#                                                                                                                                       #
#  Last modified: 05.12.2012                                                                                                            #
#########################################################################################################################################

#Algorithm:

### 1. Read in (interactively) the options given by the user: an alignment in fasta format (F) or IMGT format (I);  name of the alignment file;  name of the table containing the alleles for each individual; name of the file where final alignment will be written;  threshhold value.
### 2. Read in the table and create an array with all the exclusive allele names in the table. Create another array with all the occurrences of the alleles in the table (with repetition). 
### (subroutine POP_TABLE).
### 3. If the initial alignment is in fasta format (option F):
### 3.1. create a hash with the sequences from the input fasta file that have the first two sets of digits equal to those of any allele present in the table.
### 4. If the initial alignment is in IMGT format (option I):
### 4.1. convert this alignment to fasta format (subroutine PARSING_IMGT);
### 4.2. execute item 3.1.
### 5. For each allele from the table which represents a different allelic group (first 4 digits are unique):
### 5.1. open an output file names group_X.fasta, where "X" are the first 4 digits of the allele name, as in the table;
### 5.2. write (in fasta format) all the alleles pertaining to each of these groups into the corresponding file.
### 6. For each of the files created in 5.1, open and create a consensus sequence for that group of alleles.
### 7. Create an alignment in fasta, using the ouput file name informed by the user (item 1), with all the consensus sequences.
### 8. Once this fasta of consensus sequences is created, use the threshhold value informed by the user to edit this alignment.


#####################
##Prerequisites:

##1. A tab delimited table with four collumns, where the first one informs the population from which the sample was taken, the second one informs and id for each individual, and the 3rd and 4th collumns inform the names of the alleles based on sequencing of exons 2 and 3.
##2. A multiple alignment fasta file, or a file in IMGT "nuc" format (i.e., complete coding sequence). In the 2nd case, the user informs this is the case and the script converts this into fasta format.
## see ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/
##3. Also, the user should inform, interactively, a name for the final fasta output file and the threshhold value.

#######################

use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Clustalw;

##declare variables:

my $al; 
my $gt;
my $out;
my $name;
my $name2;
my $seq;
my %seqs_imgt;
my @fields;
my @als;
my @als2;
my $subkey;
my $cons;
my $als;
my $key;
my $alignout;
my %cons;
my $all_seqs_obj;
my $out_obj;
my $seq_obj;
my $seq_obj2;
my @params;
my $factory;
my $aln;
my $answ1;
my $answ2;
my $in;
my $out_imgt;
my @seq;
my $allele;
my $sequence;
my %alleles;
my @nome;
my @refnucs;
my $popal;
my $thr;

##########################################################################################################################################

### INPUT #####

#Read in (interactively) the options given by the user: an alignment in fasta format (F) or IMGT format (I);  
#name of the alignment file;  name of the table containing the alleles for each individual; name of the file 
#where final alignment will be written;  threshhold value.

print "Do you have an IMGT CDS alignment (I) or an aligned fasta file (F)? (I/F)\n";

$answ1=<STDIN>;
chomp ($answ1);

unless($answ1 eq "I" || $answ1 eq "F" ||$answ1 eq "i" || $answ1 eq "f"){
	die "Sorry, wrong answer! There are only two options!\n";
	}
		
print "Type the name of the file:\n";

$al=<STDIN>;
$al="B_nuc.txt";

print "Type in the name of the table with the alleles, etc:\n";

$gt=<STDIN>;
chomp ($gt);

print "What is the desired cutoff threshold for '?'?:\n";
$thr=<STDIN>;
chomp ($thr);

print "Type in the name you want for your final fasta file, after the threshold is applied:\n";
$out=<STDIN>;
chomp($out);

##########################################################################################################################################
#######SCRIPT##############
### Read in the table and create an array with all the exclusive allele names in the table. Create another array with all the occurrences of the alleles in the table (with repetition). 
### Subroutines at the bottom of the script.

POP_TABLE();

### 3. If the initial alignment is in fasta format (option F): 
### 3.1. create a hash with the sequences from the input fasta file that have the first two sets of digits equal to those of any allele present in the table.
### 4. If the initial alignment is in IMGT format (option I):
### 4.1. convert this alignment to fasta format (subroutine PARSING_IMGT);
### 4.2. execute item 3.1.

if ($answ1 eq "F" || $answ1 eq "f" ){

	SUB_FASTA($al);
}

if ($answ1 eq "I" ||$answ1 eq "i" ){

	PARSING_IMGT();   
	SUB_FASTA($out_imgt);  
}

##now we have the hash seqs_imgt with all the sequences corresponding to the alleles represented in the table.


### 5. For each allele from the table which represents a different allelic group (first 4 digits are unique):
### 5.1. open an output file names group_X.fasta, where "X" are the first 4 digits of the allele name, as in the table;
### 5.2. write (in fasta format) all the alleles pertaining to each of these groups into the corresponding file.
### 6. For each of the files created in 5.1, open and create a consensus sequence for that group of alleles.
### 7. Create an alignment in fasta, using the ouput file name informed by the user, with all the consensus sequences.


foreach $als (@als){
	
    my $consensus = "group".$als.".fasta";
    
    open(CONS, ">$consensus") or die "can't open $consensus!\n";

    
	foreach $key (keys(%seqs_imgt)){
        
        if($key=~/\w+\*(\d+)\:(\d+)/){
            $subkey=$1.$2;

        }

		if ($als eq $subkey){
			print CONS ">$key\n";
			print CONS "$seqs_imgt{$key}\n";
            
		}
	}
    
    close CONS;
    
    
###consensus for each group
    
    open(CONS2,"$consensus") or die "can't open $consensus!\n";
    
    my $fasta = Bio::AlignIO->new(-fh => \*CONS2, -format => 'fasta');
    
	my $ali = $fasta->next_aln();
    
    $cons{$als} = $ali->consensus_iupac();
    
    close CONS2;

}

open(POPIU, ">pop_consensus_iupac_align.fasta") or die "can't open pop_consensus_iupac_align.fasta!\n";

$popal="pop_consensus_align.fasta";
open(POP, ">$popal") or die "can't open $popal\n";

foreach $als (@als2){
    print POPIU ">$als\n$cons{$als}\n";
    my $seqgaps = $cons{$als};
    $seqgaps =~ tr/[a-z]/[A-Z]/;
    $seqgaps =~ tr/(R|Y|K|M|S|W|B|D|H|V|N|)/\?/;
    print POP ">$als\n$seqgaps\n";

}

close POPIU;
close POP;

### 8. Once this fasta of consensus sequences is created, use the threshhold value informed by the user to edit this alignment.

APPLY_THR($popal, $thr);

exit;


########################################################SUBROUTINES#################################################
#######################################################
#FIRST SUBROUTINE: Create alignment fasta file (only used if $answ1 eq "I").

sub PARSING_IMGT{


	open (in, "$al") or die "can't open $al\n";
	$out_imgt = "imgt".$out;
	open (OUT, ">$out_imgt") or die "can't open $out!\n";

		while (<in>){
  
		  chomp;
			  if ($_ =~ m/(^\s[A-Z]\*)/){ #if the line begins with A, B or C
  	
			    	$_ =~ s/^\s//;
				    $_ =~ s/\|//g;
				    $_ =~ s/[\s]+/K/;
				    $_ =~ s/K/\t/g;
				    $_ =~ s/ //g; 
    
    
    				($allele, $sequence) = split (/\t/, $_);
    
				    if (exists($alleles{$allele})) {
				      $alleles{$allele} .= $sequence;
				    }
				    else{
				      $alleles{$allele}  = $sequence;
    				}
    
  				}
		}



	foreach $allele(sort(keys(%alleles))){
			if($allele !~ m/(N|S|Q|L)/g){
			  push (@nome,$allele);
			  push (@seq,$alleles{$allele});
			}
	}

	$seq[0]=~ s/\./-/g;  #Replace . by - in the reference sequence
	print OUT ">$nome[0]\n$seq[0]\n";        #printing reference sequence

	my $y = @nome;
	@refnucs = split (//,$seq[0]);     # Splitting all nucleotides in the reference sequence

	my $x = @refnucs;                    
	for (my $counter=1; $counter<$y ; $counter++){
	  for (my $counter1=0; $counter1<$x; $counter1++) {
    	if (substr($seq[$counter],$counter1,1) eq "-"){
	      substr($seq[$counter],$counter1,1) = $refnucs[$counter1];
    		}
	    elsif (substr($seq[$counter],$counter1,1) eq "*" or substr($seq[$counter],$counter1,1) eq '.'){
    	  substr($seq[$counter],$counter1,1) = '-';
    		}
  		}
	print OUT ">$nome[$counter]\n$seq[$counter]\n";
	}
	return ($out_imgt);
}   # End of subroutine

############################################################################################
###############################################################
# SECOND SUBROUTINE: Read in the table and create an array with all the exclusive allele names in the table. Create another array with all the occurrences of the alleles in the table (with repetition). 

sub POP_TABLE{


	open(GT, "<$gt");
    
	readline(GT); 
	
    while(<GT>){
        
		chomp;  
		@fields=split(/\s/, $_);
		push(@als2, $fields[2]);    # @als2 contains all occurences of the alleles in the table.
		push(@als2, $fields[3]);
		
        if(!grep /$fields[2]/, @als){   #@als contains all the exclusive allele names in the table.
            push(@als,$fields[2]);
        }
		
        if(!grep /$fields[3]/, @als){
            push(@als,$fields[3]);
		}
        
	
	}
    
	close GT;
	return(@als);
	return(@als2);
} # End of subroutine

############################################################################################
###############################################################
# THIRD SUBROUTINE: create a hash with the sequences from the input fasta file that have the first two sets of digits equal to those of any allele present in the table.

sub SUB_FASTA{
	my $a=@_[0];
	$all_seqs_obj = Bio::SeqIO->new(-file => "<$a",
								   	-format => 'Fasta'); 

	

		while ($seq_obj = $all_seqs_obj->next_seq()) {
    
	    $name= $seq_obj-> id();

	    $seq = $seq_obj-> seq();
		
    	# Compare sequences from inicial fasta file with those in the table
		    foreach $als (@als){
    
                # Select first 4 digits from allele name from fasta file
        	    if($name=~/^\w+\*(\d+)\:(\d+)/){                
            	    $subkey=$1.$2;
            		}
            		
        
	            if ($als eq $subkey){
    	            # Store in %seqs_imgt all IMGT alleles compatible with each group contemplated in the user table
					$seqs_imgt{$name}=$seq;
            		}
        		}
			
			}
			
    		return %seqs_imgt;
    		
} # End of subroutine

################################################################################################
#########
########
# FOURTH SUBROUTINE:  Once this fasta of consensus sequences is created, use the threshhold value informed by the user to edit this alignment.


sub APPLY_THR{
    
    my ($alfile, $thr) = @_;
    my @doubts;
    my $length;
    my @keep;
    my $nseq=0;
    my @gaps;
    print "The alignment positions with percentage of '?' above the given cutoff threshold are: ";
    
    open(AL, "<$alfile");
  
    while(<AL>){
        
        # Read only sequence lines
        if($_!~/\>/){
            
            $nseq++;
            
            my @nts=split("",$_);
            $length=@nts;
            
            # Count positions with "?" or gaps ("-")
            for(my$i=0;$i<$length; $i++){

                if ($nts[$i] eq "?"){
                    $doubts[$i]++;  # @doubts has the alignment lenght. Each position in @doubts stores the number of occurences of "?" for that position in the alignment. 
                }
                
                if ($nts[$i] eq "-"){
                    $gaps[$i]++;    # @gaps has the alignment lenght. Each position in @gaps stores the number of occurences of "-" for that position in the alignment. 
                }
                    
            }
        }
    }
    
    print "\n";
    
    close AL;

    # Compare number of "?" and "-" in each position of the alignment with threshold
    for (my $i=0; $i<$length; $i++){
       
        my $percent=($doubts[$i]/$nseq);
        
        # if the percentage of "?" in a given position ($i) is lower or equal to the given threshold...
        if($percent<=$thr){
            # ...and if $i is not a gap in all sequences
            if($gaps[$i]<$nseq){
                push @keep, $i; #keep $i position in the alignment
            }
        }
        
        else{   # this prints for the user the alignment positions that were excluded based on the threshold
        	print $i+1;
        	print "\t";
    	}
    }
    
    open(AL, "<$alfile");
    
    open(OUT, ">$out");
    
    while(<AL>){
        
        # Prints sequence name in output file
        if($_=~/\>/){
            print OUT $_;
        }
        
        # Prints in output file only the positions that were not excluded
        if($_!~/\>/){
                        
            my @nts=split("",$_);
            
            foreach my $pos(@keep){

                print OUT $nts[$pos];
            }
        }
    }
    close OUT;
    close AL;   
}# End of subroutine
