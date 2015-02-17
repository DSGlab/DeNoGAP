#### MODULE TO PREDICT CORE GENES FROM THE CLUSTER #####

#!/usr/bin/perl
package CoreGenome;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Cluster;
use Bio::Seq;
use Bio::SeqIO;
use SQLiteDB;

    my($cluster_file)=(shift);
    my($db_dir)=(shift);
    my($db_name)=(shift);    
    my($out_dir_seq)=(shift);
    my($out_dir_aln)=(shift);
    my($core_alignment_file)=(shift);
    my($threshold)=(shift);
    my($seq_type)=(shift);
    my($include_outgroup)=(shift); 

    if(!defined($include_outgroup)){
      $include_outgroup=0;
    }

    mkdir($out_dir_seq);
    mkdir($out_dir_aln);

    ##### Read Cluster File #### 
    my($ortholog_cluster,$status)=Cluster::read_cluster($cluster_file,'xxx');
    
    ##### Create List of Strain in Cluster ###
    my %genome=();
    my %outgroup=();
   
    my($get_outgroup)="SELECT DISTINCT(abbreviation),outgroup from OrganismInfo";
    my ($outgroup_record)=SQLiteDB::get_record($db_dir,$db_name,$get_outgroup);

    if(scalar($outgroup_record)>1){           
       foreach my $row(@{$outgroup_record}){           
          my $taxa=shift(@{$row});
          my $seq_type=shift(@{$row});
          chomp($seq_type);
      
          if($seq_type eq "Yes" and $include_outgroup=~/NO/i){
            $outgroup{$taxa}=$taxa;
          }else{
            $genome{$taxa}=$taxa;
          }
       }
    }
   
    my $count_genome=keys(%genome);

    my $threshold_genome=int((($threshold/100)*$count_genome));

    print "Find Core Gene\n\n";
    
    my($core_alignment)=find_core_gene($ortholog_cluster,$db_dir,$db_name,\%genome,$out_dir_seq,$out_dir_aln,$threshold_genome,$seq_type);

    my %core_genome_alignment=%{$core_alignment};

    open(CORE_ALN,">$core_alignment_file");

    print "Printing out Core Genome Alignment\n\n";

    foreach my $genome_name(keys %core_genome_alignment){
       print CORE_ALN ">$genome_name\n$core_genome_alignment{$genome_name}\n";
    }
    close CORE_ALN;


#### Find Core Gene ####

sub find_core_gene{

  my(%ortholog_cluster)=%{(shift)};
  my($db_dir)=(shift);
  my($db_name)=(shift);
  my(%genome)=%{(shift)};
  my($out_dir_seq)=(shift);
  my($out_dir_aln)=(shift);
  my($threshold)=(shift);
  my($seq_type)=(shift);
  
  my $group_concatenated=0;
  my %core_genome_alignment=();

  foreach my $group_id(keys %ortholog_cluster){

        my $group_gene=$ortholog_cluster{$group_id};
        my @gene_id=split("\t",$group_gene);

        if(scalar(@gene_id)==1){
          next;
        }
     
        my $list_gene_id='';
        
        my $list_protein_id='';

        foreach my $gene_id(@gene_id){

           $gene_id=~/(\w+)(\|)(.+)/;

           my $id=$3;
           my $genome=$1;

           if($gene_id=~/\~/){
             next;
           }

           if(!defined($genome{$genome})){
             next;
           }
           
           if($seq_type eq "nucl"){
             $id=~s/\_1$//g;
           }
           
           $list_gene_id=$list_gene_id.",'".$id."'";  
           
           $list_protein_id=$list_protein_id.",'".$gene_id."'";            
        }
        
        $list_gene_id=~s/^\,//g;
        
        $list_protein_id=~s/^\,//g;

        my $get_sequence='';

        my $family_sequence=undef;

        my $homolog_sequence=undef;

        ####### search alignment in Homolog Alignment Table ###########

        if($seq_type eq "protein"){

              $get_sequence="Select * from Alignment WHERE feature_id IN ($list_protein_id)";  

              $homolog_sequence=SQLiteDB::get_record($db_dir,$db_name,$get_sequence); 
 
              if(scalar(@{$homolog_sequence})==0){

                   $get_sequence="SELECT * from ProteinSequence WHERE protein_id IN ($list_gene_id)";

                   $family_sequence=SQLiteDB::get_record($db_dir,$db_name,$get_sequence);
              }
           
        }elsif($seq_type eq "nucleotide"){

              $get_sequence="SELECT * from NucleotideSequence WHERE nucleotide_id IN ($list_gene_id)";

              $family_sequence=SQLiteDB::get_record($db_dir,$db_name,$get_sequence);
        }

        my %length_seq=();
        my %group_seq=();

        my %aligned_seq=();
        my $len_aln=0;

        if($homolog_sequence and scalar(@{$homolog_sequence})>=1){
 
            foreach my $row(@{$homolog_sequence}){  

                my @data_column=@{$row};
                
                my $seq_id=$data_column[0]; 
                my $genome_name=$data_column[1];      
                my $aln_seq=$data_column[5];

                my $seq=$aln_seq;
                   $seq=~s/\-//g;
 
                my $len=length($seq);

                $len_aln=length($aln_seq);

                if(!defined($length_seq{$genome_name})){

                    $aligned_seq{$genome_name}=$aln_seq;
                    $length_seq{$genome_name}=$len;

                }elsif($len>$length_seq{$genome_name}){

                    $aligned_seq{$genome_name}=$aln_seq;
                    $length_seq{$genome_name}=$len;
                }    
            }

        }elsif($family_sequence and scalar($family_sequence)>1){ 
 
             foreach my $row(@{$family_sequence}){           
                     my @data_column=@{$row};     
                     my $genome_name=$data_column[2];
                     my $gene_id=$genome_name."|".$data_column[1];
                     my $seq=$data_column[5];  
                     my $len=length($seq);
                        
                     if(!defined($length_seq{$genome_name})){

                        $group_seq{$genome_name}=$seq;
                        $length_seq{$genome_name}=$len;

                     }elsif($len>$length_seq{$genome_name}){

                        $group_seq{$genome_name}=$seq;
                        $length_seq{$genome_name}=$len;
                     }
             }
        }


       my $genome_count=keys(%length_seq);

       if($genome_count>=$threshold){

         print "$group_id\t$genome_count\n";

         if($family_sequence and scalar($family_sequence)>1){

            print_core_group($group_id,\%group_seq,$out_dir_seq);

            align_core_group($group_id,$out_dir_seq,$out_dir_aln);

            my $seqio_obj = Bio::SeqIO->new(-file => "$out_dir_aln/$group_id.afa", -format => "fasta");

            while(my $seq_obj= $seqio_obj->next_seq){
    
                 my $seq=$seq_obj->seq;
                 my $seq_id=$seq_obj->display_id;
                 $len_aln=length($seq);
                $aligned_seq{$seq_id}=$seq;
            }
         }

         my($core_genome_alignment)=concatenate_alignment($group_id,\%genome,\%aligned_seq,\%core_genome_alignment,$len_aln);

         %core_genome_alignment=%{$core_genome_alignment};

         $group_concatenated++;
       }

    }
   return(\%core_genome_alignment);
}

##### Print Core Gene in Fasta File #####

sub print_core_group{

   my($group_id)=(shift);
   my(%group_seq)=%{(shift)};
   my($out_dir_seq)=(shift);

   open(CORE_SEQ,">$out_dir_seq/$group_id.fasta");

   foreach my $genome_name( keys %group_seq){

          print CORE_SEQ ">$genome_name\n$group_seq{$genome_name}\n";   
   }
   close CORE_SEQ;
}

#### Align Core Group #####

sub align_core_group{

   my($group_id)=(shift);
   my($out_dir_seq)=(shift);
   my($out_dir_aln)=(shift);

   my $infile="$out_dir_seq/$group_id.fasta";
   my $outfile="$out_dir_aln/$group_id.afa";

   system("kalign -in  $infile -output $outfile -f fasta -q");
}

#### Concatenate Core Group Alignment ###

sub concatenate_alignment{

   my($group_id)=(shift);
   my(%genome)=%{(shift)};
   my(%aligned_seq)=%{(shift)};
   my(%core_genome_alignment)=%{(shift)};
   my($len_aln)=(shift);


   foreach my $genome_name(keys %genome){

        if($aligned_seq{$genome_name} and $core_genome_alignment{$genome_name}){

            $core_genome_alignment{$genome_name}=$core_genome_alignment{$genome_name}.$aligned_seq{$genome_name};

        }elsif($aligned_seq{$genome_name} and !defined($core_genome_alignment{$genome_name})){

            $core_genome_alignment{$genome_name}=$aligned_seq{$genome_name};
        
        }elsif(!defined($aligned_seq{$genome_name}) and $core_genome_alignment{$genome_name}){

             for(my $i=0;$i<$len_aln;$i++){             
                 $core_genome_alignment{$genome_name}=$core_genome_alignment{$genome_name}."-";    
             }
        }elsif(!defined($aligned_seq{$genome_name}) and !defined($core_genome_alignment{$genome_name})){

             $core_genome_alignment{$genome_name}="-";

             for(my $i=0;$i<$len_aln-1;$i++){             
                 $core_genome_alignment{$genome_name}=$core_genome_alignment{$genome_name}."-";     
             }
        } 
   }

   return(\%core_genome_alignment);
}
