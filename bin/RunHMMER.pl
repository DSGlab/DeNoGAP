#!/usr/bin/perl -w
###### ABOUT: This Script run hmmscan for protein sequences ############
###### AUTHOR:Shalabh Thakur###################################################################
###### DATE:15-MAY-2013########################################################################

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Parallel::ForkManager;
use File::Basename;
use File::Path qw(remove_tree);
use Hash::Merge qw( merge );
use List::Util qw(sum);
use Bio::SeqIO;
use SQLiteDB;

 my($genome_name)=(shift);
 my($hmm_program)=(shift);
 my($db_dir)=(shift);
 my($db_name)=(shift);
 my($DB_file)=(shift);
 my($tmp_dir)=(shift);
 my($db_size)=(shift); 
 my($chunk_size)=(shift);
 my($hmmer_opt)=(shift);
 my($cpu_core)=(shift);
   
 
###### get sequences for the query genome ###### 

my $seq_count=0;
my $chunk_count=1;
my %chunk=();

  my $get_sequence_sql="Select * from ProteinSequence where genome_abbreviation='$genome_name'";
          
  my $query_protein_sequences=SQLiteDB::get_record($db_dir,$db_name,$get_sequence_sql);
          
          if(scalar(@{$query_protein_sequences})>0){
          
                foreach my $row(@{$query_protein_sequences}){
               
                        my $protein_index=shift(@{$row});
                        my $protein_id=shift(@{$row});
                        my $genome_abbreviation=shift(@{$row});
                        my $seq_type=shift(@{$row});
                        my $seq_length=shift(@{$row});
                        my $seq=shift(@{$row});
                                                
                        if($seq_count>$chunk_size or $seq_count==0){
                        
                              open(CHUNK_SEQ,">$tmp_dir/Chunk_$chunk_count.fasta");
                            
                              $seq_count=1;
                              $chunk{$chunk_count}=$chunk_count;
                              $chunk_count++;
                         } 
                         
                         print CHUNK_SEQ ">$genome_abbreviation|$protein_id\n$seq\n";

                          $seq_count++;

                          if($seq_count>$chunk_size){
                              close CHUNK_SEQ;
                          }                                      
                }
          }


my $fork=Parallel::ForkManager->new($cpu_core);  

foreach my $chunk(keys %chunk){

    $fork->start and next;

    my $full_out=$tmp_dir."/"."full_".$hmm_program."_Chunk_".$chunk.".out";
    my $dom_out=$tmp_dir."/"."dom_".$hmm_program."_Chunk_".$chunk.".out";

    my $input_file="$tmp_dir/Chunk_$chunk.fasta";

    $hmmer_opt=$hmmer_opt." -o ".$full_out;
    $hmmer_opt=$hmmer_opt." --domtblout ".$dom_out;

    if($hmm_program=~/phmmer/){
      system("$hmm_program $hmmer_opt $input_file $DB_file");
    }elsif($hmm_program=~/hmmscan/){
      system("$hmm_program $hmmer_opt $DB_file $input_file");
    } 
   $fork->finish; 
}
$fork->wait_all_children;

foreach my $chunk(keys %chunk){

    my $full_out=$tmp_dir."/"."full_".$hmm_program."_Chunk_".$chunk.".out";
    my $dom_out=$tmp_dir."/"."dom_".$hmm_program."_Chunk_".$chunk.".out";

    my $genome_full=$tmp_dir."/"."full_".$hmm_program."_".$genome_name.".out";
    my $genome_dom=$tmp_dir."/"."dom_".$hmm_program."_".$genome_name.".out";

    system("cat $full_out >> $genome_full");
    system("cat $dom_out >> $genome_dom");

    unlink($full_out);
    unlink($dom_out);
}


