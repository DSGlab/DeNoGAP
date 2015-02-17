##### Module to get sequence ID in hash table #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package SequenceHash;
use strict;
use Exporter;
use File::Basename;
use File::Copy;
use Bio::Seq;
use Bio::SeqIO;
use SQLiteDB;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(getSequenceFeature sendSequenceDBToFile);


sub getSequenceFeature {

    my($input_seq_dir)=(shift);
   
    opendir(SEQ_DIR,$input_seq_dir);
    my(@SequenceFile)=readdir(SEQ_DIR); 

    my %SeqLenTable=();
    my %IDHashTable=();
    my %Sequence=();

    foreach(@SequenceFile){ 

       if($_=~/^\.+$/ or $_=~/\~/){next;}
 
       my $file_path=$input_seq_dir."/".$_;

       my $genome_name=$_;
          $genome_name=~s/(\.)(\w+)$//g;

        my $seqio_obj = Bio::SeqIO->new(-file => "$file_path", -format => "fasta");

        while(my $seq_obj= $seqio_obj->next_seq){

             my $seq=$seq_obj->seq;
             my $seq_id=$seq_obj->display_id;

             $seq=~s/\s+//g;

             chomp($seq_id);
             chomp($seq);

             $IDHashTable{$genome_name}->{$seq_id}=$seq_id;
             $Sequence{$genome_name}->{$seq_id}=$seq;
             $SeqLenTable{$genome_name}->{$seq_id}=length($seq);          
        }        
    }
    return(\%SeqLenTable,\%IDHashTable,\%Sequence);
}

sub sendSequenceDBToFile{

    my($db_dir)=(shift);
    my($db_name)=(shift);
    my($seq_dir)=(shift);

    my $sql_stmt="SELECT * from ProteinSequence";

    my($row_data)=SQLiteDB::get_record($db_dir,$db_name,$sql_stmt);
    
    if(scalar(@{$row_data})>0){
 
      my $prev_genome='';

      foreach my $row(sort{$a->[2] cmp $b->[2]}@{$row_data}){
          
          my @seq_feature=@{$row};
          
          my $genome_name=$seq_feature[2];
          my $seq_id=$seq_feature[1];
          my $seq=$seq_feature[5];

          if($genome_name ne $prev_genome){           
               if($prev_genome ne ''){
                 close PROT_FILE;
               }
               open(PROT_FILE,">$seq_dir/$genome_name.fasta");               
               print PROT_FILE ">$genome_name|$seq_id\n$seq\n"; 
               $prev_genome=$genome_name;
          }else{
              print PROT_FILE ">$genome_name|$seq_id\n$seq\n";
          }
      }
    }
   close PROT_FILE;
}
