##### Module to adjust sequence file #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package AdjustSequence;
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
@EXPORT_OK = qw(checkSequence);


sub checkSequence {

    my($Sequence_dir)=(shift);
    my($Adjust_sequence)=(shift);
    my(@SequenceFile)=@{(shift)}; 
    my($seq_type)=(shift);
  
    my(@array_sequence)=();
              
    if($Adjust_sequence eq "YES"){              
             
           foreach(@SequenceFile){     
                                                   
                    my $file_path1=$Sequence_dir."/".$_;
                    my($name)=basename($file_path1);
    
               if($name=~/\.fasta/ or $name=~/\.fa/ or $name=~/\.fas/){

                    my $seqio_obj = Bio::SeqIO->new(-file => "$file_path1", -format => "fasta");
                    $name=~s/\.(\w+)//g; 

                    while(my $seq_obj= $seqio_obj->next_seq){

                        my $seq=$seq_obj->seq;
                           $seq=~s/\s+//g;
                           $seq=~s/\n+//g;

                        my $seq_id=$seq_obj->display_id;
                           $seq_id=~s/\|$//g;

                        if($seq_id=~/\|/){
                           my @seq_id=split(/\|/,$seq_id);
                          
                           if($seq_id=~/^[(sp)|(tr)]/){
                              $seq_id=$seq_id[1];
                           }else{
                              $seq_id=pop(@seq_id);
                           }
                         }
                           
                        my $seq_length=length($seq);

                        my $seq_detail="$seq_id\t$name\t$seq_type\t$seq_length\t$seq";                           
                        push(@array_sequence,$seq_detail);                                          
                    }             

               }
          }           
    }

 return(\@array_sequence);
}
