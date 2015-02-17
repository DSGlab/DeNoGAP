#!/usr/bin/perl -w
###### ABOUT: This Script builds HMM model for protein families ############
###### AUTHOR:Shalabh Thakur###################################################################
###### DATE:15-MAY-2013########################################################################

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Cluster;
use Parallel::ForkManager;
use File::Basename;
use File::Path qw(remove_tree);
use Hash::Merge qw( merge );
use List::Util qw(sum);
use Bio::AlignIO;
use Bio::SeqIO;
use SQLiteDB;

my $cluster_file=(shift);
my $group_for_query_file=(shift);
my $group_org_file=(shift);
my $db_dir=(shift);
my $db_name=(shift);
my $hmm_file_dir=(shift);
my $singleton_dir=(shift);
my $hmm_dbatabase=(shift);
my $singleton_database=(shift);
my $tmp_log=(shift);
my $process=(shift);
my $tmp="../tmp";

my %out_dir=("tmp_log"=>$tmp_log,"hmm_file_dir"=>$hmm_file_dir,"singleton_group"=>$singleton_dir);
my %processed_genome=();

my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,'xxxx');
my %homologue_group=%{$homologue_group};

open(QUERY_GROUP,"$group_for_query_file");
my @query_group=<QUERY_GROUP>;
close QUERY_GROUP;

my %group_for_query=();

foreach my $group_id(@query_group){
   chomp($group_id);
   $group_for_query{$group_id}=$group_id;
}

open(GENOME_LIST,"$group_org_file");
my @genome_list=<GENOME_LIST>;
close GENOME_LIST;

my %group_gene=();

foreach my $genome(@genome_list){
   chomp($genome);
   $group_gene{$genome}=$genome;
}

 #### Initiate parallel processes ####
 my $fork=Parallel::ForkManager->new($process);  
          
 foreach my $group_id(keys %group_for_query){
                
             $fork->start and next;
                
             my $gene_list=$homologue_group{$group_id};
             my @gene_list=split("\t",$gene_list);
            
              if(scalar(@gene_list)==1){      
                 $fork->finish;
              }  
                    
                  my $ref_processed_genome={};

                  print "$group_id\n";                     
                      
                  my %family_sequence=();
                
                  my $list_gene_id='';
               
                  foreach my $gene_id(@gene_list){
                     $gene_id=~/(\w+)(\|)(.+)/;
                     my $id=$3;
                     $list_gene_id=$list_gene_id.",'".$id."'";               
                  }
                
                  $list_gene_id=~s/^\,//g;
                
                  my($get_sequence)="SELECT * from ProteinSequence WHERE protein_id IN ($list_gene_id)";
                
                  my($family_sequence)=SQLiteDB::get_record($db_dir,$db_name,$get_sequence);
                
                  if(scalar(@{$family_sequence})>1){           
                        foreach my $row(@{$family_sequence}){           
                            my @data_column=@{$row};     
                            my $genome_name=$data_column[2];
                            my $gene_id=$genome_name."|".$data_column[1];
                            my $seq=$data_column[5];
                            
                            $family_sequence{$genome_name}->{$gene_id}=$seq;
                        }
                   }   
                    
                    mkdir("$out_dir{tmp_log}/$group_id");
                    my $hmm_file="$out_dir{tmp_log}/$group_id/$group_id.hmm";
                    my $hmm_return_code=-1;                   

                    if(-s "$out_dir{hmm_file_dir}/$group_id.hmm"){                                                         
                       $hmm_return_code=0;
                       system("cp $out_dir{hmm_file_dir}/$group_id.hmm $hmm_file");                                                       
                    }

                    if(-s "$out_dir{singleton_group}/$group_id.fasta"){
                        unlink("$out_dir{singleton_group}/$group_id.fasta");
                    }                    

                    $ref_processed_genome=DivergentSequenceForGroup(\%group_gene,\%family_sequence,\@gene_list,$group_id,$hmm_file,$hmm_return_code,\%out_dir,$tmp);                  
                     	                        
                 $fork->finish;         
    }
    $fork->wait_all_children;

        ####### SINGLETON GROUP ######
        open(SINGLETON_DB,">$singleton_database");

        print "SINGLETON GROUP\n";

        foreach my $group_id(keys %group_for_query){

             my $gene_list=$homologue_group{$group_id};
             my @gene_list=split("\t",$gene_list);

             if(scalar(@gene_list) eq 1){

                    print "$group_id\n";

                    my $gene_id=shift(@gene_list);  
                    chomp($gene_id);                     
                    $gene_id=~/(\w+)(\|)(.+)/;  
                    my $id=$3;  
                    my $genome=$1; 
                    my $gene_seq=undef;

                    my($get_sequence)="SELECT * from ProteinSequence WHERE protein_id IN ('$id')";
                
                    my($family_sequence)=SQLiteDB::get_record($db_dir,$db_name,$get_sequence);
                
                    if(scalar(@{$family_sequence})==1){           
                        foreach my $row(@{$family_sequence}){           
                            my @data_column=@{$row};     
                            my $genome_name=$data_column[2];
                            my $gene_id=$genome_name."|".$data_column[1];
                            $gene_seq=$data_column[5];                           
                        }
                    }   

                    if(!defined($gene_seq)){
                        print "Error: sequence not found for $gene_id in $genome\n"; exit;
                    }
 
                    my %genome=();
                    $genome{$genome}=$genome;   
                    $gene_seq=~s/\*$//g; 
             
                    print SINGLETON_DB ">$group_id\t$gene_id\n$gene_seq\n";   

                    $processed_genome{$genome}=$genome;                                                                                                                                                               
             }
        }

        close SINGLETON_DB;    

 %group_gene=%{merge(\%group_gene,\%processed_genome)};

 open(GENOME_LIST,">$group_org_file");
 foreach my $genome_name(keys %group_gene){
    print GENOME_LIST "$genome_name\n";
 }
 close GENOME_LIST;
       
###### Subroutin to calculate Genetic Distance for whole Group ##############

sub DivergentSequenceForGroup{

    my(%group_gene)=%{(shift)};
    my(%sequence)=%{(shift)};
    my(@gene_list)=@{(shift)};
    my($group_id)=(shift);
    my($hmm_file)=(shift);
    my($hmm_return_code)=(shift);   
    my(%out_dir)=%{(shift)};
    my($tmp)=(shift);    
    my $bin=$Bin;  
    my %processed_genome=(); 
    my $build_db=1;
  
      
    ###### if more than one sequence exist for the group, get sequences for each member of the group####    

    my %model_seq=(); 
    my %query_seq=();               
      
         foreach my $gene_id(@gene_list){
             chomp($gene_id);
             $gene_id=~/(\w+)(\|)(.+)/;   
             
             my $genome_name=$1;

             my $gene_seq=$sequence{$genome_name}->{$gene_id};           
             $gene_seq=~s/\*$//g;         

             if($gene_seq eq ''){
                  print "Error message: No Sequence for gene $gene_id\n";
             }

             if($group_gene{$genome_name}){                                       
                 $model_seq{$gene_id}=$gene_seq; 
             }else{
                 $processed_genome{$genome_name}=$genome_name; 
                 $query_seq{$gene_id}=$gene_seq;
             }           
          }                                  

             ####### create blast db for model ####      
             #  system("cat $out_dir{tmp_log}/$group_id/model_gene.fasta  $out_dir{tmp_log}/$group_id/query_gene.fasta > $out_dir{tmp_log}/$group_id/model_db.fasta");
             
             #  system("makeblastdb -in $out_dir{tmp_log}/$group_id/model_db.fasta -dbtype prot -logfile $out_dir{tmp_log}/$group_id/db_log");

               ####### Blast query with model database ####
             #  system("blastp -query $out_dir{tmp_log}/$group_id/query_gene.fasta -db $out_dir{tmp_log}/$group_id/model_db.fasta -evalue 1e-5 -out $out_dir{tmp_log}/$group_id/blast_result.txt -outfmt '6  qseqid sseqid pident qlen slen length mismatch gaps qstart qend sstart send  evalue bitscore'");

               ####### Measure Weight for the Query Sequence #####
               my $count_model_gene=keys(%model_seq);
               chomp($count_model_gene); 
    
               my $count_query_gene=keys(%query_seq);
               chomp($count_query_gene);

               my $count_total_gene= $count_model_gene + $count_query_gene;
               chomp($count_total_gene); 

               my $weight_seq=weight_sequence(\%model_seq,\%query_seq,$count_total_gene);

               open(NEW_GENE,">$out_dir{tmp_log}/$group_id/new_gene.fasta");

               foreach my $new_gene_id(keys %{$weight_seq}){                
                       $new_gene_id=~/(\w+)(\|)(.+)/;   
                       my $genome_name=$1;
                       my $gene_seq=$sequence{$genome_name}->{$new_gene_id};           
                       $gene_seq=~s/\*$//g; 
                       print NEW_GENE ">$new_gene_id\n$gene_seq\n";
               }                                                   
               close NEW_GENE;   

               if($hmm_return_code!=0){
                    open(MODEL_SEQ,">$out_dir{tmp_log}/$group_id/model_gene.fasta"); 

                    foreach my $model_seq_id(keys %model_seq){
                      print MODEL_SEQ ">$model_seq_id\n$model_seq{$model_seq_id}\n";  
                    }
                   close MODEL_SEQ;
               }                    
        
               my $count_new_gene=keys %{$weight_seq};
               chomp($count_new_gene); 
  
              if($count_new_gene>=1 and $hmm_return_code==0){  
     
                  system("hmmalign -o $out_dir{tmp_log}/$group_id/$group_id.stockholm --outformat Stockholm $hmm_file $out_dir{tmp_log}/$group_id/new_gene.fasta")==0 or die "$?";
                  system("hmmbuild $out_dir{hmm_file_dir}/$group_id.hmm $out_dir{tmp_log}/$group_id/$group_id.stockholm > $out_dir{tmp_log}/hmmbuild.log")==0 or die "$?";
            
              }elsif(($count_new_gene>=1 and $count_model_gene>=1) and $hmm_return_code!=0){

                  system("cat $out_dir{tmp_log}/$group_id/model_gene.fasta >> $out_dir{tmp_log}/$group_id/new_gene.fasta")==0 or die "$?";
                  system("kalign -in $out_dir{tmp_log}/$group_id/new_gene.fasta -output $out_dir{tmp_log}/$group_id/$group_id.afa -format fasta -quiet")==0 or die "$?";
                  system("hmmbuild $out_dir{hmm_file_dir}/$group_id.hmm $out_dir{tmp_log}/$group_id/$group_id.afa > $out_dir{tmp_log}/hmmbuild.log")==0 or die "$?";       
                  unlink("$out_dir{singleton_group}/$group_id.fasta"); 
        
              }elsif(($count_new_gene>1 and $count_model_gene==0) and $hmm_return_code!=0){

                  system("cat $out_dir{tmp_log}/$group_id/model_gene.fasta >> $out_dir{tmp_log}/$group_id/new_gene.fasta")==0 or die "$?";
                  system("kalign -in $out_dir{tmp_log}/$group_id/new_gene.fasta -output $out_dir{tmp_log}/$group_id/$group_id.afa -format fasta -quiet")==0 or die "$?";
                  system("hmmbuild $out_dir{hmm_file_dir}/$group_id.hmm $out_dir{tmp_log}/$group_id/$group_id.afa > $out_dir{tmp_log}/hmmbuild.log")==0 or die "$?";       
                  unlink("$out_dir{singleton_group}/$group_id.fasta"); 

              }elsif((($count_new_gene==1 and $count_model_gene==0) or ($count_new_gene==0 and $count_model_gene==1)) and ($hmm_return_code!=0)){

                  system("cat $out_dir{tmp_log}/$group_id/model_gene.fasta >> $out_dir{tmp_log}/$group_id/new_gene.fasta")==0 or die "$?";
                  system("cp $out_dir{tmp_log}/$group_id/new_gene.fasta $out_dir{singleton_group}/$group_id.fasta")==0 or die "$?";
                  system("perl -pi -e 's/^>/>$group_id\t/g' $out_dir{singleton_group}/$group_id.fasta");
              }
    remove_tree("$out_dir{tmp_log}/$group_id"); 
    return(\%processed_genome);
}


##### Weight Sequence ######

sub weight_sequence {

   my(%model_seq)=%{(shift)};
   my(%query_seq)=%{(shift)};
   my($total_sequence)=(shift);


   my $sum_distance=0;
   my %seq_weight=();
   my $flag=0;

   foreach my $query_id(keys %query_seq){

        my $q_seq=$query_seq{$query_id};
           $q_seq=~s/\*//g;
        my $qlen=length($q_seq);
          
        my $flag=0;

        foreach my $model_seq_id(keys %model_seq){

            my $m_seq=$model_seq{$model_seq_id};
               $m_seq=~s/\*//g;
            my $m_seq_len=length($m_seq);
 
            if($q_seq eq $m_seq){
               $flag=1;
            }          
        }

        if($flag==0){
           $model_seq{$query_id}=$q_seq;
           $seq_weight{$query_id}=1;
        }
   }
   
 return(\%seq_weight);
}

