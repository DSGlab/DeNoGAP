##### Module to Map Partial Sequence#####
##### Author: Shalabh Thakur ################
##### Date: 23-AUG-2013 #####################

#!/usr/bin/perl
package PartialMap;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use List::MoreUtils qw(uniq);
use Cluster;
use MclClustering;
use Hash::Merge qw(merge);
use SQLiteDB;

use vars qw(@ISA @EXPORT @EXPORT_OK $db_name $db_dir);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(mapPartialSequence);


sub mapPartialSequence {

    my($group_file)=(shift);
    my(%out_dir)=%{(shift)};
    $db_dir=(shift);
    $db_name=(shift);

    ######## Read genes for each homologue group #########################
    my($homolog_group,$query_status)=Cluster::read_cluster($group_file,'xxx'); 

    ######## Read Group id for each gene #################################
    my($groupforgene)=Cluster::readGroupForGene($group_file); 

    ######## Parse all similarity file to find best partial and full length match #####
    my($superfamily_abc_file)=parseAllSimilarity($db_dir,$db_name,$out_dir{mcl_dir},$groupforgene,$homolog_group); 

    my $single_link_output=cluster_by_single_link($db_dir,$db_name,$out_dir{mcl_dir},$superfamily_abc_file);

    open(MCL_OUT,"$single_link_output");
    my @mcl_groups=<MCL_OUT>;

    ####### Count refine groups for each cluster line ###################
    print "Identify protein families related by fragmented sequences\n";
    my $homolog_cluster_file=mapHomologGroup(\@mcl_groups,$homolog_group,$groupforgene,\%out_dir);

    return($homolog_cluster_file);
}

##### Parse All Similarity File ######
sub parseAllSimilarity {

    my($db_dir)=(shift);
    my($db_name)=(shift);
    my($mcl_out_dir)=(shift);
    my(%groupforgene)=%{(shift)};
    my(%homolog_group)=%{(shift)};
    
    my $superfamily_abc_file="$Bin/$mcl_out_dir/superfamily_abc_file.abc";
    
    open(ABC_FILE,">$superfamily_abc_file");
    
    my %linked_groups=();
    
    ### IF Query Sequence is Potential Truncated #####
    my %truncated_query=();

    my %query_check=();
     
    my($stmt_truncated_query)="Select * from Similarity where query_length<=subject_length and pair_relation='TRUNCATED' ORDER BY query_id, bitscore, evalue";
    my($get_truncated_query)=SQLiteDB::get_record($db_dir,$db_name,$stmt_truncated_query); 
    
    if(scalar(@{$get_truncated_query})>=1){  
             
       foreach my $row(@{$get_truncated_query}){ 
                 
          my @column=@{$row};
          
          my $query_id=$column[0];
          my $target_id=$column[1];
          my $query_length=$column[2];
          my $target_length=$column[3];
          my $evalue=$column[10];
          my $genomeA='';
          my $genomeB='';
          
          $query_id=~/(\w+)(\|)(.+)/;
          $genomeA=$1;

          my $query_hmm_group=$groupforgene{$query_id};
          my $target_hmm_group=$groupforgene{$target_id};  
          
          if(!defined($target_hmm_group)){
             $target_hmm_group=$target_id;
          }else{
            $target_id=~/(\w+)(\|)(.+)/;
            $genomeB=$1;
          }
          
          if($query_hmm_group eq $target_hmm_group){
              next;
          }

          if($query_check{$query_id}){
            next;
          }elsif(!defined($query_check{$query_id}) and ($genomeA ne $genomeB)){

             $query_check{$query_id}=$query_id;
          }
          
          if(!defined($truncated_query{$query_hmm_group}->{$target_hmm_group})){  
          
              $truncated_query{$query_hmm_group}->{$target_hmm_group}=$evalue; 
              $truncated_query{$target_hmm_group}->{$query_hmm_group}=$evalue;  
              
              print ABC_FILE "$query_hmm_group\t$target_hmm_group\t$evalue\n";    
                                     
              $linked_groups{$query_hmm_group}=$query_hmm_group;
              $linked_groups{$target_hmm_group}=$target_hmm_group;              
          }                          
       }
    }

    ##### IF TARGET Sequence is Potential Truncated ######

    my %truncated_target=();
    
    my($stmt_truncated_target)="Select * from Similarity where query_length > subject_length and pair_relation='TRUNCATED' ORDER BY subject_id, bitscore, evalue";
    my($get_truncated_target)=SQLiteDB::get_record($db_dir,$db_name,$stmt_truncated_target); 
    
    if(scalar(@{$get_truncated_target})>=1){  
             
       foreach my $row(@{$get_truncated_target}){ 
                 
          my @column=@{$row};
          
          my $query_id=$column[0];
          my $target_id=$column[1];
          my $query_length=$column[2];
          my $target_length=$column[3];
          my $evalue=$column[10];
          
          my $query_hmm_group=$groupforgene{$query_id};
          my $target_hmm_group=$groupforgene{$target_id};  
          
          if(!defined($target_hmm_group)){
             $target_hmm_group=$target_id;
          }
          
          if($query_hmm_group eq $target_hmm_group){
             next;
          }
          
          if(!defined($truncated_query{$query_hmm_group}->{$target_hmm_group}) and !defined($truncated_query{$target_hmm_group}->{$query_hmm_group})){ 

            if(!defined($truncated_target{$target_hmm_group})){ 

                $truncated_query{$target_hmm_group}->{$query_hmm_group}=$evalue;  
                $truncated_query{$query_hmm_group}->{$target_hmm_group}=$evalue;  

                $truncated_target{$target_hmm_group}=$target_hmm_group;
            
                print ABC_FILE "$target_hmm_group\t$query_hmm_group\t$evalue\n"; 
            
                $linked_groups{$query_hmm_group}=$query_hmm_group;
                $linked_groups{$target_hmm_group}=$target_hmm_group; 
            }      
          }                      
       }
    }
    
    ##### IF Groups are well clustered with no truncated match in other group ########
    
    my %best_group=();

    my($stmt_best_group)="Select * from Similarity where pair_relation='BEST' or pair_relation='SELF_MATCH'";
    my($get_best_group)=SQLiteDB::get_record($db_dir,$db_name,$stmt_best_group);

    if(scalar(@{$get_best_group})>=1){  
             
       foreach my $row(@{$get_best_group}){ 
                 
           my @column=@{$row};
          
           my $query_id=$column[0];
           my $target_id=$column[1];
           my $query_length=$column[2];
           my $target_length=$column[3];
           my $evalue=$column[10];
           
           my $query_hmm_group=$groupforgene{$query_id};
           my $target_hmm_group=$groupforgene{$target_id};  
          
           if(!defined($target_hmm_group)){
             $target_hmm_group=$target_id;
           }
           
          if(!defined($truncated_query{$query_hmm_group}->{$target_hmm_group}) and !defined($truncated_query{$target_hmm_group}->{$query_hmm_group})){ 
        
             $truncated_query{$query_hmm_group}->{$target_hmm_group}=$evalue; 
             $truncated_query{$target_hmm_group}->{$query_hmm_group}=$evalue;
  
             print ABC_FILE "$query_hmm_group\t$target_hmm_group\t$evalue\n";    
            
             $linked_groups{$query_hmm_group}=$query_hmm_group;
             $linked_groups{$target_hmm_group}=$target_hmm_group; 
          }                      
       }
    }

    foreach my $hmm_group(keys %homolog_group){

         if(!defined($linked_groups{$hmm_group})){
            print ABC_FILE "$hmm_group\t$hmm_group\t0\n"; 
         }
    }
    close ABC_FILE;
    return($superfamily_abc_file);  
}

######## CLUSTER BY SINGLE LINK #######

sub cluster_by_single_link{

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($mcl_dir)=(shift);
   my($superfamily_abc_file)=(shift);

   my $single_link_output="$mcl_dir/single_link_output.txt";

   open(SINGLE_LINK,">$single_link_output");

   my $load_stmt="INSERT INTO LinkFamily(family_idA,family_idB,significance) VALUES(?,?,?)";
   SQLiteDB::load_data($db_dir,$db_name,$superfamily_abc_file,$load_stmt);

   my $all_groups_stmt="Select distinct(family_idA) from LinkFamily UNION Select distinct(family_idB) from LinkFamily";
   my($get_all_groups)=SQLiteDB::get_record($db_dir,$db_name,$all_groups_stmt);
     
      my %list_groups=();

      foreach my $row(@{$get_all_groups}){                 
          my @column=@{$row};
          $list_groups{$column[0]}=$column[0];
      }

   NEW_CLUSTER:

   print "NEW CLUSTER:\t";
   
   foreach my $group(keys %list_groups){

         my %search_groups=();

         print "$group\n";

         $search_groups{$group}=$group;

         RESEARCH_PAIR:

         my $search_groups='';

         foreach my $group(keys %search_groups){ 
             $search_groups=$search_groups.",'".$group."'";
         }
         
         $search_groups=~s/^\,//g;
         
         my $pair_stmt="Select distinct(family_idA) from LinkFamily where family_idB IN ($search_groups) UNION Select distinct(family_idB) from LinkFamily where family_idA IN ($search_groups)";
         my($get_all_pairs)=SQLiteDB::get_record($db_dir,$db_name,$pair_stmt);

         if(scalar(@{$get_all_pairs})>=1){

            foreach my $row(@{$get_all_pairs}){                 
              my @column=@{$row};
              $search_groups{$column[0]}=$column[0]; 
            }

            my $delete_pair_record="DELETE FROM LinkFamily where family_idA IN ($search_groups) or family_idB IN ($search_groups)";
            SQLiteDB::execute_sql($db_dir,$db_name,$delete_pair_record);

            goto RESEARCH_PAIR;
            
         }else{

             foreach my $group(keys %search_groups){ 
             
                print SINGLE_LINK "$group\t";
                
                if($list_groups{$group}){
                  delete($list_groups{$group});
                }  
             }
             print SINGLE_LINK "\n"; 
             
             goto NEW_CLUSTER; 
         }
   }
   close SINGLE_LINK;

  return($single_link_output);
}

sub mapHomologGroup {

    my(@mcl_group)=@{(shift)};
    my(%homolog_group)=%{(shift)};
    my(%groupforgene)=%{(shift)};
    my(%out_dir)=%{(shift)};
   
    my %map_hmm_to_super_family=();   
    my %super_family=();
 
    ###### Check information that is stored in GenetoSuperFamily table ######
    print "Check gene families in database\n";

    my $sql_get_family="SELECT DISTINCT(super_family_id),hmm_family_id from GenetoSuperFamily";
    
    my($row_family)=SQLiteDB::get_record($db_dir,$db_name,$sql_get_family);
  
    if(scalar(@{$row_family})>0){
    
       foreach my $row(@{$row_family}){   
       
          my $super_cluster_id=shift(@{$row}); 
          my $hmm_family_id=shift(@{$row});         
                  
          $map_hmm_to_super_family{$hmm_family_id}=$super_cluster_id;  
          $super_family{$super_cluster_id}=$super_cluster_id;                                               
       }        
    }
      
    my %homolog_family=();

    foreach my $homolog_group_line(@mcl_group){
     
         chomp($homolog_group_line);
         
         my @hmm_list=split("\t",$homolog_group_line);
         my %list_super_family=();
         my $homolog_gene='';
         
         #### get all the gene_ids for hmm group clustered together ###
         foreach my $hmm_group(@hmm_list){
         
            my @gene_list=$homolog_group{$hmm_group};
            
            my $super_cluster_id=$map_hmm_to_super_family{$hmm_group};
                                    
            my $gene_ids=join(" ",@gene_list);  
               
            $homolog_gene=$homolog_gene." ".$gene_ids;
                        
            if($super_cluster_id){
            
              $list_super_family{$super_cluster_id}=$super_cluster_id;        
              
            } 

         }
         
         ### count number of super family on each cluster line ###
         my $num_superfamily=keys %list_super_family;
         
         ### if only one super family is defined for the hmm cluster line ###
         if($num_superfamily==1){
         
             foreach my $superfamily(keys %list_super_family){        
             
                $homolog_family{$superfamily}=$homolog_gene;
                
             }
         }elsif($num_superfamily > 1){ 
            
             ### if 2 or more super family are clustered together, assign new super family id for all hmm groups ####   
             my $count_family=(keys %super_family) + 1;
             my $family_id="Cluster_".$count_family;
             $homolog_family{$family_id}=$homolog_gene;
             $super_family{$family_id}=$family_id;
             
             #### delete records for old super family ids from alignment table and homolog pair tables ####
            # foreach my $superfamily(keys %list_super_family){ 
              
            #    print "Deleteing alignment for $family_id\n";
             
               # my $delete_aln_record="DELETE FROM Alignment where alignment_id='$superfamily'";
               # SQLiteDB::execute_sql($db_dir,$db_name,$delete_aln_record);
             
            #    print "Deleting homolog pair for $family_id\n";
             
              #  my $delete_ortholog_record="DELETE FROM Ortholog where homolog_cluster_id='$superfamily'";
              #  SQLiteDB::execute_sql($db_dir,$db_name,$delete_ortholog_record);
             
              #  my $delete_inparalog_record="DELETE FROM Inparalog where homolog_cluster_id='$superfamily'";
              #  SQLiteDB::execute_sql($db_dir,$db_name,$delete_inparalog_record);
             
              #  my $delete_incongruent_record="DELETE FROM Incongruent where homolog_cluster_id='$superfamily'";
              #  SQLiteDB::execute_sql($db_dir,$db_name,$delete_incongruent_record);
             
              #  my $delete_distance_record="DELETE FROM DistancePair where homolog_cluster_id='$superfamily'";
              #  SQLiteDB::execute_sql($db_dir,$db_name,$delete_distance_record);
             
            # }
             
         }else{
             ##### if new super family cluster is there, assign new family id for hmm groups ####
             my $count_family=(keys %super_family) + 1;
             
             my $family_id="Cluster_".$count_family;
             
             $homolog_family{$family_id}=$homolog_gene;
             
             $super_family{$family_id}=$family_id;
         }    
    }
  
  
    ##### LOAD Super family information in database #####
    print "Deleting old records for Gene2Superfamily Table\n";
    
    my $delete_old_record="DELETE FROM GenetoSuperFamily";
    SQLiteDB::execute_sql($db_dir,$db_name,$delete_old_record);


    #### SUPER HOMOLOG FAMILY FILE #####    
    my $homolog_cluster_file="$out_dir{result_dir}/homolog_cluster.txt";
    open(MAP_GENE,">$homolog_cluster_file"); 

 
    ##### MAPING TABLE FOR HMM FAMILY TO SUPER HOMOLOG FAMILY ####
    my $GenetoSuperfamily_data="$out_dir{result_dir}/GenetoSuperfamily.txt";
    open(G2S,">$GenetoSuperfamily_data");
    
    my @new_records_to_load=();
    
    foreach my $super_family(keys %homolog_family){
    
        my $seq_id_line=$homolog_family{$super_family};
        
        print MAP_GENE "$super_family:\t$seq_id_line\n";
        
        my @array_seq_id=split(" ",$seq_id_line);
        
        foreach my $seq_id(@array_seq_id){
        
             my $hmm_family=$groupforgene{$seq_id};
             
             $seq_id=~/(\w+)(\|)(.+)/;
             
             my $id=$3;
             
             my $genome=$1;
             
             my $record="$id\t$genome\t$hmm_family\t$super_family";
                          
             print G2S "$record\n";
        }
    }
    
    close MAP_GENE;
    close G2S;
    
    print "Load  New Homolog Gene Family Information\n";
    
    my $load_stmt="INSERT INTO GenetoSuperFamily(gene_id,genome_name,hmm_family_id,super_family_id) VALUES(?,?,?,?)";
    
    SQLiteDB::load_data($db_dir,$db_name,$GenetoSuperfamily_data,$load_stmt);
  
    return($homolog_cluster_file); 
}  

