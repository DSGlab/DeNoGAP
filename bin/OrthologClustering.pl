##### Module to Cluster Ortholog Pairs #######
##### Author: Shalabh Thakur #################
##### Date: 28-SEP-2013 ######################

#!/usr/bin/perl
package OrthologCluster;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Cluster;
use SQLiteDB;
use Hash::Merge qw(merge );
use File::Path qw(remove_tree);
use Parallel::ForkManager;


my($homolog_cluster_file)=(shift);
my($cluster_dir)=(shift);
my($db_dir)=(shift);
my($db_name)=(shift);
my($tmp_dir)=(shift); 
my($process)=(shift);
my($transitivity_threshold)=(shift);


my($homolog_cluster,$query_status)=Cluster::read_cluster($homolog_cluster_file,'xxx');

my %homolog_cluster=%{($homolog_cluster)};

my $fork=Parallel::ForkManager->new($process);
    
foreach my $cluster_id(keys %homolog_cluster){

        $fork->start and next;
        
        print "$cluster_id\n";

        my $cluster_gene=$homolog_cluster{$cluster_id};

        my @cluster_gene=split(' ',$cluster_gene);

        my $size_cluster=scalar(@cluster_gene);

        if($size_cluster eq 1){
        
          open(SINGLE_LINK,">$cluster_dir/orthoCluster_$cluster_id.txt");
          
              print SINGLE_LINK "$cluster_id.1:\t";
              
              print SINGLE_LINK "$cluster_gene[0]\n";
              
          close SINGLE_LINK;
          
          $fork->finish();
        }

        #### CREATE TEMP DATABASE ######

        chdir($tmp_dir);   
    
        my $dbh=DBI->connect("dbi:SQLite:dbname=$cluster_id","","") or die $DBI::errstr;
        my $stmt=$dbh->prepare("SELECT SQLITE_VERSION()");
        $stmt->execute();   
        $stmt->finish();

        $dbh->do("CREATE TABLE IF NOT EXISTS $cluster_id(taxaA TEXT, idA TEXT , taxaB TEXT, idB TEXT, cluster_id TEXT)");

        chdir($Bin);

        ##### READ ORTHOLOG PAIR ###
        my %list_protein_id=();
        
        #open(ORTHOLOG_FILE,"$ortholog_pair_dir/$cluster_id.txt");
        #open(INPARALOG_FILE,"$inparalog_pair_dir/$cluster_id.txt");

        #my @ortholog_pair=<ORTHOLOG_FILE>;
        #my @inparalog_pair=<INPARALOG_FILE>;

        my @pair=();

        my $get_ortholog_pair="Select distinct(taxonA),idA,taxonB,idB,homolog_cluster_id from Ortholog where homolog_cluster_id='$cluster_id'";

        my $ortho_pair=SQLiteDB::get_record($db_dir,$db_name,$get_ortholog_pair);

        if(scalar(@{$ortho_pair})>0){
          
            foreach my $row(@{$ortho_pair}){
      
               my $pair_line=join("\t",@{$row});

               push(@pair,$pair_line);  

            }
        }

        my $get_inparalog_pair="Select distinct(taxonA),idA,taxonB,idB,homolog_cluster_id from InParalog where homolog_cluster_id='$cluster_id'";

        my $inpara_pair=SQLiteDB::get_record($db_dir,$db_name,$get_inparalog_pair);

        if(scalar(@{$inpara_pair})>0){
          
            foreach my $row(@{$inpara_pair}){
      
               my $pair_line=join("\t",@{$row});

               push(@pair,$pair_line);  
            }
        }
        

        my $sql_stmt="INSERT into $cluster_id (taxaA,idA,taxaB,idB,cluster_id) VALUES (?,?,?,?,?)";

        SQLiteDB::load_from_array($tmp_dir,$cluster_id,\@pair,$sql_stmt);


        my $all_ids_stmt="Select distinct(idA),taxaA from $cluster_id UNION Select distinct(idB),taxaB from $cluster_id";
        my($get_all_ids)=SQLiteDB::get_record($tmp_dir,$cluster_id,$all_ids_stmt);
     
        my %list_ids=();

        foreach my $row(@{$get_all_ids}){                 
          my @column=@{$row};
          $list_ids{$column[0]}=$column[1];
       }

       open(SINGLE_LINK,">$cluster_dir/orthoCluster_$cluster_id.txt");

       my $cluster_count=1;

       my %search_ids=();

       NEW_CLUSTER:
   
       foreach my $id(keys %list_ids){

         %search_ids=();

         my %search_genome=();

         my $start_genome=$list_ids{$id};

        # print "Start:$id\t$start_genome\n";

         $search_ids{$id}=$start_genome;
         $search_genome{$start_genome}=$start_genome;
        
         my $pair_stmt="Select distinct(idB),taxaB from $cluster_id where idA='$id' and taxaA='$start_genome'";
         my($get_all_pairs)=SQLiteDB::get_record($tmp_dir,$cluster_id,$pair_stmt);

         my %reciprocal_pair=();
         my $search_id='';
         my $search_genome='';

         if(scalar(@{$get_all_pairs})>=1){

            foreach my $row(@{$get_all_pairs}){
     
              my @column=@{$row};
              my $pair_id=$column[0];
              my $pair_genome=$column[1];
         
              $reciprocal_pair{$pair_id}=$pair_genome;
              $search_id=$search_id.",'".$pair_id."'";
              $search_genome=$search_genome.",'".$pair_genome."'";
           }
         }
         $search_genome=~s/^\,//g;            
         $search_id=~s/^\,//g;
           
             
         foreach my $reciprocal_id(keys %reciprocal_pair){

                 my $reciprocal_genome=$reciprocal_pair{$reciprocal_id};

                 ###### Find reciprocal gene in species N for all gene (i..n-1) from species (j..jn-1) ####

                 my $best_ortho_stmt="Select distinct(idB),taxaB,taxaA from $cluster_id where idA IN ($search_id) and taxaA IN ($search_genome) and taxaB='$reciprocal_genome'";
                 my($get_best_ortholog)=SQLiteDB::get_record($tmp_dir,$cluster_id,$best_ortho_stmt);

                 my %found_ortholog=();   #### record 1 if gene C from species C is reciprocal best match for all gene a..n else 0 ###

                 if(scalar(@{$get_best_ortholog})>=1){

                      foreach my $row(@{$get_best_ortholog}){
                       
                             my @column=@{$row};

                             my $rec_ortholog=$column[0]; 
                             my $target_genome=$column[2];  

                             #print "$target_genome\t$rec_ortholog\t$reciprocal_id\n";

                             $found_ortholog{$target_genome}->{$rec_ortholog}=1;
                      }
                 }

                 #### Check for true ortholog pair ###
                 my $num_genome=keys(%found_ortholog);
                 my %count_ortholog_genome=();
                   
                 foreach my $genome(keys %found_ortholog){

                      my %ortholog_id=%{$found_ortholog{$genome}};

                      foreach my $ortho_id(keys %ortholog_id){

                           if(!defined($count_ortholog_genome{$ortho_id})){
                              $count_ortholog_genome{$ortho_id}=1;
                           }else{
                              $count_ortholog_genome{$ortho_id}++;
                           }
                      }
                 }
                 
                 ##### If true ortholog than add for printing to the file #####
                 foreach my $ortho_id(keys %count_ortholog_genome){

                      my $count_genome=$count_ortholog_genome{$ortho_id};

                      my $is_true_ortholog=($count_genome/$num_genome)*100;

                    #  print "$reciprocal_genome\t$ortho_id\t$is_true_ortholog\n";

                      if($is_true_ortholog>=$transitivity_threshold){
                        $search_ids{$ortho_id}=$reciprocal_genome;
                        $search_genome{$reciprocal_genome}=$reciprocal_genome;
                      }
                 }  
                       
         }
           
        print SINGLE_LINK "$cluster_id".".$cluster_count:\t";

        my $remove_id_record='';

        foreach my $protein_id(keys %search_ids){ 
             
             my $ortholog_genome=$search_ids{$protein_id};
             print SINGLE_LINK "$ortholog_genome|$protein_id\t";
                
             delete($list_ids{$protein_id});

             $remove_id_record=$remove_id_record.",'".$protein_id."'";   
        }
        $remove_id_record=~s/^\,//g;

        #### delete pair record for all clustered true ortholog ids from database ####
        my $delete_pair_record="DELETE FROM $cluster_id where idA IN ($remove_id_record) or idB IN ($remove_id_record)";
        SQLiteDB::execute_sql($tmp_dir,$cluster_id,$delete_pair_record);

        print SINGLE_LINK "\n"; 
           
        $cluster_count++;
        goto NEW_CLUSTER; 

   }
   close SINGLE_LINK;

system("rm $tmp_dir/$cluster_id");

$fork->finish();

}
$fork->wait_all_children;     
