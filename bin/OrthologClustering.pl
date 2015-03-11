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
       my %reciprocal_pair=();
       my %pair_total_gene=();

       NEW_CLUSTER:
   
       foreach my $id(keys %list_ids){

         %search_ids=();
         %reciprocal_pair=();
         %pair_total_gene=();

         print "Start:$id\n";

         my $search_id='';

         $search_id=$search_id."'".$id."'";

         $pair_total_gene{$id}=$list_ids{$id};

         RESEARCH_PAIR:
        
        my $pair_stmt="Select distinct(idA),taxaA,idB,taxaB,divergence from $cluster_id where idB IN ($search_id) UNION Select distinct(idB),taxaB,idA,taxaA,divergence from $cluster_id where idA IN ($search_id)";
        my($get_all_pairs)=SQLiteDB::get_record($db_dir,$cluster_id,$pair_stmt);

        my $is_new_genome=0;

        if(scalar(@{$get_all_pairs})>=1){          

           foreach my $row(@{$get_all_pairs}){
     
              my @column=@{$row};
              my $pair_idA=$column[0];
              my $pair_genomeA=$column[1];
              my $pair_idB=$column[2];
              my $pair_genomeB=$column[3];
              my $divergence=$column[4];

              $reciprocal_pair{$pair_idB}->{$pair_idA}=$divergence;
              $reciprocal_pair{$pair_idA}->{$pair_idB}=$divergence;

              $pair_total_gene{$pair_idA}=$pair_genomeA;
              $pair_total_gene{$pair_idB}=$pair_genomeB;

             if($search_id=~/$pair_idA/){
                 next;
             }else{
                 $search_id=$search_id.",'".$pair_idA."'";                                              
                 $is_new_genome=1; 
             }        
           }
        }

         if($is_new_genome==1){
             $search_id=~s/^\,//g;
             goto RESEARCH_PAIR;
         }
          
         ####### check weight for each pair ##############     

         my $remove_id_record='';

         print SINGLE_LINK "$cluster_id".".$cluster_count:\t";

         my $total_count=keys(%pair_total_gene);
    
        foreach my $pair_geneA(keys %pair_total_gene){

                my $similar_count=0;

                my %similar_pair=();

                foreach my $pair_geneB(keys %pair_total_gene){

                      unless($pair_geneA eq $pair_geneB){

                          if($reciprocal_pair{$pair_geneA}->{$pair_geneB}){

                              my $distance=$reciprocal_pair{$pair_geneA}->{$pair_geneB};

                              $similar_pair{$pair_geneA}->{$pair_geneB}=$distance;

                              $similar_count++;
                          } 
                      }
                }

                my $percent_transitivity=($similar_count/$total_count)*100;

                if($percent_transitivity>=$transitivity_threshold){

                     my $genome=$pair_total_gene{$pair_geneA};

                     print SINGLE_LINK "$genome|$pair_geneA\t";

                     delete($list_ids{$pair_geneA});

                     $remove_id_record=$remove_id_record.",'".$pair_geneA."'"; 
                }
        }
        
        $remove_id_record=~s/^\,//g;

        if($remove_id_record eq ''){

               my $genome=$pair_total_gene{$id};

               print SINGLE_LINK "$genome|$id\t";

               #print "$genome|$id\n";

               delete($list_ids{$id});

               $remove_id_record=$remove_id_record."'".$id."'"; 

               $remove_id_record=~s/^\,//g; 
        }

        #### delete pair record for all clustered true ortholog ids from database ####
        my $delete_pair_record="DELETE FROM $cluster_id where idA IN ($remove_id_record) or idB IN ($remove_id_record)";
        SQLiteDB::execute_sql($db_dir,$cluster_id,$delete_pair_record);

        print SINGLE_LINK "\n"; 
           
        $cluster_count++;

        goto NEW_CLUSTER; 
   }
   close SINGLE_LINK;

system("rm $tmp_dir/$cluster_id");

$fork->finish();

}
$fork->wait_all_children;     
