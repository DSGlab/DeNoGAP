##### Module to Initialize the global variable for denovo protein family prediction#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

##### Version 2.0: Do not need Distance File to be inputed ###

##### Version 2.1: Do not need alignment file and hmm model file and singleton sequence fasta file to be inputed ####

#!/usr/bin/perl
package Initialize;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Cluster;
use File::Basename;
use Bio::AlignIO;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(check_input);


sub check_input {

    my %config_param=%{(shift)};
    my %out_dir=%{(shift)};
    my $cluster_file='';
    my %group_gene=();
    my %homologue_group=();
    my $cluster_num=0;

    ######### CLUSTER FILE #######
    if($config_param{GROUP}->{HMM_CLUSTER_FILE} and -s "$config_param{GROUP}->{HMM_CLUSTER_FILE}"){
    
       $cluster_file=$config_param{GROUP}->{HMM_CLUSTER_FILE};
       
       my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,'xxx'); 
       
       $cluster_num=keys %{$homologue_group};
       
       %homologue_group=%{$homologue_group};
       
       ###### Initialize group_gene hash table #############
       foreach my $group_id(keys %homologue_group){    
             
            my $gene_list=$homologue_group{$group_id};
            
            my @gene_list=split(" ",$gene_list);  
                      
            foreach(@gene_list){   
                $_=~/(\w+)(\|)(.+)/;
                $group_gene{$1}=$1;               
            } 
       }  
          
    }else{
    
      print STDERR "CLUSTER FILE [$config_param{GROUP}->{HMM_CLUSTER_FILE}] cannot be found, please check again\n"; exit;
      
    } 

   return($cluster_file,\%group_gene,\%homologue_group);
}
