##### Module to read protein family cluster#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package Cluster;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";


use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(read_cluster readGroupForGene getGroupLength);

####### Genes for each groups ######

sub read_cluster {
    my($cluster_file)=(shift);
    my($genome_name)=(shift);
    my(%homologue_group)=();
    my $genome_status=0;

    open(GROUP_OUT,"$cluster_file") or die "$0: open file: $!";

    foreach(<GROUP_OUT>){

          chomp($_);
          if($_=~/$genome_name/){
             $genome_status=1;
             last;
          }      
          my @group=split("\t",$_);
          my $group_id=shift(@group);
          my $group_line=join("\t",@group);
          $group_id=~s/\://g;
          $group_line=~s/^\t//g;          
          $group_line=~s/\s+$//g;   
          chomp($group_line);
          $homologue_group{$group_id}=$group_line;            
    }

    return(\%homologue_group,$genome_status);
}

##### Group for each gene ####

sub readGroupForGene {

    my($cluster_file)=(shift);
    my %groupforgene=();
 
    open(GROUP_OUT,"$cluster_file") or die "$0: open file: $!";

    foreach(<GROUP_OUT>){
          chomp($_);
          my @group=split("\t",$_);
          my $group_id=shift(@group);
          my $group_line=join("\t",@group);
          $group_id=~s/\://g;
          $group_line=~s/^\t//g;
          $group_line=~s/\s+$//g;
          chomp($group_line);
          my @genes=split("\t",$group_line);

          foreach my $gene(@genes){
             chomp($gene);
             $groupforgene{$gene}=$group_id;                      
          }                            
    }

   return(\%groupforgene);
}

#### Group Length #####
sub getGroupLength {

    my(%homolog_group)=%{(shift)};
    my(%seq_len_table)=%{(shift)};
    my(%grouplength)=();  
     
   
    foreach my $group_id(keys %homolog_group){

          my @gene_list=split("\t",$homolog_group{$group_id});
         
          foreach(@gene_list){  
               chomp($_);
               $_=~s/\s+$//g;            
               $_=~/(\w+)(\|)(.+)/;
               my $genome_name=$1;
               my $len=$seq_len_table{$genome_name}->{$_};            
               
               if((!defined($grouplength{$group_id})) or (($len) and ($len>=$grouplength{$group_id}))){
                 $grouplength{$group_id}=$len;                               
               }
          }
                
    }
   return(\%grouplength);
}
