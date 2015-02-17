##### Module to Predict Annotation such as Domain using InterProScan #######
##### Author: Shalabh Thakur #################
##### Date: 20-JAN-2014 ######################

#!/usr/bin/perl

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use Hash::Merge qw( merge );

my($genome_name)=(shift);
my($ipr_result_xml_file)=(shift); 
my($ipr_result_tsv_file)=(shift); 
my($ipr_domain_dir)=(shift);
my($ipr_interpro_dir)=(shift);
my($ipr_go_dir)=(shift);
my($ipr_pathway_dir)=(shift);
my($ipr_signalp_dir)=(shift);
my($ipr_tmhmm_dir)=(shift);
my($ipr_phobius_dir)=(shift);

my($xml_feature)=parse_iprscan_xml_output($ipr_result_xml_file);

parse_iprscan_result($genome_name,$ipr_result_tsv_file,$xml_feature,$ipr_domain_dir,$ipr_interpro_dir,$ipr_go_dir,$ipr_pathway_dir,$ipr_signalp_dir,$ipr_tmhmm_dir,$ipr_phobius_dir);


sub parse_iprscan_xml_output{

  my($ipr_result_xml_file)=(shift);

  open(XML_FILE,"$ipr_result_xml_file");
  my @xml_file=<XML_FILE>;
  close XML_FILE;

  my %xml_feature=();

  foreach my $xml_line(@xml_file){

     if($xml_line=~/\<signature\sname/){
        
           $xml_line=~/(name)(=)(\")(.+?)(\")/;
           my $name=$4;
   
           $xml_line=~/(desc)(=)(\")(.+?)(\")/;
           my $desc=$4;

           $xml_line=~/(ac)(=)(\")(.+?)(\")/; 
           my $accession_id=$4;

           $xml_feature{$accession_id}->{name}=$name;
           $xml_feature{$accession_id}->{desc}=$desc;
     
      }elsif($xml_line=~/(\<go\-xref)/){

           $xml_line=~/(category)(=)(\")(.+?)(\")/; 
           my $category=$4;

           $xml_line=~/(name)(=)(\")(.+?)(\")/; 
           my $name=$4;

           $xml_line=~/(id)(=)(\")(.+?)(\")/; 
           my $go_id=$4;

           $xml_feature{$go_id}->{$category}=$name;

      }elsif($xml_line=~/(\<pathway\-xref)/){

           $xml_line=~/(name)(=)(\")(.+?)(\")/; 
           my $name=$4;

           $xml_line=~/(id)(=)(\")(.+?)(\")/; 
           my $id=$4;

           $xml_line=~/(db)(=)(\")(.+?)(\")/; 
           my $db=$4;

           my $pathway_id=$db.": ".$id;

           $xml_feature{$pathway_id}=$name;
      }    
   }

  return(\%xml_feature);
}

####### Parse Iprscan TSV Result ######

sub parse_iprscan_result{

    my($genome)=(shift);
    my($iprscan_out)=(shift);
    my(%xml_feature)=%{(shift)};
    my($ipr_domain_dir)=(shift);
    my($ipr_interpro_dir)=(shift);
    my($ipr_go_dir)=(shift);
    my($ipr_pathway_dir)=(shift);
    my($ipr_signalp_dir)=(shift);
    my($ipr_tmhmm_dir)=(shift);
    my($ipr_phobius_dir)=(shift);
   
    my %iprscan_annotation=();
    
    ###### Open raw output file ####
    open(IPRSCAN_OUT,"$iprscan_out");
    my @iprscan_output=<IPRSCAN_OUT>;
    close IPRSCAN_OUT;

    open(D,">$ipr_domain_dir/DOMAIN_$genome.txt");
    open(I,">$ipr_interpro_dir/INTERPRO_$genome.txt");
    open(G,">$ipr_go_dir/GO_$genome.txt");
    open(P,">$ipr_pathway_dir/PATHWAY_$genome.txt");
    open(S,">$ipr_signalp_dir/SIGNALP_$genome.txt");
    open(T,">$ipr_tmhmm_dir/TMHMM_$genome.txt");
    open(PO,">$ipr_phobius_dir/PHOBIUS_$genome.txt");
   
    foreach my $row_line(@iprscan_output) {
        my ($interpro_annotation)=parse_iprscan_line($row_line);           
        %iprscan_annotation=%{merge(\%iprscan_annotation,$interpro_annotation)};               
    } 

    foreach my $locus_tag(keys %iprscan_annotation){

        $locus_tag=~/(\w+)(\|)(.+)/;
        my $locus_id=$3;
        my $genome_name=$1;

        my %signature_type=%{$iprscan_annotation{$locus_tag}};

       foreach my $signature_db(keys %signature_type){

         my %pfam=%{$iprscan_annotation{$locus_tag}->{$signature_db}};

         if($signature_db eq "INTERPRO" or $signature_db eq "GO" or $signature_db eq "PATHWAY"){
           next;
         }

         foreach my $domain_id(keys %pfam){
 
             my $domain_id=$pfam{$domain_id}->{id};
             my $db_source=$signature_db;
             my $domain_name=$xml_feature{$domain_id}->{name};
             my $seq_length=$pfam{$domain_id}->{locus_length};
             my $domain_desc=$xml_feature{$domain_id}->{desc}; 
             my $domain_start=$pfam{$domain_id}->{start};
             my $domain_end=$pfam{$domain_id}->{end};
             my $significance=$pfam{$domain_id}->{significance};

             print D "$locus_id\t$genome_name\t$seq_length\t$db_source\t$domain_id\t$domain_name\t$domain_start\t$domain_end\t$significance\t$domain_desc\n";
         }
       }

       if($iprscan_annotation{$locus_tag}->{INTERPRO}){

          my %interpro=%{$iprscan_annotation{$locus_tag}->{INTERPRO}};

          foreach my $ipr_domain(keys %interpro){
             my $ipr_id=$interpro{$ipr_domain}->{id};
             my $ipr_desc=$interpro{$ipr_domain}->{name};

             print I "$locus_id\t$genome_name\t$ipr_id\t$ipr_desc\n";
          }
       }

       if($iprscan_annotation{$locus_tag}->{GO}){

          my %go=%{$iprscan_annotation{$locus_tag}->{GO}};

          foreach my $go_anno(keys %go){
             my $go_id=$go{$go_anno}->{id};
             my %go_category=%{$xml_feature{$go_id}};

             foreach my $category(keys %go_category){
               my $go_desc=$go_category{$category};

               print G "$locus_id\t$genome_name\t$go_id\t$category\t$go_desc\n";
             }
          }
        }

        if($iprscan_annotation{$locus_tag}->{PATHWAY}){
          my %pathway=%{$iprscan_annotation{$locus_tag}->{PATHWAY}};
          foreach my $pathway_anno(keys %pathway){
             my $pathway_id=$pathway{$pathway_anno}->{id};
             my $pathway_name=$xml_feature{$pathway_id};

             print P "$locus_id\t$genome_name\t$pathway_id\t$pathway_name\n";
          }
        }
    }

    close D;
    close P;
    close I;
    close G;
    close S;
    close T;
    close PO;
}

#### Parse HMMPFAM ####

sub parse_iprscan_line{

   my $row_line=(shift);

   chomp($row_line);

   my $locus_tag='NULL';
   my $locus_length='NULL';
   my $domain_type='NULL';
   my $domain_id='NULL';
   my $domain_name='NULL';
   my $domain_start='NULL';
   my $domain_end='NULL';
   my $domain_evalue='NULL';
   my $interpro_id='NULL';
   my $interpro_name='NULL';
   my $go_annotation="NULL";
   my $kegg_annotation="NULL";
   my $dbxref_annotation='NULL';

    my @column=split("\t",$row_line);

       $locus_tag=$column[0];       
       $locus_length=$column[2];
       $domain_type=$column[3];  
       $domain_id=$column[4];
       $domain_name=$column[5];
       $domain_start=$column[6];
       $domain_end=$column[7];
       $domain_evalue=$column[8];
       $interpro_id=$column[11];
       $interpro_name=$column[12];
       $go_annotation=$column[13];
       $kegg_annotation=$column[14];
      
       if($go_annotation and $go_annotation ne "NULL"){
          $dbxref_annotation=$go_annotation;
       }

       if($dbxref_annotation eq "NULL" and $kegg_annotation and $kegg_annotation ne "NULL"){
          $dbxref_annotation=$kegg_annotation; 
       }elsif($dbxref_annotation ne "NULL" and $kegg_annotation and $kegg_annotation ne "NULL"){
          $dbxref_annotation=$dbxref_annotation."\t".$kegg_annotation; 
       }

    my %iprscan_annotation=();

    my %signalp_annotation=();

    my %phobius_annotation=();

    my %tmhmm_annotation=();
	
    ######## GO and PATHWAY annotation ######	

    if(($dbxref_annotation) and ($dbxref_annotation ne 'NULL' or $dbxref_annotation ne '')){  
      
         my @dbxref_annotation=split("\t",$dbxref_annotation);
	  
         foreach my $xref_annotation(@dbxref_annotation) { 
             
                 if($xref_annotation=~/\|/){
                     my @annotation=split(/\|/,$xref_annotation);
                  
                     foreach my $xref_id(@annotation){ 
                          if($xref_id=~/^GO\:/){
                             $iprscan_annotation{$locus_tag}->{GO}->{$xref_id}->{id}=$xref_id; 
                          }else{
                             $iprscan_annotation{$locus_tag}->{PATHWAY}->{$xref_id}->{id}=$xref_id; 
                          }
                     }
                }else{
                    if($xref_annotation=~/^GO\:/){
                         $iprscan_annotation{$locus_tag}->{GO}->{$xref_annotation}->{id}=$xref_annotation; 
                        
                    }elsif($xref_annotation ne "NULL"){
                         $iprscan_annotation{$locus_tag}->{PATHWAY}->{$xref_annotation}->{id}=$xref_annotation; 
                    }
                }   
         }                       	   
    } 

     if($domain_type=~/SignalP/){

         my @row_column=split("\t",$row_line);

         $row_column[0]=~/(\w+)(\|)(.+)/;
         my $protein_id=$3;
         my $genome_name=$1;
         my $start=$row_column[5];
         my $end=$row_column[6];

         print S "$protein_id\t$genome_name\t$start\t$end\n";

     }elsif($domain_type eq "Phobius"){

         my @row_column=split("\t",$row_line);

         $row_column[0]=~/(\w+)(\|)(.+)/;
         my $protein_id=$3;
         my $genome_name=$1;
         my $domain_name=$row_column[4];
         my $domain_description=$row_column[5];
         my $start=$row_column[6];
         my $end=$row_column[7];

         print PO "$protein_id\t$genome_name\t$domain_name\t$domain_description\t$start\t$end\n";

    }elsif($domain_type eq "TMHMM"){

         my @row_column=split("\t",$row_line);

         $row_column[0]=~/(\w+)(\|)(.+)/;
         my $protein_id=$3;
         my $genome_name=$1;
         my $start=$row_column[5];
         my $end=$row_column[6];

         print T "$protein_id\t$genome_name\t$start\t$end\n";

    }else{   
       ####### InterPro Domain Information #####
       if(($interpro_id) and ($interpro_id ne 'NULL' or $interpro_id ne '')){  
          $iprscan_annotation{$locus_tag}->{INTERPRO}->{$interpro_id}->{id}=$interpro_id; 
          $iprscan_annotation{$locus_tag}->{INTERPRO}->{$interpro_id}->{name}=$interpro_name;
       }

       ###### DOMAIN SIGNATURES #####
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{id}=$domain_id;
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{name}=$domain_name;
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{locus_length}=$locus_length;
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{start}=$domain_start;
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{end}=$domain_end;
       $iprscan_annotation{$locus_tag}->{$domain_type}->{$domain_id}->{significance}=$domain_evalue;
    }
   
    return(\%iprscan_annotation);   
}

