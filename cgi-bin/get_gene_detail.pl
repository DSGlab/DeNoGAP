#!/usr/bin/perl -w
###### ABOUT: This Script find gene information for selected gene id ############
###### AUTHOR:Shalabh Thakur###################################################################


use strict;
use warnings;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use Env;
use FindBin qw($Bin);
use DBI;

print "Content-type:text/html\n\n";

my $project_name=param('project_name');
my $db_dir=param('db_dir');
my $db_name=param('db_name');
my $report_dir=param('report_dir');
my $gene_id=param('gene_detail_radio');
my $genome_code=param('genome_name');

my %gene_annotation=();

###### get feature details #####
$gene_annotation{locus_tag}=$gene_id;

### get full name of the species /strain ###
my $species_stmt="SELECT * from OrganismInfo WHERE abbreviation='$genome_code'";
my ($species_data)=get_record($db_dir,$db_name,$species_stmt);

 if(scalar(@{$species_data})>0){
      foreach my $row(@{$species_data}){
          shift(@{$row});          
          $gene_annotation{genome_full_name}=shift(@{$row});
          $gene_annotation{species}=shift(@{$row});
          $gene_annotation{species_abbreviation}=shift(@{$row});
          shift(@{$row});
      }
 }

### get Taxonomy data ####
my $taxonomy_stmt="SELECT * from Taxonomy WHERE species='$gene_annotation{species}'";
my ($taxonomy_data)=get_record($db_dir,$db_name,$taxonomy_stmt);

### get Genomic Feature ###
my $feature_stmt="SELECT * from GeneFeature WHERE feature_id='$gene_id' and genome_name='$gene_annotation{genome_full_name}'";
my ($feature_data)=get_record($db_dir,$db_name,$feature_stmt);

 if(scalar(@{$feature_data})>0){
      foreach my $row(@{$feature_data}){
          shift(@{$row}); 
          shift(@{$row});            
          $gene_annotation{feature_type}=shift(@{$row});
          $gene_annotation{genome_id}=shift(@{$row});
          $gene_annotation{genome_type}=shift(@{$row});
          shift(@{$row}); 
          $gene_annotation{genome_length}=shift(@{$row});
          $gene_annotation{feature_start}=shift(@{$row});
          $gene_annotation{feature_end}=shift(@{$row});
          $gene_annotation{nuc_len}=shift(@{$row});
          $gene_annotation{aa_len}=shift(@{$row}); 
          $gene_annotation{strand}=shift(@{$row});         
          $gene_annotation{index_on_genome}=shift(@{$row}); 
          $gene_annotation{frame}=shift(@{$row});          
          $gene_annotation{product_description}=shift(@{$row});  
          $gene_annotation{gene_name}="Not available";        
      }
 }

##### get comparative genomics details #####
my $hmm_group_stmt="SELECT * from GenetoSuperFamily WHERE gene_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($hmm_group_data)=get_record($db_dir,$db_name,$hmm_group_stmt);

 if(scalar(@{$hmm_group_data})>0){
      foreach my $row(@{$hmm_group_data}){
         shift(@{$row});
         shift(@{$row});
         shift(@{$row});
         $gene_annotation{hmm_group}=shift(@{$row});
         $gene_annotation{homolog_group}=shift(@{$row});
      }
 }

my $ortholog_group_stmt="SELECT * from MapGeneIdtoGeneFamily WHERE gene_id='$gene_id' and species_abbreviation='$gene_annotation{species_abbreviation}'";
my ($ortho_group_data)=get_record($db_dir,$db_name,$ortholog_group_stmt);

if(scalar(@{$ortho_group_data})>0){
      foreach my $row(@{$ortho_group_data}){
         shift(@{$row});
         $gene_annotation{ortholog_group}=shift(@{$row});
         $gene_annotation{ortholog_group}=~s/\://g;
      }
}

#### get annotation details #######
my $pfam_stmt="SELECT * from DomainAnnotation where protein_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($pfam_data)=get_record($db_dir,$db_name,$pfam_stmt);

if(scalar(@{$pfam_data})>0){     
  $gene_annotation{pfam_domain}=$pfam_data;
}else{
  $gene_annotation{pfam_domain}="No annotation";
}

my $go_stmt="SELECT * from GOAnnotation where protein_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($go_data)=get_record($db_dir,$db_name,$go_stmt);

if(scalar(@{$go_data})>0){     
  $gene_annotation{go_annotation}=$go_data;
}else{
  $gene_annotation{go_annotation}="No annotation";
}

my $ipr_stmt="SELECT * from InterProAnnotation where protein_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($ipr_data)=get_record($db_dir,$db_name,$ipr_stmt);

if(scalar(@{$ipr_data})>0){     
  $gene_annotation{interpro_annotation}=$ipr_data;
}else{
  $gene_annotation{interpro_annotation}="No annotation";
}

my $pathway_stmt="SELECT * from PathwayAnnotation where protein_id='$gene_id' and genome_name='$gene_annotation{species_abbreviation}'";
my ($pathway_data)=get_record($db_dir,$db_name,$pathway_stmt);

if(scalar(@{$pathway_data})>0){     
  $gene_annotation{pathway_annotation}=$pathway_data;
}else{
  $gene_annotation{pathway_annotation}="No annotation";
}

### get sequence details ######


#### get record from database ####
sub get_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   my @row_data=();

   chdir($db_dir);
   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1}) or die $DBI::errstr;   
   my $stmt=$dbh->prepare("$sql_stmt");
   $stmt->execute();

   while(my @row=$stmt->fetchrow_array()){       
      push(@row_data,\@row);
   }

   $stmt->finish;
   $dbh->disconnect();
   chdir($Bin);
   return(\@row_data);
}



print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title></title>
</head>
<body>
<form method="post" name="Gene_Information"> <br>
<table style="text-align: left; width: 70%; height: 70%; margin-left: auto; margin-right: auto;" border="1" cellpadding="2" cellspacing="2">
<tbody>
<tr>
<td colspan="2" rowspan="1"
style="vertical-align: top; width: 305px;">
<div style="text-align: center;"><big style="font-weight: bold;"><big><big>GENE
: $gene_id</big></big></big><br>
<br>
</div>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="width: 305px; background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Genomic Feature:</big></big><br>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1" style="width: 305px;"></td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">Locus Tag:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_id<br>
</td>
</tr>
<tr>
<td style="width: 305px; font-weight: bold;">Species / Strain
Name:<br>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_full_name}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Species Abbreviation:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{species_abbreviation}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome ID: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_id}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome Type: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_type}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genome Length: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{genome_length}<br>
</td>
</tr>
<tr>
<tr>
<td style="font-weight: bold;">Index on Genome: <br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{index_on_genome}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Genomic Location:<br>
</td>
<td style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{feature_start} : $gene_annotation{feature_end} ($gene_annotation{strand})<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Gene Name:<br>
</td>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{gene_name}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Feature Type:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{feature_type}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Protein Length:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{aa_len}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Nucleotide Length:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{nuc_len}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Product Description:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{product_description}<br>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Comparative Genomics Information</big>:</big><br>
</td>
</tr>
<tr>
</tr>
<tr>
<td style="font-weight: bold;">Homolog Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{homolog_group}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">Ortholog Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{ortholog_group}<br>
</td>
</tr>
<tr>
<td style="font-weight: bold;">HMM Model Group:<br>
</td>
<td
style="vertical-align: center; width: 925px; text-align: left;">$gene_annotation{hmm_group}<br>
</td>
</tr>
<tr>
<td colspan="2" rowspan="1"
style="background-color: rgb(204, 204, 204);"><big
style="font-weight: bold;"><big>Annotation</big>:</big><br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">GO Annotation:<br>
</td>
<td>
<br>
<br>
<table style="vertical-align: center; width: 100%; height: 15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="2" cellspacing="2">
<tbody>
<tr><th> GO ID </th><th>GO Category</th><th>GO Description</th></tr>
*;

if($gene_annotation{go_annotation} ne "No annotation"){
      foreach my $row(@{$gene_annotation{go_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $go_id=shift(@{$row});
         my $go_category=shift(@{$row});
         my $go_description=shift(@{$row});

          print qq*
            <tr align="center">
              <td>$go_id</td> <td>$go_category</td> <td> $go_description</td>
            </tr>   
            <tr><td></td><tr>      
          *;
      }
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">PFam:<vr>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th> PFam ID</th><th>Pfam Name</th><th>Description</th><th>Start</th><th>End</th></tr>
*;

if($gene_annotation{pfam_domain} ne "No annotation"){
      foreach my $row(@{$gene_annotation{pfam_domain}}){
         shift(@{$row});
         shift(@{$row});
         shift(@{$row});
         my $pfam_id=shift(@{$row});
         my $pfam_name=shift(@{$row});
         my $pfam_description=shift(@{$row});
         my $pfam_start=shift(@{$row});
         my $pfam_end=shift(@{$row});
         shift(@{$row});
         
         print qq*
            <tr align="center">
              <td>$pfam_id</td> <td>$pfam_name</td> <td> $pfam_description</td><td> $pfam_start</td><td> $pfam_end</td>
            </tr>                    
          *;
      }
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">InterPro:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>InterPro ID</th><th>Description</th></tr>
*;

if($gene_annotation{interpro_annotation} ne "No annotation"){
      my $prev_ipr_id='';
      foreach my $row(@{$gene_annotation{interpro_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $interpro_id=shift(@{$row});
         my $interpro_name=shift(@{$row});
   
         if($interpro_id ne $prev_ipr_id){
            print qq*
              <tr align="center">
                <td>$interpro_id</td> <td>$interpro_name</td>
              </tr>                    
            *;
           $prev_ipr_id=$interpro_id
         }
      }
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
<tr valign="center">
<td style="font-weight: bold;">Pathway:<br>
</td>
<td style="vertical-align: top; width: 70%; text-align: left;"><br><br>
<table style="vertical-align: center; width: 100%; height:15%; margin-left: auto; margin-right: auto;" border="0" cellpadding="5" cellspacing="5">
<tbody>
<tr><th>Pathway ID</th><th>Description</th></tr>
*;

if($gene_annotation{pathway_annotation} ne "No annotation"){
      my $prev_path_id='';
      foreach my $row(@{$gene_annotation{pathway_annotation}}){
         shift(@{$row});
         shift(@{$row});
         my $pathway_id=shift(@{$row});
         my $pathway_name=shift(@{$row});
   
         if($pathway_id ne $prev_path_id){
            print qq*
              <tr align="center">
                <td>$pathway_id</td> <td>$pathway_name</td>
              </tr>                    
            *;
           $prev_path_id=$pathway_id
         }
      }
}else{
   print qq*
      <tr align="center">
        <td>No annotation</td>
       </tr>         
      *;
}
print qq*
</tbody>
</table>
<br>
<br>
</td>
</tr>
</tbody>
</table>
<br>
</form>
</body>
</html>
*;

