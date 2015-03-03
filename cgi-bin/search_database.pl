#!/usr/bin/perl -w
###### ABOUT: This Script Search Database ############
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
my $analysis=param('analysis');
my $core_threshold=param('core_gene_define');
my $show_result_list=param('show_result_list');
my $genome_with_homolog=param('genome_with_homolog');
my $genome_without_homolog=param('genome_without_homolog');
my $keyword=param('keyword');

#### list of genome name ####
$genome_with_homolog=~s/^\://g;
$genome_without_homolog=~s/^\://g;

my @genome_with_homolog=split(":",$genome_with_homolog);
my @genome_without_homolog=split(":",$genome_without_homolog);

my $result_id={};

##### Get Result for the Gene #######
if($analysis eq "core_gene" or $analysis eq "variable_gene"){

   $result_id=get_gene($db_dir,$db_name,$analysis);

}
elsif($analysis eq "unique_gene"){
   ######## get column name #######
   my $sql_get_column_name="PRAGMA table_info(Profile)";

   my($column_data)=get_record($db_dir,$db_name,$sql_get_column_name);

   shift(@{$column_data});
   
   my @result_id=();
   my %unique_gene_count=();

   ### For every genome selected ####
   foreach my $genome_name(@genome_with_homolog){

       my $unique_gene_sql=$genome_name."=1";

       ##### gene_count=1 and for other genomes gene_count=0
       foreach my $row(@{$column_data}){               
              my @column=@{$row};
              my $column_name=$column[1];   

              if($column_name ne $genome_name){
                $unique_gene_sql=$unique_gene_sql." AND ".$column_name."=0";
              } 
       }

       ##### Get unique gene list for the genome and add to main result_id array #####
       ##### also count number of unique gene in each genome and add to unique_gene_count hash ###
       #print $unique_gene_sql;
       my $sql_stmt_unique="SELECT id ,$genome_name FROM Profile WHERE ($unique_gene_sql)";
       my($unique_row_data)=get_record($db_dir,$db_name,$sql_stmt_unique);      

        if(scalar(@{$unique_row_data})>0){
             my $gene_count=1;
             foreach my $row(@{$unique_row_data}){
               my @column=@{$row};
               my $count_genome=0;
               my $id=shift(@column);                            
               push(@result_id,$id);               
               $gene_count++;
             }           
           $unique_gene_count{$genome_name}=$gene_count;
       }
  }
$result_id=\@result_id;  
}

my $number_group=scalar(@{$result_id});

my($result_gene)=get_gene_info($db_dir,$db_name,$result_id,\@genome_with_homolog,$keyword);


######### FIND GENES ################
sub get_gene{

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($analysis)=(shift);

   my $homolog_in_genome='';
   my $homolog_notin_genome='';
   my $show_column='';
   my $sql_stmt='';
   my $define_core_genome=0;
   my @result_id=();

   if(scalar(@genome_with_homolog)>=1) {
     $homolog_in_genome=join("=1 OR ",@genome_with_homolog);
     $homolog_in_genome=$homolog_in_genome."=1";
     $show_column=join(", ",@genome_with_homolog);
     $define_core_genome=int((($core_threshold/100)*scalar(@genome_with_homolog)));
   }

   if(scalar(@genome_without_homolog)>0){
      $homolog_notin_genome=join("=0 AND ",@genome_without_homolog);
      $homolog_notin_genome=$homolog_notin_genome."=0";
   } 

   if($homolog_notin_genome eq ''){     
     $sql_stmt="SELECT id ,$show_column FROM Profile WHERE ($homolog_in_genome)";
   }else{
     $sql_stmt="SELECT id ,$show_column FROM Profile WHERE ($homolog_in_genome) AND ($homolog_notin_genome)";
   }

   my($row_data)=get_record($db_dir,$db_name,$sql_stmt); 
 
   if(scalar(@{$row_data})>0){
      foreach my $row(@{$row_data}){
          my @column=@{$row};
          my $count_genome=0;
          my $id=shift(@column);
        
          foreach(@column){
            $count_genome=$count_genome + $_;
          }

        if($analysis eq "core_gene"){ 
          if($count_genome>=$define_core_genome){
              push(@result_id,$id);
          }
        }elsif($analysis eq "variable_gene"){          
          if($count_genome>1 and $count_genome<$define_core_genome){              
              push(@result_id,$id);
          }elsif($define_core_genome==2){
              push(@result_id,$id);
          }
        }
      }
   }
   return(\@result_id);
}

#### get gene ids in each family #####
sub get_gene_info {

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my(@result_id)=@{(shift)};
   my(@genome_with_homolog)=@{(shift)};
   my($keyword)=(shift);
   my $sql_stmt='';
   my $genome='';

   if(scalar(@genome_with_homolog)>1){
     $genome=join("' , '",@genome_with_homolog);
   }else{
     $genome=$genome_with_homolog[0];
   }
   
   my $group_ids=join("' , '",@result_id);
        
   $sql_stmt="SELECT * from MapGeneIdtoGeneFamily WHERE genefamily_id IN ('$group_ids') and species_abbreviation IN ('$genome')";

   my($result_gene)=get_record($db_dir,$db_name,$sql_stmt);

   #### If gene list by keyword #####
   if($keyword ne ''){
      my @row_data=@{$result_gene};
      my @row_data2=@{$result_gene};
      my $list_gene_id='';      
           
      foreach my $row(sort{$a->[1] cmp $b->[1]}@row_data){
         my @row=@{$row}; 
         my $index_id=$row[0];
         my $group_id=$row[1];
         my $gene_id=$row[2];
         my $genome_abbrv=$row[3];     
         $list_gene_id=$list_gene_id.",'".$gene_id."'";         
      }
      $list_gene_id=~s/^\,//g;
 
      my $keyword_sql_stmt="Select Distinct(feature_id),description from GeneFeature where feature_id IN ($list_gene_id) and description LIKE '%$keyword%'";
  
      my($gene_with_keyword)=get_record($db_dir,$db_name,$keyword_sql_stmt); 
      my %gene_id_keyword=();
      my %group_id_keyword=();
      my @keyword_result_gene=();     

      if(scalar(@{$gene_with_keyword})>=1){           
           foreach my $row(@{$gene_with_keyword}){           
             my $feature_id=shift(@{$row}); 
             my $description=shift(@{$row});                                      
             $gene_id_keyword{$feature_id}=$description;             
           }
      } 

      my $prev_group_id='';
      my @group_With_keyword=();
      my $has_keyword=0;
      
      #### Search for group with atleast one gene having assigned keyword #### 
      foreach my $row(sort{$a->[1] cmp $b->[1]}@row_data){
         my @row=@{$row}; 
         my $index_id=$row[0];
         my $group_id=$row[1];
         my $gene_id=$row[2];
         my $genome_abbrv=$row[3];   
       
         if($group_id ne $prev_group_id){

            if($has_keyword eq 1){
              @keyword_result_gene=(@keyword_result_gene , @group_With_keyword);
            }
            $prev_group_id=$group_id;
            $has_keyword=0;
            @group_With_keyword=();

             if($gene_id_keyword{$gene_id}){                          
                $has_keyword=1;
             }
             push(@group_With_keyword,$row);

         }else{
            if($gene_id_keyword{$gene_id}){                          
               $has_keyword=1;
            }
            push(@group_With_keyword,$row);
         }
      } 
   
      $result_gene=\@keyword_result_gene;
   } 
  
   return($result_gene);
}

sub get_genome_name {
 
   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($genome_abbrv)=(shift);
   my $genome_name='';

   my $sql_stmt="SELECT Distinct(genome_name),abbreviation from OrganismInfo WHERE abbreviation IN ('$genome_abbrv')";
   my ($row_data)=get_record($db_dir,$db_name,$sql_stmt);

  return($row_data);
}

sub get_gene_description{

   my($db_dir)=(shift);
   my($db_name)=(shift);
   my($gene_id)=(shift);
   my($genome_abbrv)=(shift);
   my($keyword)=(shift);
   my $genome_description='';
   my $sql_stmt='';

   if($keyword eq ''){
     $sql_stmt="SELECT description from GeneFeature WHERE feature_id='$gene_id'";
   }else{     
     $sql_stmt="SELECT description from GeneFeature WHERE feature_id='$gene_id' and description LIKE '%$keyword%'";
   }

   my ($row_data)=get_record($db_dir,$db_name,$sql_stmt);

    if(scalar(@{$row_data})>0){
      foreach my $row(@{$row_data}){
          $genome_description=shift(@{$row});
      }
    }else{
        $genome_description="Unannotated protein sequence";
    }
  return($genome_description);
}

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


##################################### HTML CODE STARTS HERE FOR GENE FAMILY WISE SORTING TABLE ###############################
if($show_result_list eq "Gene Family ID" and $analysis ne "variable_gene" and $analysis ne "unique_gene"){

print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title>result_analysis</title>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery-latest.js"></script>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery.tablesorter.js"></script>

<script type="text/javascript">

    function change_genome(genome_name_id){

       var genome_name=document.getElementById(genome_name_id).options[document.getElementById(genome_name_id).selectedIndex].value;

       var group_id_list=document.getElementsByClassName('gene_list');

       for(var j=0;j<group_id_list.length;j++){

         var group_id=group_id_list[j].id;

         var feature_class_name=genome_name+'_'+group_id;
         var id_desc_td='desc_'+group_id;
         var id_radio_detail='detail_'+group_id;       
         var gene_list_id=group_id;
         var option_count=0;
         var selected_option=1;

       var gene_feature_elements=document.getElementsByClassName(feature_class_name);

       var gene_select=document.getElementById(gene_list_id);

       for(var i=gene_select.length;i>=0;i--){
           gene_select.remove(i);
       }

       for(var i=0;i<gene_feature_elements.length;i++){

          var element_id=gene_feature_elements[i].id; 
                     
          var gene_id=document.getElementById(element_id).value;

          if(element_id=='detail_'+gene_id){

              var description_id='description_'+gene_id;

              var gene_description=document.getElementById(description_id).value;
              
              if(selected_option==1){                 
                 var gene_option=document.createElement("option");
                 gene_option.text=gene_id;
                 gene_select.add(gene_option,gene_select[option_count]);
                 gene_select.selectedIndex=option_count;

                 document.getElementById(id_desc_td).innerHTML=gene_description;
                 document.getElementById(id_radio_detail).value=gene_id;
                 selected_option++;

              }else{              
                 var gene_option=document.createElement("option");
                 gene_option.text=gene_id;
                 gene_select.add(gene_option,gene_select[option_count]);
  
                 document.getElementById(id_desc_td).innerHTML=gene_description;
                 document.getElementById(id_radio_detail).value=gene_id;                    
              }  

              if(selected_option==1){
                 document.getElementById(id_desc_td).innerHTML=gene_description;
                 document.getElementById(id_radio_detail).value=gene_id;
              }         
              option_count++;              
           }
        } 
      }    
    }

    function change_gene(id){

       var listgene_id='listgene_'+id;

       var gene_id=document.getElementById(listgene_id).options[document.getElementById(listgene_id).selectedIndex].text;

       var id_hidden_gene_desc='description_' + gene_id;    
      
       var gene_desc=document.getElementById(id_hidden_gene_desc).value;
  
       var id_desc_td='desc_'+id;
       var id_radio_detail='detail_'+id;

       document.getElementById(id_desc_td).innerHTML=gene_desc;
       document.getElementById(id_radio_detail).value=gene_id;       
    }

</script>

<script>
\$(document).ready(function() { 
  \$('#filter_word').bind('keyup', function() {
    var s = new RegExp(this.value,'i');
    \$('#result_table tbody tr').each(function() {        
        if(s.test(this.innerHTML)) \$(this).show();
        else \$(this).hide();
     });
  });
});
</script>

</head>

<body>
<div style="text-align: center;"><big><big><big>Result of Analysis</big></big></big><br>
$number_group shown in result
<br>
<form method="post" name="result_set" action="get_gene_detail.pl">
<div style="text-align: center;"></div>
<div style="text-align: left;">
Filter row with term: </span><input maxlength="100" size="25" name="filter_word" id="filter_word"><br> <br>
</div>

<div>
<table id="genome" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="0" cellpadding="1" cellspacing="1">
<tbody>
<tr>
<td colspan="4" rowspan="1" style="vertical-align: top;">
<span style="font-weight: bold;"> Genome Name:
<select name="genome_name" id="genome_name" onchange="change_genome(this.id)">
*;
##### selection list of genome names for faimlywise representation #####
my $list_genome=join("' , '",@genome_with_homolog);   #### Comman-sperated vales of genome names in which homolog is present #####
my $row_genome_name=get_genome_name($db_dir,$db_name,$list_genome); 

if(scalar(@{$row_genome_name})>0){
     my $count_genome=1;
     foreach my $row(sort{$a->[1] cmp $b->[1]} @{$row_genome_name}){
            my $genome_name=shift(@{$row});
            my $genome_abbreviation=shift(@{$row});

            if($count_genome==1){
               print qq*
                  <option selected="selected" value="$genome_abbreviation">$genome_name</option>				 
               *;
               $count_genome++;
            }else{
               print qq*
                  <option value="$genome_abbreviation">$genome_name</option>
               *; 
            }
     }
} 
print qq*
</select>
</span> 
</td>
<td rowspan="1" style="vertical-align: top; text-align: centre;"><input name="Export_Result_Table" value="Export Result Table" type="button"><br></td>
</td>
<td rowspan="1" style="vertical-align: top; text-align: right;"><input name="Submit_Gene_ID" value="Show Gene Information" type="submit"><br></td>
</td>
</tr>
</tbody>
</table>
<br>
</div>

<table id="result_table" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="1" cellpadding="1" cellspacing="1">
<thead>
<tr>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">Group
ID<br>
</th>
<th
style="vertical-align: top; font-weight: bold; width: 335px; text-align: center;">Gene
ID<br>
</th>
<th
style="vertical-align: middle; font-weight: bold; width: 632px; text-align: center; white-space: nowrap;">&nbsp;Gene
Description<br>
</th>
<th
style="vertical-align: top; width: 171px; font-weight: bold; text-align: center;">Genes per Genome<br>
</th>
<th
style="vertical-align: top; width: 48px; font-weight: bold; text-align: center;">Details<br>
</th>
</tr>
</thead>
<tbody>
*;

###### load result set in the html table #######
my @row_data=@{$result_gene};

if(scalar(@row_data)>0){
 
  my $prev_group_id='';
  my $gene_count=1;
  my $genome_name='';
  my $gene_description='';
  my $selected_genome_abbrv='';
  my $detail_gene='';  
  my @gene_feature=();
  my %genome_abbrv=();

  open(EXPORT_RESULT,">$Bin/../tmp/exported_result.txt");

  system("chmod 755 $Bin/../tmp/exported_result.txt");

  my @export_result=();

  #### sort by group id and genome name #####
  foreach my $row(sort{$a->[1] cmp $b->[1] || $a->[3] cmp $b->[3]}@row_data){

         shift(@{$row});
         my $group_id=shift(@{$row});
         my $gene_id=shift(@{$row});
         my $genome_abbrv=shift(@{$row});

         my @store_result=();

         	 
         if($group_id ne $prev_group_id){

              if($prev_group_id eq ''){

                ###### print group id in first result column ######
                print qq*
                  <tr>
                  <td style="vertical-align: top; width: 162px; text-align: left;" id="group_$group_id" class="group_id">$group_id<br>
                  </td>
                *;
                $prev_group_id=$group_id;
                $gene_count=1;
                @gene_feature=();
                %genome_abbrv=();

                print EXPORT_RESULT "$group_id\t";
                push(@store_result,$group_id);

              }else{ 

                ######## drop-down list of gene id for each gene family in selected genome ######        
                print qq*
                  </td>
                  <td style="vertical-align: top; width: 335px;" name="gene_$prev_group_id">
                  <select name="listgene_$prev_group_id" id="$prev_group_id" class="gene_list" onchange="change_gene($prev_group_id)">
                *;
                 
                    my $gene_per_genome=0;

                    for(my $i=0;$i<scalar(@gene_feature);$i=$i+4){

                            my $feature_id=$gene_feature[$i];
                            my $feature_genome=$gene_feature[$i+1];
                            my $feature_abbrv=$gene_feature[$i+2];
                            my $feature_desc=$gene_feature[$i+3];
                      
                          if($feature_abbrv eq $selected_genome_abbrv){

                               if($feature_id eq $detail_gene){
                                 print qq*
                                   <option selected="selected">$feature_id</option>
                                 *;
                               }else{
                                 print qq*
                                   <option>$feature_id</option>
                                 *;
                               }
                             $gene_per_genome++;
                           }
                     }      
                 
                print qq*
                         </select>                     
                     </td>
                     <td style="vertical-align: top; width: 632px;" id="desc_$prev_group_id">$gene_description<br>
                     </td>
                     <td style="vertical-align: top; width: 171px; font-weight: bold; text-align: center;">$gene_per_genome<br>
                     </td>
                     <td style="vertical-align: top; width: 48px; text-align: center;">
                     <input checked="checked" name="gene_detail_radio" id="detail_$prev_group_id" type="radio" value="$detail_gene">
                     <input name="project_name" id="project_name" type="hidden" value="$project_name">
                     <input name="db_name" id="db_name" type="hidden" value="$db_name">
                     <input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
                     <input name="report_dir" id="report_dir" type="hidden" value="$report_dir">
                *;
                
                ####### details for radio button, hidden values to pass for gene details #####
                for(my $i=0;$i<scalar(@gene_feature);$i=$i+4){

                    my $feature_id=$gene_feature[$i];
                    my $feature_genome=$gene_feature[$i+1];
                    my $feature_abbrv=$gene_feature[$i+2];
                    my $feature_desc=$gene_feature[$i+3];                  
                    my $class_name=$feature_abbrv."_".$prev_group_id;

                   

                     print qq*
                        <input class="$class_name" name="hidden_genome_name" id="genome_$feature_id" value="$feature_genome" type="hidden">
                        <input class="$class_name" name="hidden_gene_description" id="description_$feature_id"  value="$feature_desc" type="hidden">
                        <input class="$class_name" name="hidden_gene_detail" id="detail_$feature_id" value="$feature_id" type="hidden">
                    *; 
                 }               

                print qq*
                    </td>
                 </tr> 
                *;
                ###### print next group id in result column on next row of result table and re-initialize the variables ######
                print qq*
                  <tr>
                    <td style="vertical-align: top; width: 162px; text-align: left;" id="group_$group_id" class="group_id">$group_id<br>
                  </td>
                *;
                $prev_group_id=$group_id;
                $gene_count=1;
                @gene_feature=();
                %genome_abbrv=();
                print EXPORT_RESULT "$group_id\t";
                push(@store_result,$group_id);
            }
        }

        ##### get description for each gene in a selected genome and add to gene feature array#####
        my $row_genome=get_genome_name($db_dir,$db_name,$genome_abbrv); 
        my $description=get_gene_description($db_dir,$db_name,$gene_id,$genome_abbrv,$keyword); 
        my $detail=$gene_id; 
        my $genome='';

        if(scalar(@{$row_genome})>0){
            foreach my $row(@{$row_genome}){
              $genome=shift(@{$row});              
            }
        }  
        push(@gene_feature,$detail);
        push(@gene_feature,$genome);
        push(@gene_feature,$genome_abbrv);
        push(@gene_feature,$description);         

        ###### Selection list of Genome ####
        if(!defined($genome_abbrv{$genome_abbrv})){
            ##### first gene to be displayed ###
            if($gene_count==1){
               $genome_name=$genome;
               $gene_description=$description;
               $detail_gene=$detail;
               $selected_genome_abbrv=$genome_abbrv;
               $gene_count++; 
               print EXPORT_RESULT "$gene_description\n";
               push(@store_result,$gene_description);
            }
            $genome_abbrv{$genome_abbrv}=1;
        }   

    push(@export_result,\@store_result);    
   }

print qq*
</td>
<td style="vertical-align: top; width: 335px;" name="gene_$prev_group_id">
<select name="listgene_$prev_group_id" id="$prev_group_id" class="gene_list" onchange="change_gene($prev_group_id)">
*;
my $gene_per_genome=0;

for(my $i=0;$i<scalar(@gene_feature);$i=$i+4){

    my $feature_id=$gene_feature[$i];
    my $feature_genome=$gene_feature[$i+1];
    my $feature_abbrv=$gene_feature[$i+2];
    my $feature_desc=$gene_feature[$i+3];
                      
    if($feature_abbrv eq $selected_genome_abbrv){
        if($feature_id eq $detail_gene){
           print qq*
           <option selected="selected">$feature_id</option>
           *;
        }else{
           print qq*
           <option>$feature_id</option>
           *;
        }
       $gene_per_genome++;
    }
}

print qq*
</select>                     
</td>
<td style="vertical-align: top; width: 632px;" id="desc_$prev_group_id">$gene_description<br>
</td>
<td
style="vertical-align: top; width: 171px; font-weight: bold; text-align: center;">$gene_per_genome<br>
</td>
<td style="vertical-align: top; width: 48px; text-align: center;">
<input checked="checked" name="gene_detail_radio" id="detail_$prev_group_id" type="radio" value="$detail_gene">
*;

for(my $i=0;$i<scalar(@gene_feature);$i=$i+4){

    my $feature_id=$gene_feature[$i];
    my $feature_genome=$gene_feature[$i+1];
    my $feature_abbrv=$gene_feature[$i+2];
    my $feature_desc=$gene_feature[$i+3];
    my $class_name=$feature_abbrv."_".$prev_group_id;
    
    print qq*
        <input class="$class_name" name="hidden_genome_name" id="genome_$feature_id" value="$feature_genome" type="hidden">
        <input class="$class_name" name="hidden_gene_description" id="description_$feature_id"  value="$feature_desc" type="hidden">
        <input class="$class_name" name="hidden_gene_detail" id="detail_$feature_id" value="$feature_id" type="hidden">
    *; 
} 
print qq*
  </td>
  </tr>
*;

}          
print qq*
</tbody>
</table>
<br>
</form>
</div>
</body>
</html>
*;
}

##########################################################################
########## Show result sorted by gene and strains (for variable genes) #####

if($show_result_list eq "Gene ID" or $analysis eq "variable_gene" or $analysis eq "unique_gene"){
my $number_row=scalar(@{$result_gene});
print qq*
<html>
<head>
<meta content="text/html; charset=ISO-8859-1"
http-equiv="content-type">
<title>result_analysis</title>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery-latest.js"></script>
<script type="text/javascript" src="jquery/jquery_tablesorter/jquery.tablesorter.js"></script>

<script>
\$(document).ready(function() { 
    // call the tablesorter plugin 
    \$("#result_table").tablesorter(); 
});
</script>

<script>
\$(document).ready(function() { 
  \$('#filter_word').bind('keyup', function() {
    var s = new RegExp(this.value,'i');
    \$('#result_table tbody tr').each(function() {
        if(s.test(this.innerHTML)) \$(this).show();
        else \$(this).hide();
    });
  });
});
</script>

</head>
<body>
<div style="text-align: center;"><big><big><big>Result of Analysis</big></big></big><br>
 $number_row shown in result
<br>
<form method="post" name="result_set" action="get_gene_detail.pl">
<div style="text-align: right;">
 <input name="Submit_Gene_ID" value="Show Gene Information" type="submit">
</div>
<div style="text-align: left;">
Filter row with term: </span><input maxlength="100" size="25" name="filter_word" id="filter_word"><br>
</div>
<table id="result_table" class="tablesorter" style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;"
border="1" cellpadding="1" cellspacing="1">
<thead>
<tr>
<th
style="vertical-align: top; width: 162px; font-weight: bold; text-align: center;">Group
ID<br>
</th>
<th
style="vertical-align: top; width: 171px; font-weight: bold; text-align: center;">Gene
ID<br>
</th>
<th
style="vertical-align: top; font-weight: bold; width: 335px; text-align: center;">Genome
Name<br>
</th>
<th
style="vertical-align: middle; font-weight: bold; width: 632px; text-align: center; white-space: nowrap;">&nbsp;Gene
Description<br>
</th>
<th
style="vertical-align: top; width: 48px; font-weight: bold; text-align: center;">Details<br>
</th>
</tr>
</thead>
<tbody>
*;

###### load result set in the html table #######
my @row_data=@{$result_gene};

if(scalar(@row_data)>0){

  if($analysis eq "variable_gene" or $analysis eq "unique_gene"){
    @row_data=sort{$a->[3] cmp $b->[3] || $a->[1] cmp $b->[1]}@row_data
  }else{
    @row_data=sort{$a->[1] cmp $b->[1] || $a->[3] cmp $b->[3]}@row_data
  }
 
  foreach my $row(@row_data){
         shift(@{$row});
         my $group_id=shift(@{$row});
         my $gene_id=shift(@{$row});
         my $genome_abbrv=shift(@{$row});
         my $genome='';
       
        my $row_genome=get_genome_name($db_dir,$db_name,$genome_abbrv); 
        my $description=get_gene_description($db_dir,$db_name,$gene_id,$genome_abbrv,$keyword);

        if(scalar(@{$row_genome})>0){
           foreach my $row(@{$row_genome}){
              $genome=shift(@{$row});              
           }
        }
               print qq*
                  <tr>
                    <td style="vertical-align: top; width: 162px; text-align: left;" id="group_$group_id" class="group_id">$group_id<br>
                    </td>
                    <td style="vertical-align: top; width: 171px; text-align: left;" class="gene_id">$gene_id<br>
                    </td>
                    <td style="vertical-align: top; width: 335px;" id="genome_$group_id">$genome<br>
                        <input name="genome_name" id="genome_name" type="hidden" value="$genome_abbrv">
                    </td>
                    <td style="vertical-align: top; width: 632px;" id="desc_$group_id">$description<br>
                    </td>
                    <td style="vertical-align: top; width: 48px; text-align: center;">
                    <input checked="checked" name="gene_detail_radio" id="detail_$group_id" value="$gene_id" type="radio">
                     <input name="project_name" id="project_name" type="hidden" value="$project_name">
                     <input name="db_name" id="db_name" type="hidden" value="$db_name">
                     <input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
                     <input name="report_dir" id="report_dir" type="hidden" value="$report_dir">
                    </td>
                  </tr>
                 *;     
   }
}
print qq*
</tbody>
</table>
<br>
</form>
</div>
</body>
</html>
*;
}


