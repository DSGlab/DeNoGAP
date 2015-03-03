#!/usr/bin/perl -w
###### ABOUT: This Script Create the Index page for pipeline report in html format ############
###### AUTHOR:Shalabh Thakur###################################################################


use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use File::Basename;
use DBI;
use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);

my $project_name=param('project_name');
my $database_file_path=param('database_file');
my $web_root=param('web_root');


my $db_dir=dirname($database_file_path);
my $db_name=basename($database_file_path);
my $report_dir="$web_root/$project_name";


##### html index file ####

print "Content-type:text/html\n\n";

###### start print html code in the file #######

print qq* 

<!DOCTYPE HTML>
<html>
  <head>
     <meta content="text/html; charset=ISO-8859-1"
         http-equiv="content-type">
         <title>Report_main_page</title>
     <link rel="stylesheet" type="text/css" href="jquery/jquery_datatable/jquery.dataTables.css">
     <script type="text/javascript">

      <!-- Function to select all check box for list of genomes -->

           function select_all(id, id2){
                var id_value=id.value;
                var id2_value=id2.value;               
                
                if(id.checked==true){
                     id.checked=true;
                     check_all(id_value);
                     uncheck_all(id2_value);
                     id2.checked=false;
                }
                else if(id.checked==false){
                     id.checked=false;
                     uncheck_all(id_value);                  
                }
            }

        <!-- Function to check all check box for list of genomes -->

            function check_all(id_value){
              var list_checkbox=document.getElementsByClassName(id_value);
              for(var i=0;i<list_checkbox.length;i++){                      
                      list_checkbox[i].checked=true;              
              } 
            }

        <!-- Function to uncheck all check box for list of genomes -->

           function uncheck_all(id_value){
             var list_checkbox=document.getElementsByClassName(id_value);
             for(var i=0;i<list_checkbox.length;i++){
                 list_checkbox[i].checked=false;
             }
           }
  
        <!-- Function to deselect check box if another option is checked -->

           function uncheck_set(deselect_id){
              deselect_id.checked=false;
           }

        <!-- Function to compile list of genome with homolog and without homolog after submit form button is clicked -->

           function GenomeList(){

               var list_with_homolog='';
               var list_without_homolog='';

               var list_checkbox_set1=document.getElementsByClassName('WithHomolog');
 
               for(var i=0;i<list_checkbox_set1.length;i++){
                  if(list_checkbox_set1[i].checked==true){
                     list_with_homolog=list_with_homolog+':'+ list_checkbox_set1[i].value;
                  }
               }

               var list_checkbox_set2=document.getElementsByClassName('WithoutHomolog');
 
               for(var i=0;i<list_checkbox_set2.length;i++){
                  if(list_checkbox_set2[i].checked==true){
                     list_without_homolog=list_without_homolog + ':' + list_checkbox_set2[i].value;
                  }
               }

               document.getElementById('genome_with_homolog').value=list_with_homolog;
               document.getElementById('genome_without_homolog').value=list_without_homolog;       
           }

</script>

<script type="text/javascript" src="jquery/jquery_tablesorter/jquery-latest.js"></script>

<script type="text/javascript" src="jquery/jquery_tablesorter/jquery.tablesorter.js"></script>

<script type="text/javascript" src="jquery/jquery_datatable/jquery.dataTables.js"></script>

<script type="text/javascript" src="jquery/jquery_datatable/jquery.dataTables.min.js"></script>

<script>
\$(document).ready(function() { 
    // call the tablesorter plugin 
    \$("#genome_name").tablesorter(); 
}); 
</script>

<script>
\$(document).ready(function() {
    var rowCount = \$('#genome_name tr').length;
    \$('#genome_name').dataTable({"iDisplayLength": rowCount});    
} );
</script>

</head>

  <body>
    <h1 style="text-align: center;"><small><span style="font-style: italic; font-weight: bold;">Analysis Report of $project_name</span></small></h1>
   <div style="text-align: center;">
   <form name="report_species_compare" id="report_species_compare" method="post" action="search_database.pl" onsubmit="GenomeList()">
   <table style="text-align: left; width: 100%; margin-left: auto; margin-right: auto;" border="0" cellpadding="2" cellspacing="2">
     <tbody>
       <tr>
         <td style="vertical-align: top; width: 451px;"><big><span style="font-weight: bold;">Find in selected genomes:</span></big><br>
         </td>
       </tr>
       <tr style="font-weight: bold;">
           <td style="vertical-align: top; width: 451px;"><input checked="checked" name="analysis" value="core_gene" type="radio">Core Genes<br>
          </td>
       </tr>
    <tr style="font-weight: bold;">
        <td style="vertical-align: top; width: 451px;"><input name="analysis" value="variable_gene" type="radio">Variable Genes<br>
        </td>
    </tr>
    <tr style="font-weight: bold;">
      <td style="vertical-align: top; width: 451px;"><input name="analysis" value="unique_gene" type="radio">Unique Genes<br>
    </td>
    </tr>
    <tr>
<td style="vertical-align: top;"><span
style="font-weight: bold;">&nbsp;Define core gene as present in % of
genomes:</span>
<select name="core_gene_define">
<option selected="selected">100</option>
<option>99</option>
<option>98</option>
<option>97</option>
<option>96</option>
<option>95</option>
<option>94</option>
<option>93</option>
<option>92</option>
<option>91</option>
<option>90</option>
<option>89</option>
<option>88</option>
<option>87</option>
<option>86</option>
<option>85</option>
<option>84</option>
<option>83</option>
<option>82</option>
<option>81</option>
<option>80</option>
<option>79</option>
<option>78</option>
<option>77</option>
<option>76</option>
<option>75</option>
</select>
</td>
</tr>
<tr>
<td style="vertical-align: top;">&nbsp;<span
style="font-weight: bold;">Show result sorted by:
<select name="show_result_list">
<option selected="selected">Gene Family ID</option>
<option>Gene ID</option>
</select>
</span><br>
</td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-weight: bold;">Search for gene description : </span><input maxlength="100" size="25" name="keyword" id="keyword"><br></td>
</tr>
    <tr>
      <td><input name="project_name" id="project_name" type="hidden" value="$project_name">
      </td>
    </tr>
     <tr>
      <td><input name="db_name" id="db_name" type="hidden" value="$db_name">
      </td>
    </tr>
    <tr>
      <td><input name="db_dir" id="db_dir" type="hidden" value="$db_dir">
      </td>
    </tr>
    <tr>
      <td><input name="report_dir" id="report_dir" type="hidden" value="$report_dir">
      </td>
    </tr>
    <tr>
      <td colspan="1" style="text-align: left;" ><input name="get_genes" value="Submit" type="submit">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
          <input name="reset_field" value="Reset" type="reset"><br>
      </td>

    </tr>
</tbody>
</table>
<br>
   <table id="genome_name" class="tablesorter" style="text-align: left; height: 50%; width: 100%; margin-left: auto; margin-right: auto;" border="1" cellpadding="2" cellspacing="2">
      <thead>
         <tr>
           <td style="vertical-align: top; text-align: left; width: 7%;">
               <input style="font-weight: bold;" name="set1" value="WithHomolog" id="set1" type="checkbox" onclick="select_all(this,set2)"><span style="font-weight: bold;">With Homolog <input name="genome_with_homolog" id="genome_with_homolog"
                type="hidden"></span><br>
           </td>

           <td style="vertical-align: top; text-align: left; width: 7%;">
             <input style="font-weight: bold;" name="set2" value="WithoutHomolog" id="set2" type="checkbox" onclick="select_all(this,set1)"><span style="font-weight: bold;">Without Homolog <input name="genome_without_homolog" id="genome_without_homolog"
type="hidden"></span><br>
           </td>
*;

   ######## get column name #######
   my $sql_get_column_name="PRAGMA table_info(OrganismInfo)";
   my($row_data1)=get_record($db_dir,$db_name,$sql_get_column_name);

   if(scalar(@{$row_data1}>0)){
      shift(@{$row_data1});
      foreach my $row(@{$row_data1}){               
           my @column=@{$row};
           my $column_name=$column[1];
           if($column_name eq "outgroup" or $column_name eq "genome_type"){next;}   
           print qq*
              <th style="vertical-align: top; font-weight: bold; text-align: center;">
                  $column_name<br>
              </th>
           *;      
      } 
   }
print qq*
</tr>
</thead>
<tbody>
<tr>
*;
   ######## get organism information #############
   my $sql_org_data="SELECT * from OrganismInfo ORDER BY genome_name ASC";
   my($row_data2)=get_record($db_dir,$db_name,$sql_org_data);

   if(scalar(@{$row_data2}>0)){

      foreach my $row(@{$row_data2}){ 
              
           my @column=@{$row};
           shift(@column);
           my $abbrv=$column[2];
           splice(@column,3,1);
           splice(@column,3,1);

           print qq*

                 <td style="vertical-align: top; text-align: left;">
                    <input style="font-weight: bold;" name="set1_$abbrv" value="$abbrv" id="set1_$abbrv" class="WithHomolog" type="checkbox" onclick="uncheck_set(set2_$abbrv)"><span style="font-weight: bold;"></span><br>
                 </td>

                 <td style="vertical-align: top; text-align: leftr;">
                    <input style="font-weight: bold;" name="set2_$abbrv" value="$abbrv" id="set2_$abbrv" class="WithoutHomolog" type="checkbox" onclick="uncheck_set(set1_$abbrv)"><span style="font-weight: bold;"></span><br>
                 </td>
           *;
            
           foreach my $column_value(@column){
                print qq*
                   <td style="vertical-align: top; font-weight: bold; text-align: left;">
                     $column_value<br>
                   </td>
                *;  
           }  

           print qq*
              </tr> 
           *; 
      } 
   } 

print qq*
      </tbody>
    </table>
  </form>
</body>
</html>       
*;


sub get_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   my @row_data=();

   chdir($db_dir) or die "Cannot access directory $db_dir";

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


