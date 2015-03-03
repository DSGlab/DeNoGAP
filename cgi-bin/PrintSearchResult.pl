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

open(EXPORT_RESULT,">exported_result.txt");

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
 
       print "$sql_stmt_unique\n";

       my($unique_row_data)=get_record($db_dir,$db_name,$sql_stmt_unique);      

        if(scalar(@{$unique_row_data})>0){
             my $gene_count=1;
             foreach my $row(@{$unique_row_data}){
               my @column=@{$row};
               my $count_genome=0;
               my $id=shift(@column);                            
               push(@result_id,$id);

               print EXPORT_RESULT "$id\t$genome_name\n";
               
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

     print "$sql_stmt\n";
   }else{
     $sql_stmt="SELECT id ,$show_column FROM Profile WHERE ($homolog_in_genome) AND ($homolog_notin_genome)";

     print "$sql_stmt\n";
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
 
              print EXPORT_RESULT "$id\n";
          }
        }elsif($analysis eq "variable_gene"){          
          if($count_genome>1 and $count_genome<$define_core_genome){              
              push(@result_id,$id);

              print EXPORT_RESULT "$id\n";

          }elsif($define_core_genome==2){
              push(@result_id,$id);

              print EXPORT_RESULT "$id\n";
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


my $number_row=scalar(@{$result_gene});

###### load result set in the html table #######
my @row_data=@{$result_gene};

if(scalar(@row_data)>0){

    if($analysis eq "variable_gene" or $analysis eq "unique_gene"){
        @row_data=sort{$a->[3] cmp $b->[3] || $a->[1] cmp $b->[1]}@row_data
    }else{
        @row_data=sort{$a->[1] cmp $b->[1] || $a->[3] cmp $b->[3]}@row_data
    }
    
    system("chmod 755 exported_result.txt");

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

         print EXPORT_RESULT "$group_id\t$genome_abbrv\t$genome\t$gene_id\t$description\n";
    }

    close EXPORT_RESULT;
}



