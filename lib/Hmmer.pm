##### Module to access Hmmer package #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package Hmmer;
use strict;
use Exporter;
use File::Basename;
use Tie::File;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(Run Sort_table Read_domain_table Read_aligned_sequence);


sub Run {
    
    my($hmm_program)=(shift);
    my($hmm_file)=(shift);
    my($DB_file)=(shift);
    my($HMMFULL_OUT)=(shift);
    my($HMMTBL_OUT)=(shift);
    my($HMMDOM_OUT)=(shift);    
    my(%HMMER_Param)=%{(shift)};
    my($db_size)=(shift);

    my $Param='';

    my $outfile=basename($hmm_file);
       $outfile=~s/(\.)(\w+)/\.out/g;   

    foreach(keys %HMMER_Param){

       if($_=~/^(o)$/){
          $Param=$Param." -".$_." ".$HMMFULL_OUT."/"."full_".$hmm_program."_".$outfile;
       }elsif($_=~/^(tblout)$/){
          $Param=$Param." --".$_." ".$HMMTBL_OUT."/"."tbl_".$hmm_program."_".$outfile;
       }elsif($_=~/^(domtblout)$/){
          $Param=$Param." --".$_." ".$HMMDOM_OUT."/"."dom_".$hmm_program."_".$outfile;
       }elsif($_=~/^([E|A|T|Z])$/){
          $Param=$Param." -".$_." ".$HMMER_Param{$_};
       }else{
          $Param=$Param." --".$_." ".$HMMER_Param{$_};
       }
    }  

    $Param=$Param." -Z ".$db_size;

    #$Param=$Param." --cpu 5";
 
    if($hmm_program=~/phmmer/){
      system("$hmm_program $Param $hmm_file $DB_file");
    }elsif($hmm_program=~/hmmscan/){
      system("$hmm_program $Param $DB_file $hmm_file");
    }     
}

sub Sort_table {

    my($query_file)=(shift);
    my($HMMDOM_OUT)=(shift);
     
    my $file_line=0;
    my $file=basename($query_file);
       $file=~s/(\.)(\w+)/\.out/g;  
       $file="dom_".$file; 

    my @dom_result=();
    tie @dom_result, 'Tie::File', "$HMMDOM_OUT/$file";

    my $title_line='';
    my @result_row=();

    foreach(@dom_result){

         if($_=~/\#/ and $file_line==0) {
             $title_line=$title_line."\n".$_;                      
             next;
          }
          elsif($_=~/\#/){
            next;
          }elsif($_=~/^\s+$/){
              next;
          }else{         
             my @dom_line=split(' ',$_);             
             push(@result_row,\@dom_line);
             $file_line++;
          }
    }  

    my @sorted_result=sort {$a->[3] cmp $b->[3] ||
                            $a->[6] <=> $b->[6]                            
                           } @result_row;

   open(DOM_TAB,">$HMMDOM_OUT/$file") or die "Cannot open tabular output file\n";

   if($title_line ne ''){
      print DOM_TAB "$title_line\n";
   }

   foreach my $row(@sorted_result){
      foreach my $column(@{$row}){
         print DOM_TAB "$column\t";
      }      
      print DOM_TAB "\n";
  }
   close DOM_TAB;
}

###### READ HMMER DOMAIN OUTPUT FILE ########

sub Read_domain_table {

    my($hmm_program)=(shift);
    my($genome_name)=(shift);
    my($HMMDOM_OUT)=(shift); 
    my(%DomainTable)=();
      
    my $file_line=1;
    my $file="dom_".$hmm_program."_".$genome_name.".out";

    #### Read Tabular HMM Dom Output ###

    my @dom_result=();
    tie @dom_result, 'Tie::File', "$HMMDOM_OUT/$file";   

     my $prev_target_id='';
     my $prev_query_id='';
     my $domain_num=0;

     foreach my $domain_line(@dom_result){

          if($domain_line=~/^\s+/ or $domain_line eq ''){
             next; 
          }

         if($domain_line=~/\#/) {  
             $file_line++; 
             next; 
         }else{         
             my @dom_line=split(" ",$domain_line);

         # target name(0) accession(1) tlen(2) query_name(3) accession(4) qlen(5) E-value(6) score(7) bias(8) #(9) of(10) c-Evalue(11) i-Evalue(12) score(13) bias(14) h_from(15) h_to(16) a_from(17) a_to(18) e_from(19) e_to(20) acc(21) description of target(22)
         #### Read Parameters in scalar variables #### 

             my $query_id=$dom_line[3];
             my $target_id=$dom_line[0];

            if($target_id ne $prev_target_id or $query_id ne $prev_query_id){
               $domain_num=1;
               $prev_target_id=$target_id;
               $prev_query_id=$query_id;
            }
            
            if($hmm_program eq "phmmer"){
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_name}=$dom_line[0];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_accession}=$dom_line[1];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_length}=$dom_line[2]; 
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_name}=$dom_line[3];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_accession}=$dom_line[4];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_length}=$dom_line[5];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_start}=$dom_line[15];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_end}=$dom_line[16];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_start}=$dom_line[17];
                 $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_end}=$dom_line[18];                
             }
             elsif($hmm_program eq "hmmscan"){
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_name}=$dom_line[0];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_accession}=$dom_line[1];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_length}=$dom_line[2]; 
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_name}=$dom_line[3];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_accession}=$dom_line[4];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_length}=$dom_line[5];                
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_start}=$dom_line[17];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_end}=$dom_line[18];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_start}=$dom_line[15];
                $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_end}=$dom_line[16];
             } 
            
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{evalue}=$dom_line[6];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{bit_score}=$dom_line[7];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{bias}=$dom_line[8];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{num_domain}=$dom_line[9];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{total_domain}=$dom_line[10];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{c_evalue}=$dom_line[11];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{i_evalue}=$dom_line[12];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{domain_score}=$dom_line[13];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{domain_bias}=$dom_line[14];                   
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{accuracy}=$dom_line[21];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_description}=$dom_line[22];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{envelope_start}=$dom_line[19];
             $DomainTable{$query_id}->{$target_id}->{$domain_num}->{envelope_end}=$dom_line[20];

             $domain_num++;         
         }
      }            
      return(\%DomainTable);
}

##### Read Full Output for HMMER Programs ####

sub Read_aligned_sequence {

    my($hmm_program)=(shift);
    my($genome_name)=(shift);
    my($HMMFULL_OUT)=(shift);
    my(%DomainTable)=%{(shift)};
     
    my %SequenceAlignment=();

    my $file_line=0;
    my $start_aln=0;
    my $query_seq='';
    my $target_seq='';
    my $sim_line='';
    my $prob_line='';
    my $query_id='';
    my $target_id='';
    my $domain_num=0;

    my $file="full_".$hmm_program."_".$genome_name.".out";

    my @full_result=();
    tie @full_result, 'Tie::File', "$HMMFULL_OUT/$file";   

    my $query_range=undef;
    my $target_range=undef; 

    foreach my $full_line(@full_result){

          $file_line++;

          if($full_line=~/\#/){
             next;
          }
    
          if($full_line=~/(\=\=)(\s+)(domain)(\s+)(\d+)/){

            if($query_seq ne '' and $target_seq ne '' and $sim_line ne ''){

                 if($query_id eq '' or $target_id eq ''){
                     next;
                 }
                 
                 print "$query_id\t$target_id\n";

                  $query_seq=~s/\.*//g;
                  $query_seq=~s/\-*//g;
                  $target_seq=~s/\.*//g;
                  $target_seq=~s/\-*//g;

                  my $hsp_len=length($query_seq);
                  my $query_region=length($query_seq);
                  my $target_region=length($target_seq);

                  my $identical_region=$sim_line;
                     $identical_region=~s/[\s|\+]//g;
                  my $num_identical=length($identical_region);
                  
                  my $similar_region=$sim_line;
                     $similar_region=~s/\s//g;
                  my $num_similar=length($similar_region);

                  my $q_start=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_start};
                  my $q_end=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_end};

                  my $t_start=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_start};
                  my $t_end=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_end};
                  
                  my $new_query_range = Bio::Range->new(-start=>$q_start, -end=>$q_end, -strand=>+1);
                  my $new_target_range= Bio::Range->new(-start=>$t_start, -end=>$t_end, -strand=>+1); 
                                    
                  if($domain_num==1){
                      
                      $query_range=undef;
                      $target_range=undef;
                      $query_range=$new_query_range;
                      $target_range=$new_target_range;

                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;
                  
                  }else{
                     my($q_range,$t_range,$q_region,$t_region,$overlap)=overlap_domain($query_range,$target_range,$new_query_range,$new_target_range,$query_region,$target_region,$hsp_len);

                     if($overlap==0){

                      $query_range=$q_range;
                      $target_range=$t_range;
                     
                      $query_region=$q_region;
                      $target_region=$t_region;
                      $hsp_len=$q_region;

                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;

                    }elsif($overlap==1){
                    
                      $query_range=$q_range;
                      $target_range=$t_range;
                      
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;
                    }
                  }
            }

            $full_line=~/(\=\=)(\s+)(domain)(\s+)(\d+)/;

            $start_aln=1;
            $query_seq='';
            $target_seq='';
            $sim_line='';
            $prob_line='';
            $query_id='';
            $target_id='';
            $domain_num=$5;
            next;    
        }

         if($full_line=~/^(\s+)$/ and $start_aln!=2){
            next;
         }elsif($full_line eq ''){
            next;
         }elsif($full_line=~/(x*)(\.*)(\s+)(RF)$/ or $full_line=~/(x*)(\.*)(\s+)(CS)$/){
            next;
         }

        if($full_line=~/>>/ or $full_line=~/(Internal pipeline statistics summary:)/){
            $start_aln=0;
            next;
        }

        if($start_aln==1){
            my $aln_line=$full_line;
            $aln_line=~s/^(\s+)//g;
 
            my @aln_line=split(' ',$aln_line);           

               if($hmm_program eq "phmmer"){
                   $query_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $query_seq=$query_seq.$seq;                
               }
               elsif($hmm_program eq "hmmscan"){
                   $target_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $target_seq=$target_seq.$seq;
               }
            $start_aln=2;
            next;
         }

         if($start_aln==2){
            my $aln_line=$full_line;
            $aln_line=~s/^(\s+)//g;                    
            $sim_line=$sim_line.$aln_line;
            $start_aln=3;
            next;
          }

          if($start_aln==3){
            my $aln_line=$full_line;         
            $aln_line=~s/^(\s+)//g;
             
            my @aln_line=split(' ',$aln_line);            

              if($hmm_program eq "phmmer"){
                 $target_id=shift(@aln_line);
                 my $seq=$aln_line[1];            
                 $target_seq=$target_seq.$seq;                   
              }
              elsif($hmm_program eq "hmmscan"){
                 $query_id=shift(@aln_line);
                 my $seq=$aln_line[1];            
                 $query_seq=$query_seq.$seq;
              }
            $start_aln=4;
            next;
          }
 
          if($start_aln==4){
            my $aln_line=$full_line;           
            $aln_line=~s/^(\s+)//g;                    
            $prob_line=$prob_line.$aln_line;
            $start_aln=1;
            next;
          }
    }

    ###### IF last hit in full alignment file ##########
    if(scalar(@full_result) eq $file_line){

         if($query_seq ne '' and $target_seq ne '' and $sim_line ne ''){

                 $query_seq=~s/\.*//g;
                 $query_seq=~s/\-*//g;
                 $target_seq=~s/\.*//g;
                 $target_seq=~s/\-*//g;

                 my $identical_region=$sim_line;
                    $identical_region=~s/[\s|\+]//g;                 
                 my $num_identical=length($identical_region); 
             
                 my $similar_region=$sim_line;
                    $similar_region=~s/\s//g;
                 my $num_similar=length($similar_region);    

                 my $hsp_len=length($query_seq);
                 my $query_region=length($query_seq);
                 my $target_region=length($target_seq);                     

                  my $q_start=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_start};
                  my $q_end=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{query_end};

                  my $t_start=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_start};
                  my $t_end=$DomainTable{$query_id}->{$target_id}->{$domain_num}->{target_end};

                  my $new_query_range = Bio::Range->new(-start=>$q_start, -end=>$q_end, -strand=>+1);
                  my $new_target_range= Bio::Range->new(-start=>$t_start, -end=>$t_end, -strand=>+1); 
           
                  
                  if($domain_num==1){
                      
                      $query_range=undef;
                      $target_range=undef;
                      $query_range=$new_query_range;
                      $target_range=$new_target_range;

                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;
                  
                  }else{
                     my($q_range,$t_range,$q_region,$t_region,$overlap)=overlap_domain($query_range,$target_range,$new_query_range,$new_target_range,$query_region,$target_region,$hsp_len);

                     if($overlap==0){

                      $query_range=$q_range;
                      $target_range=$t_range;
                     
                      $query_region=$q_region;
                      $target_region=$t_region;
                      $hsp_len=$q_region;

                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;

                    }elsif($overlap==1){
                    
                      $query_range=$q_range;
                      $target_range=$t_range;
                      
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_region}=$query_region; 
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_region}=$target_region;
                      $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{hsp_len}=$hsp_len;                    
                    }
                }          
           }
    }  
    return(\%SequenceAlignment);
}

##### Check overlapping domains ###

sub overlap_domain{

    my($query_range)=(shift);
    my($target_range)=(shift);
    my($new_query_range)=(shift);
    my($new_target_range)=(shift);
    my $query_region=(shift);
    my $target_region=(shift);
    my($hsp_len)=(shift);
    my $overlap=-1;

    my $overlap_query_region=0;
    my $overlap_target_region=0;
    my $union_query_region=1;
    my $union_target_region=1;
    
    my $percent_query_overlap=0.0;
    my $percent_target_overlap=0.0;

    ####### parsing overlapping domain regions in target and query sequences #####
    my @overlap_query_range=$query_range->intersection($new_query_range); 
    my @overlap_target_range=$target_range->intersection($new_target_range);

    my @union_query_range=$query_range->union($new_query_range); 
    my @union_target_range=$target_range->union($new_target_range);

    #### if next query domain overlaps with previous domain region ####
     if(scalar(@overlap_query_range)>0){
                      
          $union_query_region=($union_query_range[1] - $union_query_range[0]) + 1;
          $overlap_query_region=($overlap_query_range[1]-$overlap_query_range[0]) + 1;
                         
          $union_target_region=($union_target_range[1] - $union_target_range[0]) + 1;

          if(scalar(@overlap_target_range)>0){
                $overlap_target_region=($overlap_target_range[1]-$overlap_target_range[0]) + 1;
          } 
                                        
          $percent_query_overlap=($overlap_query_region/$hsp_len)*100;                         
          $percent_target_overlap=($overlap_target_region/$hsp_len)*100;

          if($percent_query_overlap <25 or $percent_target_overlap <25 ){
               $query_region=$query_region-$overlap_query_region;     
               $query_range=$new_query_range;  
               
               $target_region=$target_region-$overlap_target_region;     
               $target_range=$new_target_range;
               
               $overlap=0;                                     
          }                                         
   }else{
        $overlap=1;
        $query_range=$new_query_range; 
        $target_range=$new_target_range;
   }

  return($query_range,$target_range,$query_region,$target_region,$overlap);
}

