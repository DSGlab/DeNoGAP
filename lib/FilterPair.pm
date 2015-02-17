##### Module to Filter BEST/PARTIAL/CHIMERIC hits for Query sequences based on threshold #####
##### Author: Shalabh Thakur ################
##### Date: 16-Aug-2013 #####################

#!/usr/bin/perl
package FilterPair;
use strict;
use Exporter;
use File::Basename;
use Bio::Range;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(getHitForQuery);

sub getHitForQuery {

    my($genome_name)=(shift);
    my(%PARSE_HMMER)=%{(shift)};
    my($best_output_dir)=(shift);
    my($all_output_dir)=(shift);
    my($chimeric_output_dir)=(shift);
    my(%DomainTable)=%{(shift)};
    my(%SequenceAlignment)=%{(shift)};
    
    my @prev_hsp=();
    my @domain=();    
    my @query_coordinate=();
    my @subject_coordinate=();
    my @start_coord=();
    my @end_coord=();
    my $counter=1;    
    my %PairRelation=();
    my %MatchedQuery=();
    my %QueryEvalue=();    
    my @chimera=();
    my %chimeratargetperquery=();
    my %chimeratarget=();
    my %chimeraquery=();
    my $totalidentity=0.0;
    my $totalsimilarity=0.0;
    my $total_alignment_length=0;
    my $file_size=keys(%DomainTable);

    my $best_output_file=$best_output_dir."/besthit_$genome_name.txt";
    my $all_output_file=$all_output_dir."/allhit_$genome_name.txt";
    my $chimeric_output_file=$chimeric_output_dir."/chimeric_$genome_name.txt";
    my $raw_output=$all_output_dir."/RAW_$genome_name.txt";
    
    my %hmm_homolog=();   ### Homolog to HMM / SINGLETON GROUP ###
    my %hmm_homolog_significance=();
    my %self_homolog=();  ### Homolog to other gene in same taxa ####
    my %self_homolog_significance=();

    open(RAW_OUT,">$raw_output");

    ##### Loop through query - target pairs #######

    foreach my $query_id(keys %DomainTable){  #### QUERY LOOP ####

          my %target=%{$DomainTable{$query_id}};

          foreach my $target_id(keys %target){  ##### TARGET LOOP #####

              if(!defined($SequenceAlignment{$query_id}->{$target_id})){
                 print "ERROR: Undefined value in hash SequenceAlignment for $query_id and $target_id\n";
                 next;
               }
    
               my %domain_align_feature=%{$SequenceAlignment{$query_id}->{$target_id}};
 
               my %domain_feature=%{$target{$target_id}};

               my $aligned_query_region=0;
               my $aligned_target_region=0;

               my @similarity=();
               my $hit_significance='';

               ###### Calculate total identical, similar positions in alignment #######
               ###### Calculate total sequence region covered in query and target #####

               foreach my $domain_num(keys %domain_align_feature){  
            
                    my %align_feature=%{$domain_align_feature{$domain_num}};               
                       
                    my $identical_position=$align_feature{identical};
                    my $similar_position=$align_feature{similar}; 
                    my $alignment_len=$align_feature{hsp_len};
                    my $query_hsp_len=$align_feature{query_region};
                    my $target_hsp_len=$align_feature{target_region};
                
                    $aligned_query_region=$aligned_query_region + $query_hsp_len;
                    $aligned_target_region=$aligned_target_region + $target_hsp_len;
               }

               ######### Loop through each domain feature of query - target #######
              
               foreach my $domain_num(keys %domain_align_feature){

                      if(!defined($domain_feature{$domain_num})){
                       print "ERROR: Undefined value in hash domain_feature for $query_id -> $target_id -> $domain_num\n";
                       next;
                      } 
                      my %feature=%{$domain_feature{$domain_num}};

                      my %align_feature=%{$domain_align_feature{$domain_num}};

                      my $query_id=$feature{query_name};
                      my $target_id=$feature{target_name};
                      my $accuracy=$feature{accuracy};
                      my $seq_evalue=$feature{evalue};
                      my $c_evalue=$feature{c_evalue};
                      my $i_evalue=$feature{i_evalue};
                      my $bit_score=$feature{bit_score};
                      my $query_len=$feature{query_length};
                      my $target_len=$feature{target_length}; 
                      my $qstart=$feature{query_start};
                      my $qend=$feature{query_end};
                      my $tstart=$feature{target_start};
                      my $tend=$feature{target_end};
                      my $num_domain=$feature{num_domain};
                      my $total_domain=$feature{total_domain};
                      my $target_description=$feature{target_description};

                      my $identical_position=$align_feature{identical};
                      my $similar_position=$align_feature{similar}; 
                      my $alignment_len=($qend-$qstart) + 1;
                      my $query_domain_len=($qend-$qstart) + 1;
                      my $target_domain_len=($tend-$tstart) + 1;

                      if($identical_position==0 and $similar_position==0 and $alignment_len==0){
                          next;
                      }

                      my $relation='';
                      $hit_significance=$seq_evalue;
                      
                      my $q_coverage=sprintf "%.2f", ($aligned_query_region/$query_len)*100;
                      my $t_coverage=sprintf "%.2f", ($aligned_target_region/$target_len)*100;
 
                      my $percent_identity=sprintf "%.2f", ($identical_position/$alignment_len)*100;
                      my $percent_similarity=sprintf "%.2f", ($similar_position/$alignment_len)*100;

                      #my $domain_identity=sprintf "%.2f", ($identical_position/$alignment_len)*100;
                      #my $domain_similarity=sprintf "%.2f", ($similar_position/$alignment_len)*100;
   
                      if($query_id ne $target_id){
                         $relation=find_QH_Relation(\%PARSE_HMMER,$query_id,$target_id,$seq_evalue,$qstart,$qend,$tstart,$tend,$total_domain,$query_len,$target_len,$percent_identity,$percent_similarity,$q_coverage,$t_coverage,$accuracy);
                      }else{
                         $relation="SELF_MATCH";
                      }

                      my $similarity_line="$query_id\t$target_id\t$query_len\t$target_len\t$total_domain\t$num_domain\t$qstart\t$qend\t$tstart\t$tend\t$seq_evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$q_coverage\t$t_coverage\t$relation\t$target_description";                         
                       
                      my $print_similarity_line="$query_id\t$target_id\t$query_len\t$target_len\t$total_domain\t$num_domain\t$identical_position\t$similar_position\t$alignment_len\t$qstart\t$qend\t$tstart\t$tend\t$seq_evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$q_coverage\t$t_coverage\t$relation\t$target_description";                         

                      print RAW_OUT "$print_similarity_line\n";                      
                      push(@similarity,$similarity_line);                                            
                 } 

                  if($hit_significance eq ''){
                    next;
                  }

                  if($target_id=~/(Group)(\d+)/){    
                     $hmm_homolog{$query_id}->{$target_id}=\@similarity;
                     $hmm_homolog_significance{$query_id}->{$target_id}=$hit_significance;
                  }else{
                     $self_homolog{$query_id}->{$target_id}=\@similarity;
                     $self_homolog_significance{$query_id}->{$target_id}=$hit_significance;
                  }     
          }
    }

    close RAW_OUT;

    PrintOutPairSimilairty(\%hmm_homolog,\%self_homolog,\%hmm_homolog_significance,\%self_homolog_significance,$all_output_file,$best_output_file,$chimeric_output_file);
}


##### Finds Relation between Query and Hit Sequence #######

sub find_QH_Relation {
    
    my(%PARSE_HMMER)=%{(shift)}; 
    my($query_id)=(shift);
    my($subject_id)=(shift);
    my($seq_evalue)=(shift);
    my($qstart)=(shift);
    my($qend)=(shift);
    my($tstart)=(shift);
    my($tend)=(shift);
    my($total_domain)=(shift);
    my($query_len)=(shift);
    my($subject_len)=(shift);
    my($percent_identity)=(shift);
    my($percent_similarity)=(shift);
    my($query_coverage)=(shift);
    my($subject_coverage)=(shift);
    my($accuracy)=(shift);
    
    my $relation='';

    ###### Check for Query-Subject Relation ######
    if($query_coverage>=$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{HMM_COVERAGE} and (($percent_identity>=$PARSE_HMMER{IDENTITY}) and ($percent_similarity>=$PARSE_HMMER{SIMILARITY}))){ 
          $relation="BEST";
    }elsif((($query_coverage>=$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage<$PARSE_HMMER{HMM_COVERAGE}) or ($query_coverage<$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{HMM_COVERAGE})) and (($percent_identity>=$PARSE_HMMER{IDENTITY}) and ($percent_similarity>=$PARSE_HMMER{SIMILARITY}))){
       $relation="TRUNCATED";    
    }elsif($total_domain==1 and ($query_coverage<$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage<$PARSE_HMMER{HMM_COVERAGE}) and ($query_coverage>=$PARSE_HMMER{MIN_CHIMERA_QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{MIN_CHIMERA_HMM_COVERAGE}) and ($percent_identity>=$PARSE_HMMER{MIN_CHIMERA_IDENTITY} or $percent_similarity>=$PARSE_HMMER{MIN_CHIMERA_SIMILARITY})){ 
       ##### Check for Chimera #### 
       my $isChimera=checkForChimericGene($query_id,$subject_id,$seq_evalue,$qstart,$qend,$tstart,$tend,$query_len,$subject_len);      
       if($isChimera and $accuracy>=$PARSE_HMMER{CHIMERA_ACCURACY}){            
           $relation=$isChimera;        
       }    
    }else{
      $relation="INSIGNIFICANT";
    }    
    return($relation); 
}

##### Check for Chimera Gene ########

sub checkForChimericGene {

    my($query_id)=(shift);
    my($subject_id)=(shift);
    my($evalue)=(shift);   
    my($qstart)=(shift);
    my($qend)=(shift);
    my($tstart)=(shift);
    my($tend)=(shift);
    my($query_len)=(shift);
    my($subject_len)=(shift);
    my($total_domain)=(shift);  
    my $isChimera='INSIGNIFICANT';
       
    if(($qstart<=30) and ($qend<$query_len)){
        if(($tstart<=30) and ($tend<$subject_len)){
            $isChimera="N-Chimera";                 
        }
    }elsif(($qstart>30) and  ((($query_len-$qend)+1)<=15)){
        if(($tstart>30) and ((($subject_len-$tend)+1)<=15)){
            $isChimera="C-Chimera";               
        }
    }      
   return($isChimera);    
}

##### Print pairwise similarity information ######

sub PrintOutPairSimilairty {

    my(%hmm_homolog)=%{(shift)};
    my(%self_homolog)=%{(shift)};
    my %hmm_homolog_significance=%{(shift)};
    my %self_homolog_significance=%{(shift)};
    my($all_output_file)=(shift);
    my($best_output_file)=(shift);
    my($chimeric_output_file)=(shift);

    open(BESTHITFILE,">$best_output_file");
    open(ALLHITFILE,">$all_output_file");
    open(CHIMERAFILE,">$chimeric_output_file");

    my %best_match=();
    ##### GROUP HOMOLOG SEARCH ######
    foreach my $query_id(keys %hmm_homolog_significance){

            my %target_group=%{$hmm_homolog_significance{$query_id}}; 
 
            my $hit_count=1;
            
            foreach my $target_id(sort{$target_group{$a}<=>$target_group{$b}} keys %target_group){

                 my @similarity=@{$hmm_homolog{$query_id}->{$target_id}};
                 
                 foreach my $pair_similarity(@similarity){
                    #### check best match ####
                    if($hit_count==1){ 
                         ###### If top hit is best match #####
                         if($pair_similarity=~/(BEST)/){                                   
                             print ALLHITFILE "$pair_similarity\n";
                             print BESTHITFILE "$pair_similarity\n";  
                             $best_match{$query_id}="best"; 
                             $hit_count++;
                         }elsif($pair_similarity=~/(TRUNCATED)/){                      
                             $best_match{$query_id}="truncated";     
                             ##### Do not print here, wait till get the best match ###### 
                              $hit_count++;                                                   
                         }elsif($pair_similarity=~/(Chimera)/){ 
                              ### if top hit is chimera match ####
                             print ALLHITFILE "$pair_similarity\n";
                             print CHIMERAFILE "$pair_similarity\n";
                             $best_match{$query_id}="chimera"; 
                             $hit_count++;  
                         }
                        
                    }elsif($hit_count>1){
                         if(($best_match{$query_id} eq "truncated")  and ($pair_similarity=~/(BEST)/)){
                             #### if best match found for truncated query id print it here ###
                             print ALLHITFILE "$pair_similarity\n";
                             print BESTHITFILE "$pair_similarity\n";  
                             $best_match{$query_id}="best"; 
                             $hit_count++;
                                                                             
                         }elsif(($best_match{$query_id} eq "best")  and ($pair_similarity=~/(BEST)/)){
                             ### Print other best hit ####
                             print ALLHITFILE "$pair_similarity\n";
                             $best_match{$query_id}="best"; 
                             $hit_count++;
                             
                         }elsif(($best_match{$query_id} eq "chimera") and ($pair_similarity=~/(Chimera)/)){
                             #### print other chimera hit if previous hits were also chimera ####
                             print ALLHITFILE "$pair_similarity\n";
                             print CHIMERAFILE "$pair_similarity\n";
                             $hit_count++;
                         }    
                                
                    } 
                 }
            }   
    } 
       
    ####### PRINT SELF MATCH FOR TRUNCATED GROUP AND CHIMERA ######
    
    foreach my $query_id(keys %best_match){
    
         if($best_match{$query_id} eq "truncated" or $best_match{$query_id} eq "chimera" ){
         
             ##### Truncated query id self match ####
             print "get self match for truncated query $query_id\n";
             my @similarity=@{$self_homolog{$query_id}->{$query_id}}; 
             
             foreach my $pair_similarity(@similarity){
                  print BESTHITFILE "$pair_similarity\n"; 
             }
          }

         #### Print all truncated matches for query_id #####           
         if($best_match{$query_id} eq "truncated"){  
             ### Truncated query id partial matches ######
              my %target_group=%{$hmm_homolog_significance{$query_id}}; 
              foreach my $target_id(sort{$target_group{$a}<=>$target_group{$b}} keys %target_group){

                 my @similarity=@{$hmm_homolog{$query_id}->{$target_id}};
                 
                 foreach my $pair_similarity(@similarity){               
                      if($pair_similarity=~/(TRUNCATED)/){  
                          print ALLHITFILE "$pair_similarity\n";
                      }         
                 }
                 
              }   
         }
    }
    
    ####### SELF HOMOLOG SEARCH FOR QUERY ID WITH INSIGNIFICANT MATCH or NO MATCH  #######
    my %self_homolog_match=();
    
    foreach my $query_id(keys %self_homolog_significance){
    
            my %target_group=%{$self_homolog_significance{$query_id}}; 
            
            foreach my $target_id(sort{$target_group{$a}<=>$target_group{$b}} keys %target_group){
            
                my @similarity=@{$self_homolog{$query_id}->{$target_id}};
                     
                if(!defined($best_match{$query_id}) and !defined($best_match{$target_id})){
                    foreach my $pair_similarity(@similarity){
                       
                         print "$query_id\t $target_id\n";

                         $self_homolog_match{$query_id}="";

                         if(($pair_similarity=~/(BEST)/ or $pair_similarity=~/(SELF_MATCH)/)){   

                             if($self_homolog{$query_id}){
                                if($self_homolog_match{$query_id} eq "chimera"){
                                    next;
                                }
                             }                              
                             print ALLHITFILE "$pair_similarity\n";
                             print BESTHITFILE "$pair_similarity\n"; 
                             
                             if($query_id ne $target_id){
                                $self_homolog_match{$query_id}="best";
                             }
                                                   
                         }elsif($pair_similarity=~/(TRUNCATED)/){

                             
                             if($self_homolog{$query_id}){
                                if($self_homolog_match{$query_id} eq "chimera"){
                                    next;
                                }
                             }  
                         
                             print ALLHITFILE "$pair_similarity\n";
                             $self_homolog_match{$query_id}="truncated";
                             
                         }elsif($pair_similarity=~/(Chimera)/){
                            
                             if($self_homolog{$query_id}){
                                if($self_homolog_match{$query_id} ne "chimera" and $self_homolog_match{$query_id} ne ''){
                                    next;
                                }
                             } 
                            print ALLHITFILE "$pair_similarity\n"; 
                            print CHIMERAFILE "$pair_similarity\n"; 
                            $self_homolog_match{$query_id}="chimera";
                         }   
                    } 
                }   
            }
     }

     close BESTHITFILE;
     close ALLHITFILE;
     close CHIMERAFILE; 
}



##############################################################################################################
