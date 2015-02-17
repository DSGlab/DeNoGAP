##### Module to Perform Gene Predic:tion #####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################
#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use Env;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Range;
use File::Basename;

 my($genome_name)=(shift);
 my($full_name)=(shift);
 my($genome_file)=(shift);
 my($project_dir)=(shift);
 my($glimmer_dir)=(shift);
 my($genemark_dir)=(shift);
 my($prodigal_dir)=(shift);
 my($fragscan_dir)=(shift);
 my($genome_dir)=(shift);
 my($cds_dir)=(shift);
 my($protein_dir)=(shift);
 my($feature_dir)=(shift);
 my($translation)=(shift);
 my($overlap_base)=(shift);
 my($glimmer3_opt)=(shift);
 my($long_orf_opt)=(shift);
 my($multi_extract_opt)=(shift);
 my($build_icm_opt)=(shift);
 my($genemark_opt)=(shift);
 my($prodigal_opt)=(shift);
 my($fragscan_opt)=(shift);
 
    $full_name=~s/\#/ /g;
    
    run_glimmer($genome_name,$genome_file,$glimmer_dir,$glimmer3_opt,$long_orf_opt,$multi_extract_opt,$build_icm_opt);
    run_genemark($genome_name,$genome_file,$genemark_dir,$genemark_opt);
    run_prodigal($genome_name,$genome_file,$prodigal_dir,$prodigal_opt);
    run_fragscan($genome_name,$genome_file,$fragscan_dir,$fragscan_opt);

    read_gene_prediction($genome_name,$full_name,$genome_file,$project_dir,$glimmer_dir,$genemark_dir,$prodigal_dir,$fragscan_dir,$cds_dir,$protein_dir,$feature_dir,$translation,$overlap_base);
    
    unlink("$project_dir/all_$genome_name.txt"); 

 
######## GLIMMER #####

sub run_glimmer{

    my($genome_name)=(shift);
    my($genome_file)=(shift);
    my($glimmer_dir)=(shift);
    my($glimmer3_opt)=(shift);
    my($long_orf_opt)=(shift);
    my($multi_extract_opt)=(shift);
    my($build_icm_opt)=(shift);

    print "Running Glimmer\n";

    my $out_dir="$glimmer_dir/$genome_name";
    mkdir($out_dir);

    print "$genome_file\n";

    #### Glimmer Path ####

    #### Set Glimmer options ###

    my %glimmer_predict=();
       
    my $seqio_obj = Bio::SeqIO->new(-file => "$genome_file", -format => "fasta");

    print "Finding long orfs for training\n";

    open(LONG_ORF,">$out_dir/$genome_name.longorfs");

       while(my $seq_obj= $seqio_obj->next_seq){

               my $seq=$seq_obj->seq;
               my $seq_id=$seq_obj->display_id;

               open(TMP_GENOME,">$out_dir/tmp_genome.fasta");
               print TMP_GENOME ">$seq_id\n$seq";
               close TMP_GENOME;

               #### Find long, non-overlapping orfs to use as training set ###
               my $status1=system("long-orfs $long_orf_opt $out_dir/tmp_genome.fasta $out_dir/tmp_genome.longorfs");
               if($status1!=0){
                 print "Failed to extract long-orfs from $seq_id\n"; next;
               }

               open(TMP_LONGORF,"$out_dir/tmp_genome.longorfs");
               my @long_orf=<TMP_LONGORF>;
               close TMP_LONGORF;               

               foreach my $line(@long_orf){
                    my @line=split(" ",$line);
                    print LONG_ORF "$line[0]\t$seq_id\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
               }                             
       }
 
     #### Extract Sequence for long-orfs to use as training set ####
     print "Extract long-orfs for training\n";
     my $status2=system("multi-extract $multi_extract_opt $genome_file $out_dir/$genome_name.longorfs > $out_dir/$genome_name.train");

     if($status2!=0){
       print "Failed to extract training sequence for $genome_name\n"; exit;
     }
     #### Build ICM from the training Sequence ####
     print "Building ICM\n"; 
     my $status3=system("build-icm $build_icm_opt $out_dir/$genome_name.icm < $out_dir/$genome_name.train"); 
   
     if($status3!=0){

       print "Failed to build ICM for $genome_name\n"; exit;

     }
     #### Run Glimmer ####
     print "Running Glimmer\n";
     my $status4=system("glimmer3 $glimmer3_opt $genome_file $out_dir/$genome_name.icm $out_dir/$genome_name");    

     if($status4!=0){
       print "Failed to run Glimmer for $genome_name\n"; exit;
     }

     open(GLIMMER_PREDICT,"$out_dir/$genome_name.predict");
     my @glimmer_predict=<GLIMMER_PREDICT>;
     close GLIMMER_PREDICT;

     open(GLIMMER_PREDICT1,">$out_dir/$genome_name.1.predict");
  
     my $header='';

     foreach my $line(@glimmer_predict){
           if($line=~/>/){
              $line=~s/>//g;             
              $header=$line;
              chomp($header);
              my @header=split(" ",$header); 
              $header=shift(@header);       
           }else{
              my @line=split(" ",$line);
              print GLIMMER_PREDICT1 "$line[0]\t$header\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
           }               
     }
     close GLIMMER_PREDICT1;

        my $status5=system("multi-extract $multi_extract_opt $genome_file $out_dir/$genome_name.1.predict > $out_dir/$genome_name.fna");
     
     if($status5!=0){
         print "Failed to extract gene sequence for $genome_name\n"; exit;
     }   
}


###### GeneMarkS #######

sub run_genemark {

    my($genome_name)=(shift);
    my($genome_file)=(shift);
    my($genemark_dir)=(shift);
    my($genemark_opt)=(shift);
  
    print "Running GeneMark\n";

    my $out_dir="$genemark_dir/$genome_name";

    mkdir($out_dir);


    ### temporary directory created ####
    mkdir("$genome_name");

    mkdir("$genome_name/genemark_predict");

    my %genemark_predict=();

    ##### GeneMark Path #####
     
       if($genome_file=~/^\.+$/){next;}

       print "GENOME: $genome_name\n";    
   
       chdir("$genome_name/genemark_predict");

       system("gmsn.pl --name $genome_name --species $genome_name $genemark_opt $Bin/$genome_file");

       chdir("$Bin");
   
       system("mv $genome_name/genemark_predict/* $out_dir"); 

       system("rm -r $genome_name"); 
  
    return(\%genemark_predict); 
}

##### PRODIGAL #####

sub run_prodigal {

    my($genome_name)=(shift);
    my($genome_file)=(shift);
    my($prodigal_dir)=(shift);
    my($prodigal_opt)=(shift);

    print "Running Prodigal\n";

    my $out_dir="$prodigal_dir/$genome_name";
    mkdir($out_dir); 

    my %prodigal_predict=();

    #### Prodigal Path ####      
    if($genome_file=~/^\.+$/){next;}

    print "GENOME:$genome_name\n";

    system("prodigal -g 11 -m -t $out_dir/$genome_name.training -i $genome_file");

    system("prodigal $prodigal_opt -a $out_dir/$genome_name.faa -d $out_dir/$genome_name.fnn -o $out_dir/$genome_name.gff -s $out_dir/$genome_name.score  -t $out_dir/$genome_name.training -i $genome_file");       
}

##### FRAGSCAN ####

sub run_fragscan{

    my($genome_name)=(shift);
    my($genome_file)=(shift);
    my($fragscan_dir)=(shift);
    my($fragscan_opt)=(shift);
    
    print "Running FragScan\n";

    my $out_dir="$fragscan_dir/$genome_name";

    mkdir($out_dir); 

    my %fragscan_predict=();

    #### Fragscan Path ####
      
    if($genome_file=~/^\.+$/){next;}
 
    print "GENOME: $genome_name\n";
       
    system("run_FragGeneScan.pl -genome=$genome_file -out=$out_dir/$genome_name $fragscan_opt");      
}

#### Compare Gene Prediction to find best gene sequence ####

sub read_gene_prediction {

    my($genome_name)=(shift);
    my($full_genome_name)=(shift);
    my($genome_file)=(shift);
    my($project_dir)=(shift);
    my($glimmer_dir)=(shift);
    my($genemark_dir)=(shift);
    my($prodigal_dir)=(shift);
    my($fragscan_dir)=(shift);
    my($cds_dir)=(shift);
    my($protein_dir)=(shift);
    my($feature_dir)=(shift);  
    my($translation_table)=(shift);  
    my($overlap_base)=(shift);
   
    my %glimmer_predict=();
    my %genemark_predict=();
    my %prodigal_predict=();
    my %fragscan_predict=();    
    my %sequence=();  
    
    if($genome_file=~/^\.+$/){next;}

    my $genome_file_basename=basename($genome_file);

    my $seqio_obj = Bio::SeqIO->new(-file => "$genome_file", -format => "fasta");
    
    while(my $seq_obj= $seqio_obj->next_seq){
        my $seq=$seq_obj->seq;
        my $seq_id=$seq_obj->display_id;               
        $sequence{$genome_name}->{$seq_id}=$seq;            
    }           
   
        ###### READ GLIMMER OUTPUT ####
        open(GLIMMER,"$glimmer_dir/$genome_name/$genome_name.1.predict");
        my @glimmer_predict=<GLIMMER>;          
        $glimmer_predict{$genome_name}=\@glimmer_predict;
        close GLIMMER;

        ###### READ GENEMARK OUTPUT #####
        open(GENEMARK,"$genemark_dir/$genome_name/$genome_file_basename.gff");
        my @genemark_predict=<GENEMARK>; 
        $genemark_predict{$genome_name}=\@genemark_predict;
        close GENEMARK;

        ##### READ PRODIGAL OUTPUT #####
        open(PRODIGAL,"$prodigal_dir/$genome_name/$genome_name.gff");
        my @prodigal_predict=<PRODIGAL>;
        $prodigal_predict{$genome_name}=\@prodigal_predict;
        close PRODIGAL; 

        ##### READ FRAGSCAN OUTPUT ######
        open(FRAGSCAN,"$fragscan_dir/$genome_name/$genome_name.out");
        my @fragscan_predict=<FRAGSCAN>;
        $fragscan_predict{$genome_name}=\@fragscan_predict;
        close FRAGSCAN;   

        ##### READ PREDICTED CDS SEQUENCE ########
        my %predicted_sequence=();

        ##### GLIMMER SEQ ######
        my $glimmer_seqobj = Bio::SeqIO->new(-file => "$glimmer_dir/$genome_name/$genome_name.fna", -format => "fasta");
        my %glimmer_seq=();

        while(my $glimmer_seqobj= $glimmer_seqobj->next_seq){

             my $seq_id=$glimmer_seqobj->display_id;
             my $desc=$glimmer_seqobj->desc;
             my $seq=$glimmer_seqobj->seq;
   
             my @desc=split(" ",$desc);
             my $contig=shift(@desc); 
             $glimmer_seq{$seq_id."_".$contig}=$seq;         
        }

        ###### GeneMark Seq ######
        my $genemark_seqobj = Bio::SeqIO->new(-file => "$genemark_dir/$genome_name/$genome_file_basename.fnn", -format => "fasta");
        my %genemark_seq=();
  
        while(my $genemark_seqobj= $genemark_seqobj->next_seq){
             my $seq_id=$genemark_seqobj->display_id;             
             my $seq=$genemark_seqobj->seq;
             $genemark_seq{$seq_id}=$seq;         
        }

        ####### Prodigal Seq ####

        my $prodigal_seqobj = Bio::SeqIO->new(-file => "$prodigal_dir/$genome_name/$genome_name.fnn", -format => "fasta");
        my %prodigal_seq=();
  
        while(my $prodigal_seqobj= $prodigal_seqobj->next_seq){
             my $seq_id=$prodigal_seqobj->display_id;             
             my $seq=$prodigal_seqobj->seq;
             $prodigal_seq{$seq_id}=$seq;         
        }
 
        ####### FragScan Seq #####
 
        my $fragscan_seqobj = Bio::SeqIO->new(-file => "$fragscan_dir/$genome_name/$genome_name.ffn", -format => "fasta");
        my %fragscan_seq=();
  
        while(my $fragscan_seqobj= $fragscan_seqobj->next_seq){
             my $seq_id=$fragscan_seqobj->display_id;             
             my $seq=$fragscan_seqobj->seq;
             $fragscan_seq{$seq_id}=$seq;         
        }

        $predicted_sequence{glimmer}=\%glimmer_seq;
        $predicted_sequence{genemark}=\%genemark_seq;
        $predicted_sequence{prodigal}=\%prodigal_seq;
        $predicted_sequence{fragscan}=\%fragscan_seq;

        my %glimmer_gene=();
        my %genemark_gene=();
        my %prodigal_gene=();
        my %fragscan_gene=(); 

        my $combine_prediction="$project_dir/all_$genome_name.txt";

        open(ALL_GENE_PRED,">$combine_prediction");       
         
        ###### Parse Result from Glimmer ######     
        if(%glimmer_predict){

            my @glimmer_predict=@{$glimmer_predict{$genome_name}};  
     
            my $gene_id_count="00001";         
              
            foreach my $line(@glimmer_predict){

                 my @glimmer_column=split(" ",$line);

                 my $gene_id=$genome_name."_".$gene_id_count;

                 $glimmer_gene{$gene_id}->{gene_id}=$gene_id;
                 $glimmer_gene{$gene_id}->{tool_id}=$glimmer_column[0]."_".$glimmer_column[1];                
                 $glimmer_gene{$gene_id}->{contig_id}=$glimmer_column[1];
                 $glimmer_gene{$gene_id}->{start}=$glimmer_column[2];
                 $glimmer_gene{$gene_id}->{end}=$glimmer_column[3];
                 
                 my $coordinate_length=0;
                 
                 if($glimmer_column[4]=~/\-/){
                    $glimmer_gene{$gene_id}->{strand}="-";                    
                    if($glimmer_gene{$gene_id}->{start}< $glimmer_gene{$gene_id}->{end}){
                     next;
                    }else{
                      $coordinate_length=($glimmer_gene{$gene_id}->{start} - $glimmer_gene{$gene_id}->{end}) + 1;
                    } 
                 }else{
                    $glimmer_gene{$gene_id}->{strand}="+";
                    if($glimmer_gene{$gene_id}->{start} > $glimmer_gene{$gene_id}->{end}){
                       next;
                    }else{
                      $coordinate_length=($glimmer_gene{$gene_id}->{end} - $glimmer_gene{$gene_id}->{start}) + 1;
                    } 
                 }                 
                 $glimmer_column[4]=~/(\d)/;
                 $glimmer_gene{$gene_id}->{frame}=$1;               
                 $glimmer_gene{$gene_id}->{tool}="glimmer";

                 my $tool_id=$glimmer_gene{$gene_id}->{tool_id};
                 my $seq=$predicted_sequence{glimmer}->{$tool_id};
                 
                 if(!defined($seq)){
                   next;
                 }

                 my $predicted_seq_len=length($seq);

                 if($predicted_seq_len ne $coordinate_length){
                     next;
                 }                
               
                 my $start_gene=$glimmer_gene{$gene_id}->{start};
                 my $end_gene=$glimmer_gene{$gene_id}->{end};
                 
                 if($glimmer_gene{$gene_id}->{strand} eq "-"){
                     $start_gene=$glimmer_gene{$gene_id}->{end};
                     $end_gene=$glimmer_gene{$gene_id}->{start};
                 }
                 
                 print ALL_GENE_PRED $glimmer_gene{$gene_id}->{gene_id},"\t",                       
                       $glimmer_gene{$gene_id}->{contig_id},"\t",
                       $start_gene,"\t",
                       $end_gene,"\t",
                       $glimmer_gene{$gene_id}->{strand},"\t",
                       $glimmer_gene{$gene_id}->{frame},"\t",
                       $glimmer_gene{$gene_id}->{tool},"\t",
                       $glimmer_gene{$gene_id}->{tool_id},"\n";   

                 $gene_id_count++;                         
            } 
        }
        
        ##### Parse Result from GeneMark ######
        if(%genemark_predict){
           
            my @genemark_predict=@{$genemark_predict{$genome_name}};

            my $gene_id_count="00001";
            
            foreach my $line(@genemark_predict){

                   if($line=~/\#/){
                      next; 
                   }elsif($line=~/^\s+$/){
                      next;
                   }
                   my @genemark_column=split("\t",$line);
                   
                   my $coordinate_length=0;

                   my $gene_id=$genome_name."_".$gene_id_count;
                   my @header=split(" ",$genemark_column[0]);
                   my $contig_id=shift(@header);
                   my $tool_id=$genemark_column[8];
                   chomp($tool_id);
                   $tool_id=~s/\s/\_/g;               

                   $genemark_gene{$gene_id}->{gene_id}=$gene_id;
                   $genemark_gene{$gene_id}->{tool_id}=$tool_id;
                   $genemark_gene{$gene_id}->{contig_id}=$contig_id;
                   $genemark_gene{$gene_id}->{strand}=$genemark_column[6];
                   $genemark_gene{$gene_id}->{frame}=$genemark_column[7]; 
        
                   if($genemark_gene{$gene_id}->{strand} eq "-"){
                     $genemark_gene{$gene_id}->{start}=$genemark_column[4];
                     $genemark_gene{$gene_id}->{end}=$genemark_column[3];

                     if($genemark_gene{$gene_id}->{start}< $genemark_gene{$gene_id}->{end}){
                        next;
                     }else{
                        $coordinate_length=($genemark_gene{$gene_id}->{start} - $genemark_gene{$gene_id}->{end}) + 1;
                     } 

                   }else{
                     $genemark_gene{$gene_id}->{start}=$genemark_column[3];
                     $genemark_gene{$gene_id}->{end}=$genemark_column[4];

                     if($genemark_gene{$gene_id}->{start}> $genemark_gene{$gene_id}->{end}){
                        next;
                     }else{
                        $coordinate_length=($genemark_gene{$gene_id}->{end} - $genemark_gene{$gene_id}->{start}) + 1;
                     }  
                   } 
 
                   $genemark_gene{$gene_id}->{tool}="genemark";                   
                   my $seq=$predicted_sequence{genemark}->{$tool_id};

                   if(!defined($seq)){
                     next;
                   }
                   my $predicted_seq_len=length($seq);
                   
                   my $start_gene=$genemark_gene{$gene_id}->{start};
                   my $end_gene=$genemark_gene{$gene_id}->{end};
                 
                   if($genemark_gene{$gene_id}->{strand} eq "-"){
                     $start_gene=$genemark_gene{$gene_id}->{end};
                     $end_gene=$genemark_gene{$gene_id}->{start};
                   }

                   print ALL_GENE_PRED $genemark_gene{$gene_id}->{gene_id},"\t",                         
                         $genemark_gene{$gene_id}->{contig_id},"\t",
                         $start_gene,"\t",
                         $end_gene,"\t",
                         $genemark_gene{$gene_id}->{strand},"\t",
                         $genemark_gene{$gene_id}->{frame},"\t",
                         $genemark_gene{$gene_id}->{tool},"\t",
                         $genemark_gene{$gene_id}->{tool_id},"\n";     

                   $gene_id_count++;                   
            }                     
        }

        #### Parse Result from Prodigal ##### 

        if(%prodigal_predict){

             my @prodigal_predict=@{$prodigal_predict{$genome_name}};

             my $gene_id_count="00001";       
             my $tool_id_count=1;   
             my $prev_contig='';

             foreach my $line(@prodigal_predict){

                     if($line=~/\#/){
                       next;
                     }

                     my @prodigal_column=split(" ",$line);

                     my $coordinate_length=0;

                     my $gene_id=$genome_name."_".$gene_id_count;

                     if((!defined($prev_contig)) or ($prodigal_column[0] ne $prev_contig)){
                         $tool_id_count=1;  
                         $prev_contig=$prodigal_column[0];             
                     }

                     $prodigal_gene{$gene_id}->{gene_id}=$gene_id;
                     $prodigal_gene{$gene_id}->{tool_id}=$prodigal_column[0]."_".$tool_id_count;
                     $prodigal_gene{$gene_id}->{contig_id}=$prodigal_column[0];
                     $prodigal_gene{$gene_id}->{strand}=$prodigal_column[6];
                     $prodigal_gene{$gene_id}->{frame}=$prodigal_column[7];
                     $prodigal_gene{$gene_id}->{tool}="prodigal"; 

                     if($prodigal_gene{$gene_id}->{strand} eq "-"){
                        $prodigal_gene{$gene_id}->{start}=$prodigal_column[4];
                        $prodigal_gene{$gene_id}->{end}=$prodigal_column[3];

                        if($prodigal_gene{$gene_id}->{start}< $prodigal_gene{$gene_id}->{end}){
                          next;
                        }else{
                          $coordinate_length=($prodigal_gene{$gene_id}->{start} - $prodigal_gene{$gene_id}->{end}) + 1;
                        }                

                     }else{
                       $prodigal_gene{$gene_id}->{start}=$prodigal_column[3];
                       $prodigal_gene{$gene_id}->{end}=$prodigal_column[4];

                       if($prodigal_gene{$gene_id}->{start} > $prodigal_gene{$gene_id}->{end}){
                          next;
                       }else{
                           $coordinate_length=($prodigal_gene{$gene_id}->{end} - $prodigal_gene{$gene_id}->{start}) + 1;
                       }        
                     }  

                     my $tool_id=$prodigal_gene{$gene_id}->{tool_id};
                     my $seq=$predicted_sequence{prodigal}->{$tool_id};

                     if(!defined($seq)){
                       next;
                     }
                     my $predicted_seq_len=length($seq);
                
                     $gene_id_count++;
                     $tool_id_count++;
                     
                     my $start_gene=$prodigal_gene{$gene_id}->{start};
                     my $end_gene=$prodigal_gene{$gene_id}->{end};
                 
                    if($prodigal_gene{$gene_id}->{strand} eq "-"){
                        $start_gene=$prodigal_gene{$gene_id}->{end};
                        $end_gene=$prodigal_gene{$gene_id}->{start};
                    }
                   
                     print ALL_GENE_PRED $prodigal_gene{$gene_id}->{gene_id},"\t",                           
                           $prodigal_gene{$gene_id}->{contig_id},"\t",
                           $start_gene,"\t",
                           $end_gene,"\t",
                           $prodigal_gene{$gene_id}->{strand},"\t",
                           $prodigal_gene{$gene_id}->{frame},"\t",
                           $prodigal_gene{$gene_id}->{tool},"\t",
                           $prodigal_gene{$gene_id}->{tool_id},"\n";         
             }
        } 

        #### Parse Result from FragScan #####

        if(%fragscan_predict){

             my @fragscan_predict=@{$fragscan_predict{$genome_name}};

             my $gene_id_count="00001";
             my $contig_id='';

             foreach my $line(@fragscan_predict){                                

                 if($line=~/\>/){
                    $line=~/(>)(.+)/;
                    $contig_id=$2;
                    chomp($contig_id);
                    next;
                 }

                 my $coordinate_length=0;

                 my $gene_id=$genome_name."_".$gene_id_count;              

                 my @fragscan_column=split("\t",$line);                 

                 $fragscan_gene{$gene_id}->{gene_id}=$gene_id;
                 $fragscan_gene{$gene_id}->{contig_id}=$contig_id;
                 $fragscan_gene{$gene_id}->{tool_id}=$contig_id."_".$fragscan_column[0]."_".$fragscan_column[1]."_".$fragscan_column[2];
                 $fragscan_gene{$gene_id}->{strand}=$fragscan_column[2];
                 $fragscan_gene{$gene_id}->{frame}=$fragscan_column[3];
                 $fragscan_gene{$gene_id}->{tool}="fragscan";

                 if($fragscan_gene{$gene_id}->{strand} eq "-"){
                    $fragscan_gene{$gene_id}->{start}=$fragscan_column[1];
                    $fragscan_gene{$gene_id}->{end}=$fragscan_column[0];

                    if($fragscan_gene{$gene_id}->{start}< $fragscan_gene{$gene_id}->{end}){
                       next;
                    }else{
                       $coordinate_length=(($fragscan_gene{$gene_id}->{start} - $fragscan_gene{$gene_id}->{end}) - 2) - ($fragscan_column[3]-1);
                    }

                 }else{
                     $fragscan_gene{$gene_id}->{start}=$fragscan_column[0];
                     $fragscan_gene{$gene_id}->{end}=$fragscan_column[1];

                     if($fragscan_gene{$gene_id}->{start} > $fragscan_gene{$gene_id}->{end}){
                       next;
                    }else{
                       $coordinate_length=(($fragscan_gene{$gene_id}->{end} - $fragscan_gene{$gene_id}->{start}) - 2) - ($fragscan_column[3]-1);
                    }
                 } 

                     my $tool_id=$fragscan_gene{$gene_id}->{tool_id};
                     my $seq=$predicted_sequence{fragscan}->{$tool_id};
         
                     if(!defined($seq)){
                       next;
                     }
                     my $predicted_seq_len=length($seq);                 
                
                 $gene_id_count++; 
                 
                     my $start_gene=$fragscan_gene{$gene_id}->{start};
                     my $end_gene=$fragscan_gene{$gene_id}->{end};
                 
                     if($fragscan_gene{$gene_id}->{strand} eq "-"){
                        $start_gene=$fragscan_gene{$gene_id}->{end};
                        $end_gene=$fragscan_gene{$gene_id}->{start};
                    }

                 print ALL_GENE_PRED $fragscan_gene{$gene_id}->{gene_id},"\t",                       
                       $fragscan_gene{$gene_id}->{contig_id},"\t",
                       $start_gene,"\t",
                       $end_gene,"\t",
                       $fragscan_gene{$gene_id}->{strand},"\t",
                       $fragscan_gene{$gene_id}->{frame},"\t",
                       $fragscan_gene{$gene_id}->{tool},"\t",
                       $fragscan_gene{$gene_id}->{tool_id},"\n";        
             }            
        }
        close ALL_GENE_PRED;

        #######  Call Compare Gene Prediction #######
        compare_gene_prediction($combine_prediction,$genome_name,$full_genome_name,\%sequence,\%predicted_sequence,$cds_dir,$protein_dir,$feature_dir,$translation_table,$overlap_base);       
}

########Compare Gene Prediction #######

sub compare_gene_prediction{

    my($combine_prediction)=(shift);
    my($genome_name)=(shift);
    my($full_genome_name)=(shift);
    my($genome_sequence)=(shift);
    my($predicted_sequence)=(shift);
    my($cds_dir)=(shift);
    my($protein_dir)=(shift);
    my($feature_dir)=(shift);
    my($translation_table)=(shift);
    my($overlap_base)=(shift);

    open(COMBINE_PREDICTION,"$combine_prediction");
    my @combine_pred=<COMBINE_PREDICTION>;

    #### create multidimensional array for combine results ###
    my @md_combine_pred=(); 

    foreach my $line(@combine_pred){
       chomp($line);
       if($line eq ''){next;}
       my @column=split(" ",$line);  
       push(@md_combine_pred,\@column);
    }

    #### Sort multidimensinal array by first contig and than gene start columns ####    
    @md_combine_pred=sort{$a->[1] cmp $b->[1] ||
                          $a->[2] <=> $b->[2]
                         }@md_combine_pred; 

    ##### compare prediction ####
    my $prev_contig='';
    my $prev_strand='';
    my $prev_start=0;
    my $prev_end=0;  
    my $prev_frame=-1;  
    my @gene_range=();
    my $gene_index_genome=1;
    my %gene_feature=();
    my %single_feature=();
           
    my $gene_count="00001";
    my $range=undef;

    
    mkdir("$cds_dir/MULTIPLE");
    mkdir("$protein_dir/MULTIPLE");
    mkdir("$feature_dir/MULTIPLE");
    mkdir("$cds_dir/SINGLE");
    mkdir("$protein_dir/SINGLE");
    mkdir("$feature_dir/SINGLE");

    open(GENE_CLASS1,">$cds_dir/MULTIPLE/$genome_name.fasta");

    open(GENE_CLASS2,">$cds_dir/SINGLE/$genome_name.single.fasta");
   
    foreach my $row(@md_combine_pred) {
         
         my @column=@{$row};
         
         my $current_contig=$column[1];
         my $current_start=$column[2];
         my $current_end=$column[3];
         my $current_strand=$column[4];
         my $current_frame=$column[5];
         my $current_tool=$column[6];
         
         if($current_strand eq "-"){
             $current_start=$column[3];
             $current_end=$column[2];
         }
          
         my $strand='';
         my $new_range=undef;

         if($current_strand eq "+"){
            $strand="+1";
         }else{
            $strand="-1";
         }

         if($current_strand eq "+"){

            $new_range = Bio::Range->new(-start=>$current_start, -end=>$current_end, -strand=>$current_strand);

         }elsif($current_strand eq "-"){

            $new_range = Bio::Range->new(-start=>$current_end, -end=>$current_start, -strand=>$current_strand);
         }
 
         if(($current_contig ne $prev_contig) or ($current_strand ne $prev_strand)){

            $prev_contig='';
            $prev_strand='';
            $prev_start=0;
            $prev_end=0;
            $prev_frame=-1;
            $gene_index_genome=1;      
            
            goto NEW_GENE;  ### re-initialize variables for new gene ####
                                        
         }elsif(($current_contig eq $prev_contig) and ($current_strand eq $prev_strand)){

             #### overlap between predicted gene regions ###
             my @overlap_range=$range->intersection($new_range); 

             my $overlap_region=0;
             
             if(scalar(@overlap_range)>0){

                $overlap_region=$overlap_range[1]-$overlap_range[0];
             }

             if($overlap_region>=$overlap_base){

                push(@gene_range,$row);   
                $range=$new_range;
                next;

             }elsif($range->contains($new_range)){

                push(@gene_range,$row);
                next;

             }else{
             
               $prev_contig='';
               $prev_strand='';
               $prev_start=0;
               $prev_end=0;
               $prev_frame=-1;
               $gene_index_genome=1;        
          
               goto NEW_GENE;  ### re-initialize variables for new gene ####
             }
         }

         NEW_GENE:

         $prev_contig=$current_contig;
         $prev_strand=$current_strand;
         $prev_start=$current_start;
         $prev_end=$current_end;
         $prev_frame=$current_frame;
         $range=$new_range;
                 
         if(scalar(@gene_range)>=2){

            #### get longest range for gene sequence ####            
            my($gene_feature,$gene_sequence)=get_longest_range(\@gene_range, $gene_count,$gene_index_genome,$genome_name,$full_genome_name,$genome_sequence,$predicted_sequence);

            $gene_feature{$genome_name."_".$gene_count}=$gene_feature;

            print GENE_CLASS1 "$gene_sequence\n";

            $gene_count++;

            $gene_index_genome++;

         }elsif(scalar(@gene_range)==1){

            my($gene_feature,$gene_sequence)=get_longest_range(\@gene_range, $gene_count,$gene_index_genome,$genome_name,$full_genome_name,$genome_sequence,$predicted_sequence);

            $single_feature{$genome_name."_".$gene_count}=$gene_feature;

            print GENE_CLASS2 "$gene_sequence\n";

            $gene_count++;

            $gene_index_genome++;
         }

         @gene_range=();
         push(@gene_range,$row);                 
    } 
    
    close GENE_CLASS1;
    close GENE_CLASS2;

    system("transeq -sequence $cds_dir/MULTIPLE/$genome_name.fasta -outseq $protein_dir/MULTIPLE/$genome_name.aa.fasta -table $translation_table"); 

    system("transeq -sequence $cds_dir/SINGLE/$genome_name.single.fasta -outseq $protein_dir/SINGLE/$genome_name.aa.single.fasta -table $translation_table");

    ###### PRINT GENE FEATURE INFORMATION ######
    open(GENE_FEATURE_CLASS1,">$feature_dir/MULTIPLE/$genome_name.feature.txt");
    open(GENE_FEATURE_CLASS2,">$feature_dir/SINGLE/$genome_name.feature.single.txt");

    print GENE_FEATURE_CLASS1 "#feature_id\tfeature_type\tgenome_id\tgenome_type\tgenome_name\tgenome_length\tfeature_start\tfeature_end\tnuc_lenth\taa_length\tstrand\tindex_on_genome\tdescription\n";
 
    print GENE_FEATURE_CLASS2 "#feature_id\tfeature_type\tgenome_id\tgenome_type\tgenome_name\tgenome_length\tfeature_start\tfeature_end\tnuc_lenth\taa_length\tstrand\tindex_on_genome\tdescription\n";

    foreach my $feature_id(keys %gene_feature){
 
       if($gene_feature{$feature_id}){ 
          
          print GENE_FEATURE_CLASS1 $gene_feature{$feature_id},"\n";
       }
    }

    foreach my $genefeature_id(keys %single_feature){

       if($single_feature{$genefeature_id}){     
 
          print GENE_FEATURE_CLASS2 $single_feature{$genefeature_id},"\n";
        }
    }

    close GENE_FEATURE_CLASS1;

    close GENE_FEATURE_CLASS2;
}

###### Find Longest Range #####
sub get_longest_range {

   my(@gene_range)=@{(shift)}; 
   my($gene_count)=(shift);
   my($gene_index_genome)=(shift);
   my($genome_name)=(shift);
   my($full_genome_name)=(shift);
   my(%genome_sequence)=%{(shift)};
   my(%predicted_sequence)=%{(shift)};

   my $gene_feature='';
   my $gene_sequence='';
 
     my %gene=();
     $gene{tool}=''; 
     $gene{index_on_genome}=$gene_index_genome;  

     my $longest_gene_region=0;

     foreach my $row(@gene_range){
         my @column=@{$row};                        
         my $genome_seq=$genome_sequence{$genome_name}->{$column[1]};
         my $gene_seq='';  
       
         if($column[4] eq "+"){ 

             my $gene_length=($column[3] - $column[2]) + 1;             

             if($gene_length>=$longest_gene_region){

                chomp($column[7]);
                $gene{seq}=$gene_seq;
                $gene{contig}=$column[1];
                $gene{genome_length}=length($genome_seq);
                $gene{start}=$column[2];
                $gene{end}=$column[3];
                $gene{strand}=$column[4];
                $gene{frame}=$column[5];                
                $gene{nuc_length}=($gene{end}-$gene{start})+1;
                $gene{aa_length}=($gene{nuc_length}/3);             
                $gene{description}="NULL";  
                $gene{tool}=$column[6];
                $gene{tool_id}=$column[7];
                $longest_gene_region=$gene_length;
                $gene{seq}=$predicted_sequence{$column[6]}->{$column[7]};
             } 
                                                           
         }elsif($column[4] eq "-"){ 

            my $gene_length=($column[3] - $column[2]) + 1;

            if($gene_length>=$longest_gene_region){ 
               chomp($column[7]);
               $gene{seq}=$gene_seq;
               $gene{contig}=$column[1];
               $gene{genome_length}=length($genome_seq);
               $gene{start}=$column[3];
               $gene{end}=$column[2];
               $gene{strand}=$column[4];
               $gene{frame}=$column[5];               
               $gene{nuc_length}=($gene{start}-$gene{end})+1;
               $gene{aa_length}=($gene{nuc_length}/3);             
               $gene{description}="NULL"; 
               $gene{tool}=$column[6];  
               $gene{tool_id}=$column[7];
               $longest_gene_region=$gene_length;
               $gene{seq}=$predicted_sequence{$column[6]}->{$column[7]};
            }                                     
        }                              
   }   
     if($gene{seq}){

       $gene_sequence=">".$genome_name."_".$gene_count."\t".$gene{tool}."\t".$gene{tool_id}."\n".$gene{seq};

       $gene_feature=$genome_name."_".$gene_count."\t"."CDS"."\t".$gene{contig}."\t"."contig"."\t".$full_genome_name."\t".$gene{genome_length}."\t".$gene{start}."\t".$gene{end}."\t".$gene{nuc_length}."\t".$gene{aa_length}."\t".$gene{strand}."\t".$gene{index_on_genome}."\t".$gene{description};
     
     }else{    
         print "Sequence not found\n";
     }
 
   return($gene_feature,$gene_sequence);   
}

