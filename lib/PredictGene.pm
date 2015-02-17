##### Module to Perform Gene Prediction #####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package PredictGene;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use SequenceHash;
use SQLiteDB;
use lib "$Bin/../lib";
use vars qw(@ISA @EXPORT @EXPORT_OK $db_dir $db_name);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(run_gene_prediction);


sub run_gene_prediction{

    $db_dir=(shift);
    $db_name=(shift);
    my(%gene_prediction)=%{(shift)};
    my(%sequence_directory)=%{(shift)};
    my(%out_dir)=%{(shift)};
    my $exe=(shift);

    ##### Read Genome File #####

    opendir(GENOME_DIR,"$sequence_directory{GENOME_DIR}");
    my @genome_file=readdir(GENOME_DIR);

    my @tool_selected=split('',$gene_prediction{TOOL});

    my $glimmer_predict={};
    my $genemark_predict={};
    my $prodigal_predict={};
    my $fragscan_predict={};

    foreach my $selection(@tool_selected){

       if($selection eq 1){
         $glimmer_predict=run_glimmer(\%sequence_directory,\@genome_file,$out_dir{gene_prediction},$exe);
       }
       elsif($selection eq 2){
         $genemark_predict=run_genemark(\%sequence_directory,\@genome_file,$out_dir{gene_prediction},$exe);
       }elsif($selection eq 3){
         $prodigal_predict=run_prodigal(\%sequence_directory,\@genome_file,$out_dir{gene_prediction},$exe);
       }elsif($selection eq 4){
         $fragscan_predict=run_fragscan(\%sequence_directory,\@genome_file,$out_dir{gene_prediction},$exe);
       }
    }
}

######## GLIMMER #####

sub run_glimmer{

    my(%sequence_directory)=%{(shift)};
    my(@genome_file)=@{(shift)};
    my($out_dir)=(shift);
    my $exe=(shift);

    print "Running Glimmer\n";

    $out_dir=$out_dir."/glimmer_predict";
    mkdir($out_dir);
    #### Glimmer Path ####
    my $glimmer_script=$exe."/glimmer/scripts";
    my $glimmer_bin=$exe."/glimmer/bin";

    #### Set Glimmer options ###

    my $glimmeropts="-o 0 -g 100 -t 30 -A atg,ttg,gtg -z 11";

    my %glimmer_predict=();
       
    foreach my $genome(@genome_file){
      
       if($genome=~/^\.+$/){next;}

       my $genome_name=$genome;
          $genome_name=~s/\.(.+)//g;

       if((-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fasta") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fa") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.faa")){
          next;
        }
 
       print "GENOME:$genome\n";    
   
       mkdir("$Bin/$out_dir/$genome") or die "cannot create output directory for $genome at $Bin/$out_dir/$genome\n";

       my $seqio_obj = Bio::SeqIO->new(-file => "$sequence_directory{GENOME_DIR}/$genome", -format => "fasta");

       print "Finding long orfs for training\n";

       open(LONG_ORF,">$Bin/$out_dir/$genome/$genome.longorfs");

       while(my $seq_obj= $seqio_obj->next_seq){

               my $seq=$seq_obj->seq;
               my $seq_id=$seq_obj->display_id;

               open(TMP_GENOME,">$Bin/$out_dir/$genome/tmp_genome.fasta");
               print TMP_GENOME ">$seq_id\n$seq";
               close TMP_GENOME;

               #### Find long, non-overlapping orfs to use as training set ###
               my $status1=system("long-orfs -n -t 1.15 -z 11 --min_len 100 $Bin/$out_dir/$genome/tmp_genome.fasta $Bin/$out_dir/$genome/tmp_genome.longorfs");
               if($status1!=0){
                 print "Failed to extract long-orfs from $seq_id\n"; next;
               }

               open(TMP_LONGORF,"$Bin/$out_dir/$genome/tmp_genome.longorfs");
               my @long_orf=<TMP_LONGORF>;
               close TMP_LONGORF;               

               foreach my $line(@long_orf){
                    my @line=split(" ",$line);
                    print LONG_ORF "$line[0]\t$seq_id\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
               }                             
       }
 
     #### Extract Sequence for long-orfs to use as training set ####
     print "Extract long-orfs for training\n";
     my $status2=system("multi-extract --minlen 100 $sequence_directory{GENOME_DIR}/$genome $Bin/$out_dir/$genome/$genome.longorfs > $Bin/$out_dir/$genome/$genome.train");

     if($status2!=0){
       print "Failed to extract training sequence for $genome\n"; exit;
     }
     #### Build ICM from the training Sequence ####
     print "Building ICM\n"; 
     my $status3=system("build-icm -r $Bin/$out_dir/$genome/$genome.icm < $Bin/$out_dir/$genome/$genome.train"); 
   
     if($status3!=0){
       print "Failed to build ICM for $genome\n"; exit;
     }
     #### Run Glimmer ####
     print "Running Glimmer\n";
     my $status4=system("glimmer3 $glimmeropts $sequence_directory{GENOME_DIR}/$genome $Bin/$out_dir/$genome/$genome.icm $Bin/$out_dir/$genome/$genome");    

     if($status4!=0){
       print "Failed to run Glimmer for $genome\n"; exit;
     }

     open(GLIMMER_PREDICT,"$Bin/$out_dir/$genome/$genome.predict");
     my @glimmer_predict=<GLIMMER_PREDICT>;
     close GLIMMER_PREDICT;

     open(GLIMMER_PREDICT1,">$Bin/$out_dir/$genome/$genome.1.predict");
  
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

     my $status5=system("multi-extract --minlen 100 $sequence_directory{GENOME_DIR}/$genome $Bin/$out_dir/$genome/$genome.1.predict > $Bin/$out_dir/$genome/$genome.fna");
     
     if($status5!=0){
       print "Failed to extract gene sequence for $genome\n"; exit;
     }
   }    
}

###### GeneMarkS #######

sub run_genemark {

    my(%sequence_directory)=%{(shift)};
    my(@genome_file)=@{(shift)};
    my($out_dir)=(shift);
    my $exe=(shift);

    print "Running GeneMark\n";

    $out_dir=$out_dir."/genemark_predict";
    mkdir($out_dir);

    my %genemark_predict=();

    ##### GeneMark Path #####

    my $genemark_bin=$exe."/genemark/gmsuite";

     foreach my $genome(@genome_file){
      
       if($genome=~/^\.+$/){next;}

       my $genome_name=$genome;
       $genome_name=~s/\.(.+)//g;

      if((-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fasta") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fa") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.faa")){
        next;
      } 

       $genome=~/(\w+)(\.)(\w+)/;
 
       print "GENOME:$genome\n";    
   
       mkdir("$Bin/$out_dir/$genome") or die "cannot create output directory for $genome at $Bin/$out_dir/$genome\n";
       chdir("$Bin/$out_dir/$genome");
       system("gmsn.pl --name $genome --combine --gm --species $1 --prok --fnn --faa --clean --format GFF $Bin/$sequence_directory{GENOME_DIR}/$genome");
       chdir($Bin);    
    } 
}

##### PRODIGAL #####

sub run_prodigal {

    my(%sequence_directory)=%{(shift)};
    my(@genome_file)=@{(shift)};
    my($out_dir)=(shift);
    my $exe=(shift);

    print "Running Prodigal\n";

    $out_dir=$out_dir."/prodigal_predict";
    mkdir($out_dir); 

    my %prodigal_predict=();

    #### Prodigal Path ####
    my $prodigal_bin=$exe."/prodigal";

    foreach my $genome(@genome_file){
      
       if($genome=~/^\.+$/){next;}

       my $genome_name=$genome;
       $genome_name=~s/\.(.+)//g;

       if((-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fasta") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fa") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.faa")){
          next;
       } 

       $genome=~/(\w+)(\.)(\w+)/;
 
       print "GENOME:$genome\n";

       mkdir("$Bin/$out_dir/$genome") or die "cannot create output directory for $genome at $Bin/$out_dir/$genome\n";
       system("prodigal -g 11 -m -t $Bin/$out_dir/$genome/$genome.training -i $Bin/$sequence_directory{GENOME_DIR}/$genome");
       system("prodigal -a $Bin/$out_dir/$genome/$genome.faa -d $Bin/$out_dir/$genome/$genome.fnn -c -f gff -g 11 -m -o $Bin/$out_dir/$genome/$genome.gff -s $Bin/$out_dir/$genome/$genome.score  -t $Bin/$out_dir/$genome/$genome.training -i $Bin/$sequence_directory{GENOME_DIR}/$genome");         
     } 
}

##### FRAGSCAN ####

sub run_fragscan{

    my(%sequence_directory)=%{(shift)};
    my(@genome_file)=@{(shift)};
    my($out_dir)=(shift);
    my $exe=(shift);

    print "Running FragScan\n";

    $out_dir=$out_dir."/fragscan_predict";
    mkdir($out_dir); 

    my %fragscan_predict=();

    #### Fragscan Path ####
    my $fragscan_bin=$exe."/fragscan";

    foreach my $genome(@genome_file){
      
       if($genome=~/^\.+$/){next;}

       my $genome_name=$genome;
       $genome_name=~s/\.(.+)//g;

       if((-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fasta") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fa") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.faa")){
         next;
       }

       $genome=~/(\w+)(\.)(\w+)/;
 
       print "GENOME:$genome\n";
       
       mkdir("$Bin/$out_dir/$genome") or die "cannot create output directory for $genome at $Bin/$out_dir/$genome\n";
       system("run_FragGeneScan.pl -genome=$Bin/$sequence_directory{GENOME_DIR}/$genome -out=$Bin/$out_dir/$genome/$genome -complete=1 -train=complete");  
    }   
}

#### Compare Gene Prediction to find best gene sequence ####

sub read_gene_prediction {

    $db_dir=(shift);
    $db_name=(shift);
    my(%gene_prediction)=%{(shift)};
    my(%sequence_directory)=%{(shift)};
    my(%out_dir)=%{(shift)};
    my($translation_table)=(shift);
    my($annotation_database)=(shift);
    
    opendir(GENOME_DIR,"$sequence_directory{GENOME_DIR}");
    my @genome_file=readdir(GENOME_DIR);

    my $out_dir=$out_dir{gene_prediction};

    my %glimmer_predict=();
    my %genemark_predict=();
    my %prodigal_predict=();
    my %fragscan_predict=();

    my($SequenceLenHashTable,$SequenceIDHashTable,$SequenceHash)=SequenceHash::getSequenceFeature($sequence_directory{GENOME_DIR});

    my %seq_len_table=%{$SequenceLenHashTable};
    my %seq_id_table=%{$SequenceIDHashTable};
    my %sequence=%{$SequenceHash};

    mkdir($out_dir."/combine_prediction");
    mkdir($out_dir."/gene_sequence");
    mkdir($out_dir."/protein_sequence");
    mkdir($out_dir."/gene_feature");
    mkdir($out_dir."/annotation");
   
    foreach my $genome(@genome_file){ 

        if($genome=~/^\.+$/){next;}

          $genome=~/(\w+)(\.)(\w+)/;
          
          my $genome_name=$1;

           if((-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fasta") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.fa") or (-s "$sequence_directory{PROTEIN_DIR}/$genome_name.faa")){
               print "Protein sequence for $genome_name already present in Protein directory...skip\n";
               next;
           } 

          my $get_full_name="Select DISTINCT(genome_name) from OrganismInfo WHERE abbreviation='$genome_name'";
          my($record_full_name)=SQLiteDB::get_record($db_dir,$db_name,$get_full_name);
          my $full_genome_name='';

          if(scalar($record_full_name)>0){
            foreach my $row(@{$record_full_name}){  
               $full_genome_name=shift(@{$row});
            }
          }else{
             $full_genome_name="NULL";
          }

        print "$full_genome_name\n";

        ###### READ GLIMMER OUTPUT ####
        open(GLIMMER,"$Bin/$out_dir/glimmer_predict/$genome/$genome.1.predict");
        my @glimmer_predict=<GLIMMER>;          
        $glimmer_predict{$genome_name}=\@glimmer_predict;
        close GLIMMER;

        ###### READ GENEMARK OUTPUT #####
        open(GENEMARK,"$Bin/$out_dir/genemark_predict/$genome/$genome.gff");
        my @genemark_predict=<GENEMARK>; 
        $genemark_predict{$genome_name}=\@genemark_predict;
        close GENEMARK;

        ##### READ PRODIGAL OUTPUT #####
        open(PRODIGAL,"$Bin/$out_dir/prodigal_predict/$genome/$genome.gff");
        my @prodigal_predict=<PRODIGAL>;
        $prodigal_predict{$genome_name}=\@prodigal_predict;
        close PRODIGAL; 

        ##### READ FRAGSCAN OUTPUT ######
        open(FRAGSCAN,"$Bin/$out_dir/fragscan_predict/$genome/$genome.out") or die "cannot open fragscan feature file\n";
        my @fragscan_predict=<FRAGSCAN>;
        $fragscan_predict{$genome_name}=\@fragscan_predict;
        close FRAGSCAN;   

        ##### READ PREDICTED CDS SEQUENCE ########
        my %predicted_sequence=();

        ##### GLIMMER SEQ ######
        my $glimmer_seqobj = Bio::SeqIO->new(-file => "$Bin/$out_dir/glimmer_predict/$genome/$genome.fna", -format => "fasta");
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
        my $genemark_seqobj = Bio::SeqIO->new(-file => "$Bin/$out_dir/genemark_predict/$genome/$genome.fnn", -format => "fasta");
        my %genemark_seq=();
  
        while(my $genemark_seqobj= $genemark_seqobj->next_seq){
             my $seq_id=$genemark_seqobj->display_id;             
             my $seq=$genemark_seqobj->seq;
             $genemark_seq{$seq_id}=$seq;         
        }

        ####### Prodigal Seq ####

        my $prodigal_seqobj = Bio::SeqIO->new(-file => "$Bin/$out_dir/prodigal_predict/$genome/$genome.fnn", -format => "fasta");
        my %prodigal_seq=();
  
        while(my $prodigal_seqobj= $prodigal_seqobj->next_seq){
             my $seq_id=$prodigal_seqobj->display_id;             
             my $seq=$prodigal_seqobj->seq;
             $prodigal_seq{$seq_id}=$seq;         
        }
 
        ####### FragScan Seq #####
 
        my $fragscan_seqobj = Bio::SeqIO->new(-file => "$Bin/$out_dir/fragscan_predict/$genome/$genome.ffn", -format => "fasta");
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

        my $combine_prediction=$out_dir."/combine_prediction/all_$genome.txt";

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
               
                 print ALL_GENE_PRED $glimmer_gene{$gene_id}->{gene_id},"\t",                       
                       $glimmer_gene{$gene_id}->{contig_id},"\t",
                       $glimmer_gene{$gene_id}->{start},"\t",
                       $glimmer_gene{$gene_id}->{end},"\t",
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

                   print ALL_GENE_PRED $genemark_gene{$gene_id}->{gene_id},"\t",                         
                         $genemark_gene{$gene_id}->{contig_id},"\t",
                         $genemark_gene{$gene_id}->{start},"\t",
                         $genemark_gene{$gene_id}->{end},"\t",
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
                   
                     print ALL_GENE_PRED $prodigal_gene{$gene_id}->{gene_id},"\t",                           
                           $prodigal_gene{$gene_id}->{contig_id},"\t",
                           $prodigal_gene{$gene_id}->{start},"\t",
                           $prodigal_gene{$gene_id}->{end},"\t",
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

                 print ALL_GENE_PRED $fragscan_gene{$gene_id}->{gene_id},"\t",                       
                       $fragscan_gene{$gene_id}->{contig_id},"\t",
                       $fragscan_gene{$gene_id}->{start},"\t",
                       $fragscan_gene{$gene_id}->{end},"\t",
                       $fragscan_gene{$gene_id}->{strand},"\t",
                       $fragscan_gene{$gene_id}->{frame},"\t",
                       $fragscan_gene{$gene_id}->{tool},"\t",
                       $fragscan_gene{$gene_id}->{tool_id},"\n";        
             }            
        }
        #######  Call Compare Gene Prediction #######
        compare_gene_prediction($combine_prediction,$out_dir,$genome_name,$full_genome_name,\%sequence,\%predicted_sequence,\%sequence_directory,$translation_table,$annotation_database);       
    }   

}

########Compare Gene Prediction #######

sub compare_gene_prediction{

    my($combine_prediction)=(shift);
    my($out_dir)=(shift);
    my($genome_name)=(shift);
    my($full_genome_name)=(shift);
    my($genome_sequence)=(shift);
    my($predicted_sequence)=(shift);
    my(%sequence_directory)=%{(shift)};
    my($translation_table)=(shift);
    my($annotation_database)=(shift);

    open(COMBINE_PREDICTION,"$combine_prediction");
    my @combine_pred=<COMBINE_PREDICTION>;

    #### create multidimensional array for combine results ###
    my @md_combine_pred=(); 

    foreach my $line(@combine_pred){
       my @column=split("\t",$line);       
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
    

    #### Prepare Blastable SwissProt Database ####
        
    my $gene_count="00001";
    open(GENE,">$out_dir/gene_sequence/$genome_name.fasta");
    close GENE;

    foreach my $row(@md_combine_pred) {
         
         my @column=@{$row};
         
         my $current_contig=$column[1];
         my $current_start=$column[2];
         my $current_end=$column[3];
         my $current_strand=$column[4];
         my $current_frame=$column[5];
         my $current_tool=$column[6];
 
         if(($current_contig ne $prev_contig) or ($current_strand ne $prev_strand)){
            $prev_contig='';
            $prev_strand='';
            $prev_start=0;
            $prev_end=0; 
            $prev_frame=-1;
            $gene_index_genome=1;          
            goto NEW_GENE;  ### re-initialize variables for new gene ####
                                         
         }elsif(($current_contig eq $prev_contig) and ($current_strand eq $prev_strand) and ($current_strand eq "+")){

             if(($current_start > $prev_end) or ($current_end < $prev_start)){                
                $prev_contig='';
                $prev_strand='';
                $prev_start=0;
                $prev_end=0; 
                $prev_frame=-1;              
                goto NEW_GENE;
             }else{ #### gene overlap on +ve strand ####
                  if((($current_start ne $prev_start) and ($current_end ne $prev_end))){
                      $prev_contig='';
                      $prev_strand='';
                      $prev_start=0;
                      $prev_end=0; 
                      $prev_frame=-1;                      
                      goto NEW_GENE;

                  }else{               
                      push(@gene_range,$row);
                      $prev_contig=$current_contig;
                      $prev_strand=$current_strand;
                      $prev_start=$current_start;
                      $prev_end=$current_end;
                      $prev_frame=$current_frame;
                      next;
                 }
             }
                    
         }elsif(($current_contig eq $prev_contig) and ($current_strand eq $prev_strand) and ($current_strand eq "-")){
            
            if(($current_start < $prev_end) or ($current_end > $prev_start)){
                $prev_contig='';
                $prev_strand='';
                $prev_start=0;
                $prev_end=0; 
                $prev_frame=-1; 
                goto NEW_GENE;
            }else{ ### gene overlap on negative strand ####
               if((($current_start ne $prev_start) and ($current_end ne $prev_end))){
                      $prev_contig='';
                      $prev_strand='';
                      $prev_start=0;
                      $prev_end=0; 
                      $prev_frame=-1; 
                      goto NEW_GENE;

                }else{                
                      push(@gene_range,$row);
                      $prev_contig=$current_contig;
                      $prev_strand=$current_strand;
                      $prev_start=$current_start;
                      $prev_end=$current_end;
                      $prev_frame=$current_frame;
                      next;
               }
            }
         }
         NEW_GENE:
         $prev_contig=$current_contig;
         $prev_strand=$current_strand;
         $prev_start=$current_start;
         $prev_end=$current_end;
         $prev_frame=$current_frame;
         
         if(scalar(@gene_range)>=2){
            #### get longest range for gene sequence ####            
            my($gene_feature)=get_longest_range(\@gene_range, $gene_count,$gene_index_genome,$out_dir,$genome_name,$full_genome_name,$genome_sequence,$predicted_sequence,$sequence_directory{FEATURE_DIR});

            $gene_feature{$genome_name."_".$gene_count}=$gene_feature;

            $gene_count++;
            $gene_index_genome++;
         }
         @gene_range=();
         push(@gene_range,$row);                 
    }    
    system("transeq -sequence $out_dir/gene_sequence/$genome_name.fasta -outseq $out_dir/protein_sequence/$genome_name.fasta -table $translation_table"); 

   # print "Annotating $genome_name protein sequences\n";

   # system("blastp -query $out_dir/protein_sequence/$genome_name.fasta -db $annotation_database -evalue 1e-5 -max_target_seqs 1 -outfmt 0 -out $out_dir/annotation/$genome_name.blas");

   # my $report = Bio::SearchIO->new(-format => 'blast',
   #                                 -file =>"$out_dir/annotation/$genome_name.blas"); 
    
   # my %gene_annotation=();

   # while (my $result = $report->next_result) {

   #    my $query=$result->query_name;
   #       $query=~s/(_1)$//g;

   #    while( my $hit = $result->next_hit ) {   
          
   #       while( my $hsp = $hit->next_hsp ) {
          
   #           my $hit_desc=$hit->description; 
   #           my $identity=$hsp->percent_identity;
   #           my $qstart=$hsp->start('query');
   #           my $qend=$hsp->end('query');
   #           my $sstart=$hsp->start('hit'); 
   #           my $send=$hsp->end('hit');
   #           my $qlen=$result->query_length;
   #           my $slen=$hit->length;

   #           my $qcoverage=(($qend-$qstart)+1/$qlen)*100;
   #           my $scoverage=(($send-$sstart)+1/$slen)*100;
              
   #           $hit_desc=~s/(OS=)(.+)//g;

   #           if(($identity>=30) and ($qcoverage>=50 or $scoverage>=50)){
   #               $gene_annotation{$query}=$hit_desc;
   #           }else{
   #               $gene_annotation{$query}="unknown protein";
   #           }
   #       }
   #
   #    }
   #
   #    if(!defined($gene_annotation{$query})){
   #       $gene_annotation{$query}="unknown protein";
   #    }
   # }

    ###### PRINT GENE FEATURE INFORMATION ######
    open(GENE_FEATURE,">$sequence_directory{FEATURE_DIR}/$genome_name.txt");
    print GENE_FEATURE "#feature_id\tfeature_type\tgenome_id\tgenome_type\tgenome_name\tgenome_length\tfeature_start\tfeature_end\tnuc_lenth\taa_length\tstrand\tindex_on_genome\tdescription\n";

    foreach my $genefeature_id(keys %gene_feature){
 
       if($gene_feature{$genefeature_id}){      
          print GENE_FEATURE $gene_feature{$genefeature_id},"\n";
        }
    }
 
    close GENE_FEATURE;

  
    ##### Loading Predicted Gene Information ####
    print "Copying Predicted Sequences to Query Directory\n";
    system("cp  $out_dir/gene_sequence/$genome_name.fasta $sequence_directory{NUCLEOTIDE_DIR}/$genome_name.fasta")==0 or die "ERROR: Cannot copy file\n";
    system("cp  $out_dir/protein_sequence/$genome_name.fasta $sequence_directory{PROTEIN_DIR}/$genome_name.fasta")==0 or die "ERROR: Cannot copy file\n";;

    open(FILE,"$sequence_directory{FEATURE_DIR}/$genome_name.txt");
    my @file=<FILE>;
    close FILE;   

    print "Loading Gene Feature Information\n";      
    my $gene_feature_sql_stmt="INSERT into GeneFeature (feature_id, feature_type, genome_id, genome_type, genome_name, genome_length, feature_start, feature_end, nuc_length, aa_length, strand, index_in_genome, description) VALUES 
                                (?, ?, ?, ? ,? ,? ,?, ? ,? ,? ,? ,? ,?)";
    SQLiteDB::load_from_array($db_dir,$db_name,\@file,$gene_feature_sql_stmt);        
}

###### Find Longest Range #####
sub get_longest_range {

   my(@gene_range)=@{(shift)}; 
   my($gene_count)=(shift);
   my($gene_index_genome)=(shift);
   my($out_dir)=(shift);
   my($genome_name)=(shift);
   my($full_genome_name)=(shift);
   my(%genome_sequence)=%{(shift)};
   my(%predicted_sequence)=%{(shift)};
   my($feature_dir)=(shift);

   my $predictedBy=scalar(@gene_range);
   my $gene_feature='';

 if($predictedBy >=2){ 

     open(GENE,">>$out_dir/gene_sequence/$genome_name.fasta");
     #open(GENE_FEATURE,">>$feature_dir/$genome_name.txt");

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

             if($gene_length>$longest_gene_region){
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

            my $gene_length=($column[2] - $column[3]) + 1;

            if($gene_length>$longest_gene_region){ 
               chomp($column[7]);
               $gene{seq}=$gene_seq;
               $gene{contig}=$column[1];
               $gene{genome_length}=length($genome_seq);
               $gene{start}=$column[2];
               $gene{end}=$column[3];
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
       print GENE ">".$genome_name."_".$gene_count."\t".$gene{tool}."\t".$gene{tool_id},"\n",$gene{seq},"\n";
       #print GENE_FEATURE $genome_name."_".$gene_count,"\t","CDS","\t",$gene{contig},"\t","contig","\t",$full_genome_name,"\t",$gene{genome_length},"\t",$gene{start},"\t",$gene{end},"\t",$gene{nuc_length},"\t",$gene{aa_length},"\t",$gene{strand},"\t",$gene{index_on_genome},"\t",$gene{description},"\n";
       $gene_feature=$genome_name."_".$gene_count."\t"."CDS"."\t".$gene{contig}."\t"."contig"."\t".$full_genome_name."\t".$gene{genome_length}."\t".$gene{start}."\t".$gene{end}."\t".$gene{nuc_length}."\t".$gene{aa_length}."\t".$gene{strand}."\t".$gene{index_on_genome}."\t".$gene{description};
     }
 }  
 
   return($gene_feature);   
}
