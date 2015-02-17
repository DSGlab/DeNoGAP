##### This script make genebank file from feature and sequence file #######
#!/usr/bin/perl

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;


#### INPUT #####
my $genome_name=(shift);
my $genome_file=(shift);
my $gene_file=(shift);
my $protein_file=(shift);
my $feature_file=(shift);
my $genbank_file=(shift);


#### READ GENOME #####

my %genome_seq=();
        
my $seqio_obj = Bio::SeqIO->new(-file => "$genome_file", -format => "fasta");
     
   while(my $seq_obj= $seqio_obj->next_seq){

        my $seq_id=$seq_obj->display_id;
        my $seq=$seq_obj->seq;
        $seq=~s/\s+//g;
        $seq=~s/\n+//g;
        
        if($seq_id=~/\|/){
        
          my @seq_id=split(/\|/,$seq_id);
          
          if($seq_id=~/^[(sp)|(tr)]/){
             $seq_id=$seq_id[1];
          }else{         
             $seq_id=pop(@seq_id);
          }       
        }
        $genome_seq{$genome_name}->{$seq_id}=$seq;
        
    }


#### READ CDS ####

my %cds_seq=();

#opendir(GENE,"$gene_dir/$genome_name");
#my @gene_file=readdir(GENE);

#foreach my $gene_file(@gene_file){

     #  my $genome_name=$gene_file;
     #   $genome_name=~s/(\.)(.+)//g;
        
     my $seqio_obj = Bio::SeqIO->new(-file => "$gene_file", -format => "fasta");
     
     while(my $seq_obj= $seqio_obj->next_seq){

        my $seq_id=$seq_obj->display_id;
        my $seq=$seq_obj->seq;
        $seq=~s/\s+//g;
        $seq=~s/\n+//g;
        
        if($seq_id=~/\|/){
        
          my @seq_id=split(/\|/,$seq_id);
          
          if($seq_id=~/^[(sp)|(tr)]/){
             $seq_id=$seq_id[1];
          }else{         
             $seq_id=pop(@seq_id);
          }       
        }
        $cds_seq{$genome_name}->{$seq_id}=$seq;
        
    }
#}

#### PROTEIN SEQ ####

my %protein_seq=();

#opendir(PROTEIN,"$protein_dir/$genome_name");
#my @protein_file=readdir(PROTEIN);

#foreach my $protein_file(@protein_file){

#     my $genome_name=$protein_file;
#        $genome_name=~s/(\.)(.+)//g;
        
     my $seqio_obj = Bio::SeqIO->new(-file => "$protein_file", -format => "fasta");
     
     while(my $seq_obj= $seqio_obj->next_seq){

        my $seq_id=$seq_obj->display_id;
        my $seq=$seq_obj->seq;
        $seq=~s/\s+//g;
        $seq=~s/\n+//g;
        
        if($seq_id=~/\|/){
        
          my @seq_id=split(/\|/,$seq_id);
          
          if($seq_id=~/^[(sp)|(tr)]/){
             $seq_id=$seq_id[1];
          }else{         
             $seq_id=pop(@seq_id);
          }       
        }
        
        $seq_id=~s/(_1)$//g; 
        $protein_seq{$genome_name}->{$seq_id}=$seq;
        
    }
#}


##### FEATURE FILE ######

#opendir(FEATURE_DIR,"$feature_dir/$genome_name");
#my @feature_file=readdir(FEATURE_DIR);

#foreach my $feature_file(@feature_file){

#    if($feature_file=~/^\.+$/ or $feature_file=~/^\.(\w+)$/){
#        next;
#    }
#    my $genome_name=$feature_file;
#    $genome_name=~s/(\.)(.+)//g;

     print "$genome_name\n";
     
     my %feature_info=();

     open(FEATURE_FILE,"$feature_file");
     my @feature=<FEATURE_FILE>;
        
     foreach my $feature_line(@feature){
     
          if($feature_line=~/\#/){
             next;
          }
          
          my @feature_column=split("\t",$feature_line);
          
          my $contig_id=$feature_column[2];
          my $locus_id=$feature_column[0];
          
          if($contig_id=~/\|/){ 
               my @contig_id=split(/\|/,$contig_id);        
               $contig_id=pop(@contig_id);      
          }  
          
          $feature_info{$contig_id}->{$locus_id}=\@feature_column;                                        
    }  
    
     ########### Print GeneBank File ######## 
     my $prev_contig='';
     
    foreach my $contig_id(sort{$a cmp $b} keys %feature_info){
    
          my %locus_tag=%{$feature_info{$contig_id}};
          
          my %feature=();
          
          my $seq_obj=Bio::Seq->new(-display_id => $contig_id,
                                    -seq => $genome_seq{$genome_name}->{$contig_id}
                                    );    
          
          foreach my $locus_id(sort{$a cmp $b} keys %locus_tag){
          
             my @feature_column=@{$locus_tag{$locus_id}};
          
             $feature{LOCUS_ID}=$feature_column[0];
             $feature{TYPE}=$feature_column[1];
             $feature{CONTIG_ID}=$feature_column[2];
             $feature{GENOME_TYPE}=$feature_column[3];
             $feature{STRAIN}=$feature_column[4];
             $feature{GENOME_LENGTH}=$feature_column[5];
             $feature{FEATURE_START}=$feature_column[6];
             $feature{FEATURE_END}=$feature_column[7];
             $feature{NUC_LENGTH}=$feature_column[8];
             $feature{PROTEIN_LENGTH}=$feature_column[9];
             $feature{FRAME}=1;
             
             if($feature_column[10] eq "+"){
             
                $feature{STRAND}=1;
                
             }elsif($feature_column[10] eq "-"){
                $feature{STRAND}=-1;
                $feature{FEATURE_END}=$feature_column[6];
                $feature{FEATURE_START}=$feature_column[7];           
             } 
               
             $feature{INDEX}=$feature_column[11];
             $feature{PRODUCT}=$feature_column[12];  
             chomp($feature{PRODUCT});
          
             if($feature{CONTIG_ID}=~/\|/){ 
               my @contig_id=split(/\|/,$feature{CONTIG_ID});        
               $feature{CONTIG_ID}=pop(@contig_id);      
             }  
             
             ###### if protein length is not multiple of 3 than check for codon start ######
             if($feature{PROTEIN_LENGTH}=~/\.+/){
                my $len_aa=length($protein_seq{$genome_name}->{$feature{LOCUS_ID}});
                my $len_cds=($len_aa * 3);        
                my $extra_nuc=$feature{NUC_LENGTH}-$len_cds;
                my $frame=0;
                if($extra_nuc>3){
                  $frame=($extra_nuc-3) + 1;
                }else{
                  $frame=$extra_nuc + 1;
                }              
                $feature{PROTEIN_LENGTH}=$len_aa;
                $feature{NUC_LENGTH}=$len_cds;
                $feature{FRAME}=$frame;   
             }   
                 
             if($contig_id ne $prev_contig){
             
                 my $feat_source=new Bio::SeqFeature::Generic(-start=>1,
                                                              -end=> $feature{GENOME_LENGTH},
                                                              -strand=>1,
                                                              -primary_tag => "source",
                                                              -tag => {organism => $feature{STRAIN}}
                                                             );                                  
                 $seq_obj->add_SeqFeature($feat_source); 
                 
                 $prev_contig=$contig_id;              
             }   
             
                 my $feat_gene=new Bio::SeqFeature::Generic(-start=>$feature{FEATURE_START},
                                                              -end=> $feature{FEATURE_END},
                                                              -strand=>$feature{STRAND},
                                                              -primary_tag => "gene",
                                                              -tag => {locus_tag=>$feature{LOCUS_ID}}
                                                             );                                  
                 $seq_obj->add_SeqFeature($feat_gene); 
                                                  
          
             my $feat_cds=new Bio::SeqFeature::Generic(-start => $feature{FEATURE_START},
                                                       -end => $feature{FEATURE_END},
                                                       -strand => $feature{STRAND},
                                                       -primary_tag => $feature{TYPE},
                                                       -tag => {locus_tag=>$feature{LOCUS_ID},
                                                                protein_id => $feature{LOCUS_ID},
                                                                cds_length => $feature{NUC_LENGTH},
                                                                protein_length => $feature{PROTEIN_LENGTH},
                                                                trans_table => 11,
                                                                codon_start => $feature{FRAME},
                                                                product => $feature{PRODUCT},
                                                                translation => $protein_seq{$genome_name}->{$feature{LOCUS_ID}}
                                                               }
                                                       );                                
           $seq_obj->add_SeqFeature($feat_cds);     

          
       }
       
       my $io=Bio::SeqIO->new(-format => "genbank",
                              -file => ">>$genbank_file");     
       $io->write_seq($seq_obj);
    }                             
#}



