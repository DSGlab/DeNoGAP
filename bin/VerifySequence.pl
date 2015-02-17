#!/usr/bin/perl

use strict;
use warnings;
use Env;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Range;
use Tie::File;


my $blast_alignment_file=(shift);
my $feature_dir=(shift);
my $cds_dir=(shift);
my $protein_dir=(shift);
my $genome_dir=(shift);
my $project_dir=(shift);
my $verify_cds=(shift);
my $verify_protein=(shift);
my $verify_feature=(shift);
my $verify_genbank=(shift);
my $evalue_cutoff=(shift);
my $identity_cutoff=(shift);
my $query_coverage_cutoff=(shift);
my $seq_len_cutoff=(shift);

mkdir($project_dir);
mkdir($verify_cds);
mkdir($verify_protein);
mkdir($verify_feature);
mkdir($verify_genbank);

##### PARSE BLAST ALIGNMENT FILE #####
my $report = Bio::SearchIO->new(-format => 'blast', -file =>"$blast_alignment_file"); 
    
    my %gene_annotation=();

    while (my $result = $report->next_result) {

       my $query=$result->query_name;
       $query=~s/(_1)$//g;

       while( my $hit = $result->next_hit ) {   
          
          while( my $hsp = $hit->next_hsp ) {
          
              my $hit_desc=$hit->description; 
              my $evalue=$hit->significance; 
              my $identity=$hsp->percent_identity;
              my $qstart=$hsp->start('query');
              my $qend=$hsp->end('query');
              my $sstart=$hsp->start('hit'); 
              my $send=$hsp->end('hit');
              my $qlen=$result->query_length;
              my $slen=$hit->length;

              my $qcoverage=(($qend-$qstart)+1/$qlen)*100;
              my $scoverage=(($send-$sstart)+1/$slen)*100;
              
              $hit_desc=~s/(OS=)(.+)//g;

              if(($identity>=$identity_cutoff) and ($qcoverage>=$query_coverage_cutoff) and ($evalue <= $evalue_cutoff)){
                  $gene_annotation{$query}=$hit_desc;
              }
          }
       }
    }


############# FEATURE ##################

opendir(FEATURE,"$feature_dir/MULTIPLE");
my @feature_file=readdir(FEATURE);

my %verified_seq_id=();

foreach my $feature_file(@feature_file){

        my $genome_name=$feature_file;
           $genome_name=~s/(\.)(.+)//g;

        my @feature=();
        tie @feature, 'Tie::File', "$feature_dir/MULTIPLE/$feature_file";   
        
        my @feature_s=();
        tie @feature_s, 'Tie::File', "$feature_dir/SINGLE/$genome_name.feature.single.txt";  
        
        push(@feature,@feature_s);


        print "Verifying $genome_name sequences\n";

        open(V_FEATURE,">$verify_feature/$genome_name.txt");    

        foreach my $feature_line(@feature){

              if($feature_line=~/\#/){
                next;
              }
  
              my @feature_column=split("\t",$feature_line);
              
              my $seq_id=$feature_column[0];
              my $seq_len=$feature_column[9];

              if($gene_annotation{$seq_id}){
                 
                  my $gene_description=$gene_annotation{$seq_id};
 
                  pop(@feature_column);
                  push(@feature_column,$gene_description);

                  my $feature_line=join("\t",@feature_column);

                  print V_FEATURE "$feature_line\n"; 

                  $verified_seq_id{$seq_id}=$genome_name;
                
              }elsif($seq_len>=$seq_len_cutoff){
                
                  pop(@feature_column);
                  push(@feature_column,"hypothetical protein");

                  my $feature_line=join("\t",@feature_column);

                  print V_FEATURE "$feature_line\n"; 

                  $verified_seq_id{$seq_id}=$genome_name;
              }         
        }

        close V_FEATURE;
}


#### PROTEIN ####

opendir(PROTEIN,"$protein_dir/MULTIPLE");
my @protein_file=readdir(PROTEIN);

foreach my $protein_file(@protein_file){

    my $seqio_obj = Bio::SeqIO->new(-file => "$protein_dir/MULTIPLE/$protein_file", -format => "fasta");

    my $genome_name=$protein_file;
       $genome_name=~s/(\.)(.+)//g;

    open(V_PROTEIN,">$verify_protein/$genome_name.fasta");    

    while(my $seq_obj= $seqio_obj->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;

       my $gene_id=$seq_id;
       $gene_id=~s/(_1)$//g;

       if($gene_annotation{$gene_id} or length($seq)>=$seq_len_cutoff){
          print V_PROTEIN ">$seq_id\n$seq\n";
       }   
   }
   
   my $seqio_obj_s = Bio::SeqIO->new(-file => "$protein_dir/SINGLE/$genome_name.aa.single.fasta", -format => "fasta");
   
    while(my $seq_obj= $seqio_obj_s->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;

       my $gene_id=$seq_id;
       $gene_id=~s/(_1)$//g;

       if($gene_annotation{$gene_id} or length($seq)>=$seq_len_cutoff){
          print V_PROTEIN ">$seq_id\n$seq\n";
       }   
   }

   close V_PROTEIN;

}

#### CDS ###########

opendir(CDS,"$cds_dir/MULTIPLE");
my @cds_file=readdir(CDS);

foreach my $cds_file(@cds_file){

    my $seqio_obj = Bio::SeqIO->new(-file => "$cds_dir/MULTIPLE/$cds_file", -format => "fasta");

    my $genome_name=$cds_file;
       $genome_name=~s/(\.)(.+)//g;

    open(V_CDS,">$verify_cds/$genome_name.fasta");    

    while(my $seq_obj= $seqio_obj->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;

       my $aa_length=length($seq)/3;

       if($gene_annotation{$seq_id} or $aa_length>=$seq_len_cutoff){
          print V_CDS ">$seq_id\n$seq\n";
       }   
    }
    
    my $seqio_obj_s = Bio::SeqIO->new(-file => "$cds_dir/SINGLE/$genome_name.single.fasta", -format => "fasta");
    
    while(my $seq_obj= $seqio_obj_s->next_seq){
    
       my $seq=$seq_obj->seq;
       my $seq_id=$seq_obj->display_id;

       my $aa_length=length($seq)/3;

       if($gene_annotation{$seq_id} or $aa_length>=$seq_len_cutoff){
          print V_CDS ">$seq_id\n$seq\n";
       }   
    }

   close V_CDS;
}

##### GenBank Files #####

opendir(V_FEATURE,"$verify_feature");
my @verify_feature_file=readdir(V_FEATURE);

foreach my $feature_file(@verify_feature_file){

        my @feature=();
        tie @feature, 'Tie::File', "$verify_feature/$feature_file";   

        my $genome_name=$feature_file;
           $genome_name=~s/(\.)(.+)//g;
           
        my $genome_file="$genome_dir/$genome_name.fasta";   

        print "Creating GenBank file for $genome_name\n";

        system("perl CreateGBK.pl $genome_name $genome_file $verify_cds/$genome_name.fasta $verify_protein/$genome_name.fasta $verify_feature/$genome_name.txt $verify_genbank/$genome_name.gbk"); 

}


