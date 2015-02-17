use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Path qw(remove_tree);
use Parallel::ForkManager;
use Hash::Merge qw( merge );

my($genome_name)=(shift);
my($query_sequence_file)=(shift);
my($annotation_database)=(shift);
my($significance_threshold)=(shift);


print "Annotating $genome_name protein sequences using SwissProt Database\n";

system("blastp -query $query_sequence_file -db $annotation_database -evalue 1e-5 -max_target_seqs 1 -outfmt 0 -out $out_dir/annotation/$genome_name.blas");

    my $report = Bio::SearchIO->new(-format => 'blast',
                                    -file =>"$out_dir/annotation/$genome_name.blas"); 
    
    my %gene_annotation=();

    while (my $result = $report->next_result) {

       my $query=$result->query_name;
          $query=~s/(_1)$//g;

       while( my $hit = $result->next_hit ) {   
          
          while( my $hsp = $hit->next_hsp ) {
          
              my $hit_desc=$hit->description; 
              my $identity=$hsp->percent_identity;
              my $qstart=$hsp->start('query');
              my $qend=$hsp->end('query');
              my $sstart=$hsp->start('hit'); 
              my $send=$hsp->end('hit');
              my $qlen=$result->query_length;
              my $slen=$hit->length;

              my $qcoverage=(($qend-$qstart)+1/$qlen)*100;
              my $scoverage=(($send-$sstart)+1/$slen)*100;

              my @desc=split("=",$hit_desc);   
              my $product_desc=$desc[0];
              $product_desc=~s/(OS)$//g;

              if(($identity>=30) and ($qcoverage>=50 or $scoverage>=50)){
                  $gene_annotation{$query}=$product_desc;
              }else{
                  $gene_annotation{$query}="unknown protein";
              }
          }
       }
   
       if(!defined($gene_annotation{$query})){
          $gene_annotation{$query}="unknown protein";
       }
    }
