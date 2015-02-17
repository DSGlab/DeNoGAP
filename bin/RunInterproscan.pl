##### This script run clonalorigin on each alignment block on scinet ########
use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use SQLiteDB;
use Parallel::ForkManager;

my $genome_name=(shift);
my $db_dir=(shift);
my $db_name=(shift);
my $xml_dir=(shift);
my $tab_dir=(shift);
my $parsed_domain_dir=(shift);
my $parsed_interpro_dir=(shift);
my $parsed_go_dir=(shift);
my $parsed_pathway_dir=(shift);
my $interpro_path=(shift);
my $interpro_opts=(shift);
my $process=(shift);
my $tmp_dir=(shift);


  my $get_sequence_sql="Select * from ProteinSequence where genome_abbreviation='$genome_name' LIMIT 50";
          
  my $query_protein_sequences=SQLiteDB::get_record($db_dir,$db_name,$get_sequence_sql);
  
  open(CHUNK_SEQ,">$tmp_dir/$genome_name.fasta");
  
   if(scalar(@{$query_protein_sequences})>0){
          
        foreach my $row(@{$query_protein_sequences}){
              
             my $protein_index=shift(@{$row});
             my $protein_id=shift(@{$row});
             my $genome_abbreviation=shift(@{$row});
             my $seq_type=shift(@{$row});
             my $seq_length=shift(@{$row});
             my $seq=shift(@{$row});

             $seq=~s/\*$//g;
                                            
             print CHUNK_SEQ ">$genome_abbreviation|$protein_id\n$seq\n";                           
        }
  }
  close CHUNK_SEQ;


   ##### Pass the command to execute at each core #####
   system("bash $interpro_path/interproscan.sh -i $tmp_dir/$genome_name.fasta $interpro_opts -f XML -o $xml_dir/$genome_name.xml");

   system("bash $interpro_path/interproscan.sh -mode convert -i $xml_dir/$genome_name.xml -f tsv -o $tab_dir/$genome_name.tab");
   
