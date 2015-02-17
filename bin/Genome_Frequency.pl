########################## PROGRAM TO CREATE PHYLOGENETIC PROFILE OF GROUP ABSENCE / PRESENCE IN STRAINS ################

#!/usr/bin/perl
use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use List::Util 'sum';
use Cluster;



my %options= (ortholog_cluster=> undef,
              seq_dir=> undef,
              list_genome=>undef,
              out_dir=>undef);                     ###### Read Input cluster ortholog file name ####);

GetOptions('infile=s'=> \$options{ortholog_cluster},
           'seq_dir=s'=>\$options{seq_dir},
           'list_genome=s'=>\$options{list_genome},
           'out_dir=s'=>\$options{out_dir});


 if(!(defined($options{ortholog_cluster}))){
   print "unable to open the cluster file: input valide cluster file using -infile <cluster file name> option\n";
 }

### Create output directory for group profile ####
mkdir($options{out_dir});

### read sequence directory
opendir(SEQFILE,$options{seq_dir});
my @readfilename=readdir(SEQFILE);

### create strain-list from sequence files ####
open(GENOME,">$options{list_genome}");
my %strain_name=();

foreach my $strain_name(@readfilename){     
   if($strain_name=~/^\.+/ or $strain_name=~/\~/){
      next;
   }
$strain_name=~s/(\.)(\w+)//g;
$strain_name{$strain_name}=$strain_name;
print GENOME "$strain_name\n";
}
close GENOME; 

#### read genome list ####
open(GENOME,$options{list_genome});
my @genome=<GENOME>;

### read ortholog cluster #####

my($ortholog_cluster,$status)=Cluster::read_cluster($options{ortholog_cluster},'xxx');
my %ortholog_cluster=%{$ortholog_cluster};

##### create profile #####
  
CreateProfile(\%ortholog_cluster,\@genome,$options{out_dir});

####subroutine to create profile ####

sub CreateProfile {

    my(%OrthologCluster)=%{(shift)};
    my(@genome)=@{(shift)};
    my $output_dir=(shift);

    open(GROUP_PROF,">$output_dir/GroupProfile.txt");
    open(GROUP_ID,">$output_dir/ListGroupLocus.txt");
    open(GENE_FREQ,">$output_dir/GenomeFrequency.txt");

    print GROUP_PROF "GROUP\t";
    print GENE_FREQ "cluster\t";
    print GENE_FREQ "genomes\t";
    print GENE_FREQ "total\t";

    #### print out the genome name on horizontal axis in profile ###		
    foreach my $genome(@genome) { 
       $genome=~s/\s//g;
       print GROUP_PROF "$genome\t";
       print GENE_FREQ "$genome\t";
    }
    print GROUP_PROF "\n";
    print GENE_FREQ "\n";
		
    #### print profile for each ortholog group ####
    my $cluster_number=1;
    my $num_genome=0;
    my $total_gene=0;

    while(my($key,$value)=each(%OrthologCluster)) {
		
         my %profile=();
         my %frequency=();
         my @gene_id=split("\t",$value);
			  
         foreach my $genome(@genome) {	
             $genome=~s/\s//g;			  
             $profile{$genome}=0;
             $frequency{$genome}=0;				   
	 }
	   
         foreach my $ortho_seq(@gene_id) {
            $ortho_seq=~/(\w+)(\|)(.+)/;
	    $profile{$1}=1;
            $frequency{$1}=$frequency{$1}+1; 		  
	    print GROUP_ID "$key\t$3\t$1\t".scalar(@gene_id)."\n";
         }
         
         print GROUP_PROF "$key\t";
         print GENE_FREQ "$cluster_number\t";

         $num_genome=sum values(%profile);
         $total_gene=scalar(@gene_id);

         print GENE_FREQ "$num_genome\t";
         print GENE_FREQ "$total_gene\t";
         
         $cluster_number++;         	 
			
	 foreach my $genome(@genome) {
	     chomp($genome);
             $genome=~s/\s//g;				 
				 
             print GROUP_PROF "$profile{$genome}\t";
             print GENE_FREQ "$frequency{$genome}\t";

	     print "$profile{$genome}\t";
	 }			
	 print GROUP_PROF "\n";
         print GENE_FREQ "\n";        	 
   }

   close GROUP_PROF;
   close GENE_FREQ;
   close GROUP_ID;
}
	

