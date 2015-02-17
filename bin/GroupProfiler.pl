########################## PROGRAM TO CREATE PHYLOGENETIC PROFILE OF GROUP ABSENCE / PRESENCE IN STRAINS ################

#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util 'sum';


my %options= (ortholog_cluster=> undef,
              list_strain=> undef,
              out_dir=>undef,
              profile_file=>undef,
              locus_map_file=>undef);                     ###### Read Input cluster ortholog file name ####);

GetOptions('infile=s'=> \$options{ortholog_cluster},
           'list_strain=s'=> \$options{list_strain},
           'out_dir=s'=>\$options{out_dir},
           'profile_file=s'=>\$options{profile},
           'locus_map_file=s'=>\$options{locus_map});


 if(!(defined($options{ortholog_cluster}))){
   print "unable to open the cluster file: input valide cluster file using -infile <cluster file name> option\n";
 }

mkdir($options{out_dir});
 
open(GENOME,$options{list_strain});
my @genome=<GENOME>;

open(GROUP_PROF,">$options{profile}");

open(GROUP_ID,">$options{locus_map}");


open(IOG,"$options{ortholog_cluster}") or die "Err: $!\n";
my @read_iog_file=<IOG>;


my %ortholog_cluster=ReadOrthologGroups(\@read_iog_file);
  
CreateProfile(\%ortholog_cluster);

sub ReadOrthologGroups {

     my(@ortho_group)=@{(shift)};
     my %ortho_clust=();

     foreach my $clust(@ortho_group) {

	  my @clust_members=split(' ',$clust);
	  my $clust_id=shift(@clust_members);
	  #$clust_id=~s/\://g;
          if(scalar(@clust_members)>=1){
	         $ortho_clust{$clust_id}=\@clust_members;
          }
      }
      return %ortho_clust;
}
		
 sub CreateProfile {

       my(%OrthologCluster)=%{(shift)};
       my $title=0;

        print GROUP_PROF "GROUP\t";
		
        foreach my $genome(@genome) {
	       $genome=~s/\s//g;
	       print GROUP_PROF "$genome\t";
        }
	    print GROUP_PROF "\n";
       
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
	   
           	     foreach my $ortho_seq(@{$value}) {
		       
		               $ortho_seq=~/(\w+)(\|)(\w+)/;

                       if(!defined($profile{$1})){
                         next;
                       }

		               $profile{$1}=1;
                       $frequency{$1}=$frequency{$1}+1;   
		               print GROUP_ID "$key\t$3\t$1\n";
                  }
             
                  $num_genome=sum values(%profile);
                  $total_gene=scalar(@gene_id);

                  if($num_genome==0){
                    next;
                  }
                  
                  print GROUP_PROF "$key\t";
	              print "$key\t";

                  $cluster_number++;
			
	          foreach my $genome(@genome) {
		       chomp($genome);
		       $genome=~s/\s//g;				 
				 
		       print GROUP_PROF "$profile{$genome}\t";
		       print "$profile{$genome}\t";
		  }	
	           print GROUP_PROF "\n";
        }
     close GROUP_PROF;
     close GROUP_ID;
 }
	

