##### Module to Perform MCL CLUSTERING#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package MclClustering;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(run_mcl cluster_hmm_super_family);


sub run_mcl {

    my($similarity_file)=(shift);
    my($mcl_dir)=(shift);
    my($mcl_output)='';

    print STDOUT "MCL CLUSTERING\n";
    ##### Prepare *.abc file ########
    system("cut -f 1,2,11 $mcl_dir/$similarity_file > $mcl_dir/mcxload_input.abc");

    ##### Replace 0 evalue with minimum non-zero evalue plus 1. i.e if min evalue is 1e-300 than replace 0 evalue with 1e-301 ###
    open(MCL_ABC,"$mcl_dir/mcxload_input.abc");
    my @mcl_abc=<MCL_ABC>;  
    my @arrayofarray=();
    foreach my $abc_line(@mcl_abc){
         chomp($abc_line);
         my @abc_column=split("\t",$abc_line);
         push(@arrayofarray,\@abc_column);
    } 
    my @sorted_arrayofarray=sort {$a->[2] <=> $b->[2]} @arrayofarray;
    my %evalue=();
    my @unique_evalue=grep{!$evalue{$_->[2]}++}@sorted_arrayofarray;
    $unique_evalue[1]->[2]=~/(e)(-)(\d+)/;   
    my $min_evalue=$3+1;        
    
    undef @arrayofarray;
    undef @mcl_abc;
    undef @sorted_arrayofarray;
    undef @unique_evalue;

    #### mcl clustering ######
    #chdir($mcl);
    system("mcxload -abc $Bin/$mcl_dir/mcxload_input.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil($min_evalue)' -o $Bin/$mcl_dir/mcl_input.mci -write-tab $Bin/$mcl_dir/mcl_input.tab");
    system("mcl $Bin/$mcl_dir/mcl_input.mci -I 1.5 -use-tab $Bin/$mcl_dir/mcl_input.tab -o $Bin/$mcl_dir/mcl_output");
    $mcl_output="mcl_output";
    #chdir($Bin);   
    
    return($mcl_output);
}


sub cluster_hmm_super_family{

    my($abc_file)=(shift);
    my($mcl)=(shift);
    my($mcl_dir)=(shift);
    my($mcl_output)='';

    system("mcl $abc_file --abc --abc-neg-log10 -abc-tf 'ceil(400)' -I 1.5  -o $Bin/$mcl_dir/mcl_output");
    $mcl_output="mcl_output";
    #chdir($Bin);   
    return($mcl_output);
}
