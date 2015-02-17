##### Module to read configuration file #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package Configuration;
use strict;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(getConfig);


sub getConfig {

    my $config_file=(shift);
    my %ConfigParam=();

    open(CONFIG,$config_file) or die "Error: Cannot Open Configuration file\n";

    my $Section_key='';
    my $Param_key='';
    my $Param_value='';

    while(<CONFIG>) {
        
         if($_=~/^(\#)/ or $_=~/^\s+$/){         
            next;         
         }elsif($_=~/^(\[)(\w+)(\])/){
            $Section_key=$2;
            next;         
         }elsif($_=~/^(\-)(\w+)/){
              $_=~/^(\-)(\w+)(\:)(.+)/;
              $Param_key=$2;
              $Param_value=$4;
              $Param_value=~s/^\s+//g;                

              $ConfigParam{$Section_key}->{$Param_key}=$Param_value;              
         }
    }
    return(%ConfigParam);
}
