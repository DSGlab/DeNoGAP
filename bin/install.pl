#!/usr/bin/perl -w
###### ABOUT: This Script install the require perl modules and software to run the genomics pipeline ####
###### AUTHOR:Shalabh Thakur###################################################################

#### INSTRUCTIONS TO INSTALL THE PACKAGE #####

### Required Input: <path to home directory>, directory in which all required softwares should be installed #####
### Permissions: Installation of Perl module and some external program would need root permission ########
### USAGE: perl install.pl <complete path of home directory> ####

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use Getopt::Long;

my $HOME=undef;

GetOptions('install_dir=s'=>\$HOME);

if(!defined($HOME)){
   $HOME="../exe";
}

my @perl_std_module=("FindBin",
                     "Env",
                     "Exporter",
                     "Getopt::Long",
                     "File::Basename",
                     "File::Copy",
                     "Parallel::ForkManager",
                     "List::MoreUtils",
                     "List::Util",
                     "File::Path",
                     "Hash::Merge",
                     "DBI");

my @bioperl_module= ("Bio::Perl", "Bio::SeqIO", "Bio::Seq", "Bio::SearchIO", "Bio::Tools::Phylo::Phylip::ProtDist","Bio::AlignIO");

my @interproscan_module=("CGI",
                         "English",
                         "File::Spec::Functions",
                         "FileHandle",
                         "IO::Scalar",
                         "IO::String",
                         "Mail::Send",
                         "Sys::Hostname",
                         "URI::Escape",
                         "XML::Parser",
                         "XML::Quote");

my @iprscan_file=("ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/RELEASE/4.8/iprscan_v4.8.tar.gz",
                  "ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/BIN/4.x/iprscan_bin4.x_Linux64.tar.gz",
                  "ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/iprscan_DATA_43.1.tar.gz",
                  "ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/iprscan_MATCH_DATA_43.1.tar.gz",
                  "ftp://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/DATA/iprscan_PTHR_DATA_38.0.tar.gz");
 
my $iprscan_file=\@iprscan_file;

my %external_prog=('muscle'=>"http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz",
                   'kalign'=>"http://msa.sbc.su.se/downloads/kalign/current.tar.gz", 
                   'mcl'=>"http://micans.org/mcl/src/mcl-12-135.tar.gz", 
                   'hmmer'=>"http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz", 
                   'phylip'=>"http://evolution.gs.washington.edu/phylip/download/phylip-3.695.tar.gz", 
                   'glimmer'=>"http://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz", 
                   'prodigal'=>"http://prodigal.googlecode.com/files/prodigal.v2_60.linux", 
                   'fragscan'=>"http://omics.informatics.indiana.edu/mg/get.php?software=FragGeneScan1.16.tar.gz", 
                   'genemark'=>"http://opal.biology.gatech.edu",
                   'sqlite'=>"https://sqlite.org/2014/sqlite-autoconf-3080401.tar.gz", 
                   'emboss'=>"sudo apt-get install emboss"#,
                   'blast'=>"ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz",
                   'iprscan'=>$iprscan_file                   
                   );

my %uniprot_db=('swissprot'=>"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz");

#my %apache_server=('source_file'=>"http://mirror.csclub.uwaterloo.ca/apache//httpd/httpd-2.4.9.tar.bz2",
#                   'install_path'=>"");



print "Checking if you have standard perl modules installed\n";

foreach my $mod(@perl_std_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       } 
}

print "Checking if you have required Bio-Perl modules installed\n";

foreach my $mod(@bioperl_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       } 
}

print "Checking if you have required modules for Interproscan installed\n";

foreach my $mod(@interproscan_module){
       eval("use $mod");
       #if error or not installed
       if($@){
         print "$mod is not found on your system or path\n\n";
         print "Do you want to Install $mod now? (Y | N):";
         my $install_now=<STDIN>;   
           if($install_now=~/Y/i){
              eval {system("sudo cpan install $mod")};             
              warn("Error:".$!) if $@;             
           }
       }else{
         print "$mod already installed on your system\n\n";
       }
}

print "Downloading current release of Uniprot/SwissProt Database\n";

foreach my $db_name(keys %uniprot_db){

       print "Do you want to download current release of Uniprot/SwissProt database? (Y | N):";
       my $download_sp=<STDIN>; 

       if($download_sp eq "Y"){
          print "Creating database directory for swissprot at $Bin/../data/SPROT_DB\n";
          mkdir("$Bin/../data/SPROT_DB");
          print "Downloading database file for SwissProt\n";
          system("wget $uniprot_db{$db_name} -P $Bin/../data/SPROT_DB");
          system("chmod -R 777 $Bin/../data/SPROT_DB");
          system("tar -xvzf $Bin/../data/SPROT_DB/*.tar.gz -C $Bin/../data/SPROT_DB");          
       }   
}



print "Checking if you have required external program installed\n";

foreach my $ext_prog(keys %external_prog){
    ### Install MUSCLE ###
    if($ext_prog eq 'muscle'){
       eval{
        mkdir("$HOME/muscle");        
        system("wget $external_prog{$ext_prog} -P $HOME/muscle");
        system("chmod -R 777 $HOME/muscle");
        system("tar -xvzf $HOME/muscle/*.tar.gz -C $HOME/muscle");
        my @program_path=split(/\//,$external_prog{$ext_prog});
        my $program_name=pop(@program_path);
           $program_name=~s/(\.)(tar)(\.)(gz)//g;
           chomp($program_name);
        system("mv $HOME/muscle/$program_name $HOME/muscle/muscle");
        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/muscle)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/muscle\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>;
        }
      };
      if($@){
         die($!);
      }

    }
    #### INSTALL KALIGN #####
    elsif($ext_prog eq 'kalign'){
      eval{
        mkdir("$HOME/kalign");        
        system("wget $external_prog{$ext_prog} -P $HOME/kalign");
        system("chmod -R 777 $HOME/kalign");
        system("tar -xvzf $HOME/kalign/*.tar.gz -C $HOME/kalign");
        chdir("$HOME/kalign");
        system("./configure");
        system("make");
        system("sudo make install");
        chdir($Bin);
      };
      if($@){
         die($!);
      }
    }
    ##### INSTALL HMMER ######
    elsif($ext_prog eq 'hmmer'){
      eval{
        mkdir("$HOME/hmmer");        
        system("wget $external_prog{$ext_prog} -P $HOME/hmmer");
        system("chmod -R 777 $HOME/hmmer");
        system("tar -xvzf $HOME/hmmer/*.tar.gz -C $HOME/hmmer");

        chdir("$HOME/hmmer");    
        my $program_name=`ls -d */`;
        $program_name=~s/\///g;        
        chomp($program_name);

        chdir("$HOME/hmmer/$program_name");
        system("./configure --prefix=$HOME/hmmer");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/hmmer/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/hmmer/bin\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }
    }
    #### INSTALL MCL #####
    elsif($ext_prog eq 'mcl'){

       eval{
        mkdir("$HOME/mcl");        
        system("wget $external_prog{$ext_prog} -P $HOME/mcl");
        system("chmod -R 777 $HOME/mcl");
        system("tar -xvzf $HOME/mcl/*.tar.gz -C $HOME/mcl");
    
        chdir("$HOME/mcl");    
        my $program_name=`ls -d */`;
        $program_name=~s/\///g;        
        chomp($program_name);

        chdir("$HOME/mcl/$program_name");
        system("./configure --prefix=$HOME/mcl");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/mcl/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/mcl/bin\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }
    }
    #### INSTALL PHYLIP ####
    elsif($ext_prog eq 'phylip'){
       eval{
        mkdir("$HOME/phylip");        
        system("wget $external_prog{$ext_prog} -P $HOME/phylip");
        system("chmod -R 777 $HOME/phylip");
        system("tar -xvzf $HOME/phylip/*.tar.gz -C $HOME/phylip");

        chdir("$HOME/phylip");
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
        chomp($program_name);

        chdir("$HOME/phylip/$program_name/src");
        system("mv Makefile.unx Makefile");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/phylip/$program_name/exe)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/phylip/$program_name/exe\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }

    }
    #### INSTALL GLIMMER ####
    elsif($ext_prog eq 'glimmer'){

        eval{
        mkdir("$HOME/glimmer");        
        system("wget $external_prog{$ext_prog} -P $HOME/glimmer");
        system("chmod -R 777 $HOME/glimmer");
        system("tar -xvzf $HOME/glimmer/*.tar.gz -C $HOME/glimmer"); 
        chdir("$HOME/glimmer");        
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
        chomp($program_name);
        chdir("$HOME/glimmer/$program_name/src"); 
        system("make");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/glimmer/$program_name/bin)\n";
           print PROFILE "set path=(\$path $HOME/glimmer/$program_name/scripts)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/glimmer/$program_name/bin\n";
           print PROFILE "export PATH=\$PATH:$HOME/glimmer/$program_name/scripts\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }

    }
    ### INSTALL PRODIGAL ####
    elsif($ext_prog eq 'prodigal'){
       eval{
        mkdir("$HOME/prodigal");        
        system("wget $external_prog{$ext_prog} -P $HOME/prodigal");
        system("chmod -R 777 $HOME/prodigal");       
        my @program_path=split(/\//,$external_prog{$ext_prog});
        my $program_name=pop(@program_path);
        chomp($program_name);
        system("mv $HOME/prodigal/$program_name $HOME/prodigal/prodigal");

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/prodigal)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/prodigal\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }

   }
   ##### INSTALL FRAGSCAN #####
   elsif($ext_prog eq 'fragscan'){
       eval{
        mkdir("$HOME/frag_gene_scan");   
   
        print "Download FragGeneScan manually from $external_prog{$ext_prog}\n";
        print "Copy dowloaded FragGeneScan tar file into $HOME/frag_gene_scan directory and press enter to continue,";
        <>;     

        system("chmod -R 777 $HOME/frag_gene_scan");
        system("tar -xvzf $HOME/frag_gene_scan/*.tar.gz -C $HOME/frag_gene_scan");

        chdir("$HOME/frag_gene_scan");
        my $program_name=`ls -d */`;
           $program_name=~s/\///g;        
           chomp($program_name);
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/frag_gene_scan/$program_name)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/frag_gene_scan/$program_name\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }
   }
   #### INSTALL GENEMARK ####
   elsif($ext_prog eq 'genemark'){
        eval{
        mkdir("$HOME/genemark"); 

        print "Download genemark manually from $external_prog{$ext_prog}\n";
        print "Copy dowloaded genemark tar file into $HOME/genemark directory and press enter to continue,";
        <>;

        system("chmod -R 777 $HOME/genemark");
        system("tar -xvzf $HOME/genemark/*.tar.gz -C $HOME/genemark");
        chdir("$HOME/genemark"); 
        my $program_name=`ls -d */`;
        chomp($program_name);    

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $HOME/genemark/$program_name/gmsuite)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$HOME/genemark/$program_name/gmsuite\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      }
   }
   ##### INSTALL SQLITE ######
   elsif($ext_prog eq 'sqlite'){

       print "Do you want to Install SQLITE on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
        }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/sqlite";

        print "The default path of installation for SQLITE is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

        eval{

        mkdir($path);
        system("wget $external_prog{$ext_prog} -P $path");
        system("chmod -R 777 $path"); 
        system("tar -xvzf $path/*.tar.gz -C $path");
        chdir("$path");      
        my $program_name=`ls -d */`;
        chomp($program_name); 
        chdir("$path/$program_name");
        system("./configure --prefix=$path");
        system("make");
        system("make install");
        chdir($Bin);

        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $path/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$path/bin\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      } 

     }
   }
   ##### INSTALL EMBOSS ####
    elsif($ext_prog eq 'emboss'){

       print "Do you want to Install EMBOSS on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
        }

     if($option_install eq "Y" or $option_install eq "y"){
         print "Installing EMBOSS requires root user permission\n\n";
         system($external_prog{$ext_prog});
     }
   }
   ##### INSTALL NCBI BLAST ####
   elsif($ext_prog eq 'blast'){

       print "Do you want to Install LOCAL NCBI BLAST on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
        }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/blast";

        print "The default path of installation for  NCBI BLAST is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

        eval{
        mkdir($path);
        system("wget $external_prog{$ext_prog} -P $path");
        system("chmod -R 777 $path"); 
        system("tar -xvzf $path/*.tar.gz -C $path");   
        my $program_name=`ls -d */`;
        chomp($program_name); 
       
        if($ENV{SHELL}=~/\bt?csh/){
           open(PROFILE, ">>$ENV{HOME}/.cshrc") or warn("\tCould not open your .cshrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "set path=(\$path $path/$program_name/bin)\n";
           close PROFILE;
        }elsif($ENV{SHELL} =~ /\bbash/ ){
           open(PROFILE, ">>$ENV{HOME}/.bashrc") or warn("\tCould not open your .bashrc file.  Your will have to add this directory to your path by yourself.\n" );
           print PROFILE "export PATH=\$PATH:$path/$program_name/bin\n";
           close PROFILE;
        }else{
           print "\tI did not recognize the shell you are using.  I'm sorry, but you will have to add this directory to your path by yourself.\n";
           print "press [ENTER] to continue.";
           <>; 
        }
      };
      if($@){
         die($!);
      } 

     }
   }
   ####### INSTALL INTERPRO_SCAN ######
   elsif($ext_prog eq 'iprscan'){

       print "Do you want to Install LOCAL InterProScan on your system? ( Y or N): ";
       my $option_install=<STDIN>;
          chomp($option_install);

       if($option_install eq "N" or $option_install eq "n"){
          next;
       }

     if($option_install eq "Y" or $option_install eq "y"){

        my $path="$HOME/interpro_scan";

        print "The default path of installation for  InterProScan is [$path], Give new path if you want to change: ";
        my $option_change_path=<STDIN>;
           chomp($option_change_path);

        if($option_change_path ne ''){            
           $path=$option_change_path; 
        }

       eval{
        mkdir($path);
        my @iprscan_file=@{$external_prog{$ext_prog}};
        
        foreach my $file(@iprscan_file){
           $file=~s/\s//g;
           chomp($file);
           print "$file\n";
           # system("wget $file -P $path");
           system("chmod -R 777 $path"); 
           chdir("$path");
           my @program_path=split(/\//,$file);
           my $program_name=pop(@program_path);
           chomp($program_name);
           system("gunzip -c $program_name | tar -xvf -");
           chdir($Bin);     
        }
        chdir("$path/iprscan");
        system("perl Config.pl");
        chdir($Bin);
      };
      if($@){
         die($!);
      }
     }
   }     
}






