#!/usr/bin/perl -w
###### ABOUT: This Script is the main controller for the execution of MGAT Package ############
###### AUTHOR:Shalabh Thakur###################################################################
###### DATE:15-MAY-2013########################################################################

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Configuration;
use AdjustSequence;
use SequenceHash;
use Hmmer;
use CompareReference;
use FilterPair;
use CreateModel;
use PredictGene;
use HomologScan;
use Getopt::Long;
use SQLiteDB;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use File::Basename;
use File::Copy;
use File::Path qw(remove_tree);
use Parallel::ForkManager;
use Hash::Merge qw( merge );
use ParseGBK;

#### GENERAL INPUT  PARAMETERS ###

#### MANDATORY PARAMTERS #########
my $genome_info_file=undef;
my $db_dir=undef;
my $db_name=undef;
my $configuration_file=undef;
my $output_dir=undef;
my $help=undef;

##### READ COMMAND LINE INPUT VARIABLES #####

GetOptions('genome_info=s'=>\$genome_info_file,
           'db_dir=s'=>\$db_dir,
           'db_name=s'=>\$db_name,
           'config=s'=>\$configuration_file,
           'output_dir=s'=>\$output_dir,
           'help=s'=>\$help
          );


##### DECLARE GLOBAL VARIABLE #####
my %config_param=();
my $tmp_dir="$output_dir/tmp";


###### Check Mandatory inputs #####
if(!defined($genome_info_file)){
   print "ERROR: Undefined Genome metadata file\n";
   exit;
}

if(!defined($db_dir) or !defined($db_name)){
   print "ERROR: Undefined database directory or database name\n";
   exit;
}

if(!defined($output_dir)){
   print "ERROR: Undefined output directory\n";
   exit;
}

if(!defined($configuration_file)){
   print "ERROR: Undefined configuration file\n";
   exit;
}

print "CREATING OUTPUT DIRECTORY $output_dir\n";

mkdir($output_dir);
mkdir($tmp_dir);


##### BUILD SQLITE DATABASE STRUCTURE ######
print "Create SQLITE Database\n";

SQLiteDB::create_db($db_dir,$db_name);

print "Create SQLITE Tables\n";

SQLiteDB::create_table($db_dir,$db_name);



 ####### ADD GENOME INFORMATION #########
 
if($genome_info_file){

      print "Loading Genome Information\n"; 

      unless(-s "$genome_info_file"){
         print "ERROR: genome information file is empty or not found\n";
         exit;
      }

      my($insert_column,$bind_value)=SQLiteDB::create_organism_table($db_dir,$db_name,$genome_info_file);

      open(FILE,"$genome_info_file");
      my @file=<FILE>;
      close FILE;    
        
      my $organism_sql_stmt="INSERT into OrganismInfo ($insert_column) VALUES ($bind_value)";
      SQLiteDB::load_from_array($db_dir,$db_name,\@file,$organism_sql_stmt);
}



##### READ CONFIGURATION FILE #######
print STDOUT "READ CONFIGURATION FILE\n";

if($configuration_file){

      unless(-s "$configuration_file"){
         print "ERROR: configuration file is empty or not found\n";
         exit;
      }
  
  %config_param=Configuration::getConfig($configuration_file);
}


############################ Perform Analysis #############################################

####### PARSE GENE BANK FILE ######

if($config_param{ACTIVATE_ANALYSIS}->{PARSE_GENBANK} and $config_param{ACTIVATE_ANALYSIS}->{PARSE_GENBANK}=~/^YES$/){

   if($config_param{INPUT_DIRECTORY}->{GENBANK_DIR_PATH}){

        print "Extract Sequence Information from GeneBank File\n";  
   
        opendir(GBK_DIR,"$config_param{INPUT_DIRECTORY}->{GENBANK_DIR_PATH}") or die "Cannot open directory for genebank files at $config_param{INPUT_DIRECTORY}->{GENBANK_DIR_PATH}";
        my @gbk_file=readdir(GBK_DIR);
        close GBK_DIR;

        my %gbk_info=();

        #### output directory for to store parsed genbank information ###
        print "CREATING PROJECT DIRECTORIES FOR PARSED GENBANK RESULTS\n";
        
        my %result_dir=%{$config_param{OUTPUT_DIRECTORY}};
        
        my $project_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}";
        my $genome_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{GENOME_DIR_NAME}";
        my $cds_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{CDS_DIR_NAME}";
        my $protein_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{PROTEIN_DIR_NAME}";
        my $feature_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{FEATURE_DIR_NAME}";
        
        mkdir("$project_dir");
        mkdir("$genome_dir");
        mkdir("$cds_dir");
        mkdir("$protein_dir");
        mkdir("$feature_dir");
        

        foreach my $gbk_file(@gbk_file){

              if($gbk_file=~/\~/ or $gbk_file=~/^\.+$/){
              
                next;
                
              }else{
              
                 my $genome_name=$gbk_file;
                 
                 $genome_name=~s/\.(\w+)//g;   #### remove extension from the file name ###

                 $gbk_file=$config_param{INPUT_DIRECTORY}->{GENBANK_DIR_PATH}."/".$gbk_file;

                 my($genome_seq,$coding_seq,$protein_seq,$feature)=ParseGBK::gbk_genome($gbk_file,$genome_name); 
                   
                 ParseGBK::print_genome_seq_file($genome_seq,$genome_name,$genome_dir);
                 
                 ParseGBK::print_coding_seq_file($coding_seq,$protein_seq,$genome_name,$cds_dir,$protein_dir);
                 
                 ParseGBK::print_feature_file($feature,$genome_name,$feature_dir);
              }           
         }
    }else{
        print "ERROR: genbank file directory not defined\n";
        exit;
    }

  exit;
} #### END OF PARSE GENE BANK ANALYSIS BLOCK ######


########### PREDICT GENE ##########################

if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_GENE} and $config_param{ACTIVATE_ANALYSIS}->{PREDICT_GENE}=~/YES/i){

        print "PREDICITNG GENE\n";
        
        ###### INPUT DIRECTORY ########
        
        my $genome_dir=$config_param{INPUT_DIRECTORY}->{GENOME_DIR_PATH};
        chomp($genome_dir);
        
        ###### OUTPUT_DIRECTORY #######
        
        my %result_dir=%{$config_param{OUTPUT_DIRECTORY}};
        
        my $project_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}";
        
        my $glimmer_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{GLIMMER_RESULT_DIR_NAME}";
        my $genemark_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{GENEMARK_RESULT_DIR_NAME}";
        my $prodigal_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{PRODIGAL_RESULT_DIR_NAME}";
        my $fragscan_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{FRAGSCAN_RESULT_DIR_NAME}";
        
        my $cds_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{CDS_DIR_NAME}";
        my $protein_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{PROTEIN_DIR_NAME}";
        my $feature_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}/$result_dir{FEATURE_DIR_NAME}";
        
        
        ##### OTHER PARAMETERS ######
        my $translation=$config_param{PARAMETERS}->{TRANSLATION_CODE};
        my $overlap_base=$config_param{PARAMETERS}->{OVERLAP_BASE};
        my $parallel_core=$config_param{PARAMETERS}->{PARALLEL_CPU_CORE};
        my $add_annotation=$config_param{ANNOTATION}->{ADD_ANNOTATION};
        my $significance_threshold=$config_param{ANNOTATION}->{EVALUE_THRESHOLD};
        my $annotation_database=$config_param{ANNOTATION}->{ANNOTATION_DATABASE};
        
   
        ##### CREATE OUTPUT DIRECTORIES FOR GENE PREDICTION ######     
        mkdir($project_dir);
        mkdir($glimmer_dir);
        mkdir($genemark_dir);
        mkdir($prodigal_dir);
        mkdir($fragscan_dir);
        mkdir($cds_dir);
        mkdir($protein_dir);
        mkdir($feature_dir);
        

        ##### EXTERNAL PROGRAM PARAMETERS #########################
        my $glimmer3_opt=$config_param{GLIMMER_OPTIONS}->{GLIMMER3};
        my $long_orf_opt=$config_param{GLIMMER_OPTIONS}->{LONG_ORF};
        my $multi_extract_opt=$config_param{GLIMMER_OPTIONS}->{MULTI_EXTRACT};
        my $build_icm_opt=$config_param{GLIMMER_OPTIONS}->{BUILD_ICM};
        my $genemark_opt=$config_param{GENEMARK_OPTIONS}->{GMSN};
        my $prodigal_opt=$config_param{PRODIGAL_OPTIONS}->{PRODIGAL};
        my $fragscan_opt=$config_param{FRAGSCAN_OPTIONS}->{FRAGSCAN};
        

        ##### OPEN AND READ GENOME FILES FROM INPUT DIRECTORY #####
        
        opendir(GENOME,$genome_dir) or die "cannot open genome sequence directory $genome_dir\n";
        my @genome_file=readdir(GENOME);
        
        foreach my $genome_file(@genome_file){

            if($genome_file=~/^\.+$/ or $genome_file=~/\~/){
              next;
            }

            my $genome_name=$genome_file;
            $genome_name=~s/\.(.+)//g;
            
            $genome_file=$genome_dir."/".$genome_file;

            ##### get full name for each genome in the database #####
            my $get_full_name="Select DISTINCT(genome_name) from OrganismInfo WHERE abbreviation='$genome_name'";
            my($record_full_name)=SQLiteDB::get_record($db_dir,$db_name,$get_full_name);
            
            my $full_genome_name='';

            if(scalar(@{$record_full_name})>0){
            
              foreach my $row(@{$record_full_name}){  
                 $full_genome_name=shift(@{$row});
              }
              
            }else{        
               die "Cannot find full name for $genome_name in the database\n";
            }
            
            system("perl PredictGene.pl $genome_name '$full_genome_name' $genome_file $project_dir $glimmer_dir $genemark_dir $prodigal_dir $fragscan_dir $genome_dir $cds_dir $protein_dir $feature_dir $translation $overlap_base $glimmer3_opt $long_orf_opt $multi_extract_opt $build_icm_opt $genemark_opt $prodigal_opt $fragscan_opt");       
       } 
  exit;  
}  #### END OF GENE PREDICTION BLOCK #####


############ VERIFY GENE PREDICTION USING SWISSPROT #########################

if($config_param{ACTIVATE_ANALYSIS}->{VERIFY_SEQUENCE} and $config_param{ACTIVATE_ANALYSIS}->{VERIFY_SEQUENCE}=~/^YES$/i){


      print "VERIFY SEQUENCES\n";

      my $blast_alignment_file=$config_param{INPUT}->{BLAST_ALIGNMENT_FILE};

      my $feature_dir=$config_param{INPUT}->{FEATURE_DIR};
      my $cds_dir=$config_param{INPUT}->{CDS_DIR};
      my $protein_dir=$config_param{INPUT}->{PROTEIN_DIR};
      my $genome_dir=$config_param{INPUT}->{GENOME_DIR};

      my %result_dir=%{$config_param{OUTPUT_DIRECTORY}};

      my $project_dir="$output_dir/$result_dir{PROJECT_DIR_NAME}";
      my $verify_cds="$project_dir/$result_dir{VERIFY_CDS}";
      my $verify_protein="$project_dir/$result_dir{VERIFY_PROTEIN}";
      my $verify_feature="$project_dir/$result_dir{VERIFY_GENOME_FEATURE}";
      my $verify_genbank="$project_dir/$result_dir{VERIFY_GENBANK_FILE}";

      my $evalue_cutoff=$config_param{FILTER_PARAMETERS}->{EVALUE_THRESHOLD};
      my $identity_cutoff=$config_param{FILTER_PARAMETERS}->{ALIGNMENT_IDENTITY};
      my $query_coverage_cutoff=$config_param{FILTER_PARAMETERS}->{QUERY_COVERAGE};
      my $seq_len_cutoff=$config_param{FILTER_PARAMETERS}->{MIN_PROTEIN_LENGTH};
  
      system("perl VerifySequence.pl $blast_alignment_file $feature_dir $cds_dir $protein_dir $genome_dir $project_dir $verify_cds $verify_protein $verify_feature $verify_genbank $evalue_cutoff $identity_cutoff $query_coverage_cutoff $seq_len_cutoff");

      print "VERIFICATION COMPLETE\n";
 exit;
}  ### END OF VERIFICATION BLOCK ####

############ LOAD GENOMIC DATA INTO THE DATABASE ############


if($config_param{ACTIVATE_ANALYSIS}->{LOAD_DATA} and $config_param{ACTIVATE_ANALYSIS}->{LOAD_DATA}=~/^YES$/i){

  my $protein_seq_dir=$config_param{INPUT_DIRECTORY}->{PROTEIN_DIR};
  my $coding_seq_dir=$config_param{INPUT_DIRECTORY}->{NUCLEOTIDE_DIR};
  my $genome_feature_dir=$config_param{INPUT_DIRECTORY}->{FEATURE_DIR};
 
  my $adjust_header=$config_param{PARAMETER}->{ADJUST_HEADER};


  my $get_genome_in_db="SELECT DISTINCT(abbreviation),genome_type from OrganismInfo";
  my ($genome_record)=SQLiteDB::get_record($db_dir,$db_name,$get_genome_in_db);
     
  my %genome_info=();

  if(scalar($genome_record)>1){           
     foreach my $row(@{$genome_record}){           
         my $abbrv=shift(@{$row});
         my $genome_type=shift(@{$row});                        
         $genome_info{$abbrv}=$genome_type;
     }
  }
  

  ##### get features stored in database #######

  my $get_feature_in_db="SELECT DISTINCT(feature_id),genome_name from GeneFeature";

  my ($feature_record)=SQLiteDB::get_record($db_dir,$db_name,$get_feature_in_db);

  my %feature_info=();

  if(scalar(@{$feature_record})>1){ 
          
     foreach my $row(@{$feature_record}){           
             
             my $feature_id=shift(@{$row});
             my $genome_name=shift(@{$row});  
                             
             $feature_info{$genome_name}->{$feature_id}=$genome_name;
     }
  }

  ##### READ FEATURE FILE DIRECTORY #############

   opendir(FEATURE,"$genome_feature_dir");
   my @feature_file=readdir(FEATURE);

   my @feature=();

      foreach my $feature_file(@feature_file){ 

            if($feature_file=~/^\.+$/){
               next; 
            }

            print "add $feature_file\n"; 
          
            open(FILE,"$genome_feature_dir/$feature_file");
            my @file=<FILE>;
            close FILE; 

            @feature=(@feature,@file);             
        }

   ###### REMOVE DUPLICATE FEATURE ENTRIES FROM DATABASE ##### 
    my @feature_load=();

    for(my $i=0;$i<scalar(@feature);$i++){           
        
        chomp($feature[$i]);
        my @column=split("\t",$feature[$i]);                         
        unless($feature_info{$column[4]}->{$column[0]}){               
            push(@feature_load,$feature[$i]); 
        } 
    }  
        
    print "Loading User specified Gene Feature Information\n";      
    
    my $gene_feature_sql_stmt="INSERT into GeneFeature (feature_id, feature_type, genome_id, genome_type, genome_name, genome_length, feature_start, feature_end, nuc_length, aa_length, strand, frame, index_in_genome, description) VALUES 
                               (?, ?, ?, ? ,? ,? ,?, ? ,? ,? ,? ,? ,?, ?)";
                               
    SQLiteDB::load_from_array($db_dir,$db_name,\@feature_load,$gene_feature_sql_stmt);


    ###### LOAD PROTEIN SEQUENCE IN THE DATABASE #########################################
    
    opendir(Pseq_dir,$protein_seq_dir) or die "Cannot open directory\n";
    
    my @protein_file=readdir(Pseq_dir);

    my($sequence_array)=AdjustSequence::checkSequence($protein_seq_dir,$adjust_header,\@protein_file,"protein");

    my $get_protein_seq_in_db="SELECT DISTINCT(protein_id),genome_abbreviation from ProteinSequence";
    
    my ($protein_seq_record)=SQLiteDB::get_record($db_dir,$db_name,$get_protein_seq_in_db);

    my %protein_seq_info=();

    if(scalar($protein_seq_record)>1){
           
       foreach my $row(@{$protein_seq_record}){
           
           my $feature_id=shift(@{$row});
           my $genome_abbrv=shift(@{$row});             
           $protein_seq_info{$genome_abbrv."|".$feature_id}=$feature_id;
       }
    }

    my @seq_array=@{$sequence_array};
    my @seq_array2=();

    for(my $i=0;$i<scalar(@seq_array);$i++){   
        
        chomp($seq_array[$i]);
        my @column=split("\t",$seq_array[$i]);              
        
        unless($protein_seq_info{$column[1]."|".$column[0]}){   
          
           if(!defined($genome_info{$column[1]})){

             die "Cannot find organism with abbreviation $column[1] in the organism database\n";
           }  
          
           push(@seq_array2,$seq_array[$i]);             
        }
    }   
     
    my $seq_load_stmt="INSERT INTO ProteinSequence(protein_id, genome_abbreviation, seq_type, seq_length, aminoacid_sequence) VALUES (?, ?, ?, ?, ?)";
   
    SQLiteDB::load_from_array($db_dir,$db_name,\@seq_array2,$seq_load_stmt);   
    
    
    ###### LOAD NUCLEOTIDE SEQUENCE IN THE DATABASE #########################################
    
    print STDOUT "ADDING NUCLEOTIDE SEQUENCES TO DATABASE\n";

    opendir(Nseq_dir,$coding_seq_dir) or die "Cannot open directory\n";
   
    my @cds_file=readdir(Nseq_dir);

    my($sequence_array_nuc)=AdjustSequence::checkSequence($coding_seq_dir,$adjust_header,\@cds_file,"nucleotide");

    ##### Check and Remove duplicate nucleotide sequence #######
    my $get_nuc_seq_in_db="SELECT DISTINCT(nucleotide_id),genome_abbreviation from NucleotideSequence";
    
    my ($nuc_seq_record)=SQLiteDB::get_record($db_dir,$db_name,$get_nuc_seq_in_db);
    
    my %nuc_seq_info=();
    
    my @seq_array_nuc=@{$sequence_array_nuc};
    
    my @seq_array_nuc2=();

    if(scalar($nuc_seq_record)>1){           
           
       foreach my $row(@{$nuc_seq_record}){        
          
            my $feature_id=shift(@{$row});
            my $genome_abbrv=shift(@{$row});             
            $nuc_seq_info{$genome_abbrv."|".$feature_id}=$feature_id;
       }
    }

    for(my $i=0;$i<scalar(@seq_array_nuc);$i++){  
             
            chomp($seq_array_nuc[$i]);
        
            my @column=split("\t",$seq_array_nuc[$i]);
                          
            unless($nuc_seq_info{$column[1]."|".$column[0]}){
            
                push(@seq_array_nuc2,$seq_array_nuc[$i]);             
            }
    }  
    
    my $nuc_seq_load_stmt="INSERT INTO NucleotideSequence(nucleotide_id, genome_abbreviation, seq_type, seq_length, nucleotide_sequence) VALUES (?, ?, ?, ?, ?)";
    
    SQLiteDB::load_from_array($db_dir,$db_name,\@seq_array_nuc2,$nuc_seq_load_stmt);  

  exit;    
}  ##### END OF LOAD DATA BLOCK #################


######## RUN HOMOLOG SCAN PHASES 1: COMPARE REFERENCE ########

if($config_param{ACTIVATE_ANALYSIS}->{COMPARE_REFERENCE} and $config_param{ACTIVATE_ANALYSIS}->{COMPARE_REFERENCE}=~/^YES$/){

     ######## create output directory for homolog search #####
     my $project_name=$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME};
     my $analysis="compare_reference";
     
     mkdir("$output_dir/$project_name");

     my $homolog_dir_hash =BuildOutputDirectory($project_name,$analysis,$output_dir,$tmp_dir);
          
     HomologScan::run_compare_reference($db_dir,$db_name,\%config_param,$homolog_dir_hash,$output_dir,$tmp_dir); 

  exit;               
}  ##### END OF REFERENCE COMPARISON BLOCK ####


######## RUN HOMOLOG SCAN PHASES 1: ITERATIVE PREDICTION OF HMM FAMILY ########

if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_HMM} and $config_param{ACTIVATE_ANALYSIS}->{PREDICT_HMM}=~/^YES$/i){

     ######## create output directory for homolog search #####
     my $project_name=$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME};
     my $analysis="predict_hmm";
     
     mkdir("$output_dir/$project_name");

     my $homolog_dir_hash = BuildOutputDirectory($project_name,$analysis,$output_dir,$tmp_dir);
          
     HomologScan::run_homolog_scan($db_dir,$db_name,\%config_param,$homolog_dir_hash,$output_dir,$tmp_dir); 
  
  exit;               
} #### END OF HMM FAMILY PREDICTION BLOCK #####

####### RUN HOMOLOG SCAN PHASE 3: ORTHOLOG PREDICTION ################

if(($config_param{ACTIVATE_ANALYSIS}->{PREDICT_ORTHOLOG} and $config_param{ACTIVATE_ANALYSIS}->{PREDICT_ORTHOLOG}=~/^YES$/i) or ($config_param{ACTIVATE_ANALYSIS}->{CLUSTER_ORTHOLOG} and $config_param{ACTIVATE_ANALYSIS}->{CLUSTER_ORTHOLOG}=~/^YES$/i)){

     ######## create output directory for homolog search #####
     my $project_name=$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME};
     my $analysis="predict_ortholog";
     
     mkdir("$output_dir/$project_name");

     my $homolog_dir_hash = BuildOutputDirectory($project_name,$analysis,$output_dir,$tmp_dir);
          
     HomologScan::run_ortholog_prediction($db_dir,$db_name,\%config_param,$homolog_dir_hash,$output_dir,$tmp_dir);  

  exit;              
} ###### END OF ORTHOLOG PREDICTION BLOCK #####


####################### PHYLOGENETIC PROFILE #############

if($config_param{ACTIVATE_ANALYSIS}->{PHYLOGENETIC_PROFILE} and $config_param{ACTIVATE_ANALYSIS}->{PHYLOGENETIC_PROFILE}=~/YES/i){

           my $project_name="$output_dir/$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME}";
           my $cluster_file=$config_param{GROUP}->{ORTHOLOG_CLUSTER_FILE};
           
           my $profile="$project_name/GroupProfile.txt";
           my $genetofamily="$project_name/ListGroupLocus.txt";
        
           print "Creating Phylogenetic Profile for ortholog groups\n";

           my $get_genome_in_db="SELECT DISTINCT(abbreviation)from OrganismInfo";
           my ($genome_record)=SQLiteDB::get_record($db_dir,$db_name,$get_genome_in_db);
           
           my @genome=();
           
           open(LIST_GENOME,">$tmp_dir/list_strain.txt");

           if(scalar($genome_record)>1){  
                    
             foreach my $row(@{$genome_record}){ 
                       
                my $abbrv=shift(@{$row});
                push(@genome,$abbrv);   
                
                print LIST_GENOME "$abbrv\n";                                 
             }
           }
           
           close LIST_GENOME;
          
          system("perl GroupProfiler.pl -infile $cluster_file -list_strain $tmp_dir/list_strain.txt -out_dir $project_name -profile_file $profile -locus_map_file $genetofamily");

          my($insert_species,$bind_value)=SQLiteDB::create_profile_table($db_dir,$db_name,\@genome);
  
          print "Load Profile\n";

          my $profile_insert_stmt="INSERT INTO Profile(id,$insert_species) VALUES (?,$bind_value)";

          SQLiteDB::load_data($db_dir,$db_name,$profile,$profile_insert_stmt);

          print "Load GenetoFamily Map\n";

          my $genetofamily_map_stmt="INSERT INTO MapGeneIdtoGeneFamily(genefamily_id, gene_id, species_abbreviation) VALUES (?, ?, ?)";

          SQLiteDB::load_data($db_dir,$db_name,$genetofamily,$genetofamily_map_stmt);

  exit;
}  #### END OF PHYLOGENETIC PROFILE BLOCK #####


###### CORE GENOME ALIGNMENT ######

if($config_param{ACTIVATE_ANALYSIS}->{CORE_GENOME} and $config_param{ACTIVATE_ANALYSIS}->{CORE_GENOME}=~/YES/i){
       
         print "Predict Core Genome\n"; 
       
         my $project_name=$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME};
         my $cluster_file=$config_param{GROUP}->{ORTHOLOG_CLUSTER_FILE};
         my $core_seq_dir="$output_dir/$project_name/CORE_SEQ";
         my $core_aln_dir="$output_dir/$project_name/CORE_ALN";
         my $core_alignment_file=$config_param{PARAMETERS}->{CORE_ALIGNMENT_FILE};
         my $core_threshold=$config_param{PARAMETERS}->{CORE_THRESHOLD};
         my $sequence_type=$config_param{PARAMETERS}->{SEQUENCE_TYPE};
         my $include_outgroup=$config_param{PARAMETERS}->{INCLUDE_OUTGROUP};
         
         mkdir("$output_dir/$project_name");
         mkdir($core_seq_dir);
         mkdir($core_aln_dir);
        
         system("perl CoreGenome.pl $cluster_file $db_dir $db_name $core_seq_dir $core_aln_dir $core_alignment_file $core_threshold $sequence_type $include_outgroup");
   exit;      
}  ###### END OF CORE GENOME PREDICTION BLOCK #####

##### Activate analysis for Annotation #####
if($config_param{ACTIVATE_ANALYSIS}->{PREDICT_ANNOTATION} and $config_param{ACTIVATE_ANALYSIS}->{PREDICT_ANNOTATION}=~/YES/i){  
    
         print "Predict Annotation using InterProScan\n";
         
         ##### OUTPUT DIRECTORY #####
         my $project_name=$config_param{OUTPUT_DIRECTORY}->{PROJECT_DIR_NAME};
         my $interpro_scan_out_dir="$output_dir/$project_name/INTERPRO_SCAN_OUTPUT";
         my $xml_result_dir="$interpro_scan_out_dir/XML_FILE";
         my $tabular_result_dir="$interpro_scan_out_dir/TAB_FILE";
         my $parsed_domain_dir="$interpro_scan_out_dir/DOMAIN";
         my $parsed_interpro_dir="$interpro_scan_out_dir/INTERPRO";
         my $parsed_go_dir="$interpro_scan_out_dir/GENE_ONTOLOGY";
         my $parsed_pathway_dir="$interpro_scan_out_dir/PATHWAY"; 
         my $parsed_tmhmm_dir="$interpro_scan_out_dir/TMHMM";
         my $parsed_signalp_dir="$interpro_scan_out_dir/SIGNALP";
         my $parsed_phobius_dir="$interpro_scan_out_dir/PHOBIUS";

         mkdir("$output_dir/$project_name");
         mkdir($interpro_scan_out_dir);
         mkdir($xml_result_dir);
         mkdir($tabular_result_dir);
         mkdir($parsed_domain_dir);
         mkdir($parsed_interpro_dir);
         mkdir($parsed_go_dir);
         mkdir($parsed_pathway_dir);
         mkdir($parsed_tmhmm_dir);
         mkdir($parsed_signalp_dir);
         mkdir($parsed_phobius_dir);
         
         ##### PARAMETERS #####
         
         my $cpu_core=$config_param{PARAMETER}->{PARALLEL_CPU_CORE};
         my $interpro_path=$config_param{INTERPROSCAN_PARAMETERS}->{INTERPRO_SCAN_PATH};
         my $interpro_opts=$config_param{INTERPROSCAN_PARAMETERS}->{INTERPRO_SCAN_OPTS};
         
         my $get_genome_in_db="SELECT DISTINCT(abbreviation),genome_type from OrganismInfo";
         my ($genome_record)=SQLiteDB::get_record($db_dir,$db_name,$get_genome_in_db);
     
         my @genome_info=();

         if(scalar($genome_record)>1){ 
          
            foreach my $row(@{$genome_record}){  
         
                my $abbrv=shift(@{$row});
                my $genome_type=shift(@{$row});                        
                
                unless(-s "$xml_result_dir/$abbrv.xml"){
                  push(@genome_info,$abbrv);
                }
            }
         }
         
         my $fork=Parallel::ForkManager->new($cpu_core); 
         
         foreach my $genome_name(@genome_info){
         
             $fork->start and next;
             
             print "$genome_name\n";
           
             system("perl RunInterproscan.pl $genome_name $db_dir $db_name $xml_result_dir $tabular_result_dir $parsed_domain_dir $parsed_interpro_dir $parsed_go_dir $parsed_pathway_dir $interpro_path $interpro_opts $cpu_core $tmp_dir");
         
             system("perl ParseInterproscan.pl $genome_name $xml_result_dir/$genome_name.xml $tabular_result_dir/$genome_name.tab $parsed_domain_dir $parsed_interpro_dir $parsed_go_dir $parsed_pathway_dir $parsed_signalp_dir $parsed_tmhmm_dir $parsed_phobius_dir");
   
             $fork->finish;
         }
         $fork->wait_all_children; 

         ######### LOAD INTERPRO ANNOTATION IN THE DATABASE ######

         print "Loading Domain Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_domain_dir/DOMAIN_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO DomainAnnotation(protein_id,genome_name,seq_len,domain_id,domain_name,domain_start,domain_end,significance_value,description) VALUES
                                    (?,?,?,?,?,?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }

         print "Loading InterPro Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_interpro_dir/INTERPRO_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO InterProAnnotation(protein_id, genome_name, interpro_id, interpro_name) VALUES
                                    (?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }

         print "Loading GeneOntology Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_go_dir/GO_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO GOAnnotation(protein_id, genome_name, go_id, go_category, go_description) VALUES
                                    (?,?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }

         print "Loading Pathway Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_pathway_dir/PATHWAY_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO PathwayAnnotation(protein_id, genome_name, pathway_id, pathway_name) VALUES
                                    (?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }

         print "Loading Phobius Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_phobius_dir/PHOBIUS_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO PhobiusAnnotation(protein_id, genome_name, domain_name, domain_description, domain_start, domain_end) VALUES
                                    (?,?,?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }
         
         print "Loading SignalP Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_signalp_dir/SIGNALP_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO SignalPAnnotation(protein_id, genome_name, domain_start, domain_end) VALUES
                                    (?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);

         }

         print "Loading TMHMM Annotation\n";

         foreach my $genome_name(@genome_info){

            my $domain_file="$parsed_tmhmm_dir/TMHMM_$genome_name.txt";

            my $domain_insert_stmt="INSERT INTO TMHMMAnnotation(protein_id, genome_name, domain_start, domain_end) VALUES
                                    (?,?,?,?)";

            SQLiteDB::load_data($db_dir,$db_name,$domain_file,$domain_insert_stmt);
         }    
 
  exit;     
} ##### END OF ANNOTATION BLOCK


##### SHOW HELP #####
if(defined($help)){
   showHelp();
   exit;
}

sub showHelp {

   print STDOUT "Program: Microbial Genome analysis Toolkit\n",
                "Author:  Shalabh Thakur\n",
                "Version: 1.0\n\n";

   print STDOUT "[Mandatory Parameters]\n",
                "-config: This parameters take MGAT configuration file name as input.",
                "By Default configuration file would be searched in config folder of Home directory. ",
                "If configuration file is locted in another directory, than change the url of config dir",
                "in the script.\n\n";

   print STDOUT "[Optional Parameters]\n",
                "-project: Name of existing project directory\n",
                "-help: Give information on available options\n";  
}

##### Subroutin To Build Output Directory Structure #####

sub BuildOutputDirectory {

    my($project_name)=(shift);
    my($analysis)=(shift);
    my($output)=(shift);
    my($tmp_dir)=(shift);

     my %out_dir=();

       $out_dir{project_name}="$output/$project_name";
       $out_dir{tmp_log}=$tmp_dir;
       

     if($analysis eq "compare_reference" or $analysis eq "predict_hmm" or $analysis eq "predict_ortholog"){

       my $homolog_dir="$output/$project_name/HOMOLOG_SCAN";            $out_dir{homolog_dir}=$homolog_dir;
       my $result_dir="$output/$project_name/RESULT";                   $out_dir{result_dir}=$result_dir;
       my $ortholog_dir="$homolog_dir/ORTHOLOG";                        $out_dir{ortholog_dir}=$ortholog_dir;
       my $pair_distance="$ortholog_dir/PAIR_DISTANCE";                 $out_dir{pair_distance}=$pair_distance;
       my $ortholog_cluster="$ortholog_dir/ORTHO_CLUSTER";              $out_dir{ortholog_cluster}=$ortholog_cluster;
       my $ortholog_pair="$ortholog_dir/PAIR_ORTHOLOG";                 $out_dir{pair_ortholog}=$ortholog_pair;
       my $inparalog_pair="$ortholog_dir/PAIR_INPARALOG";               $out_dir{pair_inparalog}=$inparalog_pair;
       #my $incongruent_pair="$ortholog_dir/PAIR_INCONGRUENT";           $out_dir{pair_incongruent}=$incongruent_pair;
       my $homolog_alignment="$ortholog_dir/ALIGNMENT";                 $out_dir{homolog_alignment}=$homolog_alignment;
       my $similarity_dir="$homolog_dir/BEST_PAIR";                     $out_dir{similarity_dir}=$similarity_dir;
       my $all_similarity_dir="$homolog_dir/ALL_PAIR";                  $out_dir{all_similarity_dir}=$all_similarity_dir;
       my $chimeric_similarity_dir="$homolog_dir/CHIMERA";              $out_dir{chimeric_similarity_dir}=$chimeric_similarity_dir;
       my $mcl_dir="$homolog_dir/MCL";                                  $out_dir{mcl_dir}=$mcl_dir;   
       my $hmm_dir="$homolog_dir/HMM";                                  $out_dir{hmm_dir}=$hmm_dir;
       my $cluster_dir="$homolog_dir/CLUSTER";                          $out_dir{cluster_dir}=$cluster_dir;
       my $hmm_file_dir="$hmm_dir/MODEL";                               $out_dir{hmm_file_dir}=$hmm_file_dir;
       my $singleton_group="$hmm_dir/SINGLETON";                        $out_dir{singleton_group}=$singleton_group;
       my $hmm_db_dir="$homolog_dir/HMM_DB";                            $out_dir{hmm_db_dir}=$hmm_db_dir;
       my $hmm_out_dir="$homolog_dir/HMMER_OUT";                        $out_dir{hmm_out_dir}=$hmm_out_dir;
       my $hmm_fullout_dir="$hmm_out_dir/HMM_FULL";                     $out_dir{hmm_fullout_dir}=$hmm_fullout_dir;
       my $hmm_domout_dir="$hmm_out_dir/HMM_DOM";                       $out_dir{hmm_domout_dir}=$hmm_domout_dir; 

       mkdir($homolog_dir);
       mkdir($result_dir);
       mkdir($ortholog_dir);
       mkdir($pair_distance);
       mkdir($ortholog_cluster);
       mkdir($ortholog_pair);
       mkdir($inparalog_pair);
       #mkdir($incongruent_pair);
       mkdir($homolog_alignment);
       mkdir($similarity_dir);
       mkdir($all_similarity_dir);
       mkdir($chimeric_similarity_dir);
       mkdir($mcl_dir);
       mkdir($hmm_dir);
       mkdir($cluster_dir);
       mkdir($hmm_file_dir);
       mkdir($singleton_group);
       mkdir($hmm_db_dir);  
       mkdir($hmm_out_dir);
       mkdir($hmm_fullout_dir);
       mkdir($hmm_domout_dir);   
    }

    return(\%out_dir);
}
   


