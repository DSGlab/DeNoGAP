##### Module to Predict Ortholog Genes #######
##### Author: Shalabh Thakur #################
##### Date: 28-SEP-2013 ######################

#!/usr/bin/perl
package SQLiteDB;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use DBI;


use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(create_db create_table create_organism_table create_profile_table load_data load_from_array fetch_similarity get_record execute_sql);


sub create_db{

    my $db_dir=(shift);
    my $db_name=(shift);

    #### CREATE NEW SQLite DATABASE ####
    chdir($db_dir);   
    print "$db_dir\n";
    my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $DBI::errstr;

    my $stmt=$dbh->prepare("SELECT SQLITE_VERSION()");
    $stmt->execute();

    my $ver=$stmt->fetch();
    print "$db_name Datbase Created in SQLite with version ",@$ver,"\n";    
    $stmt->finish();
    $dbh->disconnect();    
    chdir($Bin);
}

###### Create Tables in Database ####

sub create_table{

    my $db_dir=(shift);
    my $db_name=(shift);

    chdir($db_dir);
    my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $DBI::errstr;

    #### Create Taxonomy Table ####    
    $dbh->do("CREATE TABLE IF NOT EXISTS Taxonomy(taxa_index_id INTEGER PRIMARY KEY AUTOINCREMENT, taxonomy_id INTEGER UNIQUE, kingdom TEXT, phylum TEXT, class TEXT, order_name TEXT, family TEXT, genus TEXT, species TEXT UNIQUE)");

    #### Create OrganismInfo Table ####
    #$dbh->do("CREATE TABLE OrganismInfo(genome_index_id INTEGER PRIMARY KEY AUTOINCREMENT, genome_name TEXT UNIQUE, species TEXT, abbreviation TEXT, use_outgroup TEXT, additional_info TEXT)");

    #### Create GeneFeature Table ####
    $dbh->do("CREATE TABLE IF NOT EXISTS GeneFeature(feature_index_id INTEGER PRIMARY KEY AUTOINCREMENT, feature_id TEXT UNIQUE, feature_type TEXT, genome_id TEXT, genome_type TEXT, genome_name TEXT, genome_length INT, feature_start INT, feature_end INT, nuc_length INT, aa_length INT, strand TEXT, frame INT, index_in_genome INT, description TEXT)");
  
    #### Create Protein Sequence Table #####
    $dbh->do("CREATE TABLE IF NOT EXISTS ProteinSequence(prot_index_id INTEGER PRIMARY KEY AUTOINCREMENT, protein_id TEXT, genome_abbreviation TEXT, seq_type TEXT, seq_length INT, aminoacid_sequence TEXT, UNIQUE(protein_id,genome_abbreviation,seq_type))");
    
    #### Create Nucleotide Sequence Table ####
    #$dbh->do("DROP TABLE IF EXISTS NucleotideSequence");
    $dbh->do("CREATE TABLE IF NOT EXISTS NucleotideSequence(nuc_index_id INTEGER PRIMARY KEY AUTOINCREMENT, nucleotide_id TEXT, genome_abbreviation TEXT, seq_type TEXT, seq_length INT, nucleotide_sequence TEXT, UNIQUE(nucleotide_id,genome_abbreviation,seq_type))");

    #### Create Homolog Group Sequence Alignment Table ####
    #$dbh->do("DROP TABLE IF EXISTS Alignment");
    $dbh->do("CREATE TABLE IF NOT EXISTS Alignment(feature_id TEXT, genome_abbreviation TEXT, seq_type TEXT, alignment_id TEXT, alignment_length INT, alignment_sequence TEXT, UNIQUE(feature_id,genome_abbreviation,alignment_id,seq_type))");

    #### Create GeneFamilyRelation Table ####
    #$dbh->do("DROP TABLE IF EXISTS Similarity");
    $dbh->do("CREATE TABLE IF NOT EXISTS Similarity(query_id TEXT, subject_id TEXT, query_length INT, subject_length INT, total_domain INT, high_scoring_domain INT, qstart INT, qend INT, sstart INT, send INT, evalue REAL, bitscore REAL, percent_identity REAL, percent_similarity REAL, query_coverage REAL, subject_coverage REAL, pair_relation TEXT, note TEXT, UNIQUE(query_id, subject_id, query_length, subject_length, total_domain, high_scoring_domain, qstart, qend, sstart, send, evalue, bitscore, percent_identity, percent_similarity, query_coverage, subject_coverage, pair_relation))");

    ####### Create Parent-Child protein family relationship Table ####
    $dbh->do("DROP TABLE IF EXISTS LinkFamily");
    $dbh->do("CREATE TABLE IF NOT EXISTS LinkFamily(family_idA TEXT, family_idB TEXT, significance REAL, UNIQUE(family_idA, family_idB))");

    ####### Create Distance Pair Table #######  
    #$dbh->do("DROP TABLE IF EXISTS DistancePair");  
    $dbh->do("CREATE TABLE IF NOT EXISTS DistancePair(taxonA TEXT, idA TEXT , taxonB TEXT, idB TEXT, divergence REAL, homolog_cluster_id TEXT , UNIQUE(taxonA,idA,taxonB,idB,homolog_cluster_id))");

    ####### Create Ortholog Pair Table #######
    $dbh->do("CREATE TABLE IF NOT EXISTS Ortholog(ortho_index_id INTEGER PRIMARY KEY AUTOINCREMENT, taxonA TEXT, idA TEXT , taxonB TEXT, idB TEXT, divergence REAL, homolog_cluster_id TEXT, UNIQUE(taxonA,idA,taxonB,idB,homolog_cluster_id))");

    ####### Create Inparalog Pair Table ######
    $dbh->do("CREATE TABLE IF NOT EXISTS Inparalog(inpar_index_id INTEGER PRIMARY KEY AUTOINCREMENT, taxonA TEXT, idA TEXT , taxonB TEXT, idB TEXT, divergence REAL, min_ortholog_divergence REAL, homolog_cluster_id TEXT, UNIQUE(taxonA,idA,taxonB,idB,homolog_cluster_id))");   

    ####### Create Inparalog Pair Table ######
    #$dbh->do("CREATE TABLE IF NOT EXISTS Incongruent(inpar_index_id INTEGER PRIMARY KEY AUTOINCREMENT, taxonA TEXT, idA TEXT , taxonB TEXT, idB TEXT, divergence REAL,min_ortholog_divergence REAL,max_ortholog_divergence REAL, homolog_cluster_id TEXT, UNIQUE(taxonA,idA,taxonB,idB,homolog_cluster_id))");   

    #### Create GenetoSuperFamily Map Table ######
    $dbh->do("CREATE TABLE IF NOT EXISTS GenetoSuperFamily(gene_superfamily_index_id INTEGER PRIMARY KEY AUTOINCREMENT, gene_id TEXT, genome_name TEXT, hmm_family_id TEXT, super_family_id TEXT, UNIQUE(gene_id,genome_name,hmm_family_id,super_family_id))");

    #### Create MapGeneToOrthologGeneFamily #####
    #$dbh->do("DROP TABLE IF EXISTS MapGeneIdtoGeneFamily");
    $dbh->do("CREATE TABLE IF NOT EXISTS MapGeneIdtoGeneFamily(familymap_index_id INTEGER PRIMARY KEY AUTOINCREMENT, genefamily_id TEXT, gene_id TEXT, species_abbreviation TEXT, UNIQUE(gene_id, genefamily_id, species_abbreviation))");
 
    ########### ANNOTATION TABLES #########

    #### Create Gene Ontology Table ####
    #$dbh->do("DROP TABLE IF EXISTS GOAnnotation");
    $dbh->do("CREATE TABLE IF NOT EXISTS GOAnnotation(protein_id TEXT, genome_name TEXT, go_id TEXT, go_category TEXT, go_description TEXT)");

    #### Create InterPro Table #########
    #$dbh->do("DROP TABLE IF EXISTS InterProAnnotation");
    $dbh->do("CREATE TABLE IF NOT EXISTS InterProAnnotation(protein_id TEXT, genome_name TEXT, interpro_id TEXT, interpro_name TEXT)");

    #### Create Pfam Table #############
    #$dbh->do("DROP TABLE IF EXISTS PfamAnnotation");
    $dbh->do("CREATE TABLE IF NOT EXISTS DomainAnnotation(protein_id TEXT, genome_name TEXT, seq_len INT, domain_id TEXT, domain_name TEXT, domain_start TEXT, domain_end TEXT, significance_value REAL, description TEXT)");
 
    $dbh->do("CREATE TABLE IF NOT EXISTS PathwayAnnotation(protein_id TEXT, genome_name TEXT, pathway_id TEXT, pathway_name TEXT)");

    $dbh->do("CREATE TABLE IF NOT EXISTS PhobiusAnnotation(protein_id TEXT, genome_name TEXT, domain_name TEXT, domain_description TEXT, domain_start INT, domain_end INT)");

    $dbh->do("CREATE TABLE IF NOT EXISTS SignalPAnnotation(protein_id TEXT, genome_name TEXT, domain_start INT, domain_end INT)");

    $dbh->do("CREATE TABLE IF NOT EXISTS TMHMMAnnotation(protein_id TEXT, genome_name TEXT, domain_start INT, domain_end INT)");

    
    $dbh->disconnect();
    chdir($Bin);
}

sub create_organism_table {
    my $db_dir=(shift);
    my $db_name=(shift);
    my $load_file=(shift);

    my $insert_column='';
    my $bind_value='';

    open(FILE,"$load_file");
    my @input_file=<FILE>;
    close FILE;

    chdir($db_dir);
    my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $DBI::errstr;

    foreach my $line(@input_file){

      chomp($line);

      if($line=~/\#/){

        $line=~s/\#//g;
        $line=~s/\n//g;

        my @table_column=split("\t",$line);  
        my $table_column=join(" TEXT,",@table_column);
           $table_column=$table_column." TEXT"; 
           $table_column=~s/(genome_name TEXT)/genome_name TEXT UNIQUE/g;

        $dbh->do("DROP TABLE IF EXISTS OrganismInfo");
        $dbh->do("CREATE TABLE IF NOT EXISTS OrganismInfo(genome_index_id INTEGER PRIMARY KEY AUTOINCREMENT, $table_column)"); 

        $insert_column=join(",",@table_column);
        $bind_value=$insert_column;
        $bind_value=~s/(\w+)/\?/g;  
        last;     
      }
   }
    
   $dbh->disconnect();
   chdir($Bin); 
   
   return($insert_column,$bind_value);   
}

###### Create Table for Profile Matrix ######
sub create_profile_table {
    my $db_dir=(shift);
    my $db_name=(shift);
    my @genome=@{(shift)};
     
    my $species_column=join(" INT,",@genome);
       $species_column=$species_column." INT";
       $species_column=~s/\n//g; 

    my $create_sql_stmt="CREATE TABLE Profile(id TEXT, $species_column)";

    chdir($db_dir);
    my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","") or die $DBI::errstr;
    
    $dbh->do("DROP TABLE IF EXISTS Profile");
    $dbh->do($create_sql_stmt);

    my $insert_genome=join(",",@genome);
       $insert_genome=~s/\n//g;
    my $bind_value=$insert_genome;
       $bind_value=~s/(\w+)/\?/g;

    $dbh->disconnect();
    chdir($Bin);
  
    return($insert_genome,$bind_value);
}

###### Load data from file #####

sub load_data{

    my $db_dir=(shift);
    my $db_name=(shift);
    my $load_file=(shift);
    my $load_sql_stmt=(shift);

    open(LOAD_FILE,"$load_file");
    my @load_data=<LOAD_FILE>;
    close LOAD_FILE;

   chdir($db_dir);
   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1, AutoCommit =>0}) or die $DBI::errstr;

   my $stmt=$dbh->prepare($load_sql_stmt);

   foreach my $line(@load_data){
      chomp($line);      
      if($line=~/\#/){next;}
      $line=~s/\s+$//g;   
      print "$line\n";   
      my @load_column=split("\t",$line);

    $stmt->execute(@load_column) or die $DBI::errstr," @load_column\n";
   }

   $stmt->finish;
   $dbh->commit;
   $dbh->disconnect();
   chdir($Bin);
}
######## Load Sequence  in Database ######

sub load_from_array{

    my $db_dir=(shift);
    my $db_name=(shift);
    my @seq_array=@{(shift)};
    my $load_seq_stmt=(shift);

    chdir($db_dir);
    my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1, AutoCommit =>0}) or die $DBI::errstr;

    my $stmt=$dbh->prepare($load_seq_stmt);

   foreach my $line(@seq_array){
      chomp($line);
      if($line=~/\#/){next;}
      print "$line\n";
      my @load_column=split("\t",$line); 
      $stmt->execute(@load_column) or die $DBI::errstr," @load_column\n";
   }

   $stmt->finish;
   $dbh->commit;
   $dbh->disconnect();
   chdir($Bin);
}
##### Fetch record from table #####
sub get_record {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   my @row_data=();

   chdir($db_dir);
   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1}) or die $DBI::errstr; 

   my $stmt=$dbh->prepare("$sql_stmt");
   $stmt->execute();

   while(my @row=$stmt->fetchrow_array()){
      push(@row_data,\@row);
   } 

   $stmt->finish;
   $dbh->disconnect();
   chdir($Bin);
   return(\@row_data);   
}

#### insert record in linked family #####

sub execute_sql {

   my $db_dir=(shift);
   my $db_name=(shift);
   my $sql_stmt=(shift);

   chdir($db_dir);
   my $dbh=DBI->connect("dbi:SQLite:dbname=$db_name","","",{RaiseError =>1}) or die $DBI::errstr; 

   my $stmt=$dbh->prepare("$sql_stmt");
   $stmt->execute();

   $stmt->finish;
   $dbh->disconnect();
   chdir($Bin);
}
