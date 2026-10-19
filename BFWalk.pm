=head1 LICENSE

Copyright [2026] Jędrzej Kubica, Nicolas Thierry-Mieg

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  BFWalk <https://github.com/jedrzejkubica/BFWalk>

=cut

=head1 NAME

BFWalk - Add BFWalk score to the output

=head1 SYNOPSIS

 python BFWalk.py --network interactome.sif --seeds seeds.txt 1>scores.tsv
 ./vep -i variants.vcf --plugin BFWalk,file=scores.tsv
 ./vep -i variants.vcf --plugin BFWalk,file=MMAF:scores_MMAF.tsv,file=NOA:scores_NOA.tsv

=head1 DESCRIPTION


 An Ensembl VEP plugin that adds BFWalk scores to VEP output.
 BFWalk (https://github.com/jedrzejkubica/BFWalk) is a network-propagation
 algorithm based on non-backtracking walks. It scores every protein in an
 interactome by its proximity to a set of seeds (e.g. known phenotype-associated
 genes/proteins). Higher scores mean more likely to be associated with the phenotype.

 Scores depend on the seed set. Run BFWalk once per phenotype, then give
 the plugin one file= per phenotype. Each phenotype gets its own column.

 BFWalk scores UniProt proteins when run on the interactome. The plugin
 matches them to transcripts using the UniProt/Swiss-Prot accessions that
 Ensembl stores for each transcript's protein product in the VEP cache.
 Accession versions (P04637.293) and isoform suffixes
 (P04637-2) are ignored when matching.

 Scores are reported per transcript: only transcripts whose
 protein product has a Swiss-Prot accession can receive a score.
 Non-coding transcripts, and coding transcripts that Ensembl maps only to
 TrEMBL, get no value. If a transcript has several Swiss-Prot accessions,
 the highest score is reported.

 Parameters:
   file=PHENOTYPE:FILE  BFWalk scores with header NODE\tSCORE for phenotype (PHENOTYPE).
   format=FMT      sprintf format for scores, default %.3g

 Output columns:
   BFWalk_NAME_score   BFWalk score for the transcript's protein, one column per phenotype
   BFWalk_score        the same, when a single file is given without a name

 Examples:
   --plugin BFWalk,file=scores.tsv
   --plugin BFWalk,file=MMAF:scores_MMAF.tsv,file=NOA:scores_NOA.tsv

 Requires a VEP cache (or database) that includes UniProt annotation,
 which the standard Ensembl human caches do.


=cut
 
package BFWalk;
 
use strict;
use warnings;
 
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
 
my $COLUMN_DESCRIPTION = "BFWalk score of the transcript's Swiss-Prot protein (higher means closer to the seed proteins)";
my $DEFAULT_SCORE_FORMAT = '%.3g';
 
 
# new() is called once, when VEP starts.
# It reads all score files into memory, so that run() only has to look up values.
sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
 
  # VEP's params_to_hash() keeps only the last value of a repeated key,
  # but file= is given once per phenotype, so we read the parameters ourselves.
  my @file_parameters = ();
  my $score_format = $DEFAULT_SCORE_FORMAT;
 
  foreach my $parameter (@{ $self->params() }) {
    my ($key, $value) = split(/=/, $parameter, 2);
 
    if (!defined($value) || $value eq '') {
      die("ERROR: BFWalk: parameters must be key=value, got '$parameter'\n");
    }
 
    if ($key eq 'file') {
      push(@file_parameters, $value);
    }
    elsif ($key eq 'format') {
      $score_format = $value;
    }
    else {
      die("ERROR: BFWalk: unknown parameter '$key' (expected file, format)\n");
    }
  }
 
  if (scalar(@file_parameters) == 0) {
    die("ERROR: BFWalk requires at least one file=PHENOTYPE:FILE\n");
  }
 
  $self->{score_format} = $score_format;
 
  # Scores from all phenotypes are stored in a single table, keyed by protein:
  #   $self->{scores_by_protein}{UniProt accession}{output column name} = score
  # so that run() can find every phenotype's score for a protein with one lookup.
  $self->{scores_by_protein} = {};
  $self->{column_names} = [];
 
  foreach my $file_parameter (@file_parameters) {
    my $phenotype;
    my $scores_file;
 
    if ($file_parameter =~ /^(\w+):(.+)$/) {
      $phenotype = $1;
      $scores_file = $2;
    }
    else {
      $scores_file = $file_parameter;
    }
 
    # Without a phenotype the column is just "BFWalk_score".
    # With several files, every file needs its own phenotype,
    # otherwise they would all write into the same column.
    my $column_name;
    if (defined($phenotype)) {
      $column_name = "BFWalk_" . $phenotype . "_score";
    }
    elsif (scalar(@file_parameters) == 1) {
      $column_name = "BFWalk_score";
    }
    else {
      die("ERROR: BFWalk: with several files, each needs a phenotype: file=PHENOTYPE:FILE ('$file_parameter')\n");
    }
 
    foreach my $existing_column_name (@{ $self->{column_names} }) {
      if ($existing_column_name eq $column_name) {
        die("ERROR: BFWalk: phenotype '$phenotype' used twice\n");
      }
    }
 
    if (!-e $scores_file) {
      die("ERROR: BFWalk file '$scores_file' not found\n");
    }
 
    my $score_by_protein = read_scores_file($scores_file);
 
    foreach my $protein (keys(%{$score_by_protein})) {
      $self->{scores_by_protein}{$protein}{$column_name} = $score_by_protein->{$protein};
    }
 
    push(@{ $self->{column_names} }, $column_name);
  }
 
  return $self;
}
 
 
# Reads one BFWalk scores file and returns a reference to a hash:
#   { UniProt accession => score }
# The file must have a header line with the columns NODE and SCORE.
sub read_scores_file {
  my ($scores_file) = @_;
 
  open(my $file_handle, '<', $scores_file) or die("ERROR: cannot open $scores_file: $!\n");
 
  my $header_line = <$file_handle>;
  if (!defined($header_line)) {
    die("ERROR: BFWalk: $scores_file is empty\n");
  }
  chomp($header_line);
 
  # Find the columns by name rather than by position,
  # so the file still works if the columns are in another order.
  my @header_columns = split(/\t/, $header_line);
  my $node_column_index;
  my $score_column_index;
 
  for (my $index = 0; $index < scalar(@header_columns); $index++) {
    my $column_title = uc($header_columns[$index]);
    if ($column_title eq 'NODE') {
      $node_column_index = $index;
    }
    elsif ($column_title eq 'SCORE') {
      $score_column_index = $index;
    }
  }
 
  if (!defined($node_column_index)) {
    die("ERROR: BFWalk: column 'NODE' not found in header of $scores_file\n");
  }
  if (!defined($score_column_index)) {
    die("ERROR: BFWalk: column 'SCORE' not found in header of $scores_file\n");
  }
 
  my %score_by_protein = ();
  my $number_of_malformed_lines = 0;
 
  while (my $line = <$file_handle>) {
    chomp($line);
 
    if ($line eq '') {
      next;
    }
 
    my @fields = split(/\t/, $line);
    my $protein = $fields[$node_column_index];
    my $score = $fields[$score_column_index];
 
    # Accept plain numbers and scientific notation (e.g. 1.2e-05).
    # Anything else is counted and reported once, instead of being silently ignored.
    if (!defined($protein) || !defined($score) || $score !~ /^-?[\d.]+(e[-+]?\d+)?$/i) {
      $number_of_malformed_lines++;
      next;
    }
 
    $score_by_protein{ normalise_accession($protein) } = $score;
  }
 
  close($file_handle);
 
  if ($number_of_malformed_lines > 0) {
    warn("WARNING: BFWalk: skipped $number_of_malformed_lines malformed lines in $scores_file\n");
  }
  if (scalar(keys(%score_by_protein)) == 0) {
    die("ERROR: no scores read from $scores_file\n");
  }
 
  return \%score_by_protein;
}
 
 
# The VEP cache stores versioned accessions (P04637.293) and BFWalk nodes may
# be isoforms (P04637-2), while the interactome is built on plain accessions.
# Removing both suffixes lets the two sources match.
sub normalise_accession {
  my ($accession) = @_;
 
  $accession = uc($accession);
  $accession =~ s/\.\d+$//;
  $accession =~ s/-\d+$//;
 
  return $accession;
}
 
 
# Tells VEP to call run() only for variants overlapping transcripts,
# since only transcripts have a protein product to look up.
sub feature_types {
  return ['Transcript'];
}
 
 
# Tells VEP which columns this plugin adds, so it can describe them
# in the output header (and in the CSQ field description for VCF output).
sub get_header_info {
  my $self = shift;
 
  my %column_descriptions = ();
  foreach my $column_name (@{ $self->{column_names} }) {
    $column_descriptions{$column_name} = $COLUMN_DESCRIPTION;
  }
 
  return \%column_descriptions;
}
 
 
# run() is called once per variant allele and overlapping transcript.
# It returns a hash of { column name => value }, which VEP adds to that output row.
# An empty hash means: add nothing to this row.
sub run {
  my ($self, $transcript_variation_allele) = @_;
 
  my $transcript = $transcript_variation_allele->transcript();
  if (!defined($transcript)) {
    return {};
  }
 
  # The VEP cache stores the Swiss-Prot accession(s) of the transcript's protein
  # as a comma-separated string. Transcripts without one have no value or '-'.
  my $swissprot_accessions = $transcript->{_swissprot};
  if (!defined($swissprot_accessions) || $swissprot_accessions eq '' || $swissprot_accessions eq '-') {
    return {};
  }
 
  # A transcript can have several accessions; for each phenotype keep the highest score.
  my %highest_score_by_column = ();
 
  foreach my $accession (split(/,/, $swissprot_accessions)) {
    my $scores_by_column = $self->{scores_by_protein}{ normalise_accession($accession) };
    if (!defined($scores_by_column)) {
      next;
    }
 
    foreach my $column_name (keys(%{$scores_by_column})) {
      my $score = $scores_by_column->{$column_name};
      if (!exists($highest_score_by_column{$column_name}) || $score > $highest_score_by_column{$column_name}) {
        $highest_score_by_column{$column_name} = $score;
      }
    }
  }
 
  my %output = ();
  foreach my $column_name (keys(%highest_score_by_column)) {
    $output{$column_name} = sprintf($self->{score_format}, $highest_score_by_column{$column_name});
  }
 
  return \%output;
}
 
1;