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
 ./vep -i variants.vcf --plugin BFWalk,file=scores.tsv,map=uniprot_parsed.tsv
 ./vep -i variants.vcf --plugin BFWalk,file=scores.tsv,map=uniprot_parsed.tsv,label=MMAF

=head1 DESCRIPTION
 
 
 An Ensembl VEP plugin that adds BFWalk scores to transcript-level output.
 BFWalk (https://github.com/jedrzejkubica/BFWalk) is a network-propagation
 algorithm based on non-backtracking walks. It scores every protein in an
 interactome by its proximity to a set of seeds (e.g. known phenotype-associated
 genes/proteins). Higher scores mean more likely to be associated with the phenotype.

 Scores depend on the seed set. Re-run BFWalk for each phenotype, and use
 label= to run the plugin several times in one VEP command.

 BFWalk scores UniProt proteins when run on the interactome.
 The plugin maps UniProt ACs to genes with the
 uniprot_parsed.tsv file produced by BFWalk's Interactome/uniprot_parser.py
 (columns PrimaryAC, SecondaryACs, TaxID, GeneName, Synonyms). Use the same
 file that was used to build the interactome.

 Parameters (key=value, comma-separated):
   file=FILE      BFWalk scores with header NODE\tSCORE (required)
   map=FILE       uniprot_parsed.tsv (required)
   label=NAME     column prefix becomes BFWalk_NAME_ (optional)
   format=FMT     sprintf format for scores, default %.3g

 Output columns:
   BFWalk_score    BFWalk score for the gene


=cut

package BFWalk;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $SCORE_DESC = 'BFWalk score (interactome-based, higher means closer more likely to be associated with the phenotype)';

sub new {
  my $class = shift;
  my $self  = $class->SUPER::new(@_);

  my $p = $self->params_to_hash();
  my $file    = $p->{file} or die "ERROR: BFWalk requires file=<BFWalk scores.tsv>\n";
  my $map     = $p->{map}  or die "ERROR: BFWalk requires map=<uniprot_parsed.tsv>\n";
  my $format  = $p->{format} // '%.3g';

  for my $f (grep { defined } $file, $map) {
    die "ERROR: BFWalk file '$f' not found\n" unless -e $f;
  }

  $self->{column} = (defined $p->{label} ? "BFWalk_$p->{label}" : 'BFWalk') . '_score';


  # BFWalk scores: parse UniProt AC to score
  my %prot_score;
  open my $fh, '<', $file or die "ERROR: cannot open $file: $!\n";

  my $hdr = <$fh>;
  die "ERROR: BFWalk: $file is empty\n" unless defined $hdr;
  chomp $hdr;
  my %scol; my $i = 0;
  $scol{uc $_} = $i++ for split /\t/, $hdr;
  for (qw(NODE SCORE)) {
    die "ERROR: BFWalk: column '$_' not found in header of $file\n" unless exists $scol{$_};
  }

  my $n_bad = 0;
  while (my $line = <$fh>) {
    chomp $line;
    next unless length $line;
    my @f = split /\t/, $line;
    my ($ac, $score) = @f[$scol{NODE}, $scol{SCORE}];
    unless (defined $ac && defined $score && $score =~ /^-?[\d.]+(?:e[-+]?\d+)?$/i) {
      $n_bad++;
      next;
    }
    $ac =~ s/-\d+$//;
    $prot_score{uc $ac} = $score;
  }
  close $fh;

  warn "WARNING: BFWalk: skipped $n_bad malformed lines in $file\n" if $n_bad;
  die "ERROR: no scores read from $file\n" unless %prot_score;

  # uniprot_parsed.tsv: parse the UniProt AC to gene name mapping
  my %ac2gene;
  open my $mh, '<', $map or die "ERROR: cannot open $map: $!\n";
  my $mhdr = <$mh>;
  die "ERROR: BFWalk: $map is empty\n" unless defined $mhdr;
  chomp $mhdr;
  my %col; $i = 0;
  $col{$_} = $i++ for split /\t/, $mhdr;
  for (qw(PrimaryAC TaxID GeneName Synonyms)) {
    die "ERROR: BFWalk: column '$_' not found in $map header\n" unless exists $col{$_};
  }
  while (my $line = <$mh>) {
    chomp $line;
    my @f = split /\t/, $line, -1;
    next unless ($f[$col{TaxID}] // '') =~ /\b9606\b/;
    my $gene = lc($f[$col{GeneName}] // '');
    next unless length $gene;
    my $ac = uc($f[$col{PrimaryAC}] // '');
    $ac2gene{$ac} = $gene if exists $prot_score{$ac};
  }
  close $mh;

  # map UniProt ACs to genes
  my %raw;
  for my $ac (keys %prot_score) {
    my $g = $ac2gene{$ac} // next;
    $raw{$g} = $prot_score{$ac}
      if !exists $raw{$g} || $prot_score{$ac} > $raw{$g};
  }
  my $n_unmapped = grep { !exists $ac2gene{$_} } keys %prot_score;
  warn sprintf("WARNING: BFWalk: %d of %d proteins have no human gene in %s\n",
               $n_unmapped, scalar keys %prot_score, $map) if $n_unmapped;
  die "ERROR: BFWalk: no proteins could be mapped to genes using $map\n" unless %raw;


  $self->{scores} = { map { $_ => sprintf($format, $raw{$_}) } keys %raw };

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;
  return { $self->{column} => $SCORE_DESC };
}

sub run {
  my ($self, $tva) = @_;
  my $tr  = $tva->transcript or return {};
  my $key = $tr->{_gene_symbol} || $tr->{_gene_hgnc} or return {};
  $key = lc $key;
  return {} unless exists $self->{scores}{$key};
  return { $self->{column} => $self->{scores}{$key} };
}

1;
