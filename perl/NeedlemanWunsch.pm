package NeedlemanWunsch;

# perl/NeedlemanWunsch.pm
# url: https://github.com/DanielTillett/nucleics-align
# original author: Isaac Turner <turner.isaac@gmail.com>
# license: BSD 3-Clause, see LICENSE
# date: Oct 2012
#
# Runs bin/needleman_wunsch and reads its alignments back:
#
#   my $nw = NeedlemanWunsch->new(gapopen => -5, linear => 1);
#   my $aln = $nw->do_alignment('ACGTTGCA', 'ACGTGCA');
#   print "$aln->{align1}\n$aln->{sep}\n$aln->{align2}\n$aln->{score}\n";
#   my $score = $nw->do_score('ACGTTGCA', 'ACGTGCA'); # score only
#
# Options (names ignore case):
#   match, mismatch, gapopen, gapextend      scores (a gap of n letters
#                                            scores gapopen + n*gapextend)
#   scoring                                  PAM30, PAM70, BLOSUM80, BLOSUM62
#                                            or DNA_HYBRIDIZATION
#   substitution_matrix, substitution_pairs  a scoring file (see scoring/)
#   wildcard                                 'N 0', or { N => 0, X => -1 }
#   case_sensitive, freestartgap, freeendgap, nogaps, nogapsin1, nogapsin2,
#   nogapsinlongest, nomismatches            flags (true to set)
#   linear                                   align in linear memory (same
#                                            score; may pick another
#                                            equal-scoring alignment)
#   linearfallback                           linear memory only where too big
#                                            for the score matrices
#   linearabove                              linear memory for pairs whose
#                                            longer sequence has at least
#                                            this many letters
#   wfa, wfabudget                           the wavefront algorithm first on
#                                            similar pairs, within wfabudget
#                                            diagonal updates (same score;
#                                            may pick another equal-scoring
#                                            alignment)
#   threads                                  threads within an alignment
#   timeout                                  seconds to wait for each line of
#                                            output (default 0: no limit)
#   cmd                                      the needleman_wunsch to run

use strict;
use warnings;

use Carp;
use File::Basename ();
use lib File::Basename::dirname(__FILE__); # NucleicsAlign/Tool.pm beside it
use parent 'NucleicsAlign::Tool';

my %SPEC = (
  name => 'needleman_wunsch',
  cmd => File::Basename::dirname(__FILE__)."/../bin/needleman_wunsch",
  flags => [qw(case_sensitive freestartgap freeendgap linear linearfallback
               wfa nogaps nogapsin1 nogapsin2 nogapsinlongest nomismatches)],
  values => [qw(match mismatch gapopen gapextend substitution_matrix
                substitution_pairs linearabove wfabudget threads)],
);

sub new
{
  my ($class, %given) = @_;
  my $self = $class->_new(\%SPEC, %given);
  $self->{_align_number} = 0;
  return $self;
}

# The tool for mode ('align' or 'score'), started the first time
sub _proc
{
  my ($self, $mode) = @_;
  if(!$self->{_procs}{$mode}) {
    my @extra = $mode eq 'score' ? ('--scoreonly')
                                 : ('--pretty', '--printscores');
    $self->_start($mode, @extra, '--stdin');
  }
  return $mode;
}

sub _score_of
{
  my ($line) = @_;
  if(!defined($line) || $line !~ /^score: (-?\d+)$/) {
    croak("Cannot find the score in '".($line // '')."'");
  }
  return $1;
}

# Best global alignment of seq1 and seq2: a hash of seq1, seq2, align1, sep
# (| match, * mismatch, space gap), align2, score and number (0, 1, ...)
sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;
  $self->_check_seqs($seq1, $seq2);
  my %result = ('seq1' => $seq1, 'seq2' => $seq2,
                'number' => $self->{_align_number}++);
  my @lines;

  if($seq1 eq '' || $seq2 eq '') {
    # The pair comes after "--", so a sequence may start with '-'
    @lines = $self->_run_once('--pretty', '--printscores', '--', $seq1, $seq2);
  }
  else {
    my $mode = $self->_proc('align');
    $self->_send($mode, "$seq1\n$seq2\n");
    @lines = map { $self->_read_line($mode) } 1..5;
  }

  @result{'align1', 'sep', 'align2'} = @lines[0..2];
  $result{'score'} = _score_of($lines[3]);
  return \%result;
}

# The best global alignment's score alone (one pass, linear memory)
sub do_score
{
  my ($self, $seq1, $seq2) = @_;
  $self->_check_seqs($seq1, $seq2);
  if($seq1 eq '' || $seq2 eq '') {
    return _score_of(($self->_run_once('--scoreonly', '--', $seq1, $seq2))[0]);
  }
  my $mode = $self->_proc('score');
  $self->_send($mode, "$seq1\n$seq2\n");
  my $score = _score_of($self->_read_line($mode));
  $self->_read_line($mode);
  return $score;
}

sub print_alignment
{
  my ($self, $hit, $out) = @_;
  $out //= \*STDOUT;
  print $out $hit->{'align1'}."\n";
  print $out $hit->{'sep'}."\n";
  print $out $hit->{'align2'}."\n";
  print $out "score: $hit->{'score'}\n\n";
}

1;
