package SmithWaterman;

# perl/SmithWaterman.pm
# url: https://github.com/DanielTillett/nucleics-align
# original author: Isaac Turner <turner.isaac@gmail.com>
# license: BSD 3-Clause, see LICENSE
# date: Oct 2012
#
# Runs bin/smith_waterman and reads its local alignments back, best first:
#
#   my $sw = SmithWaterman->new(minscore => 5, maxhits => 3);
#   $sw->do_alignment('ttacgtacgtaa', 'acgtacgt');
#   while(defined(my $hit = $sw->get_next_hit())) { $sw->print_hit($hit); }
#   my @hits = $sw->get_all_hits(); # or all at once, after do_alignment()
#
# Options (names ignore case):
#   match, mismatch, gapopen, gapextend      scores (a gap of n letters
#                                            scores gapopen + n*gapextend)
#   scoring                                  PAM30, PAM70, BLOSUM80, BLOSUM62
#                                            or DNA_HYBRIDIZATION
#   substitution_matrix, substitution_pairs  a scoring file (see scoring/)
#   wildcard                                 'N 0', or { N => 0, X => -1 }
#   minscore                                 hits scoring at least this
#                                            (default 0 here: every hit)
#   maxhits                                  at most this many hits per pair
#   fastlocal                                with maxhits 1: the best hit of
#                                            a long pair in little memory
#                                            (same score; may pick another
#                                            equal-scoring hit)
#   case_sensitive, nogaps, nogapsin1, nogapsin2, nogapsinlongest,
#   nomismatches                             flags (true to set)
#   threads                                  threads within an alignment
#   timeout                                  seconds to wait for each line of
#                                            output (default 0: no limit)
#   cmd                                      the smith_waterman to run

use strict;
use warnings;

use Carp;
use File::Basename ();
use lib File::Basename::dirname(__FILE__); # NucleicsAlign/Tool.pm beside it
use parent 'NucleicsAlign::Tool';

# Printed before each hit, or the '==' after the last, when asked for one
use constant {PROMPT_LINE => 'next [h]it or [a]lignment: '};

my %SPEC = (
  name => 'smith_waterman',
  cmd => File::Basename::dirname(__FILE__)."/../bin/smith_waterman",
  flags => [qw(case_sensitive nogaps nogapsin1 nogapsin2 nogapsinlongest
               nomismatches fastlocal)],
  values => [qw(match mismatch gapopen gapextend substitution_matrix
                substitution_pairs minscore maxhits threads)],
);

sub new
{
  my ($class, %given) = @_;
  my $self = $class->_new(\%SPEC, %given);
  $self->{_align_number} = -1; # of the alignment being read
  $self->{_seq1} = $self->{_seq2} = undef;
  $self->{_waiting} = 1; # the tool is waiting for a pair, not for h or a
  return $self;
}

# The tool, started the first time it is needed and again after an error
# stopped it (its stdin is then out of step with its prompts)
sub _proc
{
  my ($self) = @_;
  if(!$self->{_procs}{tool}) {
    $self->_start('tool', '--stdin', '--pretty');
    $self->{_waiting} = 1;
    $self->{_align_number} = -1; # a new tool numbers its alignments from 0
  }
  return 'tool';
}

# A stopped tool waits for nothing
sub _stop
{
  my ($self, @rest) = @_;
  $self->SUPER::_stop(@rest);
  $self->{_waiting} = 1;
}

# The pair the tool failed on, with its error
sub _tool_error
{
  my ($self) = @_;
  print STDERR "ErrSeq1: '".($self->{_seq1} // '')."'\n";
  print STDERR "ErrSeq2: '".($self->{_seq2} // '')."'\n";
}

sub _expect
{
  my ($self, $pattern) = @_;
  my $line = $self->_read_line('tool');
  croak("smith_waterman said '$line' (expected $pattern)")
    if $line !~ /$pattern/;
  return $line;
}

# Align seq1 against seq2; read the hits with get_next_hit(). Returns false
# (with a warning) if either sequence is empty: there are no hits then
sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;
  my $tool = $self->_proc();

  if(!$self->{_waiting})
  {
    # Skip the hits left from the previous alignment
    $self->_send($tool, "a\n");
    $self->_expect('^'.quotemeta(PROMPT_LINE).'==$');
    $self->{_waiting} = 1;
  }

  $self->_check_seqs($seq1, $seq2);
  if($seq1 eq '' || $seq2 eq '') {
    carp("Cannot align lengths of zero");
    return 0;
  }

  $self->{_align_number}++;
  $self->{_seq1} = $seq1;
  $self->{_seq2} = $seq2;
  my $expected = $self->{_align_number};

  $self->_send($tool, "$seq1\n$seq2\n");
  $self->{_waiting} = 0;

  $self->_expect("^== Alignment $expected ");
  $self->_expect('^$');
  return 1;
}

# The next hit of the current alignment, best first, or undef after the
# last: a hash of seq1, seq2, hit (0, 1, ...), score, align1, pos1, len1,
# sep (| match, * mismatch, space gap), align2, pos2, len2
sub get_next_hit
{
  my ($self) = @_;
  return undef if $self->{_waiting};

  $self->_send($self->_proc(), "h\n");
  my $line = $self->_read_line('tool');
  if(substr($line, 0, length(PROMPT_LINE)) ne PROMPT_LINE) {
    croak("smith_waterman said '$line' (expected its prompt)");
  }
  $line = substr($line, length(PROMPT_LINE));

  if($line eq '==') {
    $self->{_waiting} = 1;
    return undef;
  }
  if($line !~ /^hit \d+\.(\d+) score: (-?\d+)$/) {
    croak("smith_waterman said '$line' (expected a hit)");
  }

  my %result = ('seq1' => $self->{_seq1}, 'seq2' => $self->{_seq2},
                'hit' => $1, 'score' => $2);
  my $part = qr/^  (.*)  \[pos: (\d+); len: (\d+)\]$/;

  @result{'align1', 'pos1', 'len1'} = $self->_expect($part) =~ $part;
  $result{'sep'} = substr($self->_read_line('tool'), 2);
  @result{'align2', 'pos2', 'len2'} = $self->_expect($part) =~ $part;
  $self->_expect('^$');

  return \%result;
}

# Every remaining hit of the current alignment
sub get_all_hits
{
  my ($self) = @_;
  my @hits;
  while(defined(my $hit = $self->get_next_hit())) { push(@hits, $hit); }
  return @hits;
}

sub print_hit
{
  my ($self, $hit, $out) = @_;
  $out //= \*STDOUT;
  print $out "hit $self->{_align_number}.$hit->{'hit'} ".
             "score: $hit->{'score'}\n";
  print $out "  $hit->{'align1'}  [pos: $hit->{'pos1'}; len: $hit->{'len1'}]\n";
  print $out "  $hit->{'sep'}\n";
  print $out "  $hit->{'align2'}  [pos: $hit->{'pos2'}; len: $hit->{'len2'}]\n";
}

1;
