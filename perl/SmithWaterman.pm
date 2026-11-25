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

use File::Basename; # dirname
use Carp;
use IPC::Open3;
use Symbol qw(gensym);

# Printed before each hit, or the '==' after the last, when asked for one
use constant {PROMPT_LINE => 'next [h]it or [a]lignment: '};

my %FLAGS = map { $_ => 1 } qw(case_sensitive nogaps nogapsin1 nogapsin2
                               nogapsinlongest nomismatches fastlocal);
my %VALUES = map { $_ => 1 } qw(match mismatch gapopen gapextend
                                substitution_matrix substitution_pairs
                                minscore maxhits threads);
my @SCORING = qw(PAM30 PAM70 BLOSUM80 BLOSUM62 DNA_HYBRIDIZATION);

sub new
{
  my ($class, %given) = @_;
  my %options = map { lc($_) => $given{$_} } keys %given;

  my @args;
  my $cmd = dirname(__FILE__)."/../bin/smith_waterman";
  my $timeout = 0;

  if(!$options{'nomismatches'} &&
     defined($options{'match'}) != defined($options{'mismatch'})) {
    carp("Setting only one of match/mismatch");
  }

  for my $key (sort keys %options)
  {
    my $value = $options{$key};
    if($FLAGS{$key}) {
      push(@args, "--$key") if $value;
    }
    elsif($VALUES{$key}) {
      push(@args, "--$key", $value);
    }
    elsif($key eq 'scoring') {
      my ($name) = grep { uc($value) eq $_ } @SCORING;
      croak("Invalid scoring '$value' (one of @SCORING)") if !defined($name);
      push(@args, '--scoring', $name);
    }
    elsif($key eq 'wildcard') {
      my %wild = ref($value) eq 'HASH' ? %$value : split(' ', $value);
      push(@args, '--wildcard', $_, $wild{$_}) for sort keys %wild;
    }
    elsif($key eq 'timeout') { $timeout = $value; }
    elsif($key eq 'cmd') { $cmd = $value; }
    else { croak("Unknown option '$key'"); }
  }

  my $self = {
    _cmd => $cmd,
    _args => \@args,
    _proc => undef, # the running tool: { pid, to, from }, started on demand
    _timeout => $timeout,
    _align_number => -1, # of the alignment being read
    _seq1 => undef,
    _seq2 => undef,
    _waiting => 1, # the tool is waiting for a pair, not for h or a
  };

  return bless($self, $class);
}

# The tool, started the first time it is needed and again after an error
# stopped it (its stdin is then out of step with its prompts)
sub _proc
{
  my ($self) = @_;
  if(!$self->{_proc})
  {
    my ($to, $from) = (undef, gensym());
    # Standard error comes back with the output (its errors end a pair)
    my $pid = eval { open3($to, $from, undef, $self->{_cmd},
                           @{$self->{_args}}, '--stdin', '--pretty') };
    croak("Cannot run $self->{_cmd}: $@") if !$pid;
    $to->autoflush(1);
    $self->{_proc} = { pid => $pid, to => $to, from => $from };
    $self->{_waiting} = 1;
    $self->{_align_number} = -1; # a new tool numbers its alignments from 0
  }
  return $self->{_proc};
}

# A line to the tool. If it has gone, the write fails quietly (the next
# read then reports it) instead of SIGPIPE ending this process
sub _send
{
  my ($to, $text) = @_;
  local $SIG{PIPE} = 'IGNORE';
  print {$to} $text;
}

# Stops the tool. kill: it may be mid-alignment (a timeout), so it would
# not notice its input closing; terminate it before waiting
sub _stop
{
  my ($self, $kill) = @_;
  my $p = delete($self->{_proc}) or return;
  kill('TERM', $p->{pid}) if $kill;
  # The tool may already have gone: writing to it must not kill this process
  local $SIG{PIPE} = 'IGNORE';
  close($p->{to});
  close($p->{from});
  waitpid($p->{pid}, 0);
  $self->{_waiting} = 1;
}

# Lines of output with no line breaks; croaks on the tool's errors, which
# stop it (the next alignment starts it again)
sub _read_line
{
  my ($self) = @_;
  my $from = $self->_proc()->{from};
  my $line;
  eval {
    local $SIG{ALRM} = sub { die("timeout\n") };
    alarm($self->{_timeout});
    $line = <$from>;
    alarm(0);
  };
  alarm(0);
  if($@) {
    $self->_stop(1);
    croak($@ eq "timeout\n" ? "Timeout reading smith_waterman output"
                            : "Cannot read smith_waterman output: $@");
  }
  if(!defined($line)) {
    $self->_stop();
    croak("smith_waterman stopped (have you compiled?)");
  }
  $line =~ s/\r?\n$//;
  if($line =~ /^(Alignment )?(Error|Warning):/) {
    print STDERR "ErrSeq1: '".($self->{_seq1} // '')."'\n";
    print STDERR "ErrSeq2: '".($self->{_seq2} // '')."'\n";
    $self->_stop();
    croak($line);
  }
  return $line;
}

sub _expect
{
  my ($self, $pattern) = @_;
  my $line = $self->_read_line();
  croak("smith_waterman said '$line' (expected $pattern)")
    if $line !~ /$pattern/;
  return $line;
}

# Align seq1 against seq2; read the hits with get_next_hit(). Returns false
# (with a warning) if either sequence is empty: there are no hits then
sub do_alignment
{
  my ($self, $seq1, $seq2) = @_;
  my $to = $self->_proc()->{to};

  if(!$self->{_waiting})
  {
    # Skip the hits left from the previous alignment
    _send($to, "a\n");
    $self->_expect('^'.quotemeta(PROMPT_LINE).'==$');
    $self->{_waiting} = 1;
  }

  for my $seq ($seq1, $seq2) {
    croak("Sequences must be defined") if !defined($seq);
    croak("Whitespace or control characters not allowed in sequences")
      if $seq =~ /[\x00-\x20\x7f]/;
    # The tool would read these as a FASTA or FASTQ header and wait for
    # the rest of the record
    croak("Sequences cannot start with '>' or '@'") if $seq =~ /^[>@]/;
  }
  if($seq1 eq '' || $seq2 eq '') {
    carp("Cannot align lengths of zero");
    return 0;
  }

  $self->{_align_number}++;
  $self->{_seq1} = $seq1;
  $self->{_seq2} = $seq2;
  my $expected = $self->{_align_number};

  _send($to, "$seq1\n$seq2\n");
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

  _send($self->_proc()->{to}, "h\n");
  my $line = $self->_read_line();
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
  $result{'sep'} = substr($self->_read_line(), 2);
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

# Stop the tool (also done when the object goes)
sub destructor
{
  my ($self) = @_;
  $self->_stop();
}

sub DESTROY
{
  my ($self) = @_;
  local ($?, $@);
  $self->destructor();
}

1;
