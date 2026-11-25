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
#   threads                                  threads within an alignment
#   timeout                                  seconds to wait for each line of
#                                            output (default 0: no limit)
#   cmd                                      the needleman_wunsch to run

use strict;
use warnings;

use File::Basename; # dirname
use Carp;
use IPC::Open3;
use Symbol qw(gensym);

my %FLAGS = map { $_ => 1 } qw(case_sensitive freestartgap freeendgap linear
                               linearfallback nogaps nogapsin1 nogapsin2
                               nogapsinlongest nomismatches);
my %VALUES = map { $_ => 1 } qw(match mismatch gapopen gapextend
                                substitution_matrix substitution_pairs
                                threads);
my @SCORING = qw(PAM30 PAM70 BLOSUM80 BLOSUM62 DNA_HYBRIDIZATION);

sub new
{
  my ($class, %given) = @_;
  my %options = map { lc($_) => $given{$_} } keys %given;

  my @args;
  my $cmd = dirname(__FILE__)."/../bin/needleman_wunsch";
  my $timeout = 0;

  if(defined($options{'match'}) != defined($options{'mismatch'})) {
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
    _timeout => $timeout,
    _align_number => 0,
    _procs => {}, # running tools, by mode: 'align', 'score'
  };

  return bless($self, $class);
}

# The tool for mode ('align' or 'score'), started the first time
sub _proc
{
  my ($self, $mode) = @_;
  if(!$self->{_procs}{$mode})
  {
    my @extra = $mode eq 'score' ? ('--scoreonly')
                                 : ('--pretty', '--printscores');
    my ($to, $from) = (undef, gensym());
    # Standard error comes back with the output (its errors end a pair)
    my $pid = eval { open3($to, $from, undef, $self->{_cmd},
                           @{$self->{_args}}, @extra, '--stdin') };
    croak("Cannot run $self->{_cmd}: $@") if !$pid;
    $to->autoflush(1);
    $self->{_procs}{$mode} = { pid => $pid, to => $to, from => $from };
  }
  return $self->{_procs}{$mode};
}

sub _stop
{
  my ($self, $mode) = @_;
  my $p = delete($self->{_procs}{$mode}) or return;
  close($p->{to});
  close($p->{from});
  waitpid($p->{pid}, 0);
}

# Lines of output with no line breaks; croaks on the tool's errors
sub _read_line
{
  my ($self, $from, $mode) = @_;
  my $line;
  eval {
    local $SIG{ALRM} = sub { die("timeout\n") };
    alarm($self->{_timeout});
    $line = <$from>;
    alarm(0);
  };
  alarm(0);
  if($@) {
    $self->_stop($mode) if defined($mode);
    croak($@ eq "timeout\n" ? "Timeout reading needleman_wunsch output"
                            : "Cannot read needleman_wunsch output: $@");
  }
  if(!defined($line)) {
    $self->_stop($mode) if defined($mode);
    croak("needleman_wunsch stopped (have you compiled?)");
  }
  $line =~ s/\r?\n$//;
  if($line =~ /^(Alignment )?(Error|Warning):/) {
    $self->_stop($mode) if defined($mode);
    croak($line);
  }
  return $line;
}

sub _check
{
  my ($seq1, $seq2) = @_;
  for my $seq ($seq1, $seq2) {
    croak("Sequences must be defined") if !defined($seq);
    croak("Whitespace or control characters not allowed in sequences")
      if $seq =~ /[\x00-\x20\x7f]/;
    # The tool would read these as a FASTA or FASTQ header and wait for
    # the rest of the record
    croak("Sequences cannot start with '>' or '@'") if $seq =~ /^[>@]/;
  }
}

# Lines of the tool's output for one pair given on its command line: for
# pairs with an empty sequence, which its line-by-line input cannot carry.
# The pair comes after "--", so a sequence may start with '-'
sub _run_once
{
  my ($self, @extra) = @_;
  my ($to, $from) = (undef, gensym());
  my $pid = eval { open3($to, $from, undef, $self->{_cmd},
                         @{$self->{_args}}, @extra) };
  croak("Cannot run $self->{_cmd}: $@") if !$pid;
  close($to);
  my @lines;
  my $error;
  eval {
    local $SIG{ALRM} = sub { die("timeout\n") };
    alarm($self->{_timeout});
    while(defined(my $line = <$from>)) {
      $line =~ s/\r?\n$//;
      if($line =~ /^(Alignment )?(Error|Warning):/) { $error = $line; last; }
      push(@lines, $line);
    }
    alarm(0);
  };
  alarm(0);
  my $failed = $@;
  close($from);
  kill('TERM', $pid) if $failed || defined($error);
  waitpid($pid, 0);
  croak($failed eq "timeout\n" ? "Timeout reading $self->{_cmd} output"
                               : "Cannot read $self->{_cmd} output: $failed")
    if $failed;
  croak($error) if defined($error);
  return @lines;
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
  _check($seq1, $seq2);
  my %result = ('seq1' => $seq1, 'seq2' => $seq2,
                'number' => $self->{_align_number}++);
  my @lines;

  if($seq1 eq '' || $seq2 eq '') {
    @lines = $self->_run_once('--pretty', '--printscores', '--', $seq1, $seq2);
  }
  else {
    my $p = $self->_proc('align');
    print {$p->{to}} "$seq1\n$seq2\n";
    @lines = map { $self->_read_line($p->{from}, 'align') } 1..5;
  }

  @result{'align1', 'sep', 'align2'} = @lines[0..2];
  $result{'score'} = _score_of($lines[3]);
  return \%result;
}

# The best global alignment's score alone (one pass, linear memory)
sub do_score
{
  my ($self, $seq1, $seq2) = @_;
  _check($seq1, $seq2);
  if($seq1 eq '' || $seq2 eq '') {
    return _score_of(($self->_run_once('--scoreonly', '--', $seq1, $seq2))[0]);
  }
  my $p = $self->_proc('score');
  print {$p->{to}} "$seq1\n$seq2\n";
  my $score = _score_of($self->_read_line($p->{from}, 'score'));
  $self->_read_line($p->{from}, 'score');
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

# Stop the tool (also done when the object goes)
sub destructor
{
  my ($self) = @_;
  $self->_stop($_) for keys %{$self->{_procs}};
}

sub DESTROY
{
  my ($self) = @_;
  local ($?, $@);
  $self->destructor();
}

1;
