package NucleicsAlign::Tool;

# perl/NucleicsAlign/Tool.pm
# url: https://github.com/DanielTillett/nucleics-align
# license: BSD 3-Clause, see LICENSE
# date: Nov 2026
#
# What NeedlemanWunsch.pm and SmithWaterman.pm share: their options read
# into the tool's command line, and the tool run over pipes. Each module
# gives a spec (the tool's name, its default path, its flags and options
# with values) to _new() and starts its tools with _start(), under a key;
# then _send() writes a line to one, _read_line() reads one back (croaking
# on the tool's errors or a timeout, after stopping it) and _stop() stops
# it. The objects stop their tools when they go.

use strict;
use warnings;

use Carp;
use IPC::Open3;
use Symbol qw(gensym);

my @SCORING = qw(PAM30 PAM70 BLOSUM80 BLOSUM62 DNA_HYBRIDIZATION);

# An object for the tool in spec, from the options given (names ignore
# case): spec->{flags} are set by a true value, spec->{values} take one;
# scoring, wildcard, timeout and cmd are understood for every tool
sub _new
{
  my ($class, $spec, %given) = @_;
  my %options = map { lc($_) => $given{$_} } keys %given;
  my %flags = map { $_ => 1 } @{$spec->{flags}};
  my %values = map { $_ => 1 } @{$spec->{values}};

  my @args;
  my $cmd = $spec->{cmd};
  my $timeout = 0;

  # The tools want both or neither, unless mismatches are not allowed
  if(!$options{'nomismatches'} &&
     defined($options{'match'}) != defined($options{'mismatch'})) {
    carp("Setting only one of match/mismatch");
  }

  for my $key (sort keys %options)
  {
    my $value = $options{$key};
    if($flags{$key}) {
      push(@args, "--$key") if $value;
    }
    elsif($values{$key}) {
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
    _name => $spec->{name},
    _cmd => $cmd,
    _args => \@args,
    _timeout => $timeout,
    _procs => {}, # running tools, by key: { pid, to, from }
  };
  return bless($self, $class);
}

# Start the tool under key, with its options then @extra. Standard error
# comes back with the output (its errors end a pair)
sub _start
{
  my ($self, $key, @extra) = @_;
  my ($to, $from) = (undef, gensym());
  my $pid = eval { open3($to, $from, undef, $self->{_cmd},
                         @{$self->{_args}}, @extra) };
  croak("Cannot run $self->{_cmd}: $@") if !$pid;
  $to->autoflush(1);
  return $self->{_procs}{$key} = { pid => $pid, to => $to, from => $from };
}

# Text to the tool under key. If it has gone, the write fails quietly (the
# next read then reports it) instead of SIGPIPE ending this process
sub _send
{
  my ($self, $key, $text) = @_;
  local $SIG{PIPE} = 'IGNORE';
  print {$self->{_procs}{$key}{to}} $text;
}

# Stop the tool under key. kill: it may be mid-alignment (a timeout), so it
# would not notice its input closing; terminate it before waiting
sub _stop
{
  my ($self, $key, $kill) = @_;
  my $p = delete($self->{_procs}{$key}) or return;
  kill('TERM', $p->{pid}) if $kill;
  # The tool may already have gone: writing to it must not kill this process
  local $SIG{PIPE} = 'IGNORE';
  close($p->{to});
  close($p->{from});
  waitpid($p->{pid}, 0);
}

# Called with a tool's error line before it is stopped and the error
# croaked; for a module to add what it knows
sub _tool_error { }

# A line of the tool's output under key, with no line break. Croaks, after
# stopping the tool, on a timeout, the end of its output or its errors
sub _read_line
{
  my ($self, $key) = @_;
  my $from = $self->{_procs}{$key}{from};
  my $line;
  eval {
    local $SIG{ALRM} = sub { die("timeout\n") };
    alarm($self->{_timeout});
    $line = <$from>;
    alarm(0);
  };
  alarm(0);
  if($@) {
    $self->_stop($key, 1);
    croak($@ eq "timeout\n" ? "Timeout reading $self->{_name} output"
                            : "Cannot read $self->{_name} output: $@");
  }
  if(!defined($line)) {
    $self->_stop($key);
    croak("$self->{_name} stopped (have you compiled?)");
  }
  $line =~ s/\r?\n$//;
  if($line =~ /^(Alignment )?(Error|Warning):/) {
    $self->_tool_error($line);
    $self->_stop($key);
    croak($line);
  }
  return $line;
}

# Croaks on sequences the tool cannot be given a line at a time
sub _check_seqs
{
  my ($self, @seqs) = @_;
  for my $seq (@seqs) {
    croak("Sequences must be defined") if !defined($seq);
    croak("Whitespace or control characters not allowed in sequences")
      if $seq =~ /[\x00-\x20\x7f]/;
    # The tool would read these as a FASTA or FASTQ header and wait for
    # the rest of the record
    croak("Sequences cannot start with '>' or '@'") if $seq =~ /^[>@]/;
  }
}

# Lines of the tool's output for one pair given on its command line (in
# @extra): for pairs with an empty sequence, which its line-by-line input
# cannot carry
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

# Stop the tools (also done when the object goes)
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
