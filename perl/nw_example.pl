#!/usr/bin/env perl

# perl/nw_example.pl
# url: https://github.com/DanielTillett/nucleics-align
# original author: Isaac Turner <turner.isaac@gmail.com>
# license: BSD 3-Clause, see LICENSE
# date: Oct 2012

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use NeedlemanWunsch;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./nw_example.pl [seq1 seq2]\n";

  exit(1);
}

print_usage() if @ARGV != 0 && @ARGV != 2;

my $nw = NeedlemanWunsch->new();

my @aligns = @ARGV ? ([@ARGV])
                   : (['asdfasdf','asdfasdf'],
                      ['dogg', 'ggod'],
                      ['asdfas', 'fadsas'],
                      ['', ''],
                      ['asdf', 'lkj']);

for my $align (@aligns)
{
  print "Doing '$align->[0]' vs '$align->[1]'\n";
  my $result = $nw->do_alignment(@$align);
  $nw->print_alignment($result);
}

# The score alone, in linear memory; and a linear-memory alignment
my ($s1, $s2) = @{$aligns[-1]};
print "score only: ".$nw->do_score($s1, $s2)."\n";
print "linear memory: ".
      NeedlemanWunsch->new(linear => 1)->do_alignment($s1, $s2)->{score}."\n";

$nw->destructor();
