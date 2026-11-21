#!/usr/bin/env perl

# perl/sw_example.pl
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

use SmithWaterman;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./sw_example.pl [seq1 seq2]\n";

  exit(1);
}

print_usage() if @ARGV != 0 && @ARGV != 2;

my $sw = SmithWaterman->new('minscore' => 2);

my @aligns = @ARGV ? ([@ARGV])
                   : (['asdfasdf','asdfasdf'],
                      ['dogg', 'ggod'],
                      ['asdfas', 'fadsas'],
                      ['hello_henry', 'wellohenry'],
                      ['asdf', 'lkj']);

for my $align (@aligns)
{
  print "Doing '$align->[0]' vs '$align->[1]'\n";
  $sw->do_alignment(@$align);

  # The best hit, then how many others there were
  if(defined(my $hit = $sw->get_next_hit()))
  {
    print "hit:\n";
    print "  ".$hit->{'align1'}."\n";
    print "  ".$hit->{'sep'}."\n";
    print "  ".$hit->{'align2'}."\n";
    print "  (and ".scalar($sw->get_all_hits())." more)\n";
  }
}

$sw->destructor();
