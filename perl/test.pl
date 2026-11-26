#!/usr/bin/env perl

# perl/test.pl: tests of NeedlemanWunsch.pm and SmithWaterman.pm against
# the built tools. Run after make: perl perl/test.pl
# url: https://github.com/DanielTillett/nucleics-align
# license: BSD 3-Clause, see LICENSE

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;
use Test::More;

use NeedlemanWunsch;
use SmithWaterman;

my $input = "$FindBin::Bin/../tests/cli/input";

# Needleman-Wunsch
{
  my $nw = NeedlemanWunsch->new();
  my $aln = $nw->do_alignment('ACGT', 'AGT');
  is_deeply([@$aln{qw(align1 sep align2 score number)}],
            ['ACGT', '| ||', 'A-GT', -2, 0], 'one pair');
  is($nw->do_alignment('ACGT', 'ACGT')->{number}, 1, 'numbered from 0');
  is($nw->do_score('ACGT', 'AGT'), -2, 'score only');

  # Empty sequences, with the options in force (free end gaps here)
  is($nw->do_alignment('', '')->{score}, 0, 'both empty');
  my $empty = $nw->do_alignment('ACGT', '');
  is_deeply([@$empty{qw(align1 align2 score)}], ['ACGT', '----', -8],
            'one empty');
  is(NeedlemanWunsch->new(freestartgap => 1)->do_alignment('', 'ACGT')
       ->{score}, 0, 'empty with a free start gap');
  is($nw->do_score('', 'AG'), -6, 'score only, empty');

  # Options: mixed case once, linear memory, scoring by name, wildcards
  is(NeedlemanWunsch->new(GapOpen => -10)->do_alignment('ACGT', 'AGT')
       ->{score}, -8, 'mixed-case option');
  is(NeedlemanWunsch->new(linear => 1, threads => 2)
       ->do_alignment('ACGTTGCA', 'ACGTGCA')->{score},
     $nw->do_alignment('ACGTTGCA', 'ACGTGCA')->{score}, 'linear memory');
  my ($p, $q) = ('MKTAYIAKQRQISFVKSHFSRQ', 'MKTAHIAKQRQISFVKSHFSRQ');
  is(NeedlemanWunsch->new(scoring => 'blosum62')->do_alignment($p, $q)
       ->{score},
     NeedlemanWunsch->new(substitution_matrix =>
       "$FindBin::Bin/../scoring/BLOSUM62.txt")->do_alignment($p, $q)
       ->{score}, 'scoring by name and from a file');
  is(NeedlemanWunsch->new(wildcard => {N => 1})->do_alignment('ACNT', 'ACGT')
       ->{score}, 4, 'wildcard');
  is(NeedlemanWunsch->new(freestartgap => 1, freeendgap => 1)
       ->do_alignment('TTTTACGTTTTT', 'ACGT')->{score}, 4, 'free end gaps');

  # Errors
  eval { NeedlemanWunsch->new(scoring => 'NOPE') };
  like($@, qr/Invalid scoring/, 'bad scoring name');
  eval { NeedlemanWunsch->new(nonsense => 1) };
  like($@, qr/Unknown option/, 'unknown option');
  eval { $nw->do_alignment("AC GT", 'ACGT') };
  like($@, qr/Whitespace/, 'whitespace in a sequence');
  my $pairs = NeedlemanWunsch->new(substitution_pairs => "$input/ac_pairs.txt");
  is($pairs->do_alignment('ACCA', 'ACA')->{score}, -2, 'pairs file');
  eval { $pairs->do_alignment('ACGT', 'AGT') };
  like($@, qr/Unknown character pair/, 'tool error');
  is($pairs->do_alignment('CAAC', 'CAC')->{score}, -2, 'restarts after one');
  # The tool reads a leading '>' or '@' as a FASTA/FASTQ header and would
  # wait for the rest of the record: refused before it is sent
  my $quick = NeedlemanWunsch->new(timeout => 5);
  eval { $quick->do_alignment('>ACGT', 'ACGT') };
  like($@, qr/cannot start with/, 'leading > refused');
  eval { $quick->do_alignment('ACGT', '@ACGT') };
  like($@, qr/cannot start with/, 'leading @ refused');
  is($quick->do_alignment('-A', '')->{align1}, '-A', 'leading dash, empty');
  is($quick->do_alignment('-A', 'A')->{align1}, '-A', 'leading dash');
  $quick->destructor();
  # The tool's other ways to the same score
  my $plain = $nw->do_alignment('ACGTTGCAACGT', 'ACGTGCAACGT')->{score};
  is(NeedlemanWunsch->new(wfa => 1)->do_alignment('ACGTTGCAACGT',
       'ACGTGCAACGT')->{score}, $plain, 'wfa');
  is(NeedlemanWunsch->new(wfa => 1, wfabudget => 100000)->do_alignment(
       'ACGTTGCAACGT', 'ACGTGCAACGT')->{score}, $plain, 'wfabudget');
  is(NeedlemanWunsch->new(linearabove => 5)->do_alignment('ACGTTGCAACGT',
       'ACGTGCAACGT')->{score}, $plain, 'linearabove');
  $nw->destructor();
}

# Smith-Waterman
{
  my $sw = SmithWaterman->new(maxhits => 2);
  ok($sw->do_alignment('acgtacgtacgt', 'ttacgtacgtaa'), 'aligned');
  my $hit = $sw->get_next_hit();
  is_deeply([@$hit{qw(hit score align1 sep align2 pos1 len1 pos2 len2)}],
            [0, 18, 'acgtacgta', '|||||||||', 'acgtacgta', 0, 9, 2, 9],
            'best hit');
  my $printed = '';
  open(my $out, '>', \$printed) or die;
  $sw->print_hit($hit, $out);
  close($out);
  is($printed, "hit 0.0 score: 18\n  acgtacgta  [pos: 0; len: 9]\n".
               "  |||||||||\n  acgtacgta  [pos: 2; len: 9]\n", 'print_hit');
  is(scalar($sw->get_all_hits()), 1, 'maxhits');
  ok(!defined($sw->get_next_hit()), 'no more hits');

  # A new alignment skips hits left unread
  $sw->do_alignment('acgtacgtacgt', 'ttacgtacgtaa');
  $sw->get_next_hit();
  $sw->do_alignment('AAAACCCCGGGG', 'CCCC');
  my @hits = $sw->get_all_hits();
  is($hits[0]{score}, 8, 'next alignment');

  {
    local $SIG{__WARN__} = sub {};
    ok(!$sw->do_alignment('', 'ACGT'), 'empty: no alignment');
  }
  ok(!defined($sw->get_next_hit()), 'empty: no hits');
  is(scalar(SmithWaterman->new(minscore => 9)->do_alignment('acgtacgt',
       'acgtacgt') && 1), 1, 'minscore accepted');
  eval { SmithWaterman->new(Scoring => 'nope') };
  like($@, qr/Invalid scoring/, 'bad scoring name');
  $sw->destructor();

  # After a tool error the next alignment starts a fresh tool
  my $pairs = SmithWaterman->new(substitution_pairs => "$input/ac_pairs.txt",
                                 minscore => 1, timeout => 5);
  eval { $pairs->do_alignment('ACGT', 'AGT') };
  like($@, qr/Unknown character pair/, 'sw tool error');
  ok($pairs->do_alignment('ACCA', 'ACA'), 'sw restarts after one');
  is($pairs->get_next_hit()->{score}, 2, 'sw hit after restart');
  eval { $pairs->do_alignment('>ACGT', 'ACGT') };
  like($@, qr/cannot start with/, 'sw leading > refused');
  $pairs->destructor();
  ok(SmithWaterman->new()->do_alignment('-ACCA', 'ACA'), 'sw leading dash');
  my $fast = SmithWaterman->new(maxhits => 1, fastlocal => 1);
  ok($fast->do_alignment('acgtacgtacgt', 'ttacgtacgtaa'), 'fastlocal');
  is($fast->get_next_hit()->{score}, 18, 'fastlocal best hit');
  $fast->destructor();
}

# A tool that stops answering: the read times out and the tool is killed
# rather than waited for (a tool mid-alignment would not notice its input
# closing)
SKIP: {
  skip('needs a POSIX shell', 6) if $^O eq 'MSWin32';
  require File::Temp;
  # The tool ignores TERM too, so only the escalation to KILL ends it
  my ($fh, $path) = File::Temp::tempfile(UNLINK => 1);
  print $fh "#!/bin/sh\ntrap '' TERM\nexec sleep 30\n";
  close($fh);
  chmod(0755, $path);
  for my $class ('NeedlemanWunsch', 'SmithWaterman') {
    my $tool = $class->new(cmd => $path, timeout => 1);
    my $start = time();
    eval { $tool->do_alignment('ACGT', 'AGT'); $tool->get_next_hit() };
    like($@, qr/Timeout reading/, "$class: timeout");
    $tool->destructor();
    cmp_ok(time() - $start, '<', 10, "$class: the stalled tool was killed");
  }
  # The same for a pair with an empty sequence (run once, on the command
  # line)
  my $start = time();
  eval { NeedlemanWunsch->new(cmd => $path, timeout => 1)->do_alignment('', 'ACGT') };
  like($@, qr/Timeout reading/, 'empty pair: timeout');
  cmp_ok(time() - $start, '<', 10, 'empty pair: the stalled tool was killed');
}

done_testing();
