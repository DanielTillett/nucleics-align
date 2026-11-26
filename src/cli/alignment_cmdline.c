/*
 alignment_cmdline.c
 url: https://github.com/DanielTillett/nucleics-align
 original author: Isaac Turner <turner.isaac@gmail.com>
 license: BSD 3-Clause, see LICENSE
 date: Nov 2013
 */

// Turn on debugging output by defining NUCLEICS_ALIGN_VERBOSE
//#define NUCLEICS_ALIGN_VERBOSE

// request decent POSIX version (and, on macOS, which hides names beyond it
// once one is asked for, sysconf(_SC_NPROCESSORS_ONLN) and the like)
#define _XOPEN_SOURCE 700
#define _DEFAULT_SOURCE
#define _DARWIN_C_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h> // strcasecmp
#include <ctype.h> // isspace
#include <errno.h> // errno, ERANGE
#include <limits.h> // INT_MIN
#include <stdarg.h> // for va_list
#include <pthread.h>
#include <unistd.h> // sysconf

#include "cli_seqio.h"

#include "alignment.h"
#include "alignment_cmdline.h"
#include "alignment_scoring_load.h"
#include "alignment_macros.h"

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);
  if(len == 0) return 0; // strtol would read "" as 0

  errno = 0;
  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(errno == ERANGE || tmp > INT_MAX || tmp < INT_MIN ||
     strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);
  if(len == 0) return 0;

  const char *p = str;
  while(isspace((unsigned char)*p)) p++;
  if(*p == '-' || *p == '\0') return 0;

  errno = 0;
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(errno == ERANGE || tmp > UINT_MAX || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (unsigned int)tmp;
    return 1;
  }
}

// Larger scores could overflow once added together (see SCORE_ABS_MAX)
static bool score_in_range(int score)
{
  return score >= -SCORE_ABS_MAX && score <= SCORE_ABS_MAX;
}

static void vprint_usage(enum NucleicsAlignCmdType cmd_type,
                         const score_t defaults[4], const char *cmdstr,
                         const char *errfmt, va_list argptr)
  __attribute__((noreturn));

// The usage text, after the error given (NULL: none), then exit
static void vprint_usage(enum NucleicsAlignCmdType cmd_type,
                         const score_t defaults[4], const char *cmdstr,
                         const char *errfmt, va_list argptr)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    vfprintf(stderr, errfmt, argptr);
    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
  }

  fprintf(stderr, "usage: %s [OPTIONS] [--] [seq1 seq2]\n", cmdstr);

  fprintf(stderr,
"  %s optimal %s alignment (maximises score).  \n"
"  Takes a pair of sequences on the command line, or can read from a\n"
"  file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.\n\n",
          cmd_type == NUCLEICS_ALIGN_SW_CMD ? "Smith-Waterman" : "Needleman-Wunsch",
          cmd_type == NUCLEICS_ALIGN_SW_CMD ? "local" : "global");

  fprintf(stderr,
"  OPTIONS:\n"
"    --file <file>        Sequence file reading with gzip support - read two\n"
"                         sequences at a time and align them\n"
"    --files <f1> <f2>    Read one sequence from each file to align at one time\n"
"    --stdin              Read pairs from STDIN one at a time, interactively\n"
"                         (for a batch, use '--file -')\n"
"    --                   End of options: the two sequences follow, even\n"
"                         if one starts with '-'\n"
"\n"
"    --case_sensitive     Use case sensitive character comparison [default: off]\n"
"\n"
"    --match <score>      [default: %i]\n"
"    --mismatch <score>   [default: %i]\n"
"    --gapopen <score>    [default: %i]\n"
"    --gapextend <score>  [default: %i]\n"
"\n"
"    --scoring <%s>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n"
"    --wildcard <w> <s>   Character <w> matches all characters with score <s>\n"
"\n"
"    --threads <n>        Threads to use, within an alignment and for aligning\n"
"                         pairs from a file at once [default: half the CPUs]\n\n",
          defaults[0], defaults[1],
          defaults[2], defaults[3], scoring_system_names());

  if(cmd_type == NUCLEICS_ALIGN_SW_CMD)
  {
    // SW specific
    fprintf(stderr,
"    --minscore <score>   Minimum required score\n"
"                         [default: match * MAX(0.2 * length, 2)]\n"
"    --maxhits <hits>     Maximum number of results per alignment\n"
"                         [default: no limit]\n"
"    --fastlocal          Opt-in scan + BiWFA/linear sub-rectangle traceback\n"
"                         for single-hit runs (use with --maxhits 1; exact\n"
"                         optimal score, may pick another equal-scoring hit)\n"
"\n"
"    --context <n>        Print <n> bases of context\n"
"    --printseq           Print sequences before local alignments\n");
  }
  else
  {
    // NW specific
    fprintf(stderr,
"\n"
"    --freestartgap       No penalty for gap at start of alignment\n"
"    --freeendgap         No penalty for gap at end of alignment\n"
"\n"
"    --printscores        Print optimal alignment scores\n"
"    --zam                A funky type of output\n"
"    --linear             Align in linear memory (Hirschberg, Myers-Miller)\n"
"    --scoreonly          Print only the score (fast, linear memory)\n"
"    --linearfallback     Align in linear memory where too big for the matrices\n"
"                         (same score, but may pick another equal-scoring\n"
"                         alignment)\n"
"    --linearabove <n>    Align in linear memory where either sequence has at\n"
"                         least <n> letters (as --linearfallback)\n"
"    --wfa                Align similar pairs by the wavefront algorithm, much\n"
"                         faster for near-identical ones (same score, but may\n"
"                         pick another equal-scoring alignment)\n"
"    --wfabudget <n>      Diagonal-update cap for --wfa (default: automatic)\n");
  }

  fprintf(stderr,
"    --printmatrices      Print dynamic programming matrices\n"
"    --printfasta         Print fasta header lines\n"
"    --pretty             Print with a descriptor line\n"
"    --colour             Print with colour\n"
"\n"
"  Experimental Options:\n"
"    --nogapsin1          No gaps allowed within the first sequence\n"
"    --nogapsin2          No gaps allowed within the second sequence\n"
"    --nogaps             No gaps allowed in either sequence\n"
"    --nogapsinlongest    No gaps allowed within the longer sequence\n");

  fprintf(stderr,
"    --nomismatches       No mismatches allowed%s\n",
          cmd_type == NUCLEICS_ALIGN_SW_CMD ? "" : " (cannot be used with --nogaps..)");

  fprintf(stderr,
"\n"
" DETAILS:\n"
"  * For help choosing scoring, see docs/scoring.md.\n"
"  * Gap (of length N) penalty is: (open+N*extend)\n"
"  * To do alignment without affine gap penalty, set '--gapopen 0'.\n"
"  * Scoring files should be matrices, with entries separated by a single\n"
"    character or whitespace. See files in the 'scoring' directory for examples.\n"
"\n"
"  https://github.com/DanielTillett/nucleics-align  (compiled: %s %s)\n", __DATE__, __TIME__);

  exit(EXIT_FAILURE);
}

void cmdline_free(cmdline_t *cmd)
{
  free(cmd->file_paths1);
  free(cmd->file_paths2);
  free(cmd);
}

// Everything an option's handler needs while the arguments are read
typedef struct {
  cmdline_t *cmd;
  scoring_t *scoring;
  enum NucleicsAlignCmdType cmd_type;
  score_t defaults[4];
  const char *cmdstr;
  bool scoring_set, substitutions_set, match_set, mismatch_set;
} parse_ctx_t;

static void ctx_usage(const parse_ctx_t *c, const char *errfmt, ...)
  __attribute__((format(printf, 2, 3)))
  __attribute__((noreturn));

static void ctx_usage(const parse_ctx_t *c, const char *errfmt, ...)
{
  va_list argptr;
  va_start(argptr, errfmt);
  vprint_usage(c->cmd_type, c->defaults, c->cmdstr, errfmt, argptr);
}

// An option: its name, how many values follow it, which tool takes it
// (TOOL_ANY: both), whether the first pass over the arguments applies it
// (--help, and the scoring system and case sensitivity, which the other
// options build on), and what it does with its values v[0], v[1]
enum { TOOL_ANY, TOOL_NW, TOOL_SW };

typedef struct {
  const char *name;
  int values;
  int tool;
  bool early;
  void (*apply)(parse_ctx_t *c, char *const *v);
} option_t;

// --- Handlers: the scoring's rules ---

static void opt_freestartgap(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_start_gap_penalty = true; }
static void opt_freeendgap(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_end_gap_penalty = true; }
static void opt_nogaps(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_gaps_in_a = c->scoring->no_gaps_in_b = true; }
static void opt_nogapsin1(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_gaps_in_a = true; }
static void opt_nogapsin2(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_gaps_in_b = true; }
static void opt_nogapsinlongest(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_gaps_in_longest = true; }
static void opt_nomismatches(parse_ctx_t *c, char *const *v)
{ (void)v; c->scoring->no_mismatches = true; }

// --- Handlers: the first pass's options ---

static void opt_help(parse_ctx_t *c, char *const *v)
{ (void)v; ctx_usage(c, NULL); }
static void opt_case_sensitive(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->case_sensitive = 1; }
static void opt_scoring(parse_ctx_t *c, char *const *v)
{
  if(c->scoring_set)
    ctx_usage(c, "More than one scoring system specified - not permitted");
  if(!scoring_system_by_name(c->scoring, v[0]))
    ctx_usage(c, "Unknown --scoring choice, not one of %s",
              scoring_system_names());
  c->scoring_set = true;
}

// --- Handlers: output and modes ---

static void opt_printseq(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_seq = true; }
static void opt_printmatrices(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_matrices = true; }
static void opt_printscores(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_scores = true; }
static void opt_printfasta(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_fasta = true; }
static void opt_pretty(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_pretty = true; }
static void opt_colour(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->print_colour = true; }
static void opt_zam(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->zam_stle_output = true; }
static void opt_scoreonly(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->score_only = true; }
static void opt_linear(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->linear_space = true; }
static void opt_linearfallback(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->linear_fallback = true; }
static void opt_wfa(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->wfa_align = true; }
static void opt_fastlocal(parse_ctx_t *c, char *const *v)
{ (void)v; c->cmd->fast_local = true; }
static void opt_stdin(parse_ctx_t *c, char *const *v)
{
  (void)v;
  // Similar to --file below
  cmdline_add_files(c->cmd, "", NULL);
  c->cmd->interactive = true;
}

// --- Handlers: options with values ---

static void opt_substitution_matrix(parse_ctx_t *c, char *const *v)
{
  gzFile file = gzopen(v[0], "r");
  if(file == NULL) ctx_usage(c, "Couldn't read: %s", v[0]);
  align_scoring_load_matrix(file, v[0], c->scoring);
  gzclose(file);
  c->substitutions_set = true;
}
static void opt_substitution_pairs(parse_ctx_t *c, char *const *v)
{
  gzFile file = gzopen(v[0], "r");
  if(file == NULL) ctx_usage(c, "Couldn't read: %s", v[0]);
  align_scoring_load_pairwise(file, v[0], c->scoring);
  gzclose(file);
  c->substitutions_set = true;
}
static void opt_minscore(parse_ctx_t *c, char *const *v)
{
  // 0 and 1 alike: every hit with a positive score
  if(!parse_entire_int(v[0], &c->cmd->min_score) || c->cmd->min_score < 0)
    ctx_usage(c, "Invalid --minscore <score> argument (must be an int >= 0)");
  c->cmd->min_score_set = true;
}
static void opt_maxhits(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_uint(v[0], &c->cmd->max_hits_per_alignment) ||
     c->cmd->max_hits_per_alignment == 0)
    ctx_usage(c, "Invalid --maxhits <hits> argument (must be at least 1)");
  c->cmd->max_hits_per_alignment_set = true;
}
static void opt_linearabove(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_uint(v[0], &c->cmd->linear_above) ||
     c->cmd->linear_above == 0)
    ctx_usage(c, "Invalid --linearabove <n> argument (must be at least 1)");
}
static void opt_wfabudget(parse_ctx_t *c, char *const *v)
{
  const char *arg = v[0];
  if(arg[0] == '\0' || !isdigit((unsigned char)arg[0]))
    ctx_usage(c, "Invalid --wfabudget <n> argument (must be a positive integer)");
  char *end = NULL;
  errno = 0;
  unsigned long long n = strtoull(arg, &end, 10);
  if(errno == ERANGE || end == NULL || *end != '\0' || n == 0 ||
     n > (unsigned long long)SIZE_MAX)
    ctx_usage(c, "Invalid --wfabudget <n> argument (must be a positive integer)");
  c->cmd->wfa_budget = (size_t)n;
}
static void opt_threads(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_uint(v[0], &c->cmd->threads) ||
     c->cmd->threads == 0 || c->cmd->threads > 4096)
    ctx_usage(c, "Invalid --threads <n> argument (must be 1 to 4096)");
}
static void opt_context(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_uint(v[0], &c->cmd->print_context))
    ctx_usage(c, "Invalid --context <c> argument (must be >= 0)");
}
static void opt_match(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_int(v[0], &c->scoring->match))
    ctx_usage(c, "Invalid --match argument ('%s') must be an int", v[0]);
  c->match_set = true;
}
static void opt_mismatch(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_int(v[0], &c->scoring->mismatch))
    ctx_usage(c, "Invalid --mismatch argument ('%s') must be an int", v[0]);
  c->mismatch_set = true;
}
static void opt_gapopen(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_int(v[0], &c->scoring->gap_open))
    ctx_usage(c, "Invalid --gapopen argument ('%s') must be an int", v[0]);
}
static void opt_gapextend(parse_ctx_t *c, char *const *v)
{
  if(!parse_entire_int(v[0], &c->scoring->gap_extend))
    ctx_usage(c, "Invalid --gapextend argument ('%s') must be an int", v[0]);
}
static void opt_file(parse_ctx_t *c, char *const *v)
{ cmdline_add_files(c->cmd, v[0], NULL); }
static void opt_files(parse_ctx_t *c, char *const *v)
{
  // "- -": both from stdin
  if(strcmp(v[0], "-") == 0 && strcmp(v[1], "-") == 0)
    cmdline_add_files(c->cmd, v[0], NULL);
  else
    cmdline_add_files(c->cmd, v[0], v[1]);
}
static void opt_wildcard(parse_ctx_t *c, char *const *v)
{
  int wildscore = 0;
  if(strlen(v[0]) != 1 || !parse_entire_int(v[1], &wildscore))
    ctx_usage(c, "--wildcard <w> <s> takes a single character and a number");
  if(!score_in_range(wildscore))
    ctx_usage(c, "Scores must be between %i and %i", -SCORE_ABS_MAX,
              SCORE_ABS_MAX);
  scoring_add_wildcard(c->scoring, v[0][0], wildscore);
}

// The options, matched ignoring case. An unknown option is an error, as is
// a value-taking option without its values
static const option_t options[] = {
  {"--help",                0, TOOL_ANY, true,  opt_help},
  {"-help",                 0, TOOL_ANY, true,  opt_help},
  {"-h",                    0, TOOL_ANY, true,  opt_help},
  {"--case_sensitive",      0, TOOL_ANY, true,  opt_case_sensitive},
  {"--scoring",             1, TOOL_ANY, true,  opt_scoring},
  {"--freestartgap",        0, TOOL_NW,  false, opt_freestartgap},
  {"--freeendgap",          0, TOOL_NW,  false, opt_freeendgap},
  {"--nogaps",              0, TOOL_ANY, false, opt_nogaps},
  {"--nogapsin1",           0, TOOL_ANY, false, opt_nogapsin1},
  {"--nogapsin2",           0, TOOL_ANY, false, opt_nogapsin2},
  {"--nogapsinlongest",     0, TOOL_ANY, false, opt_nogapsinlongest},
  {"--nomismatches",        0, TOOL_ANY, false, opt_nomismatches},
  {"--printseq",            0, TOOL_SW,  false, opt_printseq},
  {"--printmatrices",       0, TOOL_ANY, false, opt_printmatrices},
  {"--printscores",         0, TOOL_NW,  false, opt_printscores},
  {"--printfasta",          0, TOOL_ANY, false, opt_printfasta},
  {"--pretty",              0, TOOL_ANY, false, opt_pretty},
  {"--colour",              0, TOOL_ANY, false, opt_colour},
  {"--zam",                 0, TOOL_NW,  false, opt_zam},
  {"--scoreonly",           0, TOOL_NW,  false, opt_scoreonly},
  {"--linear",              0, TOOL_NW,  false, opt_linear},
  {"--linearfallback",      0, TOOL_NW,  false, opt_linearfallback},
  {"--wfa",                 0, TOOL_NW,  false, opt_wfa},
  {"--fastlocal",           0, TOOL_SW,  false, opt_fastlocal},
  {"--stdin",               0, TOOL_ANY, false, opt_stdin},
  {"--substitution_matrix", 1, TOOL_ANY, false, opt_substitution_matrix},
  {"--substitution_pairs",  1, TOOL_ANY, false, opt_substitution_pairs},
  {"--minscore",            1, TOOL_SW,  false, opt_minscore},
  {"--maxhits",             1, TOOL_SW,  false, opt_maxhits},
  {"--linearabove",         1, TOOL_NW,  false, opt_linearabove},
  {"--wfabudget",           1, TOOL_NW,  false, opt_wfabudget},
  {"--threads",             1, TOOL_ANY, false, opt_threads},
  {"--context",             1, TOOL_SW,  false, opt_context},
  {"--match",               1, TOOL_ANY, false, opt_match},
  {"--mismatch",            1, TOOL_ANY, false, opt_mismatch},
  {"--gapopen",             1, TOOL_ANY, false, opt_gapopen},
  {"--gapextend",           1, TOOL_ANY, false, opt_gapextend},
  {"--file",                1, TOOL_ANY, false, opt_file},
  {"--files",               2, TOOL_ANY, false, opt_files},
  {"--wildcard",            2, TOOL_ANY, false, opt_wildcard},
};

static const option_t *option_find(const char *arg)
{
  size_t i;
  for(i = 0; i < sizeof(options) / sizeof(options[0]); i++)
    if(strcasecmp(arg, options[i].name) == 0) return &options[i];
  return NULL;
}

// The options the other options build on: --help, --case_sensitive and
// --scoring, applied before the rest (case sensitivity shapes the tables
// the scoring and substitution files fill), stepping over each option's
// values so a value that happens to look like one of these (a file called
// --help) is not taken for it. "--" ends the options. Also the first check
// that every value-taking option has its values
static void parse_early(parse_ctx_t *c, int argc, char **argv)
{
  int argi;
  for(argi = 1; argi < argc && strcmp(argv[argi], "--") != 0; argi++)
  {
    const option_t *opt = option_find(argv[argi]);
    const int values = opt != NULL ? opt->values : 0;
    if(argi + values >= argc)
      ctx_usage(c, "%s requires %s", argv[argi],
                values == 2 ? "two arguments" : "an argument");
    if(opt != NULL && opt->early) opt->apply(c, &argv[argi+1]);
    argi += values;
  }
}

// The checks across options, once every option is read
static void validate(parse_ctx_t *c, bool options_ended)
{
  cmdline_t *cmd = c->cmd;
  scoring_t *scoring = c->scoring;

  if(!score_in_range(scoring->match) || !score_in_range(scoring->mismatch) ||
     !score_in_range(scoring->gap_open) || !score_in_range(scoring->gap_extend))
  {
    ctx_usage(c, "Scores must be between %i and %i", -SCORE_ABS_MAX,
              SCORE_ABS_MAX);
  }

  if((c->match_set && !c->mismatch_set && !scoring->no_mismatches) ||
     (!c->match_set && c->mismatch_set))
  {
    ctx_usage(c, "--match --mismatch must both be set or neither set");
  }
  else if(c->substitutions_set && !c->match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  if(scoring->use_match_mismatch && scoring->match < scoring->mismatch) {
    ctx_usage(c, "Match value should not be less than mismatch penalty");
  }

  // A gap scoring above 0 would grow out of Smith-Waterman's empty cells
  if(c->cmd_type == NUCLEICS_ALIGN_SW_CMD &&
     alignment_sw_gaps_unsupported(scoring))
    ctx_usage(c, "Smith-Waterman gaps must not score above 0: --gapextend, "
                 "and --gapopen plus --gapextend, must be at most 0");

  // Cannot guarantee that we can perform a global alignment if nomismatches
  // and nogaps is true
  if(c->cmd_type == NUCLEICS_ALIGN_NW_CMD && scoring->no_mismatches &&
     (scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
      scoring->no_gaps_in_longest))
  {
    ctx_usage(c, "--nogaps.. --nomismatches cannot be used at together");
  }

  // An option after the first sequence is a mistake, not a sequence
  if(cmd->seq2 != NULL && !options_ended && cmd->seq2[0] == '-')
    ctx_usage(c, "Unknown argument '%s' (put -- before sequences that start "
                 "with -)", cmd->seq2);

  if(cmd->seq1 == NULL && cmd->file_list_length == 0)
  {
    ctx_usage(c, "No input specified");
  }

  if(cmd->zam_stle_output &&
     (cmd->print_pretty || cmd->print_scores ||
      cmd->print_colour || cmd->print_fasta))
  {
    ctx_usage(c, "Cannot use --printscore, --printfasta, --pretty or --colour "
                 "with --zam");
  }

  if(cmd->score_only &&
     (cmd->zam_stle_output || cmd->print_pretty || cmd->print_colour ||
      cmd->print_matrices))
  {
    ctx_usage(c, "--scoreonly cannot be used with --zam, --pretty, --colour or "
                 "--printmatrices");
  }

  if(cmd->linear_above > 0 && cmd->print_matrices)
    ctx_usage(c, "--linearabove cannot be used with --printmatrices");
  if(cmd->linear_fallback && cmd->print_matrices)
    ctx_usage(c, "--linearfallback cannot be used with --printmatrices");
  if(cmd->wfa_align && cmd->print_matrices)
    ctx_usage(c, "--wfa cannot be used with --printmatrices");

  if(cmd->linear_space &&
     (scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
      scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
      scoring->no_gaps_in_longest || scoring->no_mismatches ||
      cmd->print_matrices))
  {
    ctx_usage(c, "--linear cannot be used with --freestartgap, --freeendgap, "
                 "--nogaps.., --nomismatches or --printmatrices");
  }
}

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum NucleicsAlignCmdType cmd_type)
{
  cmdline_t* cmd = calloc(1, sizeof(cmdline_t));
  if(cmd == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  cmd->file_list_length = 0;
  cmd->file_list_capacity = 256;
  cmd->file_paths1 = malloc(sizeof(char*) * cmd->file_list_capacity);
  cmd->file_paths2 = malloc(sizeof(char*) * cmd->file_list_capacity);
  if(cmd->file_paths1 == NULL || cmd->file_paths2 == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  cmd->seq1 = cmd->seq2 = NULL;
  // All values initially 0

  // Store defaults
  parse_ctx_t ctx = {cmd, scoring, cmd_type,
                     {scoring->match, scoring->mismatch, scoring->gap_open,
                      scoring->gap_extend},
                     argv[0], false, false, false, false};
  parse_ctx_t *c = &ctx;

  if(argc == 1) ctx_usage(c, NULL);

  parse_early(c, argc, argv);

  // Before wildcards and substitution files are added: they keep their case
  if(cmd->case_sensitive) scoring_make_case_sensitive(scoring);

  const int tool = cmd_type == NUCLEICS_ALIGN_NW_CMD ? TOOL_NW : TOOL_SW;
  bool options_ended = false; // after "--": sequences only, any first byte
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcmp(argv[argi], "--") == 0)
    {
      options_ended = true;
      argi++;
      if(argc - argi != 2)
        ctx_usage(c, "-- must be followed by the two sequences");
      break;
    }
    else if(argv[argi][0] == '-')
    {
      const option_t *opt = option_find(argv[argi]);
      if(opt == NULL) {
        if(argi == argc-1)
          ctx_usage(c, "Unknown argument without parameter: %s", argv[argi]);
        ctx_usage(c, "Unknown argument '%s'", argv[argi]);
      }
      if(opt->tool != TOOL_ANY && opt->tool != tool)
        ctx_usage(c, "%s only valid with %s", opt->name,
                  opt->tool == TOOL_NW ? "Needleman-Wunsch" : "Smith-Waterman");
      // The first pass applied the early ones (and checked every option
      // has its values)
      if(!opt->early) opt->apply(c, &argv[argi+1]);
      argi += opt->values;
    }
    else
    {
      if(argc - argi != 2) ctx_usage(c, "Unknown options: '%s'", argv[argi]);
      break;
    }
  }

  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    cmd->seq1 = argv[argi];
    cmd->seq2 = argv[argi+1];
  }

  validate(c, options_ended);

  return cmd;
}


void cmdline_add_files(cmdline_t *cmd, char* p1, char* p2)
{
  if(cmd->file_list_length == cmd->file_list_capacity)
  {
    cmd->file_list_capacity *= 2;
    size_t mem = sizeof(char*) * cmd->file_list_capacity;
    cmd->file_paths1 = realloc(cmd->file_paths1, mem);
    cmd->file_paths2 = realloc(cmd->file_paths2, mem);

    if(cmd->file_paths1 == NULL || cmd->file_paths2 == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  cmd->file_paths1[cmd->file_list_length] = p1;
  cmd->file_paths2[cmd->file_list_length] = p2;
  cmd->file_list_length++;
}

size_t cmdline_get_num_of_file_pairs(cmdline_t *cmd)
{
  return cmd->file_list_length;
}

char* cmdline_get_file1(cmdline_t *cmd, size_t i)
{
  return cmd->file_paths1[i];
}

char* cmdline_get_file2(cmdline_t *cmd, size_t i)
{
  return cmd->file_paths2[i];
}

static seq_file_t* open_seq_file(const char *path, bool use_zlib)
{
  return (strcmp(path,"-") != 0 || use_zlib) ? seq_open(path)
                                             : seq_dopen(fileno(stdin), false);
}

// Open path1, and path2 unless it is NULL (pairs are then read from path1
// alone). Reports and returns false if either cannot be opened
static bool open_seq_files(const char *path1, const char *path2,
                           bool use_zlib, seq_file_t **sf1, seq_file_t **sf2)
{
  if((*sf1 = open_seq_file(path1, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return false;
  }

  if(path2 == NULL)
  {
    *sf2 = *sf1;
  }
  else if((*sf2 = open_seq_file(path2, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path2);
    fflush(stderr);
    seq_close(*sf1);
    return false;
  }
  return true;
}

static void close_seq_files(seq_file_t *sf1, seq_file_t *sf2)
{
  seq_close(sf1);
  if(sf2 != sf1) seq_close(sf2);
}

// How the input ended: cleanly, with a sequence left without a pair, or
// with malformed input (the reader's message in *msg)
typedef enum { INPUT_OK, INPUT_ODD, INPUT_BAD } input_end_t;

// Once the input ends: a sequence left without a pair, or malformed input,
// is an error (the pairs before it are already written), no pairs a
// warning. two_files: the pairs came from two files, one sequence from
// each. Returns false on an error
static bool report_input_end(input_end_t how, const char *msg,
                             bool two_files, size_t alignments)
{
  if(how == INPUT_BAD)
  {
    fprintf(stderr, "Alignment Error: %s\n", msg);
    fflush(stderr);
    return false;
  }
  if(how == INPUT_ODD)
  {
    if(two_files)
      fprintf(stderr, "Alignment Error: the files have different numbers "
                      "of sequences\n");
    else
      fprintf(stderr, "Alignment Error: Odd number of sequences - "
                      "I read in pairs!\n");
    fflush(stderr);
    return false;
  }

  // warn if no bases read
  if(alignments == 0)
  {
    fprintf(stderr, "Alignment Warning: empty input\n");
    fflush(stderr);
  }
  return true;
}

// The pair's second sequence: read, or the input ended (how says how)
static bool read_second(seq_file_t *sf2, read_t *r2, input_end_t *how,
                        const char **msg)
{
  int rc = seq_read(sf2, r2);
  if(rc > 0) return true;
  if(rc < 0) { *how = INPUT_BAD; *msg = seq_error(sf2); }
  else *how = INPUT_ODD;
  return false;
}

// Once sf1 is exhausted (rc1: what its last read returned): with two
// files, a sequence left in sf2 is odd too
static input_end_t input_left_over(int rc1, seq_file_t *sf1, seq_file_t *sf2,
                                   read_t *spare, const char **msg)
{
  if(rc1 < 0) { *msg = seq_error(sf1); return INPUT_BAD; }
  if(sf2 == sf1) return INPUT_OK;
  int rc = seq_read(sf2, spare);
  if(rc < 0) { *msg = seq_error(sf2); return INPUT_BAD; }
  return rc > 0 ? INPUT_ODD : INPUT_OK;
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
bool align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib)
{
  seq_file_t *sf1, *sf2;
  if(!open_seq_files(path1, path2, use_zlib, &sf1, &sf2)) return false;

  read_t read1, read2;
  seq_read_alloc(&read1);
  seq_read_alloc(&read2);

  // Loop while we can read a sequence from the first file
  size_t alignments;
  input_end_t how = INPUT_OK;
  const char *msg = NULL;
  int rc1;

  for(alignments = 0; (rc1 = seq_read(sf1, &read1)) > 0; alignments++)
  {
    if(!read_second(sf2, &read2, &how, &msg)) break;
    (align)(&read1, &read2);
  }
  if(how == INPUT_OK) how = input_left_over(rc1, sf1, sf2, &read2, &msg);

  const bool ok = report_input_end(how, msg, path2 != NULL, alignments);
  close_seq_files(sf1, sf2);

  // Free memory
  seq_read_dealloc(&read1);
  seq_read_dealloc(&read2);
  return ok;
}

unsigned int cmdline_threads(const cmdline_t *cmd)
{
  if(cmd->threads > 0) return cmd->threads;
  long cpus = sysconf(_SC_NPROCESSORS_ONLN);
  return cpus > 1 ? (unsigned int)(cpus / 2) : 1;
}

//
// align_from_file_threaded(): this thread reads pairs into a ring of jobs
// (queueing), worker threads align them (running), and whichever worker
// finishes the next pair in input order writes it (ordered output)
//

// A pair read from the file, and the output aligning it gave
typedef struct
{
  read_t r1, r2;
  bool allocated, done;
  size_t number;
  char *out, *err;
  size_t out_len, err_len;
  int status;
} cmd_job_t;

typedef int (*cmd_run_t)(void *worker, read_t *r1, read_t *r2, size_t number,
                         const alignment_options_t *opts,
                         FILE *out, FILE *err);

// Pairs being aligned, in a ring of cap slots: the i-th pair read is in slot
// i % cap. Counts are of pairs since the start of the file
typedef struct
{
  pthread_mutex_t lock;
  pthread_cond_t work, space; // a pair to take (or the end); a pair written
  cmd_job_t *jobs;
  size_t cap, read, taken, written;
  bool end; // no more pairs
  int status; // a pair's non-zero status: stop after its output
  cmd_run_t run;
  const alignment_options_t *opts; // for pairs aligned beside others
} cmd_batch_t;

typedef struct
{
  cmd_batch_t *batch;
  void *worker;
  pthread_t thread;
} cmd_worker_t;

static void cmd_batch_init(cmd_batch_t *b, size_t cap, cmd_run_t run,
                           const alignment_options_t *opts)
{
  memset(b, 0, sizeof(*b));
  pthread_mutex_init(&b->lock, NULL);
  pthread_cond_init(&b->work, NULL);
  pthread_cond_init(&b->space, NULL);
  b->cap = cap;
  b->run = run;
  b->opts = opts;
  b->jobs = calloc(cap, sizeof(cmd_job_t));
  if(b->jobs == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

static void cmd_batch_free(cmd_batch_t *b)
{
  size_t i;
  for(i = 0; i < b->cap; i++) {
    // Pairs finished after one that failed are not written
    if(b->jobs[i].done) {
      free(b->jobs[i].out);
      free(b->jobs[i].err);
    }
    if(b->jobs[i].allocated) {
      seq_read_dealloc(&b->jobs[i].r1);
      seq_read_dealloc(&b->jobs[i].r2);
    }
  }
  free(b->jobs);
  pthread_mutex_destroy(&b->lock);
  pthread_cond_destroy(&b->work);
  pthread_cond_destroy(&b->space);
}

// Wait until the first n pairs are written. Returns false if a pair failed
static bool cmd_wait_written(cmd_batch_t *b, size_t n)
{
  pthread_mutex_lock(&b->lock);
  while(b->written < n && b->status == 0)
    pthread_cond_wait(&b->space, &b->lock);
  bool ok = b->status == 0;
  pthread_mutex_unlock(&b->lock);
  return ok;
}

// Write the finished pairs next in order, up to one that failed. Called
// holding the lock
static void cmd_write_done(cmd_batch_t *b)
{
  while(b->status == 0 && b->written < b->taken)
  {
    cmd_job_t *job = &b->jobs[b->written % b->cap];
    if(!job->done) break;
    fwrite(job->out, 1, job->out_len, stdout);
    fflush(stdout);
    fwrite(job->err, 1, job->err_len, stderr);
    fflush(stderr);
    free(job->out);
    free(job->err);
    job->done = false;
    b->status = job->status;
    b->written++;
    // As one pair at a time, end here: this thread may be the only one free,
    // with the reader waiting for more input (e.g. on a pipe kept open)
    if(b->status != 0) exit(b->status);
  }
  pthread_cond_broadcast(&b->space);
}

static void* cmd_worker_thread(void *ptr)
{
  cmd_worker_t *wk = ptr;
  cmd_batch_t *b = wk->batch;
  pthread_mutex_lock(&b->lock);
  while(true)
  {
    while(b->taken == b->read && !b->end && b->status == 0)
      pthread_cond_wait(&b->work, &b->lock);
    if(b->taken == b->read || b->status != 0) break;
    cmd_job_t *job = &b->jobs[b->taken++ % b->cap];
    pthread_mutex_unlock(&b->lock);

    FILE *out = open_memstream(&job->out, &job->out_len);
    FILE *err = open_memstream(&job->err, &job->err_len);
    if(out == NULL || err == NULL) {
      fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    job->status = b->run(wk->worker, &job->r1, &job->r2, job->number,
                         b->opts, out, err);
    fclose(out);
    fclose(err);

    pthread_mutex_lock(&b->lock);
    job->done = true;
    cmd_write_done(b);
  }
  pthread_mutex_unlock(&b->lock);
  return NULL;
}

// Read the next pair into its slot (once the slot's last pair is written).
// Returns 1 with a pair, 0 at the end of the input (*how and *msg say how
// it ended; *rc1 what the first file's read returned) or once a pair has
// failed
static int cmd_read_pair(cmd_batch_t *b, seq_file_t *sf1, seq_file_t *sf2,
                         int *rc1, input_end_t *how, const char **msg)
{
  if(!cmd_wait_written(b, b->read >= b->cap ? b->read - b->cap + 1 : 0))
    return 0;
  cmd_job_t *job = &b->jobs[b->read % b->cap];
  if(!job->allocated) {
    seq_read_alloc(&job->r1);
    seq_read_alloc(&job->r2);
    job->allocated = true;
  }
  if((*rc1 = seq_read(sf1, &job->r1)) <= 0) return 0;
  if(!read_second(sf2, &job->r2, how, msg)) return 0;
  return 1;
}

// Hand the pair just read to the workers
static void cmd_queue_pair(cmd_batch_t *b)
{
  pthread_mutex_lock(&b->lock);
  b->read++;
  pthread_cond_signal(&b->work);
  pthread_mutex_unlock(&b->lock);
}

// Let the workers finish the pairs taken (writing them) and stop
static void cmd_stop_workers(cmd_batch_t *b, cmd_worker_t *workers,
                             size_t nworkers, void (*free_worker)(void*))
{
  size_t i;
  pthread_mutex_lock(&b->lock);
  b->end = true;
  pthread_cond_broadcast(&b->work);
  pthread_mutex_unlock(&b->lock);
  for(i = 0; i < nworkers; i++) {
    pthread_join(workers[i].thread, NULL);
    free_worker(workers[i].worker);
  }
}

bool align_from_file_threaded(const char *path1, const char *path2,
                              bool use_zlib, unsigned int threads,
                              size_t (*number)(read_t *r1, read_t *r2),
                              const alignment_options_t *opts,
                              int (*run)(void *worker, read_t *r1, read_t *r2,
                                         size_t number,
                                         const alignment_options_t *opts,
                                         FILE *out, FILE *err),
                              void* (*new_worker)(void),
                              void (*free_worker)(void *worker))
{
  seq_file_t *sf1, *sf2;
  if(!open_seq_files(path1, path2, use_zlib, &sf1, &sf2)) return false;

  // An alignment beside others gets one thread; one on its own gets them all
  threads = MAX2(threads, 1);
  alignment_options_t one, all;
  if(opts != NULL) all = *opts;
  else alignment_options_default(&all);
  all.max_threads = all.hb_max_threads = (int)threads;
  one = all;
  one.max_threads = one.hb_max_threads = 1;

  // Pairs read ahead: enough to keep every thread busy. Workers start as
  // pairs arrive, so no more than there are pairs
  cmd_batch_t b;
  cmd_batch_init(&b, 8 * (size_t)threads, run, &one);
  cmd_worker_t *workers = calloc(threads, sizeof(cmd_worker_t));
  if(workers == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  void *main_worker = NULL;
  size_t nworkers = 0, alignments = 0;
  input_end_t how = INPUT_OK;
  const char *msg = NULL;
  int rc1 = 0;

  while(cmd_read_pair(&b, sf1, sf2, &rc1, &how, &msg))
  {
    cmd_job_t *job = &b.jobs[b.read % b.cap];
    job->number = number != NULL ? number(&job->r1, &job->r2) : alignments;
    alignments++;

    // Big pairs run on their own with all the threads
    bool alone = (double)job->r1.seq.end * (double)job->r2.seq.end >
                 (double)all.block_min_cells;

    if(!alone && nworkers < threads)
    {
      cmd_worker_t *wk = &workers[nworkers];
      wk->batch = &b;
      wk->worker = new_worker();
      if(pthread_create(&wk->thread, NULL, cmd_worker_thread, wk) == 0)
        nworkers++;
      else
        free_worker(wk->worker);
    }

    if(!alone && nworkers > 0) { cmd_queue_pair(&b); continue; }

    // Here, once every pair before is written: no worker is aligning or
    // writing then. The slot is free again afterwards
    if(!cmd_wait_written(&b, b.read)) break;
    if(main_worker == NULL) main_worker = new_worker();
    int status = run(main_worker, &job->r1, &job->r2, job->number, &all,
                     stdout, stderr);
    if(status != 0) {
      pthread_mutex_lock(&b.lock);
      b.status = status;
      pthread_mutex_unlock(&b.lock);
      break;
    }
  }

  cmd_stop_workers(&b, workers, nworkers, free_worker);
  if(main_worker != NULL) free_worker(main_worker);

  // Reported after the output of the pairs before, as one at a time
  const int status = b.status;
  bool ok = false;
  if(status == 0) {
    if(how == INPUT_OK)
      how = input_left_over(rc1, sf1, sf2, &b.jobs[b.read % b.cap].r2, &msg);
    ok = report_input_end(how, msg, path2 != NULL, alignments);
  }

  close_seq_files(sf1, sf2);
  cmd_batch_free(&b);
  free(workers);

  // As one pair at a time: after the failing pair's output
  if(status != 0) exit(status);
  return ok;
}

// The tool being run by cmdline_run_tool(), for the callbacks below
static const cmd_tool_t *run_tool;
static const cmdline_t *run_cmd;
// One pair at a time, on this thread
static tool_worker_t *main_worker;
static alignment_options_t main_opts;

static void* tool_worker_new(void)
{
  tool_worker_t *wk = malloc(sizeof(tool_worker_t));
  if(wk == NULL) {
    fprintf(stderr, "%s:%i: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  wk->aligner = run_tool->aligner_new(run_cmd);
  wk->result = alignment_create(256);
  return wk;
}

static void tool_worker_free(void *ptr)
{
  tool_worker_t *wk = ptr;
  run_tool->aligner_free(wk->aligner);
  alignment_free(wk->result);
  free(wk);
}

static size_t tool_number(read_t *read1, read_t *read2)
{
  return run_tool->number(read1->seq.b, read2->seq.b);
}

static int tool_run(void *ptr, read_t *read1, read_t *read2, size_t number,
                    const alignment_options_t *opts, FILE *out, FILE *err)
{
  tool_worker_t *wk = ptr;
  run_tool->set_options(wk->aligner, opts);
  return run_tool->align(wk, read1->seq.b, read2->seq.b,
                         (read1->name.end == 0 ? NULL : read1->name.b),
                         (read2->name.end == 0 ? NULL : read2->name.b),
                         number, out, err);
}

static void tool_pair_from_file(read_t *read1, read_t *read2)
{
  size_t number = run_tool->number != NULL ? tool_number(read1, read2) : 0;
  int status = tool_run(main_worker, read1, read2, number, &main_opts,
                        stdout, stderr);
  if(status != 0) exit(status);
}

int cmdline_run_tool(cmdline_t *cmd, const cmd_tool_t *tool)
{
  run_tool = tool;
  run_cmd = cmd;

  // Threads within an alignment, and for pairs from files at once
  unsigned int threads = cmdline_threads(cmd);
  alignment_options_default(&main_opts);
  main_opts.max_threads = main_opts.hb_max_threads = (int)threads;

  // Align!
  main_worker = tool_worker_new();
  tool->set_options(main_worker->aligner, &main_opts);

  if(cmd->seq1 != NULL)
  {
    // Align seq1 and seq2 pair passed on the command line
    size_t number = tool->number != NULL ? tool->number(cmd->seq1, cmd->seq2)
                                         : 0;
    int status = tool->align(main_worker, cmd->seq1, cmd->seq2, NULL, NULL,
                             number, stdout, stderr);
    if(status != 0) exit(status);
  }

  // Align from files. A file that cannot be opened or ends with a sequence
  // left over is reported there and then; the files after it are still
  // aligned, and the status at the end says it failed
  bool input_ok = true;
  size_t i, num_of_file_pairs = cmdline_get_num_of_file_pairs(cmd);
  for(i = 0; i < num_of_file_pairs; i++)
  {
    const char *file1 = cmdline_get_file1(cmd, i);
    const char *file2 = cmdline_get_file2(cmd, i);
    if(file1 != NULL && *file1 == '\0' && file2 == NULL) {
      file1 = "-";
      if(tool->on_stdin != NULL) tool->on_stdin();
    }
    // Interactive input, and printed matrices, one pair at a time
    bool ok;
    if(threads > 1 && !cmd->interactive && !cmd->print_matrices)
      ok = align_from_file_threaded(file1, file2, !cmd->interactive, threads,
                                    tool->number != NULL ? tool_number : NULL,
                                    &main_opts, tool_run, tool_worker_new,
                                    tool_worker_free);
    else
      ok = align_from_file(file1, file2, &tool_pair_from_file,
                           !cmd->interactive);
    input_ok = input_ok && ok;
  }

  // Free memory for storing alignment results
  tool_worker_free(main_worker);

  cmdline_free(cmd);

  return input_ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
