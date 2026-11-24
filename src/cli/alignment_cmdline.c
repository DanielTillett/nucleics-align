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

#include "seq_file/seq_file.h"

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

static void print_usage(enum NucleicsAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
  __attribute__((format(printf, 4, 5)))
  __attribute__((noreturn));

static void print_usage(enum NucleicsAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
  }

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmdstr);

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

#define usage(fmt,...) print_usage(cmd_type,defaults,argv[0],fmt, ##__VA_ARGS__)

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
  score_t defaults[4] = {scoring->match, scoring->mismatch,
                         scoring->gap_open, scoring->gap_extend};

  if(argc == 1) usage(NULL);

  // First run through arguments to set up case_sensitive and scoring system

  // case sensitive needs to be dealt with first
  // (is is used to construct hash table for swap_scores)
  char scoring_set = 0, substitutions_set = 0, match_set = 0, mismatch_set = 0;

  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "--help") == 0 ||
       strcasecmp(argv[argi], "-help") == 0 ||
       strcasecmp(argv[argi], "-h") == 0)
    {
      usage(NULL);
    }
    else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
    {
      cmd->case_sensitive = 1;
    }
    else if(strcasecmp(argv[argi], "--scoring") == 0)
    {
      if(scoring_set)
      {
        usage("More than one scoring system specified - not permitted");
      }

      if(argi + 1 >= argc)
      {
        usage("--scoring requires an argument");
      }

      if(!scoring_system_by_name(scoring, argv[argi+1]))
      {
        usage("Unknown --scoring choice, not one of %s",
              scoring_system_names());
      }

      scoring_set = 1;
      argi++; // took an argument
    }
  }

  // Before wildcards and substitution files are added: they keep their case
  if(cmd->case_sensitive) scoring_make_case_sensitive(scoring);

  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
      // strcasecmp does case insensitive comparison
      if(strcasecmp(argv[argi], "--freestartgap") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--freestartgap only valid with Needleman-Wunsch");
        scoring->no_start_gap_penalty = true;
      }
      else if(strcasecmp(argv[argi], "--freeendgap") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--freeendgap only valid with Needleman-Wunsch");
        scoring->no_end_gap_penalty = true;
      }
      else if(strcasecmp(argv[argi], "--nogaps") == 0)
      {
        scoring->no_gaps_in_a = true;
        scoring->no_gaps_in_b = true;
      }
      else if(strcasecmp(argv[argi], "--nogapsin1") == 0)
      {
        scoring->no_gaps_in_a = true;
      }
      else if(strcasecmp(argv[argi], "--nogapsin2") == 0)
      {
        scoring->no_gaps_in_b = true;
      }
      else if(strcasecmp(argv[argi], "--nogapsinlongest") == 0)
      {
        scoring->no_gaps_in_longest = true;
      }
      else if(strcasecmp(argv[argi], "--nomismatches") == 0)
      {
        scoring->no_mismatches = true;
      }
      else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
      }
      else if(strcasecmp(argv[argi], "--printseq") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_SW_CMD)
          usage("--printseq only valid with Smith-Waterman");
        cmd->print_seq = true;
      }
      else if(strcasecmp(argv[argi], "--printmatrices") == 0)
      {
        cmd->print_matrices = true;
      }
      else if(strcasecmp(argv[argi], "--printscores") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--printscores only valid with Needleman-Wunsch");
        cmd->print_scores = true;
      }
      else if(strcasecmp(argv[argi], "--printfasta") == 0)
      {
        cmd->print_fasta = true;
      }
      else if(strcasecmp(argv[argi], "--pretty") == 0)
      {
        cmd->print_pretty = true;
      }
      else if(strcasecmp(argv[argi], "--colour") == 0)
      {
        cmd->print_colour = true;
      }
      else if(strcasecmp(argv[argi], "--zam") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--zam only valid with Needleman-Wunsch");
        cmd->zam_stle_output = true;
      }
      else if(strcasecmp(argv[argi], "--scoreonly") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--scoreonly only valid with Needleman-Wunsch");
        cmd->score_only = true;
      }
      else if(strcasecmp(argv[argi], "--linear") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--linear only valid with Needleman-Wunsch");
        cmd->linear_space = true;
      }
      else if(strcasecmp(argv[argi], "--linearfallback") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--linearfallback only valid with Needleman-Wunsch");
        cmd->linear_fallback = true;
      }
      else if(strcasecmp(argv[argi], "--wfa") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--wfa only valid with Needleman-Wunsch");
        cmd->wfa_align = true;
      }
      else if(strcasecmp(argv[argi], "--fastlocal") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_SW_CMD)
          usage("--fastlocal only valid with Smith-Waterman");
        cmd->fast_local = true;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        // Similar to --file argument below
        cmdline_add_files(cmd, "", NULL);
        cmd->interactive = true;
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        usage("Unknown argument without parameter: %s", argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        // This handled above
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        gzFile sub_matrix_file = gzopen(argv[argi+1], "r");
        if(sub_matrix_file == NULL) usage("Couldn't read: %s", argv[argi+1]);

        align_scoring_load_matrix(sub_matrix_file, argv[argi+1],
                                  scoring, cmd->case_sensitive);

        gzclose(sub_matrix_file);
        substitutions_set = true;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_pairs") == 0)
      {
        gzFile sub_pairs_file = gzopen(argv[argi+1], "r");
        if(sub_pairs_file == NULL) usage("Couldn't read: %s", argv[argi+1]);

        align_scoring_load_pairwise(sub_pairs_file, argv[argi+1],
                                    scoring, cmd->case_sensitive);

        gzclose(sub_pairs_file);
        substitutions_set = true;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--minscore") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_SW_CMD)
          usage("--minscore only valid with Smith-Waterman");

        if(!parse_entire_int(argv[argi+1], &cmd->min_score))
          usage("Invalid --minscore <score> argument (must be a +ve int)");

        cmd->min_score_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--maxhits") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_SW_CMD)
          usage("--maxhits only valid with Smith-Waterman");

        if(!parse_entire_uint(argv[argi+1], &cmd->max_hits_per_alignment))
          usage("Invalid --maxhits <hits> argument (must be a +ve int)");

        cmd->max_hits_per_alignment_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--linearabove") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--linearabove only valid with Needleman-Wunsch");
        if(!parse_entire_uint(argv[argi+1], &cmd->linear_above) ||
           cmd->linear_above == 0)
          usage("Invalid --linearabove <n> argument (must be at least 1)");
        argi++;
      }
      else if(strcasecmp(argv[argi], "--wfabudget") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_NW_CMD)
          usage("--wfabudget only valid with Needleman-Wunsch");
        const char *arg = argv[argi+1];
        if(arg[0] == '\0' || !isdigit((unsigned char)arg[0]))
          usage("Invalid --wfabudget <n> argument (must be a positive integer)");
        char *end = NULL;
        errno = 0;
        unsigned long long v = strtoull(arg, &end, 10);
        if(errno == ERANGE || end == NULL || *end != '\0' || v == 0 ||
           v > (unsigned long long)SIZE_MAX)
          usage("Invalid --wfabudget <n> argument (must be a positive integer)");
        cmd->wfa_budget = (size_t)v;
        argi++;
      }
      else if(strcasecmp(argv[argi], "--threads") == 0)
      {
        if(!parse_entire_uint(argv[argi+1], &cmd->threads) ||
           cmd->threads == 0 || cmd->threads > 4096)
          usage("Invalid --threads <n> argument (must be 1 to 4096)");
        argi++;
      }
      else if(strcasecmp(argv[argi], "--context") == 0)
      {
        if(cmd_type != NUCLEICS_ALIGN_SW_CMD)
          usage("--context only valid with Smith-Waterman");

        if(!parse_entire_uint(argv[argi+1], &cmd->print_context))
          usage("Invalid --context <c> argument (must be >= 0)");

        argi++;
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->match))
        {
          usage("Invalid --match argument ('%s') must be an int", argv[argi+1]);
        }

        match_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->mismatch))
        {
          usage("Invalid --mismatch argument ('%s') must be an int", argv[argi+1]);
        }

        mismatch_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_open))
        {
          usage("Invalid --gapopen argument ('%s') must be an int", argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_extend))
        {
          usage("Invalid --gapextend argument ('%s') must be an int",
                argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        cmdline_add_files(cmd, argv[argi+1], NULL);
        argi++; // took an argument
      }
      // Remaining options take two arguments but check themselves
      else if(strcasecmp(argv[argi], "--files") == 0)
      {
        if(argi >= argc-2)
        {
          usage("--files option takes 2 arguments");
        }
        else if(strcmp(argv[argi+1], "-") == 0 && strcmp(argv[argi+2], "-") == 0)
        {
          // Read both from stdin
          cmdline_add_files(cmd, argv[argi+1], NULL);
        }
        else
        {
          cmdline_add_files(cmd, argv[argi+1], argv[argi+2]);
        }

        argi += 2; // took two arguments
      }
      else if(strcasecmp(argv[argi], "--wildcard") == 0)
      {
        int wildscore = 0;

        if(argi == argc-2 || strlen(argv[argi+1]) != 1 ||
           !parse_entire_int(argv[argi+2], &wildscore))
        {
          usage("--wildcard <w> <s> takes a single character and a number");
        }
        if(!score_in_range(wildscore))
          usage("Scores must be between %i and %i", -SCORE_ABS_MAX, SCORE_ABS_MAX);

        scoring_add_wildcard(scoring, argv[argi+1][0], wildscore);

        argi += 2; // took two arguments
      }
      else usage("Unknown argument '%s'", argv[argi]);
    }
    else
    {
      if(argc - argi != 2) usage("Unknown options: '%s'", argv[argi]);
      break;
    }
  }

  if(!score_in_range(scoring->match) || !score_in_range(scoring->mismatch) ||
     !score_in_range(scoring->gap_open) || !score_in_range(scoring->gap_extend))
  {
    usage("Scores must be between %i and %i", -SCORE_ABS_MAX, SCORE_ABS_MAX);
  }

  if((match_set && !mismatch_set && !scoring->no_mismatches) ||
     (!match_set && mismatch_set))
  {
    usage("--match --mismatch must both be set or neither set");
  }
  else if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  if(scoring->use_match_mismatch && scoring->match < scoring->mismatch) {
    usage("Match value should not be less than mismatch penalty");
  }

  // Cannot guarantee that we can perform a global alignment if nomismatches
  // and nogaps is true
  if(cmd_type == NUCLEICS_ALIGN_NW_CMD && scoring->no_mismatches &&
     (scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
      scoring->no_gaps_in_longest))
  {
    usage("--nogaps.. --nomismatches cannot be used at together");
  }

  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    cmd->seq1 = argv[argi];
    cmd->seq2 = argv[argi+1];
  }

  if(cmd->seq1 == NULL && cmd->file_list_length == 0)
  {
    usage("No input specified");
  }

  if(cmd->zam_stle_output &&
     (cmd->print_pretty || cmd->print_scores ||
      cmd->print_colour || cmd->print_fasta))
  {
    usage("Cannot use --printscore, --printfasta, --pretty or --colour with "
          "--zam");
  }

  if(cmd->score_only &&
     (cmd->zam_stle_output || cmd->print_pretty || cmd->print_colour ||
      cmd->print_matrices))
  {
    usage("--scoreonly cannot be used with --zam, --pretty, --colour or "
          "--printmatrices");
  }

  if(cmd->linear_above > 0 && cmd->print_matrices)
    usage("--linearabove cannot be used with --printmatrices");
  if(cmd->wfa_align && cmd->print_matrices)
    usage("--wfa cannot be used with --printmatrices");

  if(cmd->linear_space &&
     (scoring->no_start_gap_penalty || scoring->no_end_gap_penalty ||
      scoring->no_gaps_in_a || scoring->no_gaps_in_b ||
      scoring->no_gaps_in_longest || scoring->no_mismatches ||
      cmd->print_matrices))
  {
    usage("--linear cannot be used with --freestartgap, --freeendgap, "
          "--nogaps.., --nomismatches or --printmatrices");
  }

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
                                             : seq_dopen(fileno(stdin), false, false, 0);
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

// Warnings once the input ends: a sequence left without a pair, no pairs
static void report_input_end(bool odd, size_t alignments)
{
  if(odd)
  {
    fprintf(stderr, "Alignment Error: Odd number of sequences - "
                    "I read in pairs!\n");
    fflush(stderr);
  }

  // warn if no bases read
  if(alignments == 0)
  {
    fprintf(stderr, "Alignment Warning: empty input\n");
    fflush(stderr);
  }
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib)
{
  seq_file_t *sf1, *sf2;
  if(!open_seq_files(path1, path2, use_zlib, &sf1, &sf2)) return;

  read_t read1, read2;
  seq_read_alloc(&read1);
  seq_read_alloc(&read2);

  // Loop while we can read a sequence from the first file
  size_t alignments;
  bool odd = false;

  for(alignments = 0; seq_read(sf1, &read1) > 0; alignments++)
  {
    if(seq_read(sf2, &read2) <= 0) { odd = true; break; }
    (align)(&read1, &read2);
  }

  report_input_end(odd, alignments);
  close_seq_files(sf1, sf2);

  // Free memory
  seq_read_dealloc(&read1);
  seq_read_dealloc(&read2);
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
// Returns 1 with a pair, 0 at the end of the input (*odd: a sequence left
// over) or once a pair has failed
static int cmd_read_pair(cmd_batch_t *b, seq_file_t *sf1, seq_file_t *sf2,
                         bool *odd)
{
  if(!cmd_wait_written(b, b->read >= b->cap ? b->read - b->cap + 1 : 0))
    return 0;
  cmd_job_t *job = &b->jobs[b->read % b->cap];
  if(!job->allocated) {
    seq_read_alloc(&job->r1);
    seq_read_alloc(&job->r2);
    job->allocated = true;
  }
  if(seq_read(sf1, &job->r1) <= 0) return 0;
  if(seq_read(sf2, &job->r2) <= 0) { *odd = true; return 0; }
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

void align_from_file_threaded(const char *path1, const char *path2,
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
  if(!open_seq_files(path1, path2, use_zlib, &sf1, &sf2)) return;

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
  bool odd = false;

  while(cmd_read_pair(&b, sf1, sf2, &odd))
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
  if(status == 0) report_input_end(odd, alignments);

  close_seq_files(sf1, sf2);
  cmd_batch_free(&b);
  free(workers);

  // As one pair at a time: after the failing pair's output
  if(status != 0) exit(status);
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

  // Align from files
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
    if(threads > 1 && !cmd->interactive && !cmd->print_matrices)
      align_from_file_threaded(file1, file2, !cmd->interactive, threads,
                               tool->number != NULL ? tool_number : NULL,
                               &main_opts, tool_run, tool_worker_new,
                               tool_worker_free);
    else
      align_from_file(file1, file2, &tool_pair_from_file, !cmd->interactive);
  }

  // Free memory for storing alignment results
  tool_worker_free(main_worker);

  cmdline_free(cmd);

  return EXIT_SUCCESS;
}
