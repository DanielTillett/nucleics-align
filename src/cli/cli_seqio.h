/*
 src/cli/cli_seqio.h
 Self-contained FASTA / FASTQ / plain / gzip sequence reader for CLI tools.
 Replaces the external libs/seq_file dependency.
 */

#ifndef CLI_SEQIO_HEADER_SEEN
#define CLI_SEQIO_HEADER_SEEN

#include <stddef.h>
#include <stdbool.h>

typedef struct {
  char *b;
  size_t end;
  size_t size;
} cli_str_t;

typedef struct {
  cli_str_t name;
  cli_str_t seq;
} read_t;

typedef struct cli_seq_file seq_file_t;

/* Open path for reading ("-": standard input), gzip or not */
seq_file_t *seq_open(const char *path);
/* Read from an open descriptor. Standard input, not gzip, is read through
   stdin itself, so other reads of stdin (prompt answers) stay in step */
seq_file_t *seq_dopen(int fd, bool is_gzip);
void        seq_close(seq_file_t *sf);

void seq_read_alloc(read_t *r);
void seq_read_dealloc(read_t *r);
/* The next sequence (FASTA, FASTQ or one per line, told from the first
   line). Returns 1 with the record in r, 0 at the end of the input, -1 on
   malformed input (see seq_error(); nothing more is read) */
int  seq_read(seq_file_t *sf, read_t *r);
/* The error message after seq_read() returned -1, else NULL */
const char *seq_error(const seq_file_t *sf);

#endif /* CLI_SEQIO_HEADER_SEEN */
