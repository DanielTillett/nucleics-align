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

seq_file_t *seq_open(const char *path);
seq_file_t *seq_dopen(int fd, bool is_gzip, bool unused1, size_t unused2);
void        seq_close(seq_file_t *sf);

void seq_read_alloc(read_t *r);
void seq_read_dealloc(read_t *r);
int  seq_read(seq_file_t *sf, read_t *r);

#endif /* CLI_SEQIO_HEADER_SEEN */
