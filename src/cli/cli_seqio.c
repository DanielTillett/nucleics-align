/*
 src/cli/cli_seqio.c
 Self-contained FASTA / FASTQ / plain / gzip sequence reader for CLI tools.
 */

#ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 200809L
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>

#include "cli_seqio.h"

enum cli_seq_fmt {
  CLI_FMT_UNKNOWN = 0,
  CLI_FMT_FASTA   = 1,
  CLI_FMT_FASTQ   = 2,
  CLI_FMT_PLAIN   = 3
};

struct cli_seq_file {
  gzFile gz;
  FILE  *fp;
  enum cli_seq_fmt fmt;
  cli_str_t line;
  cli_str_t pending;
  bool   has_pending;
  bool   eof;
  bool   borrowed; // fp is stdin: not ours to close
  char   err[160]; // why seq_read() returned -1 ("" otherwise)
};

static void str_init(cli_str_t *s, size_t cap)
{
  s->size = cap > 0 ? cap : 128;
  s->b = malloc(s->size);
  if(s->b == NULL) {
    fprintf(stderr, "%s:%d: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  s->b[0] = '\0';
  s->end = 0;
}

static void str_free(cli_str_t *s)
{
  free(s->b);
  s->b = NULL;
  s->end = s->size = 0;
}

static void str_ensure(cli_str_t *s, size_t need)
{
  if(need >= SIZE_MAX - 1) {
    fprintf(stderr, "%s:%d: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if(need + 1 <= s->size) return;
  size_t nsz = s->size ? s->size : 128;
  while(nsz < need + 1) {
    if(nsz > SIZE_MAX / 2) { nsz = need + 1; break; }
    nsz *= 2;
  }
  char *nb = realloc(s->b, nsz);
  if(nb == NULL) {
    fprintf(stderr, "%s:%d: Out of memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  s->b = nb;
  s->size = nsz;
}

static void str_reset(cli_str_t *s)
{
  s->end = 0;
  if(s->b != NULL) s->b[0] = '\0';
}

static void str_append(cli_str_t *s, const char *data, size_t len)
{
  str_ensure(s, s->end + len);
  if(len > 0) memcpy(s->b + s->end, data, len);
  s->end += len;
  s->b[s->end] = '\0';
}

static void str_set(cli_str_t *s, const char *data, size_t len)
{
  str_reset(s);
  str_append(s, data, len);
}

static void str_chomp(cli_str_t *s)
{
  while(s->end > 0 && (s->b[s->end - 1] == '\n' || s->b[s->end - 1] == '\r'))
    s->b[--s->end] = '\0';
}

static bool str_all_ws(const cli_str_t *s)
{
  for(size_t i = 0; i < s->end; i++)
    if(!isspace((unsigned char)s->b[i])) return false;
  return true;
}

/* Record why the read failed; seq_read() then returns -1 */
static void seq_fail(seq_file_t *sf, const char *fmt, const char *name)
{
  if(sf->err[0] != '\0') return; // the first error stands
  snprintf(sf->err, sizeof(sf->err), fmt, name);
  sf->eof = true;
}

static bool at_eof(const seq_file_t *sf)
{
  return sf->gz != NULL ? gzeof(sf->gz) != 0 : feof(sf->fp) != 0;
}

/* Read one physical line (CRLF trimmed) into sf->line. Returns false on
   EOF, or on an error (sf->err then says which) */
static bool read_raw_line(seq_file_t *sf)
{
  if(sf->has_pending) {
    str_set(&sf->line, sf->pending.b, sf->pending.end);
    sf->has_pending = false;
    return true;
  }
  if(sf->eof) return false;

  str_reset(&sf->line);
  for(;;) {
    str_ensure(&sf->line, sf->line.end + 256);
    size_t space = sf->line.size - sf->line.end;
    int avail = space > (size_t)INT_MAX ? INT_MAX : (int)space;
    char *res = sf->gz != NULL
                ? gzgets(sf->gz, sf->line.b + sf->line.end, avail)
                : fgets(sf->line.b + sf->line.end, avail, sf->fp);
    if(res == NULL) {
      if(sf->gz != NULL) {
        int errnum = 0;
        const char *msg = gzerror(sf->gz, &errnum);
        if(errnum != Z_OK && errnum != Z_STREAM_END) {
          seq_fail(sf, "input read error: %s", msg ? msg : "zlib error");
          return false;
        }
      } else if(sf->fp != NULL && ferror(sf->fp)) {
        seq_fail(sf, "input read error%s", "");
        return false;
      }
      sf->eof = true;
      break;
    }
    // A line read in full ends with '\n', or fills the space given; one
    // that does neither stopped early at the end of the input, or at a
    // NUL byte, which the sequences (C strings) could not carry
    size_t got = strlen(sf->line.b + sf->line.end);
    sf->line.end += got;
    const bool newline = got > 0 && sf->line.b[sf->line.end - 1] == '\n';
    if(!newline && got + 1 < (size_t)avail && !at_eof(sf)) {
      seq_fail(sf, "input contains a NUL byte%s", "");
      return false;
    }
    if(newline || at_eof(sf)) break;
  }
  if(sf->line.end == 0 && sf->eof) return false;
  str_chomp(&sf->line);
  return true;
}

static void push_pending(seq_file_t *sf)
{
  str_set(&sf->pending, sf->line.b, sf->line.end);
  sf->has_pending = true;
}

static seq_file_t *sf_alloc(void)
{
  seq_file_t *sf = calloc(1, sizeof(*sf));
  if(sf == NULL) return NULL;
  str_init(&sf->line, 256);
  str_init(&sf->pending, 256);
  return sf;
}

seq_file_t *seq_open(const char *path)
{
  if(path == NULL) return NULL;
  gzFile gz = NULL;
  if(strcmp(path, "-") == 0) {
    int fd = dup(fileno(stdin));
    if(fd < 0) return NULL;
    gz = gzdopen(fd, "r");
    if(gz == NULL) { close(fd); return NULL; }
  } else {
    gz = gzopen(path, "r");
    if(gz == NULL) return NULL;
  }
  seq_file_t *sf = sf_alloc();
  if(sf == NULL) { gzclose(gz); return NULL; }
  sf->gz = gz;
  return sf;
}

seq_file_t *seq_dopen(int fd, bool is_gzip)
{
  seq_file_t *sf = sf_alloc();
  if(sf == NULL) return NULL;
  if(!is_gzip && fd == STDIN_FILENO) {
    // Read through stdin itself, not a dup() with its own buffer: the
    // interactive tools read their prompt answers with getc(stdin), which
    // must see the bytes after the sequences, not those a second buffer
    // has read ahead
    sf->fp = stdin;
    sf->borrowed = true;
    return sf;
  }
  int dfd = dup(fd);
  if(dfd < 0) { seq_close(sf); return NULL; }
  if(is_gzip) {
    sf->gz = gzdopen(dfd, "r");
    if(sf->gz == NULL) { close(dfd); seq_close(sf); return NULL; }
  } else {
    sf->fp = fdopen(dfd, "r");
    if(sf->fp == NULL) { close(dfd); seq_close(sf); return NULL; }
  }
  return sf;
}

void seq_close(seq_file_t *sf)
{
  if(sf == NULL) return;
  if(sf->gz != NULL) gzclose(sf->gz);
  if(sf->fp != NULL && !sf->borrowed) fclose(sf->fp);
  str_free(&sf->line);
  str_free(&sf->pending);
  free(sf);
}

void seq_read_alloc(read_t *r)
{
  str_init(&r->name, 128);
  str_init(&r->seq, 256);
}

void seq_read_dealloc(read_t *r)
{
  str_free(&r->name);
  str_free(&r->seq);
}

static void append_seq_nonws(cli_str_t *dst, const char *src, size_t len)
{
  str_ensure(dst, dst->end + len);
  for(size_t i = 0; i < len; i++) {
    if(!isspace((unsigned char)src[i]))
      dst->b[dst->end++] = src[i];
  }
  dst->b[dst->end] = '\0';
}

const char *seq_error(const seq_file_t *sf)
{
  return sf != NULL && sf->err[0] != '\0' ? sf->err : NULL;
}

/* 0 at the end of the input, -1 if an error was recorded on the way */
static int seq_read_end(const seq_file_t *sf)
{
  return sf->err[0] != '\0' ? -1 : 0;
}

int seq_read(seq_file_t *sf, read_t *r)
{
  if(sf == NULL || r == NULL) return 0;
  str_reset(&r->name);
  str_reset(&r->seq);

  /* Detect format on first non-blank line */
  if(sf->fmt == CLI_FMT_UNKNOWN) {
    for(;;) {
      if(!read_raw_line(sf)) return seq_read_end(sf);
      if(sf->line.end == 0 || str_all_ws(&sf->line)) continue;
      if(sf->line.b[0] == '>')      sf->fmt = CLI_FMT_FASTA;
      else if(sf->line.b[0] == '@') sf->fmt = CLI_FMT_FASTQ;
      else                          sf->fmt = CLI_FMT_PLAIN;
      push_pending(sf);
      break;
    }
  }

  if(sf->fmt == CLI_FMT_PLAIN) {
    /* One sequence per line, whitespace around it ignored, blank lines
       skipped */
    while(read_raw_line(sf)) {
      size_t s = 0, e = sf->line.end;
      while(s < e && isspace((unsigned char)sf->line.b[s])) s++;
      while(e > s && isspace((unsigned char)sf->line.b[e - 1])) e--;
      if(e == s) continue;
      str_set(&r->seq, sf->line.b + s, e - s);
      return 1;
    }
    return seq_read_end(sf);
  }

  if(sf->fmt == CLI_FMT_FASTA) {
    /* Find next '>' header */
    bool got_hdr = false;
    while(read_raw_line(sf)) {
      if(sf->line.end > 0 && sf->line.b[0] == '>') {
        str_set(&r->name, sf->line.b + 1, sf->line.end - 1);
        got_hdr = true;
        break;
      }
    }
    if(!got_hdr) return seq_read_end(sf);
    while(read_raw_line(sf)) {
      if(sf->line.end > 0 && sf->line.b[0] == '>') {
        push_pending(sf);
        break;
      }
      if(sf->line.end > 0)
        append_seq_nonws(&r->seq, sf->line.b, sf->line.end);
    }
    return sf->err[0] != '\0' ? -1 : 1;
  }

  /* CLI_FMT_FASTQ: '@name', sequence lines, '+', then quality lines with
     exactly as many characters as the sequence (they may start with '@') */
  bool got_hdr = false;
  while(read_raw_line(sf)) {
    if(sf->line.end > 0 && sf->line.b[0] == '@') {
      str_set(&r->name, sf->line.b + 1, sf->line.end - 1);
      got_hdr = true;
      break;
    }
  }
  if(!got_hdr) return seq_read_end(sf);

  /* Sequence lines until '+' */
  bool saw_plus = false;
  while(read_raw_line(sf)) {
    if(sf->line.end > 0 && sf->line.b[0] == '+') {
      saw_plus = true;
      break;
    }
    if(sf->line.end > 0 && sf->line.b[0] == '@') break; // the next record
    if(sf->line.end > 0)
      append_seq_nonws(&r->seq, sf->line.b, sf->line.end);
  }
  if(!saw_plus) {
    seq_fail(sf, "FASTQ record '%s': no '+' line", r->name.b);
    return -1;
  }

  /* Quality lines until they cover the sequence */
  size_t qlen = 0;
  while(qlen < r->seq.end) {
    if(!read_raw_line(sf)) break;
    for(size_t i = 0; i < sf->line.end; i++)
      if(!isspace((unsigned char)sf->line.b[i])) qlen++;
  }
  if(qlen != r->seq.end) {
    seq_fail(sf, "FASTQ record '%s': quality length differs from sequence "
                 "length", r->name.b);
    return -1;
  }
  return 1;
}
