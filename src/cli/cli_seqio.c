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

/* Read one physical line (CRLF trimmed) into sf->line. Returns false on EOF. */
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
          fprintf(stderr, "Alignment Error: input read error: %s\n", msg ? msg : "zlib error");
          exit(EXIT_FAILURE);
        }
      } else if(sf->fp != NULL && ferror(sf->fp)) {
        fprintf(stderr, "Alignment Error: input read error\n");
        exit(EXIT_FAILURE);
      }
      sf->eof = true;
      break;
    }
    size_t got = strlen(sf->line.b + sf->line.end);
    sf->line.end += got;
    if(got == 0 || sf->line.b[sf->line.end - 1] == '\n') break;
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

seq_file_t *seq_dopen(int fd, bool is_gzip, bool unused1, size_t unused2)
{
  (void)unused1; (void)unused2;
  seq_file_t *sf = sf_alloc();
  if(sf == NULL) return NULL;
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
  if(sf->fp != NULL) fclose(sf->fp);
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

int seq_read(seq_file_t *sf, read_t *r)
{
  if(sf == NULL || r == NULL) return 0;
  str_reset(&r->name);
  str_reset(&r->seq);

  /* Detect format on first non-blank line */
  if(sf->fmt == CLI_FMT_UNKNOWN) {
    for(;;) {
      if(!read_raw_line(sf)) return 0;
      if(sf->line.end == 0 || str_all_ws(&sf->line)) continue;
      if(sf->line.b[0] == '>')      sf->fmt = CLI_FMT_FASTA;
      else if(sf->line.b[0] == '@') sf->fmt = CLI_FMT_FASTQ;
      else                          sf->fmt = CLI_FMT_PLAIN;
      push_pending(sf);
      break;
    }
  }

  if(sf->fmt == CLI_FMT_PLAIN) {
    while(read_raw_line(sf)) {
      if(sf->line.end == 0 || isspace((unsigned char)sf->line.b[0])) continue;
      str_set(&r->seq, sf->line.b, sf->line.end);
      return 1;
    }
    return 0;
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
    if(!got_hdr) return 0;
    while(read_raw_line(sf)) {
      if(sf->line.end > 0 && sf->line.b[0] == '>') {
        push_pending(sf);
        break;
      }
      if(sf->line.end > 0)
        append_seq_nonws(&r->seq, sf->line.b, sf->line.end);
    }
    return 1;
  }

  /* CLI_FMT_FASTQ */
  bool got_hdr = false;
  while(read_raw_line(sf)) {
    if(sf->line.end > 0 && sf->line.b[0] == '@') {
      str_set(&r->name, sf->line.b + 1, sf->line.end - 1);
      got_hdr = true;
      break;
    }
  }
  if(!got_hdr) return 0;

  /* Sequence lines until '+' */
  bool saw_plus = false;
  while(read_raw_line(sf)) {
    if(sf->line.end > 0 && sf->line.b[0] == '+') {
      saw_plus = true;
      break;
    }
    if(sf->line.end > 0)
      append_seq_nonws(&r->seq, sf->line.b, sf->line.end);
  }
  if(!saw_plus) return 1;

  /* Quality lines (may wrap over multiple lines) until qual >= seq and next '@' */
  size_t qlen = 0;
  while(read_raw_line(sf)) {
    if(qlen >= r->seq.end && sf->line.end > 0 && sf->line.b[0] == '@') {
      push_pending(sf);
      break;
    }
    for(size_t i = 0; i < sf->line.end; i++)
      if(!isspace((unsigned char)sf->line.b[i])) qlen++;
  }
  return 1;
}
