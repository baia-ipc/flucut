/* (c) 2017 Giorgio Gonnella */

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include "alloc.h"
#include "strutils.h"
#include "file.h"

bool xgetline(char **line, size_t *line_size, FILE *file, Logger *logger)
{
  ssize_t getline_retval;
  getline_retval = getline(line, line_size, file);
  if (getline_retval == -1)
  {
    LCHECK(feof(file), logger, "Reading line from file failed, %s\n",
        strerror(errno));
    return false;
  }
  return true;
}

bool xgetdelim(char **line, size_t *line_size, int delim, FILE *file,
    Logger *logger)
{
  ssize_t getdelim_retval;
  getdelim_retval = getdelim(line, line_size, delim, file);
  if (getdelim_retval == -1)
  {
    LCHECK(feof(file), logger, "Reading from file until delimiter failed, %s\n",
      strerror(errno));
    return false;
  }
  return true;
}

char *compute_nth_filename(const char *basename, const char *sfx, size_t n)
{
  char *fname;
  int reqsize;
  reqsize = (sfx == NULL) ?
    snprintf(NULL, 0, "%s.%zu", basename, n) :
    snprintf(NULL, 0, "%s%s.%zu", basename, sfx, n);
  reqsize++;
  XMALLOC(fname, reqsize);
  if (sfx == NULL)
    snprintf(fname, reqsize, "%s.%zu", basename, n);
  else
    snprintf(fname, reqsize, "%s%s.%zu", basename, sfx, n);
  return fname;
}

#define MODEDESC(MODE) \
  (((MODE)[0] == 'r') ? "for reading" : "for writing")

FILE* _xfopen(const char *fname, const char *sfx, const char *mode,
            const char *description, Logger *logger,
            const char *_file, int _line) {
  FILE* file = NULL;
  char *fname_and_sfx;
  assert(fname != NULL && mode != NULL);
  fname_and_sfx = (sfx != NULL) ? strcat_alloc(fname, sfx) : (char*)fname;
  file = fopen(fname_and_sfx, mode);
  _LCHECK(file != NULL, logger, _file, _line,
      "Opening file '%s' failed, %s\n",
            fname_and_sfx, strerror(errno));
  _lprintf(logger, LOG_DEBUG, _file, _line,
      "Opened %s %s %s (%p)\n", description ?
        description : "file", fname_and_sfx, MODEDESC(mode), file);
  if (sfx != NULL)
    free(fname_and_sfx);
  return file;
}

void* _xfmap(const char *fname, const char *sfx, bool writable,
    size_t *fsize, int *fd, Logger *logger, const char *_file, int _line)
{
  char *fname_and_sfx;
  struct stat sb;
  void *map;
  assert(fname != NULL);
  fname_and_sfx = (sfx != NULL) ? strcat_alloc(fname, sfx) : (char*)fname;
  *fd = open(fname_and_sfx, writable ? O_RDWR : O_RDONLY);
  _LCHECK(*fd != -1, logger, _file, _line, "Opening file '%s' failed, %s\n",
          fname_and_sfx, strerror(errno));
  _LCHECK(fstat(*fd, &sb) == 0, logger, _file, _line,
          "fstat failed on file %s, %s\n", fname_and_sfx, strerror(errno));
  /* mmap() fails with an empty file (EINVAL), not clear as message */
  _LCHECK(sb.st_size > 0, logger, _file, _line,
          "Mapping file '%s' failed, the file is empty\n", fname_and_sfx);
  *fsize = sb.st_size;
  map = mmap(0, sb.st_size, PROT_READ | (writable ? PROT_WRITE : 0),
             MAP_SHARED, *fd, 0);
  _LCHECK(map != MAP_FAILED, logger, _file, _line,
          "Mapping file '%s' failed, %s\n", fname_and_sfx, strerror(errno));
  _lprintf(logger, LOG_DEBUG, _file, _line,
      "Mapped file %s for %s (%p)\n", fname_and_sfx,
      writable ? "writing" : "reading", map);
  if (sfx != NULL)
    free(fname_and_sfx);
  return map;
}

void xfunmap(void *map, size_t fsize, int fd, Logger *logger)
{
  LCHECK(munmap(map, fsize) == 0, logger, "Unmapping failed\n");
  LCHECK(close(fd) == 0, logger, "Closing file descriptor failed, %s\n",
         strerror(errno));
}

void _xfclose(FILE** file, Logger *logger, const char *_file, int _line)
{
  assert(file != NULL);
  assert(*file != NULL);
  _LCHECK(fclose(*file) == 0, logger, _file, _line,
           "Closing file failed, %s\n", strerror(errno));
  _lprintf(logger, LOG_DEBUG, _file, _line,
            "closed file... (%p)\n", *file);
  *file = NULL;
}

void xfrename(const char *fname_old, const char *sfx_old,
              const char *fname_new, const char *sfx_new,
              Logger *logger)
{
  char *fname_and_sfx_old, *fname_and_sfx_new;
  assert(fname_old != NULL);
  assert(fname_new != NULL);
  fname_and_sfx_old = (sfx_old != NULL) ?
    strcat_alloc(fname_old, sfx_old) : (char*)fname_old;
  fname_and_sfx_new = (sfx_new != NULL) ?
    strcat_alloc(fname_new, sfx_new) : (char*)fname_new;
  LCHECK(rename(fname_and_sfx_old, fname_and_sfx_new) == 0, logger,
    "Renaming file '%s' to '%s' failed, %s\n", fname_and_sfx_old,
    fname_and_sfx_new, strerror(errno));
  if (sfx_old != NULL)
    free(fname_and_sfx_old);
  if (sfx_new != NULL)
    free(fname_and_sfx_new);
}

void xfunlink(const char *fname, const char *sfx, Logger *logger)
{
  char *fname_and_sfx;
  assert(fname != NULL);
  fname_and_sfx = (sfx != NULL) ? strcat_alloc(fname, sfx) : (char*)fname;
  LCHECK(unlink(fname_and_sfx) == 0, logger, "Deleting file '%s' failed, %s\n",
              fname_and_sfx, strerror(errno));
  if (sfx != NULL)
    free(fname_and_sfx);
}

bool file_exists(const char *fname, const char *sfx)
{
  bool result;
  struct stat buffer;
  char *fname_and_sfx;
  assert(fname != NULL);
  fname_and_sfx = (sfx != NULL) ? strcat_alloc(fname, sfx) : (char*)fname;
  result = stat(fname_and_sfx, &buffer) == 0 ? true : false;
  if (sfx != NULL)
    free(fname_and_sfx);
  return result;
}

/* size of a seekable file, by seeking to its end; the position is then
   set to <pos> (-1: the previous one) */
static long seek_size(FILE *file, long pos, Logger *logger)
{
  long curpos = ftell(file), fsize;
  LCHECK(curpos >= 0 && fseek(file, 0, SEEK_END) == 0 &&
         (fsize = ftell(file)) >= 0 &&
         fseek(file, pos < 0 ? curpos : pos, SEEK_SET) == 0, logger,
         "Determining the size of a file failed, %s\n", strerror(errno));
  return fsize;
}

size_t file_size(FILE *file)
{
  assert(file != NULL);
  return (size_t)seek_size(file, -1, NULL);
}

void read_file(FILE *file, char **buffer, size_t *bufsize, Logger *logger)
{
  assert(file != NULL);
  *bufsize = (size_t)seek_size(file, 0, logger);
  XMALLOC(*buffer, *bufsize + 1);
  if (*bufsize > 0)
    xfread(*buffer, 1, *bufsize, file, logger);
  (*buffer)[*bufsize] = '\0';
}

#define LONG_LINE_PFXLEN 40

size_t read_lines_in_buffer(FILE *file, char *buffer, const size_t bufsize,
    char ***lines, size_t *linesbufsize, Logger *logger)
{
  size_t nlines;
  size_t nitemsread, i, lastnewline;

  assert(file != NULL);
  assert(buffer != NULL);
  LCHECK(bufsize >= READ_LINES_MIN_BUFSIZE, logger,
         "read_lines_in_buffer: the buffer size must be at least %d bytes "
         "(%zu)\n", READ_LINES_MIN_BUFSIZE, bufsize);
  if (feof(file))
    return 0;
  nitemsread = fread(buffer, 1, bufsize-2, file);
  LCHECK(nitemsread == bufsize-2 || feof(file), logger,
      "IO Error during file read access\n");
  if (nitemsread == 0)
    return 0;
  if (feof(file) && buffer[nitemsread-1] != '\n')
  {
    buffer[nitemsread] = '\n';
    nitemsread++;
  }
  lastnewline = nitemsread;
  nlines = 0;
  XCHECKSPACEINCPERC(*lines, linesbufsize, 2*sizeof(char*), 50);
  (*lines)[0] = buffer;
  for (i = 0; i < nitemsread; i++)
  {
    if (buffer[i] == '\n')
    {
      lastnewline = i;
      XCHECKSPACEINCPERC(*lines, linesbufsize, (nlines+2)*sizeof(char*),
                         50);
      (*lines)[nlines+1] = buffer + i + 1;
      nlines++;
    }
  }
  buffer[lastnewline+1] = 0;
  if (nlines == 0)
  {
    /* the bytes read are buffer[0..nitemsread-1] */
    bool prefix = nitemsread > LONG_LINE_PFXLEN;
    buffer[prefix ? LONG_LINE_PFXLEN : nitemsread] = '\0';
    lprintf(logger, LOG_ERROR, "A line of the input file is larger "
        "than the available buffer size (%zu bytes).\n", bufsize);
    lprintf(logger, LOG_ERROR, "Long line%s: %s\n",
        prefix ? " starts with" : "", buffer);
    lprintf(logger, LOG_ERROR, "Suggestion: determine the maximal length "
        " of lines in the file "
        "(e.g. awk 'length>m{m=length} END{print m}' filename)"
        "  and re-run the command with a larger value of the buffer "
        "using the buffer size options (see --help)\n");
    exit(EXIT_FAILURE);
  }
  if (!feof(file))
  {
    long rel_offset = ((long)lastnewline - ((long)nitemsread-1));
    LCHECK(fseek(file, rel_offset, SEEK_CUR) == 0, logger,
           "Seeking in a file failed, %s\n", strerror(errno));
  }
  return nlines;
}

static inline
void wrongsig(const char *fname, const char *sfx, const char *desc,
              const char *helpmsg, Logger *logger)
{
  lprintf(logger, LOG_ERROR, "%s %s%s has a wrong file signature\n",
      desc, fname, sfx != NULL ? sfx : "");
  if (helpmsg)
    lprintf(logger, LOG_ERROR, "%s\n", helpmsg);
  exit(EXIT_FAILURE);
}

void xchksig(const char *signature, const char *expected,
             const char *fname, const char *sfx,
             const char *desc, const char *helpmsg, Logger *logger)
{
  if (strncmp(signature, expected, strlen(expected)) != 0)
    wrongsig(fname, sfx, desc, helpmsg, logger);
}

FILE* _xfopen_chksig(const char *fname, const char *sfx, const char *mode,
                     const char *description,
                     const char *expected_signature, const char *helpmsg,
                     Logger *logger, const char *_file, int _line)
{
  FILE *result = _xfopen(fname, sfx, mode, description, logger, _file, _line);
  char *signature;
  size_t len = strlen(expected_signature);
  assert(mode[0] == 'r');
  XMALLOC(signature, len);
  xfread(signature, (size_t)1, len, result, logger);
  xchksig(signature, expected_signature, fname, sfx, description, helpmsg,
      logger);
  free(signature);
  return result;
}

FILE* _xfopen_putsig(const char *fname, const char *sfx, const char *mode,
                     const char *description,
                     const char *signature, Logger *logger,
                     const char *_file, int _line)
{
  FILE *result = _xfopen(fname, sfx, mode, description, logger, _file, _line);
  size_t len = strlen(signature);
  assert(mode[0] == 'w');
  xfwrite(signature, (size_t)1, len, result, logger);
  return result;
}

void *_xfmap_chksig(const char *fname, const char *sfx,
    size_t *fsize, int *fd, const char *expected_signature,
    const char *filedesc, const char *helpmsg, Logger *logger,
    const char *_file, int _line)
{
  void *result = _xfmap(fname, sfx, false, fsize, fd, logger, _file, _line);
  size_t len = strlen(expected_signature);
  if (*fsize < len)
  {
    lprintf(logger, LOG_ERROR, "%s %s%s is too short (%zu bytes) to contain "
        "the file signature (%zu bytes)\n", filedesc, fname,
        sfx != NULL ? sfx : "", *fsize, len);
    if (helpmsg)
      lprintf(logger, LOG_ERROR, "%s\n", helpmsg);
    exit(EXIT_FAILURE);
  }
  xchksig((char*)result, expected_signature, fname, sfx, filedesc, helpmsg,
      logger);
  return result;
}

bool is_directory(const char *path)
{
  struct stat st;
  return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

void xmakedirs(const char *path, Logger *logger)
{
  char *p = xstrdup(path);
  /* from the second character (the first is the root, if '/'); none for
     an empty path, for which mkdir() fails (ENOENT) */
  for (char *s = p + (p[0] != '\0'); *s != '\0'; s++)
    if (*s == '/') {
      *s = '\0';
      mkdir(p, 0777);
      *s = '/';
    }
  LCHECK(mkdir(p, 0777) == 0 || errno == EEXIST, logger,
         "Creating directory '%s' failed, %s\n", p, strerror(errno));
  free(p);
}

char *path_join(const char *dir, const char *name) {
  assert(dir != NULL);
  assert(name != NULL);
  size_t len = strlen(dir), name_len = strlen(name);
  bool slash = len > 0 && dir[len - 1] != '/';
  char *path;
  XMALLOC(path, len + slash + name_len + 1);
  memcpy(path, dir, len);
  if (slash)
    path[len] = '/';
  memcpy(path + len + slash, name, name_len + 1);
  return path;
}

char *path_basename(const char *path)
{
  size_t len = strlen(path);
  while (len > 1 && path[len - 1] == '/')
    len--;
  size_t start = len;
  while (start > 0 && path[start - 1] != '/')
    start--;
  char *result;
  XMALLOC(result, len - start + 1);
  memcpy(result, path + start, len - start);
  result[len - start] = '\0';
  return result;
}

static int cmpstrptr(const void *a, const void *b)
{
  return strcmp(*(char *const *)a, *(char *const *)b);
}

char **list_dir(const char *dir, list_dir_filter_t filter, void *data,
    size_t *n, Logger *logger)
{
  DIR *d = opendir(dir);
  LCHECK(d != NULL, logger, "Opening directory '%s' failed, %s\n", dir,
         strerror(errno));
  char **paths = NULL;
  size_t cap = 0;
  *n = 0;
  struct dirent *e;
  /* readdir() returns NULL both at the end of the directory and on error,
     setting errno only on error */
  while (errno = 0, (e = readdir(d)) != NULL) {
    if (e->d_name[0] == '.')
      continue;
    char *path = path_join(dir, e->d_name);
    if (filter == NULL || filter(e->d_name, path, data)) {
      XCHECKSPACEINCPERC(paths, &cap, (*n + 1) * sizeof(char *), 50);
      paths[(*n)++] = path;
    } else
      free(path);
  }
  LCHECK(errno == 0, logger, "Reading directory '%s' failed, %s\n", dir,
         strerror(errno));
  LCHECK(closedir(d) == 0, logger, "Closing directory '%s' failed, %s\n",
         dir, strerror(errno));
  if (*n > 1)
    qsort(paths, *n, sizeof(char *), cmpstrptr);
  return paths;
}
