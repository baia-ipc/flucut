/* (c) 2017-2021 Giorgio Gonnella */
/* File system interaction */

#ifndef FILE_H
#define FILE_H

#include <stdbool.h>
#include "logger.h"

/* Logger can be always set to NULL, if no logger available */

/*
   Open a file and exit the program on failure.

   Input parameters:
   - <fname> the filename
   - <sfx> an optional suffix; if given, it is appended to fname;
           use NULL if not necessary
   - <mode> the file opening mode
   - <description> description of the file
   - <logger> logger for the verbose msg (may be NULL)

   The function calls fopen() and the program is exited with an error
   message if it fails (returns NULL).
*/
FILE* _xfopen(const char *fname, const char *sfx, const char *mode,
            const char *description, Logger *logger, const char *_file,
            int _line);
#define xfopen(fname, sfx, mode, description, logger) \
  _xfopen(fname, sfx, mode, description, logger, __FILE__, __LINE__)

/*
   Close a file pointer and exit the program on failure.

   Input parameters:
   - <file> the file

   The function calls fclose() and the program is exited if
   fclose() does not return 0.
*/
void _xfclose(FILE** file, Logger *logger, const char *_file, int _line);
#define xfclose(file, logger) _xfclose(file, logger, __FILE__, __LINE__)

/*
   Delete the specified file.

   Input parameters:
   - <fname> the filename
   - <sfx> suffix to append to the filename;
           set to NULL if not required

  The function calls unlink() and the program is exited if
  unlink() does not return 0.
*/
void xfunlink(const char *fname, const char *sfx, Logger *logger);

/*
   Rename the specified file.

   Input parameters:
   - <fname_old>, <fname_new> the filenames
   - <sfx_old>, <sfx_new> suffix to append to the filenames;
                          set to NULL if not required

  The function calls rename() and the program is exited if
  rename() does not return 0.
*/
void xfrename(const char *fname_old, const char *sfx_old,
              const char *fname_new, const char *sfx_new,
              Logger *logger);

/*
   Allocate a buffer with the size of the given file
   and read the file content into the buffer.

   Input parameters:
   - <file> the input file (seekable)

   Output parameters:
   - <buffer> set to the allocated buffer, containing the content followed
              by a \0 (so that a text can be used as a string)
   - <bufsize> set to the size of the content (0 for an empty file)

   The program is exited with an error message if the file cannot be read.
*/
void read_file(FILE *file, char **buffer, size_t *bufsize, Logger* logger);

/*
   Read from stream using fread() and exit program with
   an error message if the number of elements read is not
   equal to nmemb.
*/
#define xfread(PTR, SIZE, NMEMB, STREAM, LOGGER) do {\
  size_t _NREAD = fread(PTR, SIZE, NMEMB, STREAM), _NMEMB = (NMEMB);\
  LCHECK(_NREAD == _NMEMB, LOGGER, "I/O error while reading a file\n");\
} while (0)

/*
   Write to stream using fwrite() and exit program with
   an error message if the number of elements written is not
   equal to nmemb.
*/
#define xfwrite(PTR, SIZE, NMEMB, STREAM, LOGGER) do { \
  size_t _NWRITTEN = fwrite(PTR, SIZE, NMEMB, STREAM), _NMEMB = (NMEMB);\
  LCHECK(_NWRITTEN == _NMEMB, LOGGER,\
      "I/O error while writing a file\n");\
} while (0)

/* minimum buffer size of read_lines_in_buffer(): one byte of the file, the
   newline appended at the end of the file, and the \0 */
#define READ_LINES_MIN_BUFSIZE 3

/*
   Read as many lines as possible from file content into the buffer
   (at most bufsize-2 bytes) and set the next available byte to \0.
   If a text file is read, then the \0 is guaranteed to be after a
   newline or the EOF. If the EOF is reached and the last char is
   not a newline, a newline is appended before the \0. Count the number
   of lines read.

   Input parameters:
   - <file> the input file
   - <buffer> the buffer
   - <bufsize> the size of the buffer, at least READ_LINES_MIN_BUFSIZE

   Return value:
   - number of lines read

   Fails if:
   - bufsize is less than READ_LINES_MIN_BUFSIZE
   - fread on the file fails
   - no newline has been found and EOF has not been reached
     (ie. a line is longer than the buffer)
*/
size_t read_lines_in_buffer(FILE *file, char *buffer, const size_t bufsize,
    char ***lines, size_t *linesbufsize, Logger *logger);

/*
  Allocate a char buffer a compute in it the concatenation of
  basename, an optional sfx (use NULL if not required), a dot
  followed by an unsigned integer n.

  Examples:
  basename = "x", sfx = ".y", n = 1, return value: "x.y.1";
  basename = "x", sfx = NULL, n = 0, return value: "x.0".

  The user must remember to free the buffer containing the returned value.
*/
char *compute_nth_filename(const char *basename, const char *sfx, size_t n);

/*
  Call getline() on a file. If the return value is -1,
  return false if end of file, otherwise exit the program
  with an error message. If the return value is not -1,
  return true.
*/
bool xgetline(char **line, size_t *line_size, FILE *file, Logger *logger);

/*
  Call getdelim() on a file. If the return value is -1,
  return false if end of file, otherwise exit the program
  with an error message. If the return value is not -1,
  return true.
*/
bool xgetdelim(char **line, size_t *line_size, int delim, FILE *file,
               Logger *logger);

/* Map file using mmap; an empty file cannot be mapped, the program is
   then exited with an error message */
void* _xfmap(const char *fname, const char *sfx, bool writable,
    size_t *fsize, int *fd, Logger *logger, const char *_file, int _line);
#define xfmap(fname, sfx, writable, fsize, fd, logger) \
   _xfmap(fname, sfx, writable, fsize, fd, logger, __FILE__, __LINE__)

/* Unmap file mapped using mmap */
void xfunmap(void *map, size_t fsize, int fd, Logger *logger);

/* check if a file exists */
bool file_exists(const char *fname, const char *sfx);

/* get file size (the file must be seekable; the program is exited with an
   error message otherwise); the position in the file is not changed */
size_t file_size(FILE *file);

/* check that a file contains the given expected signature;
   if not, the program is exited and an error message is visualized:

   signature: string containing the signature, it is not required to be
              zero-terminated
   expected: expected signature, zero-terminated
   fname/sfx: used only for error messages
   desc: a descriptive string of the file type/function
   helpmsg (can be NULL): an help message (a newline is appended to it)
*/
void xchksig(const char *signature, const char *expected,
             const char *fname, const char *sfx,
             const char *desc, const char *helpmsg,
             Logger *logger);

/* same as xfopen (see above), but also checks the signature;
   the mode string must start with 'r';
   the offset in the file is set to just after the signature */
FILE* _xfopen_chksig(const char *fname, const char *sfx, const char *mode,
                     const char *description,
                     const char *expected_signature, const char *helpmsg,
                     Logger *logger, const char *_file, int _line);
#define xfopen_chksig(fname, sfx, mode, description, \
                    expected_signature, helpmsg, logger) \
 _xfopen_chksig(fname, sfx, mode, description, \
                 expected_signature, helpmsg, logger, __FILE__, __LINE__)

/* same as xfopen (see above), but also writes the signature;
   the mode string must start with 'w';
   the offset in the file is set to just after the signature */
FILE* _xfopen_putsig(const char *fname, const char *sfx, const char *mode,
                     const char *description,
                     const char *signature, Logger *logger, const char *_file,
                     int _line);
#define xfopen_putsig(fname, sfx, mode, description, signature, logger) \
  _xfopen_putsig(fname, sfx, mode, description, signature, logger,\
      __FILE__, __LINE__)

/* Map file for reading using mmap and check the signature; a file
   shorter than the signature is an error, as a wrong signature */
void* _xfmap_chksig(const char *fname, const char *sfx,
    size_t *fsize, int *fd, const char *expected_signature,
    const char *filedesc, const char *helpmsg, Logger *logger,
    const char *_file, int _line);
#define xfmap_chksig(fname, sfx, fsize, fd, expected_signature, \
                     filedesc, helpmsg, logger) \
  _xfmap_chksig(fname, sfx, fsize, fd, expected_signature, \
                 filedesc, helpmsg, logger, __FILE__, __LINE__)

/* Whether path is an existing directory */
bool is_directory(const char *path);

/* Create a directory and its parents, if they do not exist (as mkdir -p);
   the program is exited with an error message if this fails */
void xmakedirs(const char *path, Logger *logger);

/* <dir>/<name> (allocated, to be freed by the caller), without doubling
   a slash at the end of dir; name if dir is empty */
char *path_join(const char *dir, const char *name);

/* Last component of a path, without trailing slashes (allocated, to be
   freed by the caller); e.g. "run/fastq_pass/" -> "fastq_pass" */
char *path_basename(const char *path);

/* Filter for list_dir(): name: the entry name; path: <dir>/<name>; data:
   the data passed to list_dir() */
typedef bool (*list_dir_filter_t)(const char *name, const char *path,
    void *data);

/* Paths (<dir>/<name>) of the entries of the directory dir, in
   lexicographic order, without the hidden entries (name starting with a
   dot); if filter is not NULL, only those for which it returns true;
   the number of paths is stored in n; the array and the paths are
   allocated, to be freed by the caller; the program is exited with an
   error message if the directory cannot be read */
char **list_dir(const char *dir, list_dir_filter_t filter, void *data,
    size_t *n, Logger *logger);

#endif
