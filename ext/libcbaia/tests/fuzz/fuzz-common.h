/* (c) 2026 Giorgio Gonnella */

/*
   Common functions of the fuzzing harnesses (libFuzzer).

   The parsers exit the program on invalid input; in the harnesses, exit()
   is wrapped (linker option --wrap=exit), so that it returns to the
   harness (FUZZ_RUN), which then closes the files left open; the memory
   left allocated is not freed (thus the leak detection is disabled, see
   tests/fuzz/run.sh). Crashes, memory errors and undefined behaviour are
   reported by the sanitizers.

   Also used by the harnesses of the projects using the library (with
   their directory in the include path); FUZZ_NAME (default: "cbaia") is
   used in the name of the temporary input file.
*/

#ifndef FUZZ_COMMON_H
#define FUZZ_COMMON_H

#include <dirent.h>
#include <setjmp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef FUZZ_NAME
#define FUZZ_NAME "cbaia"
#endif

static jmp_buf fuzz_env;
static bool fuzz_in_target = false;

/* file descriptors open at the start (not closed after an exit) */
static bool fuzz_initial_fds[1024];

void __real_exit(int status) __attribute__((noreturn));

void __wrap_exit(int status) {
  if (fuzz_in_target)
    longjmp(fuzz_env, 1);
  __real_exit(status);
}

static void fuzz_record_fds(void) {
  DIR *dir = opendir("/proc/self/fd");
  if (dir == NULL) {
    perror("/proc/self/fd");
    abort();
  }
  struct dirent *e;
  while ((e = readdir(dir)) != NULL) {
    int fd = atoi(e->d_name);
    if (fd >= 0 && fd < 1024)
      fuzz_initial_fds[fd] = true;
  }
  closedir(dir);
}

/* close the files left open by a run interrupted by exit() */
static void fuzz_close_new_fds(void) {
  for (int fd = 0; fd < 1024; fd++)
    if (!fuzz_initial_fds[fd])
      close(fd);
}

/* write the input to a temporary file; returns its path */
static const char *fuzz_write_input(const uint8_t *data, size_t size) {
  static char path[64] = "";
  if (path[0] == '\0')
    snprintf(path, sizeof(path), "/tmp/" FUZZ_NAME "-fuzz-%d", (int)getpid());
  FILE *f = fopen(path, "wb");
  if (f == NULL) {
    perror(path);
    abort();
  }
  fwrite(data, 1, size, f);
  fclose(f);
  return path;
}

/* run CODE; if it calls exit(), run CLEANUP (e.g. to free the objects
   which are still reachable) and continue after it */
#define FUZZ_RUN_CLEANUP(CODE, CLEANUP) \
  do { \
    fuzz_in_target = true; \
    if (setjmp(fuzz_env) == 0) { \
      CODE; \
    } else { \
      CLEANUP; \
      fuzz_close_new_fds(); \
    } \
    fuzz_in_target = false; \
  } while (0)

#define FUZZ_RUN(CODE) FUZZ_RUN_CLEANUP(CODE, (void)0)

#endif
