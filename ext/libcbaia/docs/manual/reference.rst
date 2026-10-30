.. _reference:

API reference
-------------

The declarations of the headers of the library (generated from the headers by ``docs/manual/make_reference.py``; the documentation of each function is in the header, before its declaration). The unit test functions (``<module>_test()``) and the internal functions (``_`` or ``__`` prefix, called by macros, or ``_`` suffix, called by inline functions) are not listed.

Core (``core``)
~~~~~~~~~~~~~~~~

**alloc.h**

Memory allocation macros

.. code-block:: c

   #define XMALLOC(PTR, BUFSIZE)
   #define XCALLOC(PTR, NITEMS, ITMSIZE)
   #define XREALLOC(PTR, REQSIZE)
   #define XCHECKSPACE(PTRPTR, SIZEPTR, REQSIZE)
   #define XCHECKSPACEINC(PTR, SIZEPTR, REQSIZE, INC)
   #define XCHECKSPACEINCPERC(PTR, SIZEPTR, REQSIZE, INCPERC)
   #define XFREE(PTR)

**argparse-def.h**

.. code-block:: c

   #define VERSION
   #define STR(S)
   #define XSTR(S)
   #define BOOLSTR(BOOL)
   #define PRINT_PURPOSE(STREAM)
   #define HELP_AND_EXIT
   #define PURPOSE_AND_EXIT
   #define VERSION_AND_EXIT
   #define LOGGER_OPTS_VAR_DECLARE
   #define PARSE_START
   #define PARSE_ERR(FMT, ...)
   #define OPT_ERR(C, FMT, ...)
   #define ARG_ERR(N, FMT, ...)
   static inline bool argparse_long(const char *s, long *value);
   #define ARGPARSE_STORE(DESTPTR, VALUE)
   #define PARSE_INTEGER_OPT_MIN(C, DESTPTR, MINVALUE)
   #define PARSE_INTEGER_OPT_MIN_MAX(C, DESTPTR, MINVALUE, MAXVALUE)
   #define PARSE_INTEGER_ARG_MIN(N, DESTPTR, MINVALUE)
   #define LOGFILE_OPT
   #define LOGFILE_OPT_STR
   #define LOGFILE_OPT_HELP
   #define LOGGER_LONGOPTS
   #define COMMON_LONGOPTS
   #define LONGOPTS_END
   #define LOGGER_OPTS_HELP
   #define COMMON_OPTS_HELP
   #define LOGGER_OPTS_STR
   #define COMMON_OPTS_STR
   #define LOGGER_OPTS_CHKMULTI
   #define LOGGER_OPTS_PARSE_CASES
   #define COMMON_OPTS_PARSE_CASES
   #define CHECK_NARGS(EXPECTED)
   #define CHECK_MIN_NARGS(EXPECTED)
   #define LOGGER_OPTS_APPLY(LOGGER)
   #define INIT_LOGGER(LOGGER)
   #define PARSE_LOG(FMT, ...)
   #define PARSE_LOG_N(FMT, ...)

**colors.h**

ANSI colors of the text written to a terminal

.. code-block:: c

   #define START_COLOR(STREAM, COLOR)
   #define STOP_COLOR(STREAM)
   #define COLOR_BLACK
   #define COLOR_RED
   #define COLOR_GREEN
   #define COLOR_YELLOW
   #define COLOR_BLUE
   #define COLOR_MAGENTA
   #define COLOR_CYAN
   #define COLOR_WHITE
   #define COLOR_BRIGHT_BLACK
   #define COLOR_BRIGHT_RED
   #define COLOR_BRIGHT_GREEN
   #define COLOR_BRIGHT_YELLOW
   #define COLOR_BRIGHT_BLUE
   #define COLOR_BRIGHT_MAGENTA
   #define COLOR_BRIGHT_CYAN
   #define COLOR_BRIGHT_WHITE

**file.h**

File system interaction

.. code-block:: c

   #define xfopen(fname, sfx, mode, description, logger)
   #define xfclose(file, logger)
   void xfunlink(const char *fname, const char *sfx, Logger *logger);
   void xfrename(const char *fname_old, const char *sfx_old,
                 const char *fname_new, const char *sfx_new,
                 Logger *logger);
   void read_file(FILE *file, char **buffer, size_t *bufsize, Logger* logger);
   #define xfread(PTR, SIZE, NMEMB, STREAM, LOGGER)
   #define xfwrite(PTR, SIZE, NMEMB, STREAM, LOGGER)
   #define READ_LINES_MIN_BUFSIZE
   size_t read_lines_in_buffer(FILE *file, char *buffer, const size_t bufsize,
       char ***lines, size_t *linesbufsize, Logger *logger);
   char *compute_nth_filename(const char *basename, const char *sfx, size_t n);
   bool xgetline(char **line, size_t *line_size, FILE *file, Logger *logger);
   bool xgetdelim(char **line, size_t *line_size, int delim, FILE *file,
                  Logger *logger);
   #define xfmap(fname, sfx, writable, fsize, fd, logger)
   void xfunmap(void *map, size_t fsize, int fd, Logger *logger);
   bool file_exists(const char *fname, const char *sfx);
   size_t file_size(FILE *file);
   void xchksig(const char *signature, const char *expected,
                const char *fname, const char *sfx,
                const char *desc, const char *helpmsg,
                Logger *logger);
   #define xfopen_chksig(fname, sfx, mode, description, expected_signature, helpmsg, logger)
   #define xfopen_putsig(fname, sfx, mode, description, signature, logger)
   #define xfmap_chksig(fname, sfx, fsize, fd, expected_signature, filedesc, helpmsg, logger)
   bool is_directory(const char *path);
   void xmakedirs(const char *path, Logger *logger);
   char *path_join(const char *dir, const char *name);
   char *path_basename(const char *path);
   typedef bool (*list_dir_filter_t)(const char *name, const char *path,
       void *data);
   char **list_dir(const char *dir, list_dir_filter_t filter, void *data,
       size_t *n, Logger *logger);

**logger.h**

Logging library

.. code-block:: c

   typedef struct Logger Logger;
   typedef enum LoggingLevel {
     LOG_ERROR   = 0,
     LOG_WARN    = 1,
     LOG_SUCCESS = 2,
     LOG_PHASE   = 3,
     LOG_INFO    = 4,
     LOG_PARAM   = 5,
     LOG_MEM     = 6,
     LOG_DEBUG   = 7,
     LOG_TRACE   = 8
   } LoggingLevel;
   Logger* new_logger(FILE *logfile, LoggingLevel level);
   void delete_logger(Logger* logger);
   #define lprintf(LOGGER, LEVEL, FMT, ...)
   void lnprintf(Logger *logger, LoggingLevel level, const char *fmt, ...)
        __attribute__ ((format (printf,3,4)));
   #define LCHECK(CONDITION, LOGGER, FMT, ...)
   #define fatal_printf(FMT, ...)
   #define showmem(logger)
   #define TRACE(FMT, ...)
   void NTRACE(const char *fmt, ...)
        __attribute__ ((format (printf,1,2)));
   #define DEBUG_TRACE(FMT, ...)
   #define DEBUG_NTRACE(FMT, ...)

**numutils.h**

Numeric utilities: integer division rounding towards minus infinity, the median of a histogram, and the rounding of a value to a "nice" number (for the axes of plots).

.. code-block:: c

   static inline long floordiv(long a, long b);
   size_t hist_median(const size_t *hist, size_t n_bins, size_t bin_width,
       size_t total);
   double nice_number(double v);

**progress.h**

Progress display in the style of the Python tqdm package, e.g.:

.. code-block:: c

   #define PROGRESS_MIN_INTERVAL
   typedef struct progress_t progress_t;
   progress_t* progress_new(FILE *file, const char *description, size_t total);
   void progress_update(progress_t *self, size_t n);
   void progress_finish(progress_t *self);
   size_t progress_count(const progress_t *self);
   void progress_delete(progress_t *self);

**strbuf.h**

Growable string buffer, for building output text without allocations for each piece (the memory is kept when the buffer is cleared).

.. code-block:: c

   #define STRBUF_INC
   typedef struct {
     char *s;
     size_t len;
     size_t cap;
   } strbuf_t;
   static inline void strbuf_init(strbuf_t *self);
   static inline void strbuf_free(strbuf_t *self);
   static inline void strbuf_clear(strbuf_t *self);
   static inline void strbuf_truncate(strbuf_t *self, size_t len);
   static inline void strbuf_append(strbuf_t *self, const char *s,
       size_t len);
   static inline void strbuf_append_str(strbuf_t *self, const char *s);
   static inline void strbuf_append_char(strbuf_t *self, char c);
   static inline char *strbuf_reserve(strbuf_t *self, size_t n);
   static inline void strbuf_commit(strbuf_t *self, size_t n);
   static inline void strbuf_append_long(strbuf_t *self, long value);

**strutils.h**

.. code-block:: c

   char *strcat_alloc(const char *s1, const char *s2);
   char *xstrdup(const char *s);
   bool str_has_suffix(const char *s, const char *sfx);
   void fputs_json_string(const char *s, FILE *f);
   void fputs_html_escaped(const char *s, FILE *f);
   void format_duration(double seconds, char *buf, size_t size);
   bool parse_duration(const char *s, double default_unit, long *seconds);
   static inline bool is_space(char c);

**test.h**

.. code-block:: c

   #define test_assert(TEST)
   #define test_eq(FOUND, EXP, TYPE, FMT)
   #define test_eq_size(FOUND, EXP)
   #define test_eq_int(FOUND, EXP)
   #define test_eq_long(FOUND, EXP)
   #define test_eq_str(FOUND, EXP)
   #define test_eq_char(FOUND, EXP)
   #define test_eq_mem(FOUND, LEN, EXP)

Data structures (``ds``)
~~~~~~~~~~~~~~~~~~~~~~~~~

**bitwise.h**

.. code-block:: c

   #define GET_BIT(E, BITNUM)
   #define SET_BIT(E, BITNUM)
   #define UNSET_BIT(E, BITNUM)
   static inline uint64_t bitwise_mix64(uint64_t x);

**bloom-filter.h**

Bloom filter for integer keys (two hash functions, derived from the splitmix64 hash of the key).

.. code-block:: c

   typedef struct bloom_filter_t bloom_filter_t;
   bloom_filter_t *bloom_filter_new(size_t size_exp);
   void bloom_filter_insert(bloom_filter_t *self, size_t key);
   bool bloom_filter_query(const bloom_filter_t *self, size_t key);
   void bloom_filter_delete(bloom_filter_t *self);

**heap.h**

Binary min-heap of pointers, ordered by a comparison function (the smallest element is extracted first), e.g. for priority queues or for merging sorted sequences.

.. code-block:: c

   typedef int (*heap_cmp_t)(const void *a, const void *b);
   typedef struct heap_t heap_t;
   heap_t* heap_new(heap_cmp_t cmp, size_t capacity);
   void heap_push(heap_t *self, void *element);
   void *heap_pop(heap_t *self);
   void *heap_peek(const heap_t *self);
   size_t heap_size(const heap_t *self);
   bool heap_is_empty(const heap_t *self);
   void heap_delete(heap_t *self);

**hyperloglog.h**

HyperLogLog: estimation of the number of distinct elements of a multiset, using a small fixed amount of memory: the estimator and its small-range correction (linear counting, for estimates up to 2.5 times the number of registers) of Flajolet et al. 2007; 64-bit hash values, as in Heule et al. 2013, so that no large-range correction is needed.

.. code-block:: c

   #define HLL_PRECISION
   #define HLL_REGISTERS
   typedef struct {
     uint8_t registers[HLL_REGISTERS];
   } hll_t;
   void hll_init(hll_t *self);
   uint64_t hll_hash(const void *data, size_t len);
   void hll_add(hll_t *self, uint64_t hash);
   void hll_merge(hll_t *self, const hll_t *other);
   double hll_estimate(const hll_t *self);

**set.h**

Implementation: - Sorted array S plus an unsorted buffer B - Query is by binary search on S plus linear search on B - Insertion is in B - B is sorted and merged in S when full

.. code-block:: c

   typedef size_t set_elem_t;
   typedef struct Set Set;
   Set *set_new(size_t insertion_buffer_size, size_t bloom_filter_size,
                size_t k);
   void set_insert_simple(Set *set, set_elem_t key);
   void set_insert(Set *set, set_elem_t key);
   void set_finalize_insertions(Set *set);
   bool set_query(const Set *set, set_elem_t key);
   void set_delete(Set *set);

Large files and compression (``io``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**bgzf-compressor.h**

BGZF compression (the format of BAM files): the data are compressed in blocks of at most 64 KiB, each a gzip member with an extra field (BC) containing the size of the block. As for the gzip members of gzip-compressor.h, blocks can be compressed in parallel and written in order; a BGZF file ends with an empty block (the end-of-file marker).

.. code-block:: c

   #define BGZF_EOF_LEN
   extern const uint8_t bgzf_eof[BGZF_EOF_LEN];
   typedef struct bgzf_compressor_t bgzf_compressor_t;
   bgzf_compressor_t* bgzf_compressor_new(int level);
   void bgzf_compressor_run(bgzf_compressor_t *self, const char *data,
       size_t len, strbuf_t *out);
   void bgzf_compressor_delete(bgzf_compressor_t *self);

**file-merger.h**

Merging of sorted files

.. code-block:: c

   typedef bool (*merge_getelem_f)(char **elem, size_t *elemsize, FILE *file,
       Logger *logger);
   typedef void (*merge_writeelem_f)(FILE *file, char *elem, Logger *logger);
   typedef int (*merge_cmpelems_f)(const void *elem1, const void *elem2);
   void merge_n_files(size_t nfiles, const char *basename, const char *sfx,
       bool use_stdout, merge_getelem_f getelem, merge_writeelem_f writeelem,
       merge_cmpelems_f cmpelems, Logger *logger);

**file-sorter.h**

Sorting of large files in parts which fit in memory (to be merged with merge_n_files() of file-merger.h)

.. code-block:: c

   size_t sort_file_parts(
       const char *infname,
       const char *infsfx,
       const char *outfname,
       const char *outfsfx,
       size_t bufsize,
       size_t elems_bufsize_init,
       size_t read_elems(FILE *file, char *buffer, size_t bufsize,
                         char ***elems, size_t *elemsbufsize, Logger *logger),
       void sort_elems(char **elems, size_t nelems),
       void output_sorted_elems(FILE *file, char **elems, size_t nelems,
                                Logger *logger),
       Logger *logger);

**gzip-compressor.h**

Compression of data blocks as gzip members; the concatenation of gzip members is a valid gzip file, thus blocks can be compressed in parallel and written in order.

.. code-block:: c

   typedef struct gzip_compressor_t gzip_compressor_t;
   gzip_compressor_t* gzip_compressor_new(int level);
   void gzip_compressor_run(gzip_compressor_t *self, const char *data,
       size_t len, strbuf_t *out);
   void gzip_compressor_delete(gzip_compressor_t *self);

**html-page.h**

Self-contained HTML reports: a page with a style for light and dark themes (no external resources), tiles with the main values, and plots as inline SVG (horizontal bar chart, histogram, rank plot).

.. code-block:: c

   extern const char *const html_page_css;
   #define HTML_PAGE_MIN_RELOAD
   void html_page_start(html_writer_t *w, const char *title, long refresh);
   typedef struct {
     const char *src, *src_dark, *alt;
     int height;
   } html_logo_t;
   void html_page_start_logo(html_writer_t *w, const char *title,
       const html_logo_t *logo, const char *heading, long refresh);
   void html_page_head_end(html_writer_t *w);
   void html_dashboard_start(html_writer_t *w, const char *title,
       const html_logo_t *logo, const char *heading, long refresh);
   void html_dashboard_bar_end(html_writer_t *w);
   void html_panel_start(html_writer_t *w, const char *title, const char *cls);
   void html_panel_end(html_writer_t *w);
   void html_page_end(html_writer_t *w, const char *footer);
   void html_tile(html_writer_t *w, const char *value, const char *label);
   typedef struct {
     char s[32];
   } html_num_t;
   html_num_t html_num(size_t n);
   double html_percent(size_t a, size_t b);
   void html_hbars(html_writer_t *w, size_t n, const char *const *labels,
       const size_t *values, const char *const *classes, size_t total,
       double mean);
   void html_histogram(html_writer_t *w, const size_t *hist, size_t n_bins,
       size_t bin_width, const char *what, const char *unit,
       const char *empty);
   void html_knee(html_writer_t *w, const uint64_t *counts, size_t n,
       size_t n_selected, const char *x_label, const char *y_label,
       const char *selected_label);

**html-template.h**

HTML templates: static markup (e.g. the skeleton of a page) written as a normal, indented HTML file, embedded in the program at compile time (dev/embed.sh), with named placeholders filled when it is written:

.. code-block:: c

   #define HTML_TEMPLATE_MAX_VALUES
   #define HTML_TEMPLATE_MAX_NESTING
   #define html_template(W, TMPL, ...)
   void html_template_values(html_writer_t *w, const char *tmpl,
       const char *const *kv, size_t n);

**html-writer.h**

Writing HTML (and inline SVG) element by element: the text and the attribute values are escaped in one place, and the output is indented, one element per line, so that the nesting is visible both in the code and in the page.

.. code-block:: c

   #define HTML_WRITER_MAX_DEPTH
   typedef struct {
     FILE *f;
     const char *tags[HTML_WRITER_MAX_DEPTH];
     unsigned char content[HTML_WRITER_MAX_DEPTH];
     size_t depth;
     size_t inline_from;
     size_t svg;
     const char *start;
     bool start_void;
   } html_writer_t;
   void html_writer_init(html_writer_t *w, FILE *f);
   void html_writer_end(html_writer_t *w);
   void html_open(html_writer_t *w, const char *tag);
   void html_elem(html_writer_t *w, const char *tag);
   void html_attr(html_writer_t *w, const char *name, const char *fmt, ...)
     __attribute__((format(printf, 3, 4)));
   void html_inline(html_writer_t *w);
   void html_text(html_writer_t *w, const char *text);
   void html_textf(html_writer_t *w, const char *fmt, ...)
     __attribute__((format(printf, 2, 3)));
   void html_vtextf(html_writer_t *w, const char *fmt, va_list ap);
   void html_raw(html_writer_t *w, const char *markup);
   void html_close(html_writer_t *w);
   void html_tag_text(html_writer_t *w, const char *tag, const char *text);
   void html_tag_textf(html_writer_t *w, const char *tag, const char *fmt, ...)
     __attribute__((format(printf, 3, 4)));
   size_t html_begin_lines(html_writer_t *w);

**md5.h**

MD5 message digest (RFC 1321), e.g. for the checksums of files in metadata (as the md5 of the onlists of seqspec files) or for checking downloads. Not to be used for security purposes.

.. code-block:: c

   #define MD5_H
   #define MD5_HEX_LEN
   typedef struct {
     uint32_t state[4];
     uint64_t n_bytes;
     uint8_t buffer[64];
   } md5_t;
   void md5_init(md5_t *self);
   void md5_update(md5_t *self, const void *data, size_t len);
   void md5_final(md5_t *self, char *hex);
   void md5_file(const char *path, char *hex, Logger *logger);

Sequences and reads (``seq``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**bam-format.h**

Constants and encoding helpers of the BAM format, shared by bam-reader.c and ubam-writer.c (internal: not part of the API, it can change; it is installed with the other headers, as it is self-contained).

.. code-block:: c

   #define BAM_MAGIC
   #define BAM_MAGIC_LEN
   #define BAM_CORE_LEN
   #define BAM_BASES
   #define BAM_FUNMAP
   #define BAM_FREVERSE
   #define BAM_FSECONDARY
   #define BAM_FSUPPLEMENTARY
   #define BAM_UNMAPPED_BIN
   #define BAM_QUAL_MISSING
   #define BAM_PHRED_OFFSET
   #define BAM_QUAL_MAX
   static inline char bam_qual_char(uint8_t q);
   static inline uint8_t bam_qual_value(char c);
   static inline void bam_put16(strbuf_t *b, uint16_t v);
   static inline void bam_set32(char *p, uint32_t v);
   static inline void bam_put32(strbuf_t *b, uint32_t v);

**bam-reader.h**

Reading of BAM files, e.g. the unaligned BAM files written by Dorado and MinKNOW, in batches of FASTQ-like records (see fastq-reader.h).

.. code-block:: c

   typedef struct bam_reader_t bam_reader_t;
   bool bam_reader_is_bam(const char *filename, Logger *logger);
   bam_reader_t* bam_reader_new(const char *filename, Logger *logger);
   const char *bam_reader_header(const bam_reader_t *self);
   size_t bam_reader_read_batch(bam_reader_t *self, fastq_batch_t *batch,
       size_t max_records);
   void bam_reader_delete(bam_reader_t *self);
   bool bam_reader_has_eof_marker(const char *filename);
   size_t bam_reader_count(const char *filename, Logger *logger);

**cell-caller.h**

Identification of the cell barcodes of a single-cell library: the barcodes of a whitelist (e.g. all barcodes of a platform) found exactly in the reads are counted, and the barcodes of the cells are selected from the counts, given the expected number of cells N, as in Cell Ranger (version 2): the barcodes with a count of at least 1/10 of the 99th percentile of the counts of the N most frequent barcodes; as in Cell Ranger, the threshold is rounded to the nearest integer (halves upwards; e.g. 14 for a percentile of 135..144), and is at least 1.

.. code-block:: c

   typedef struct cell_caller_t cell_caller_t;
   cell_caller_t* cell_caller_new(void);
   void cell_caller_add(cell_caller_t *self, uint64_t barcode);
   void cell_caller_merge(cell_caller_t *self, cell_caller_t *other);
   whitelist_t* cell_caller_run(const cell_caller_t *self,
       size_t expected_cells, size_t length, const whitelist_t *known,
       const char *path, Logger *logger);
   size_t cell_caller_sorted_counts(const cell_caller_t *self,
       uint64_t **counts);
   void cell_caller_delete(cell_caller_t *self);

**complexity.h**

Library complexity of a sample: estimation of the number of distinct molecules of the library, from the number of reads and of distinct molecule keys (e.g. the UMI and the first and last bases of the insert, which identify the copies of a molecule without an alignment), for deciding whether deeper sequencing of a library is worthwhile.

.. code-block:: c

   #define COMPLEXITY_MIN_READS
   #define COMPLEXITY_MIN_ERROR_FREE
   typedef struct {
     bool estimated;
     bool high_error_rate;
     double error_free;
     double molecules_seen;
     double duplication;
     bool saturated_unknown;
     double molecules;
     double molecules_2x;
   } complexity_t;
   void complexity_estimate(size_t reads, double distinct_keys,
       double mean_key_length, double error_rate, complexity_t *out);

**edit-distance.h**

Bounded edit distance of short sequences (e.g. indices, barcodes, UMIs, compared with the sequences of a list): dynamic programming in the diagonal band of width 2k + 1 around the main diagonal, stopping as soon as a row exceeds the maximum distance k; for sequences of at most SHORT_EDIT_DISTANCE_MAX bases (for longer sequences, use e.g. edlib).

.. code-block:: c

   #define SHORT_EDIT_DISTANCE_MAX
   static inline long short_edit_distance(const char *a, size_t alen,
       const char *b, size_t blen, long k);
   static inline long short_prefix_distance(const char *query, size_t qlen,
       const char *target, size_t tlen, long k, size_t *end);

**fastq-reader.h**

Reading of FASTQ files (optionally gzip-compressed) in batches of records.

.. code-block:: c

   typedef struct fastq_reader_t fastq_reader_t;
   typedef struct {
     size_t title_from, title_len;
     size_t seq_from, seq_len;
     size_t qual_from;
     bool qual_missing;
   } fastq_record_t;
   typedef struct {
     strbuf_t data;
     fastq_record_t *records;
     size_t n_records, records_cap;
   } fastq_batch_t;
   fastq_reader_t* fastq_reader_new(const char *filename, Logger *logger);
   size_t fastq_reader_read_batch(fastq_reader_t *self, fastq_batch_t *batch,
       size_t max_records);
   void fastq_reader_delete(fastq_reader_t *self);
   size_t fastq_reader_count(const char *filename, Logger *logger);
   void fastq_batch_init(fastq_batch_t *batch);
   void fastq_batch_free(fastq_batch_t *batch);

**fastq-writer.h**

Writing of FASTQ records to a buffer (e.g. for a batch of records, which is then compressed with gzip-compressor.h and written to a file).

.. code-block:: c

   void fastq_writer_record(strbuf_t *out, const char *title, size_t title_len,
       const char *seq, const char *qual, size_t len);

**kmer-index.h**

Index of the k-mers of a set of sequences (e.g. reference sequences), for finding where a query sequence may align: the k-mers of the query which occur in the indexed sequences give candidate diagonals (position in the sequence minus position in the query), around which the query is then aligned by the caller (e.g. with edlib), instead of being aligned to the whole of every sequence.

.. code-block:: c

   #define KMER_INDEX_MAX_K
   typedef struct {
     uint32_t seq;
     long diagonal;
   } kmer_index_hit_t;
   typedef struct kmer_index_t kmer_index_t;
   kmer_index_t* kmer_index_new(const char *const *seqs, const size_t *lens,
       size_t n, size_t k);
   size_t kmer_index_k(const kmer_index_t *self);
   size_t kmer_index_count(const kmer_index_t *self, const char *kmer);
   size_t kmer_index_diagonals(const kmer_index_t *self, const char *query,
       size_t len, size_t max_occ, kmer_index_hit_t **hits, size_t *cap);
   void kmer_index_delete(kmer_index_t *self);

**minknow-metadata.h**

Read metadata of MinKNOW: the FASTQ files written by MinKNOW contain the read metadata in the title as fields key=value, separated by spaces (e.g. runid=... ch=12 start_time=...), while Dorado writes SAM tags (see sam-tags.h). The fields can be converted to SAM tags: runid: ri:Z, read: rn:i, ch: ch:i, start_time: st:Z, flow_cell_id: PU:Z, protocol_group_id: pg:Z, sample_id: si:Z, barcode: ob:Z, barcode_alias: oa:Z, parent_read_id: pi:Z, basecall_model_version_id: bm:Z.

.. code-block:: c

   #define MINKNOW_BARCODE_SAM_TAG
   bool minknow_metadata_append_sam_tag(strbuf_t *out, const char *field,
       size_t len);

**outfiles.h**

Output files of a program writing the reads to several destinations, e.g. a demultiplexer (a destination for each sample), with one compressed file for each destination: - for an input file: <outdir>/<input basename>.<destination>.fastq.gz, where the suffix of the input file (.fastq, .fq, optionally followed by .gz, or .bam; see reads_reader_is_reads_filename) is removed; - for a group of input files (directory mode): <outdir>/<group>/<destination>.fastq.gz.

.. code-block:: c

   #define OUTFILES_COMPRESSLEVEL
   typedef struct outfiles_t outfiles_t;
   outfiles_t* outfiles_new(const char *outdir, const char *reads_file,
       size_t n_destinations, Logger *logger);
   outfiles_t* outfiles_new_group(const char *outdir, const char *group,
       size_t n_destinations, Logger *logger);
   void outfiles_set_append(outfiles_t *self, bool append);
   void outfiles_set_bam(outfiles_t *self, size_t destination,
       const char *header, size_t len);
   #define OUTFILES_FASTQ_SFX
   #define OUTFILES_BAM_SFX
   void outfiles_write(outfiles_t *self, size_t destination, const char *name,
       const char *data, size_t len);
   void outfiles_delete(outfiles_t *self);

**reads-reader.h**

Reader of the reads of a FASTQ file (optionally gzipped) or of an unaligned BAM file, recognized from the content, with the same interface (batches of FASTQ records, see fastq-reader.h); functions for recognizing the files of reads by their name, counting the reads, and checking whether a file is complete (e.g. when it is written by another program, such as a basecaller).

.. code-block:: c

   typedef struct reads_reader_t reads_reader_t;
   reads_reader_t* reads_reader_new(const char *filename, Logger *logger);
   bool reads_reader_is_bam(const reads_reader_t *self);
   size_t reads_reader_read_batch(reads_reader_t *self, fastq_batch_t *batch,
       size_t max_records);
   void reads_reader_delete(reads_reader_t *self);
   bool reads_reader_is_reads_filename(const char *name);
   size_t reads_reader_filename_suffix_len(const char *name);
   size_t reads_reader_count(const char *filename, Logger *logger);
   bool reads_reader_file_complete(const char *path);

**sam-header.h**

SAM header text of BAM files written from BAM input (e.g. the unaligned BAM of Dorado), keeping the provenance of the reads: - the read groups (@RG) of the inputs, so that the records keep the RG tag of their input; a program which distributes the reads to several destinations (e.g. a demultiplexer) can rename them per destination; - the programs (@PG) of the inputs, followed by the @PG line of the program, linked to the last one of the inputs (PP), as samtools does.

.. code-block:: c

   typedef struct {
     strbuf_t rg_lines, pg_lines;
     char **rg_ids;
     size_t n_rg;
     char *last_pg;
   } sam_header_inputs_t;
   void sam_header_inputs_init(sam_header_inputs_t *self);
   void sam_header_inputs_clear(sam_header_inputs_t *self);
   void sam_header_inputs_free(sam_header_inputs_t *self);
   void sam_header_inputs_add(sam_header_inputs_t *self, const char *text);
   void sam_header_inputs_add_files(sam_header_inputs_t *self,
       char *const *files, size_t n_files, Logger *logger);
   void sam_header_append_read_groups(strbuf_t *text,
       const sam_header_inputs_t *inputs, const char *name);
   void sam_header_append_programs(strbuf_t *text,
       const sam_header_inputs_t *inputs, const char *program,
       const char *version, const char *command_line);
   char *sam_header_command_line(int argc, char *const argv[]);
   const char *sam_header_line_id(const char *line, size_t line_len,
       size_t *len);

**sam-tags.h**

SAM tags in text format, in read titles: the title of a read (FASTQ header without '@', or the title returned by the BAM reader) is the read name followed by fields; in the FASTQ written by Dorado, and in the output of programs which keep the tags (e.g. for minimap2 -y or samtools import -T '*'), the fields are SAM tags in text format (e.g. ch:i:12), separated by tabs. The uBAM writer encodes them in binary form. For the read metadata of MinKNOW (fields key=value), see minknow-metadata.h.

.. code-block:: c

   bool sam_tag_is_valid(const char *field, size_t len);
   void sam_tag_start(strbuf_t *out, const char *tag, char type);
   void sam_tag_append_string(strbuf_t *out, const char *tag,
       const char *value, size_t len);
   void sam_tag_append_int(strbuf_t *out, const char *tag, long value);

**sequtil.h**

Sequence utilities: reverse complement and substrings with the semantics of the Python slices (useful for reproducing exactly the output of Python programs).

.. code-block:: c

   void reverse_complement(const char *seq, size_t len, char *dest);
   void py_slice(long len, long start, long stop, long *from, long *slice_len);

**test-bam.h**

BAM data and files for the unit tests of the modules reading or writing BAM files (internal: not part of the API, it can change).

**ubam-writer.h**

Encoding of unaligned BAM data (header and records), to be compressed with BGZF (see bgzf-compressor.h).

.. code-block:: c

   #define UBAM_MAX_READ_NAME
   #define UBAM_NO_READ_NAME
   void ubam_writer_header(strbuf_t *out, const char *text);
   void ubam_writer_record(strbuf_t *out, const char *title, size_t title_len,
       const char *seq, const char *qual, size_t len);

**whitelist.h**

Whitelist of barcodes (e.g. the cell barcodes of a single-cell platform) and correction of the barcodes extracted from the reads.

.. code-block:: c

   #define WHITELIST_MAX_BARCODE_LENGTH
   #define WHITELIST_MAX_DIST
   typedef struct whitelist_t whitelist_t;
   typedef enum {
     WHITELIST_EXACT,
     WHITELIST_CORRECTED,
     WHITELIST_AMBIGUOUS,
     WHITELIST_NOT_FOUND,
   } whitelist_status_t;
   typedef struct {
     whitelist_status_t status;
     uint64_t barcode;
     long dist;
     size_t consumed;
   } whitelist_match_t;
   whitelist_t* whitelist_new(const char *path, Logger *logger);
   whitelist_t* whitelist_new_from_codes(const uint64_t *barcodes, size_t n,
       size_t length);
   bool whitelist_encode(const whitelist_t *self, const char *seq,
       uint64_t *barcode);
   bool whitelist_contains(const whitelist_t *self, uint64_t barcode);
   void whitelist_delete(whitelist_t *self);
   size_t whitelist_barcode_length(const whitelist_t *self);
   size_t whitelist_size(const whitelist_t *self);
   size_t whitelist_codes(const whitelist_t *self, uint64_t **codes);
   void whitelist_match(const whitelist_t *self, const char *seq, size_t len,
       bool from_end, long max_dist, whitelist_match_t *match);
   void whitelist_decode(const whitelist_t *self, uint64_t barcode, char *out);
