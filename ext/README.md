# External libraries

- `libcbaia/`: libcbaia (https://github.com/baia-ipc/libcbaia, ISC license),
  vendored with `git subtree` (see `docs/dev_manual.md`); it includes zlib;
- `libyaml/`: libyaml 0.2.5 (https://github.com/yaml/libyaml, MIT license,
  see `libyaml/License`), the source files of the library from the release
  tarball (`yaml-0.2.5.tar.gz`, SHA-256
  `c642ae9b75fee120b2d96c712538bd2cf283228d2337df2cf2988e3c02678ef4`),
  unchanged.

zlib and libyaml are not used (the system libraries are used instead)
with `make SYSTEM_LIBS=1`.

Copied unchanged from their repositories:

- `edlib/`: edlib 1.2.7 (https://github.com/Martinsos/edlib, MIT license,
  see `edlib/LICENSE`);
- `libdeflate/`: libdeflate 1.26 (https://github.com/ebiggers/libdeflate,
  tag `v1.26`, commit `92e6a0db9fa848d742f9eb286c92afc60f2c3dda`, MIT
  license, see `libdeflate/COPYING`), the sources of the library only:
  `lib/`, `libdeflate.h` and `common_defs.h` (included by `lib/`); it
  compresses the output; it is used also with `make SYSTEM_LIBS=1`.

khash (klib) is included in libcbaia (`libcbaia/ext/khash`).
