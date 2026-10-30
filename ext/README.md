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

Copied unchanged from its repository:

- `edlib/`: edlib 1.2.7 (https://github.com/Martinsos/edlib, MIT license,
  see `edlib/LICENSE`).

khash (klib) is included in libcbaia (`libcbaia/ext/khash`).
