# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

xnirspec is an interactive Fortran program (PGPLOT GUI on X11) to visualise and
manipulate FITS images, originally written for Keck/NIRSPEC spectroscopic
reductions (Cardiel et al. 2003, ApJ, 584, 76). It was originally Fortran 77
and was converted to Fortran 90 in January 2023 without a full
modernisation, so most code is still F77-style logic in `.f90` files.

## Build

Autotools project (`configure.ac`, `Makefile.am`, `src/Makefile.am`). External
libraries required: PGPLOT, CFITSIO, X11.

```sh
./autogen.sh                      # autoreconf -s -i -f
./configure FC=gfortran FCFLAGS=-Wall
make                              # builds src/xnirspec and src/zzz_memory_test
sudo make install
```

- Maximum image size is fixed at configure time: `NXMAX`, `NYMAX`, `NXYMAX`
  (default 4300; `NXYMAX` must equal `max(NXMAX, NYMAX)`, checked at runtime).
  SExtractor support needs `SEXPATH=/path/to/sex`. Other `AC_ARG_VAR`s in
  `configure.ac` (`NXMAXB9`, `NYMAXB9`, `NOVERSAMPMAX`, `NZMAX`, `NMAXBUFF`) are
  marked "do not change".
- Optional GNU readline (file name completion): `--with-readline[=DIR]`,
  `--without-readline`; Homebrew and MacPorts locations are detected
  automatically on macOS. `configure` reports whether it was enabled
  (`HAVE_READLINE` in `config.h`).
- After editing `configure.ac` or any `Makefile.am`, rerun `./autogen.sh` and
  `./configure`.
- The version string comes from `AC_INIT` in `configure.ac` (via
  `src/version.inc`); bump it there, not in the Fortran source.
- All routines except the main program are compiled into the convenience
  library `src/libxnirspec.a` (not installed), linked by `xnirspec` and by the
  tests. New source files must be added to `libxnirspec_a_SOURCES` in
  `src/Makefile.am`. The `dynamic_array*.f90` modules are listed first so their
  `.mod` files exist before the files that `USE` them.
- Automake does not track Fortran `INCLUDE` dependencies: after changing a
  generated `.inc` or `button.inc`, do `make clean && make` (or add an explicit
  dependency, as done for `version.inc`).
- Builds are routinely checked against gfortran warnings (`-Wall`,
  `-Wsurprising`, `-Wcharacter-truncation`); many commits exist only to
  silence such warnings, so new code should compile warning-free.

## Tests

```sh
make check                                # build and run all the tests
make -C tests check TESTS=test_fpoly      # run a single test
```

- Unit tests live in `tests/` (one program per group of routines, linked
  against `src/libxnirspec.a`); a test passes when it exits with status 0.
  Results are in `tests/<test>.log` and `tests/test-suite.log`.
- Checks use the module `TEST_UTILS` (`tests/test_utils.f90`):
  `CHECK_REAL(name, got, expected, rtol, atol)` (same criterion as
  `numpy.isclose`), `CHECK_INT`, `CHECK_TRUE`, `CHECK_STR`, and
  `TEST_SUMMARY()` at the end (`ERROR STOP 1` if any check failed).
- To add a test: add the program to `check_PROGRAMS` and `TESTS` in
  `tests/Makefile.am`, with its `_SOURCES` and its object depending on
  `libtestutils.a` (module ordering).
- `test_readc_file.sh` drives `readc_file_driver` with redirected stdin
  (fallback path of `READC_FILE`, without readline).
- Error paths that end in `STOP` (exit status 0) cannot be tested this way.
  Avoid inputs that make a routine wait for `<CR>` (e.g. singular matrices
  in `LUDCMP`/`LUDCMPD`).
- Only routines with explicit arguments are tested so far. Routines that use
  `IMAGEN`, `COMMON` state or PGPLOT, and the interactive program itself,
  are not covered. `zzz_memory_test` only checks that the large arrays can be
  allocated with the configured dimensions.

## Architecture

- **Main program** `src/xnirspec.f90` (~2900 lines): declares global state,
  allocates the dynamic arrays, draws the button panel and runs the event loop
  (`RPGBAND` to get a click/key, `IFBUTTON` to map it to a button number `NB`,
  then a large `IF(NB.EQ.n)` dispatch). Each feature is a subroutine in its own
  file under `src/` (one routine per file, file name = routine name in lower
  case).
- **Image buffers**: images live in `IMAGEN(NXMAX, NYMAX, NMAXBUFF)`; the
  current buffer index is `NCBUFF`. Indexing is 1-based and
  `IMAGEN(J, I, NBUFF)` with `J` = column (X, NAXIS1) and `I` = row (Y,
  NAXIS2). Spectroscopic routines assume dispersion along X and the spatial
  direction along Y (images can be rotated on load via `IROTATE`). `IMAGEN_` (`NXYMAX x NXYMAX`) is a scratch image
  (also used for box-9 mosaics); `LNULL` is the matching NaN/null mask.
- **Large arrays are dynamic, not static**: each one is a module in
  `src/dynamic_array_<name>.f90` with `Initialize_Dynamic_Array_<NAME>` and
  `Deallocate_Array_<NAME>` subroutines, plus an explicit interface in
  `src/interface_<name>.inc`. A routine that needs one does
  `USE Dynamic_Array_<NAME>` and `INCLUDE 'interface_<name>.inc'`. Follow this
  pattern for any new large array. Large local arrays that stay static should
  be declared with `SAVE` (see `fitpol.f90`).
- **Before any `STOP`** in the main program, `INCLUDE 'deallocate_arrays.inc'`
  so global dynamic arrays are released; keep that file in sync when adding
  globally allocated arrays.
- **Other global state** is shared through named `COMMON` blocks (`BLK*` in the
  main program and subroutines; `BUTTON*` in `src/button.inc`). A `COMMON`
  block must be declared identically in every routine that uses it.
- **Configuration constants** come from configure-generated include files:
  `dimensions.inc`, `largest.inc`, `sexpath.inc`, `version.inc` (templates are
  the `*.inc.in` files; edit those, never the generated `.inc`).
- **GUI layer**: `button.f90`, `buttq*.f90` / `butts*.f90` (query/set button
  properties), `ifbutton.f90`, `activebut.f90`, and `rpg*.f90` wrappers around
  PGPLOT.
- **Terminal input**: `READC`, `READI`, `READF`, ... take a question, a
  default value (`'@'` = no default) and, for `READC`, the valid characters
  (`'@'` = any). File names are read with `READC_FILE` (same interface as
  `READC`), which uses readline through `freadline.c` (C interop via
  `BIND(C)`) when stdin/stdout are terminals, and a plain `READ` otherwise.
  Wildcards (`*`, `?`) in the answer are handled by the caller, which lists
  the matching files and asks again. Use `READC_FILE` for any new file name
  prompt.
- **I/O**: `sleefits.f90` reads FITS (CFITSIO; handles extensions, 3D cubes by
  slice/collapsed range, rotation), `sescrfits.f90` writes FITS;
  `try_reduceme.f90` / `sescrredu.f90` handle the legacy REDUCEME format.
- **Spectroscopic reduction**: boundary definition (`boundary.f90`,
  `refinebl.f90`, `refineba.f90`; "BL" = polynomials y=f(x) along the
  dispersion direction, "BA" = polynomials x=g(y) along the spatial direction,
  names are historical), bivariate mapping and rectification (`mapping.f90`,
  `rectify.f90`, `unrectify.f90`), sky subtraction (`skysub.f90`), extraction
  (`xextract.f90`). Box-9 (3x3 dithered mosaic) handling is in `*b9.f90`.
- `SExtractor/` holds the default SExtractor configuration files used by
  `sextractor.f90`.

## Conventions

- Source uses free form with legacy layout: 8-space indentation, upper-case
  keywords, `IMPLICIT NONE`, F77-style declarations
  (`CHARACTER*255`, `PARAMETER (...)`). Keep this style in existing files.
- Many existing comments are in Spanish; new comments, messages and commit
  messages are in English.
- Lines in comments marked `!delete` are remnants of the F77-to-F90
  conversion.
