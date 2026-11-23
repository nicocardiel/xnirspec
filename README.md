# xnirspec

This is xnirspec, a Fortran program that helps to visualise and manipulate
FITS images. The program was originally written by its author to reduce
spectroscopic data obtained with near-infrared spectrograph NIRSPEC, installed
in one of the Keck Telescopes (Hawaii). The code was initially developed
in Fortran 77, although some minor changes have been introduced (January 2023)
to compile it as Fortran 90 code (without trying to make full use of the
advantages of Fortran 90).

The program is first mentioned in Cardiel et al. (2003, ApJ, 584, 76;
see http://adsabs.harvard.edu/abs/2003ApJ...584...76C)

This software is distributed under GPL, either version 3 of the License, or
(at your option) any later version. See the file COPYING for details.

Maintainer: Nicolás Cardiel, cardiel@ucm.es  
Webpage (source): https://github.com/nicocardiel/xnirspec

## Requirements

- A Fortran compiler (e.g. `gfortran`) and a C compiler.
- GNU Autotools (`autoconf` and `automake`), to generate the `configure`
  script from the source repository.
- [PGPLOT](https://sites.astro.caltech.edu/~tjp/pgplot/), with the X11 driver.
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/).
- The X11 library.
- Optional: the GNU readline library, including its development headers
  (file name completion, see below).
- Optional: [SExtractor](https://www.astromatic.net/software/sextractor/).

The `configure` script automatically adds `/opt/local/lib` (MacPorts) on
macOS and `/usr/local/pgplot` on Linux to the library search path, when these
directories exist. Libraries installed in other locations can be indicated
with `LDFLAGS`, for example:

```sh
./configure FC=gfortran LDFLAGS="-L/path/to/pgplot -L/path/to/cfitsio/lib"
```

## Installation

Download the source code and build it:

```sh
git clone https://github.com/nicocardiel/xnirspec.git
cd xnirspec
./autogen.sh
./configure FC=gfortran FCFLAGS=-Wall
make
make check
sudo make install
```

`make check` is optional: it builds and runs the unit tests (see
[Tests](#tests)). `make install` copies the executables `xnirspec` and
`zzz_memory_test` (a check that the large image arrays can be allocated with
the configured dimensions) into `/usr/local/bin`. A different installation
directory can be selected with `--prefix`, in which case `sudo` is not
needed if the directory is writable by the user:

```sh
./configure FC=gfortran FCFLAGS=-Wall --prefix=$HOME/local
```

The name of the Fortran compiler depends on how it was installed. For
example, with MacPorts on macOS it is `gfortran-mp-<version>`:

```sh
./configure FC=gfortran-mp-14 FCFLAGS=-Wall
```

### Maximum image dimensions

By default, the maximum dimensions of the images are set to 4300 x 4300. To
allow different dimensions the code should be reinstalled defining the
parameters `NXMAX` (equivalent to `NAXIS1`), `NYMAX` (equivalent to `NAXIS2`)
and `NXYMAX` (the maximum value of `NXMAX` and `NYMAX`). For example:

```sh
./configure FC=gfortran FCFLAGS=-Wall NXMAX=6200 NYMAX=5000 NXYMAX=6200
```

The memory required by the program grows with these dimensions.

### File name completion (GNU readline)

Optionally, if the GNU readline library (including its development headers)
is available, xnirspec uses it when asking for file names: pressing TAB
completes the file name up to the longest unambiguous prefix, pressing it
again lists the possible completions, and previous entries can be recalled
with the arrow keys. If readline is not found, xnirspec is compiled without
this feature and works as usual.

On macOS, the Homebrew and MacPorts installations of GNU readline are
detected automatically. Otherwise, the location of the library can be
indicated with `--with-readline=DIR` (the files are searched in
`DIR/include` and `DIR/lib`):

```sh
./configure FC=gfortran FCFLAGS=-Wall --with-readline=$HOME/local
```

Use `--without-readline` to disable this feature. If `--with-readline` is
given explicitly and readline cannot be found, `configure` stops with an
error. The `configure` output reports whether readline has been enabled:

```
configure: file name completion with readline: yes
```

### SExtractor

To include compatibility with SExtractor it is necessary to indicate the
directory that contains its executable, which must be named `sex` (by
default, `/opt/local/bin`):

```sh
./configure FC=gfortran FCFLAGS=-Wall SEXPATH=/usr/local/bin
```

SExtractor is run in the current working directory, where the configuration
files `default.conv`, `default.nnw`, `default.param` and `default.sex` must
be present. Default versions of these files are provided in the `SExtractor`
directory of the source code.

## Tests

The unit tests are in the `tests` directory. They are built and run with:

```sh
make check
```

The results are written in `tests/test-suite.log` and in a
`tests/<test>.log` file for each test. A single test can be run with:

```sh
make -C tests check TESTS=test_fpoly
```

## Updating

After updating the source code (e.g. with `git pull`), rebuild and reinstall
the program:

```sh
./autogen.sh
./configure FC=gfortran FCFLAGS=-Wall
make clean
make
sudo make install
```

Rerunning `./autogen.sh` and `./configure` is required when `configure.ac`
or any `Makefile.am` has changed; repeat the same `configure` options used
in the previous installation. `make clean` guarantees that all the files are
recompiled after changes in the include files, whose dependencies are not
tracked automatically.

## Uninstalling

From the source directory, after running `configure` with the same options
used for the installation:

```sh
sudo make uninstall
```
