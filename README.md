# SCIAMACHY 2d satellite trace gas retrieval

Two-dimensional trace gas retrieval from
[SCIAMACHY](http://www.sciamachy.org) limb scans as described in
[Scharringhausen et al., 2008](http://www.atmos-chem-phys.net/8/1963/2008/)
([doi:10.5194/acp-8-1963-2008](http://dx.doi.org/10.5194/acp-8-1963-2008)),
[Langowski et al., 2014](http://www.atmos-meas-tech.net/7/29/2014/)
([doi:10.5194/amt-7-29-2014](http://dx.doi.org/10.5194/amt-7-29-2014)),
[Bender et al., 2013](http://www.atmos-meas-tech.net/6/2521/2013/)
([doi:10.5194/amt-6-2521-2013](http://dx.doi.org/10.5194/amt-6-2521-2013)),
[Bender et al., 2017](http://www.atmos-meas-tech.net/10/209/2017/)
([doi:10.5194/amt-10-209-2017](http://dx.doi.org/10.5194/amt-10-209-2017)).

For a more detailed history and more references,
see the [history.md](./history.md).

This version currently supports retrieving the following species:

- Mg
- Mg+
- Fe
- NO

## Prerequisites

Configuration and pre-build setup is best managed via [CMake](https://cmake.org)
version 2.6 or higher.
Other pieces needed to successfully compile the program are:

### required
- gcc and g++ version 4.9 or higher
- a linear algebra system, for example [OpenBLAS](http://www.openblas.net),
  [ATLAS](http://math-atlas.sourceforge.net),
  [Intel MKL](https://software.intel.com/en-us/mkl),
  [sunperf](https://docs.oracle.com/cd/E60778_01/html/E61993/sunperf-3p.html)
  on OpenSolaris/OpenIndiana or similar
- the netcdf library and headers (www.unidata.ucar.edu)
- the hdf5 library and headers (www.hdfgroup.org)

### optional
- the Nitric Oxide Empirical Model (NOEM), see
  [Marsh et al., 2004](http://onlinelibrary.wiley.com/doi/10.1029/2003JA010199/full)
  ([doi:10.1029/2003JA010199](http://dx.doi.org/10.1029/2003JA010199)),
  for calculating the prior densities.
  The code is not (yet) freely distributable but available on request from the author.

## Compile

Out of tree builds are supported and recommended, the simplest build steps are:

```sh
$ mkdir build
$ cd build
$ cmake [OPTIONS] ..
$ make [-j <n>]
```

This should produce an executable called `scia_retrieval` in the (current)
`build/` subdirectory.

## Install

There is no automatic install available. Simply copy the executable to a place
where it is accessible, for example `/usr/local/bin` or `~/bin`.

## Usage

### Preparation

Running the retrieval needs SCIAMACHY geo-located atmospheric level 1b spectra
(SCI\_NL\_\_1P) (for example version 8.02) provided by ESA via the
[ESA data browser](https://earth.esa.int/web/guest/data-access/browse-data-products).

The spectra need to be calibrated first which requires two steps:

1. Producing a `.child` or HDF5 (`.h5`) file;
   the first can be created for example with the
   [SciaL1C](https://earth.esa.int/web/guest/software-tools/content/-/article/scial1c-command-line-tool-4073)
   command line tool or with the free software
   [nadc\_tools](https://github.com/rmvanhees/nadc_tools).
   The latter can also output to HDF5 (`.h5`).

2. Extracting the limb spectra from the `.child` or HDF5 files to ASCII
   (or a special binary format, see [mpl\_binary.md](./mpl_binary.md)).
   For example using the `sciapy` python tools which will be available publicly soon.

#### Orbitlist

After calibration, the orbitlist file contains the pathnames of all
limb spectra belonging to one orbit. Those are typically 20--30 per orbit,
depending on the chosen category (or categories) in the steps above.
For example `/home/user/SCIA/Spectra/orbitlists/MLT/2008/orbitlist-sb34343_scia1.dat`
may contain the following:

    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070223_1_0_34343.dat
    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070319_1_0_34343.dat
    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070415_1_0_34343.dat
    ...

or using the binary format:

    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070223_1_0_34343.dat.l_mpl_binary
    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070319_1_0_34343.dat.l_mpl_binary
    /home/user/SCIA/Spectra/MLT/channel_1_v8.02/2008/20080924/SCIA_limb_20080924_070415_1_0_34343.dat.l_mpl_binary
    ...

These lists can be easily created by, for example:
```sh
$ find /home/user/SCIA/Spectra/MLT/channel_1_v8.02/<year>/<date>/SCIA_limb_*_<orbit>.dat > /home/user/SCIA/Spectra/orbitlists/MLT/<year>/orbitlist-sb<orbit>_scia1.dat
```
or:
```sh
$ find /home/user/SCIA/Spectra/MLT/channel_1_v8.02/<year>/<date>/SCIA_limb_*_<orbit>.dat.l_mpl_binary > /home/user/SCIA/Spectra/orbitlists/MLT/<year>/orbitlist-sb<orbit>_scia1.dat
```
with `<year>`, `<date>`, and `<orbit>` set appropriately.

### Command line arguments

The executable `scia_retrieval` is invoked as follows:

```sh
$ scia_retrieval <orbitlist> <output_path> <sol_scia_orbit> <sol_ref> <script1> <script2> <config>
```

- orbitlist: _filename_

  File containing the orbital spectra files one per line ([see above](#orbitlist)).<br />
  Example: `/home/user/SCIA/Spectra/orbitlists/MLT/2008/orbitlist-sb34343_scia1.dat`

- output\_path: _pathname_

  Path to save the result files (number densities, averaging kernels, etc.) to.<br />
  Example: `/home/user/SCIA/Retrieval/Species_2008`

- sol\_scia\_orbit: _filename_ or `auto`

  Filename of the orbit measured solar spectrum. When set to `auto`,
  the program tries to determine the solar spectrum location from the limb spectra
  file names (given by the entries in the orbitlist). It then assumes that the 
  solar spectrum is located in the same directory using the glob string
  `SCIA_solar_*_<orbitnumber>.dat`.<br />
  Example: `/home/user/SCIA/Spectra/2008/SCIA_solar_20080924_065146_D0_34343.dat`

- sol\_ref: _filename_

  The filename of the highly resolved solar reference spectrum, for example the
  [SAO2010 solar reference spectrum](http://www.sciencedirect.com/science/article/pii/S0022407310000610)
  ([doi:10.1016/j.jqsrt.2010.01.036](http://dx.doi.org/10.1016/j.jqsrt.2010.01.036)),
  available from
  https://www.cfa.harvard.edu/atmosphere/links/sao2010.solref.converted<br />
  Example: `/home/user/SCIA/Spectra/SolRef/sao2010.solref.dat`

- script1: _filename_

  A script for combining multiple postscript or pdf files into one large
  pdf document. A simple approach is using ghostscript, see the `multips2pdf`
  script in the source tree.<br />
  Example: `/home/user/SCIA/run/multips2pdf`

- script2: _filename_

  A script for combining multiple postscript or pdf files into one large
  postscript document. A simple approach is using ghostscript, see the `multips2ps`
  script in the source tree.<br />
  Example: `/home/user/SCIA/run/multips2ps`

- config: _filename_

  Detailed retrieval configuration, see the `SCIA2D_mlt.conf` and
  `SCIA2D_nom.conf` files in the source tree.<br />
  Example: `/home/user/SCIA/run/SCIA2D_mlt.conf`

## Licenses

This program is free software: you can redistribute it or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 2 (GPLv2), see [local copy](./COPYING.GPL2)
or [online version](http://www.gnu.org/licenses/gpl-2.0.html).

**Exemptions** from the GPLv2 are the following parts:

- The NRLMSISE-00 source files `nrlmsise-00.c`, `nrlmsise-00_data.c`,
  `nrlmsise-00_doc.txt`, and `nrlmsise-00.h` are in the public domain,
  see [COPYING.NRLMSISE-00](./COPYING.NRLMSISE-00).
- The [gzstream](http://www.cs.unc.edu/Research/compgeom/gzstream)
  source files `gzstream.cpp` and `gzstream.h` are subject
  to the GNU Lesser General Public License, version 2.1,
  see [local copy](./COPYING.gzstream) or
  [online version](http://www.gnu.org/licenses/lgpl-2.1.html).

## Acknowledgements

The initial FORTRAN version was developed at the University of Bremen, Germany,
supported by the Air Force Office of Scientific Research (AFOSR),
Air Force Material Command, USAF, under grant number FA8655-03-1-3035.
The transition to C++ was developed at the University of Bremen, Germany, also
supported by the AFOSR and the European Office of Aerospace Research and Development (EOARD),
under grant number FA8655-09-3012.
The latest version of the code was developed at the Karlsruhe Institute of Technology (KIT),
Karlsruhe, Germany, supported by the Helmholtz-society under the grant number VH-NG-624.
