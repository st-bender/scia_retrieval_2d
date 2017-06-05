# Level 1c Spectra Binary Format

This is the layout of the `.l_mpl_binary` level 1c spectra files
used as input for `scia_retrieval_2d`.
The purpose is to save disk space compared to plain ascii spectra files
and the format is basically using the binary representation of the numbers.
This format is currently not very portable and relies on
float and int being 32 bits (4 bytes) and endianness is not considered.
However, these limitations should not be important
as long as the binary files are prepared on the same machine on which
the retrieval code is run on.

Loading is handled in [`LimbNadir_IO.cpp`](./LimbNadir_IO.cpp) by:
```c++
int Load_Limb_l_mpl_binary()
```

 address | size | description
--------:|------|-------------
  0x0000 | (31 x 100) x char | text header containing metadata (number of header lines x 100 chars)<br />(Looks for a hash `#` at the beginning of each 100 characters read and adjusts the number of header lines if the first character differs.)
  0x0c1c |   int      | number of tangent points `nalt`
  0x0c20 |   int      | number of spectral points `npix`
  0x0c24 |  5 x int   | orbit state tuple: `(orbit, state_in_orbit, state_id, nr_of_profiles, actual_profile)`
  0x0c38 |  6 x int   | date tuple: `(year, month, day, hour, minute, second)`
  0x0c50 | 10 x float | centre and corner latitudes and longitudes:<br />`(centre_lat, centre_lon, corner0_lat, corner0_lon, ..., corner3_lat, corner3_lon)`
  0x0c78 |   float    | fractional phase of the orbit [0, 1) `orbit_phase`
  0x0c7c | `npix` x float | spectral wavelengths
  0x0c7c + `npix` x `sizeof(float)` | `nalt` x `limb_data` | tangent point limb data

`limb_data`:

 number | size  | description
-------:|-------|-----------------------------
   0    | float | sub satellite point latitude
   1    | float | sub satellite point longitude
   2    | float | tangent point latitude
   3    | float | tangent point longitude
   4    | float | tangent point altitude
   5    | float | tangent point solar zenith angle
   6    | float | tangent point solar azimuth angle
   7    | float | tangent point line-of-sight angle
   8    | float | top-of-atmosphere solar zenith angle
   9    | float | top-of-atmosphere solar azimuth angle
  10    | float | top-of-atmosphere line-of-sight angle
  11    | float | satellite solar zenith angle
  12    | float | satellite solar azimuth angle
  13    | float | satellite line-of-sight angle
  14    | float | satellite altitude
  15    | float | Earth radius
  16    | `npix` x float | spectral radiance values
  17    | `npix` x float | spectral radiance relative uncertainties


