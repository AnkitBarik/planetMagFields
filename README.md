# planetMagFields

Routines to plot magnetic fields of planets in our solar system.
Requires the [cartopy](https://scitools.org.uk/cartopy/docs/latest/) library.

## Usage:

```
$ ./magField <planet> <radius>`
```

This will plot the radial magnetic field of a planet (any of the names from the list
below, case insensitive) at a radius given in terms of the surface radius. The default
is the surface field. For example,

```
$ ./magField earth
```

displays some information about Earth's field and produces the surface field of Earth 

```
Radius not specified, using surface

Planet: Earth
Depth (fraction of surface radius) = 1.00
l_max = 13
Dipole tilt (degrees) = -9.410531
```

![Earth's field](/images/br_earth.png)

while

```
$ ./magField earth 0.55
```

produces the field at the core-mantle boundary (CMB)

```
Planet: Earth
Depth (fraction of surface radius) = 0.55
l_max = 13
Dipole tilt (degrees) = -9.410531
```

![Earth's CMB field](/images/br_earth_cmb.png)

```
$ ./magField all <radius>
```

would produce a plot of all magnetic field maps of different planets in a single figure
along with a table of information about dipole tilt for each.

# Data sources

Mercury : [Anderson et. al. 2012](https://doi.org/10.1029/2012JE004159)

Earth   : [IGRF 13](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)

Jupiter : [JRM09, Connerny et. al. 2018](https://doi.org/10.1002/2018GL077312)

Saturn  : [Cassini 11+, Cao et. al. 2020](https://doi.org/10.1016/j.icarus.2019.113541)

Uranus  : [Connerny et. al. 1987](https://doi.org/10.1029/JA092iA13p15329)

Neptune : [Connerny et. al. 1991](https://doi.org/10.1029/91JA01165)

Ganymede: [Kivelson et. al. 2002](https://doi.org/10.1006/icar.2002.6834)
