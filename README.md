# planetMagFields
Routines to plot magnetic fields of planets in our solar system.
Requires [cartopy](https://scitools.org.uk/cartopy/docs/latest/).

## Use
#### For a single planet

```
$ ./magField <planet name> <radius>
```

The planet name is case insensitive. The radius of the surface is 1 and is the default.
Examples:

```
$ ./magField earth
$ ./magField jupiter 0.75
```
#### All planets

```
$ ./magField all <radius>
```

would produce a plot of all magnetic field maps of different planets in a single figure
along with a table of information about dipole tilt for each.

## References

Will add soon
