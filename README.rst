=======
 map2db
=======

Extract a geospatial database from a MapsForge map

map2db is a python 3 module which creates an FDO RFC 16 conforming 
sqlite geospatial database from the contents of a mapsforge format 
map created by forgemap with the dbl option.  

This allows a map created from a derivative database of a database 
licensed under the `Open Data Commons Open Database License (ODbL) v1.0
<https://opendatacommons.org/licenses/odbl/1-0>` to satisfy the 
"Share alike" conditions of use described in section 4.6 of that license 
by simply providing a link to this tool.

If you have received a copy of a mapsforge format map with a link to 
this page, It may be have been created from a derivative database of a 
database licensed under ODbL v1.0.  If so, the creator of such a map is 
required to offer you a copy in machine readable form of the entire 
derivative database that was used to create the map.  map2db allows you 
to extract that derivative database from the map file that you already 
have, thereby satisfying this requirement.  Note that section 4.6b
allows that an algorithm may have to be applied to reconstruct the 
derivative database.  map2db is the implementation of such an algorithm.

The license statement encoded in the map file for the recreated database
is both written to the README table of the database and printed.  The 
creator of the map file is responsible for ensuring that this statement
defines the license under which the database is released, or the lack 
thereof, and that it satisfies any licencing requirements imposed by the
data used to create the database.

A TOML configuration file suitable for building a new mapsforge map with
the dbl option from the recreated database using forgemap is also
produced.

If map2db is applied to any mapsforge map that was created by any program
other than forgemap or which was created by forgemap, but without the dbl
option, then a database can usually still be created from it.  However,
this database will contain the raw content of every tile in the map
rather than a properly reconstructed database with various parts of each
feature merged together and more simplified versions of such features
discarded in favor of the less simplified versions.  The dbl option of
forgemap also explicitly encodes whether the source geometry was a line
or area feature.  Since this information is not available, those way 
features in which every coordinate block forms a closed ring with at least 
four points are assumed to be areas, while the remainder are assumed to be
lines.  This may misclassify some features as areas that should have
been lines.  Since not all closed rings of at least four points are
valid, this may produce unexpected results.  When used with such a map,
a brief message explaining this is both written to the README table of
the database and printed in place of a license statement.  Since no
database license statement is encoded in such a map file, the user is
responsible for determining what legal restrictions apply to use of the
recreated database.

map2db requires shapely version >= 2.0

Usage:

 >>> from map2db import map2db
 >>> map2db(src_pathname, dst_pathname)

or

 $ python3 map2db.py src_pathname dst_pathname

From the command line, if dst_pathname is not provided, then one is
derived from src_pathname.
