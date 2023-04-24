# map2db
Extract a geospatial database from a MapsForge map

map2db is a python 3 module which creates an FDO RFC 16 conforming 
sqlite geospatial database from the contents of a mapsforge format 
map created by the python 3 forgemap module with the odbl=true
option.  

This allows a map created from a derivative database of a database 
licensed under the `Open Data Commons Open Database License (ODbL) v1.0
<https://opendatacommons.org/licenses/odbl/1-0>` to satisfy the 
"Share alike" conditions of use described in section 4.6 of that license 
by simply providing a link to this tool.

If you have received a copy of a mapsforge format map with a link to 
this page, it was likely created by the python 3 forgemap module with
the odbl=true option and is a produced work created from a derivative 
database of a database licensed under ODbL v1.0.  The creator of such
a map is required under section 4.6 of that license to offer you a copy 
in machine readable form of the entire derivative database that was used 
to create the map.  This map2db python 3 module allows you to extract 
that derivative database from the map file that you have already 
recieved, thereby satisfying this requirement.  Note that section 4.6b
allows that an algorithm may have to be applied to reconstruct the 
derivative database.  map2db is the implementation of such an algorithm.
