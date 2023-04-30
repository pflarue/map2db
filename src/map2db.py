#!/usr/bin/python3
"""
Recreates the geospatial database used to build a mapsforge map.

Recreates the FDO RFC 16 <https://trac.osgeo.org/fdo/wiki/FDORfc16>
conforming geospatial database that was used by forgemap to create a
mapsforge map using the "dbl" option.

This allows anyone possessing such a map file to also have access to the
database from which it was created.  This satisfies the "share alike"
conditions of the Open Data Commons Open Database License (ODbL) v1.0
<https://opendatacommons.org/licenses/odbl/1-0> if the database is
derived from a database covered by this license.  It can also be used
to recreate a database not covered by that license.

The license statement encoded in the map file for the recreated database
is both written to the README table of the database and printed.  The 
creator of the map file is responsible for ensuring that this statement
defines the license under which the database is released, or lack 
thereof, and that it satisfies any licencing requirements imposed by the
data used to create the database.

A TOML configuration file suitable for building a new mapsforge map with
the "dbl" option from the recreated database using forgemap is also
produced.

If applied to any mapsforge map that was created by any program other
than forgemap or which was created by forgemap, but without the "dbl"
option, then a database can usually still be created from it.  However,
this database will contain the raw content of every tile in the map
rather than a properly reconstructed database with various parts of each
feature merged together and more simplified versions of such features
discarded in favor of the less simplified versions.  The dbl option of
forgemap also explicitly encodes whether the source geometry was a line
or area feature.  Since this information is not available, way features
in which every coordinate block forms a closed ring with at least four
points are assumed to be areas, while the remainder are assumed to be
lines.  This may misclassify some features as areas that should have
been lines.  Since not all closed rings of at least four points are
valid, this may produce unexpected results.  When used with such a map,
a brief message explaining this is both written to the README table of
the database and printed in place of a license statement.  Since no
database license statement is encoded such a map file, the user is
responsible for determining what legal restrictions apply to use of the
recreated database."

map2db requires shapely version >= 2.0

Usage:

>>> from map2db import map2db

>>> map2db(src_pathname, dst_pathname)

or

$ python3 map2db.py src_pathname dst_pathname

From the command line, if dst_pathname is not provided, then one is
derived from src_pathname.
"""

import json
import math
import os
import sqlite3
from struct import unpack
import sys

try:
    from shapely import __version__ as sversion
    from shapely import geometry as sg
    from shapely import wkb as swkb
    from shapely import validation as sv
    from shapely.ops import linemerge
except ModuleNotFoundError:
    print(
        "\nmap2db requires shapely (version >= 2.0).\n"
        "See https://pypi.org/project/Shapely for information about"
        " installing it.\n"
    )
    raise
if sversion.startswith("1."):
    raise ModuleNotFoundError(
        f"map2db requires shapely version >= 2.0, but '{sversion}' is"
        f" installed."
    )

# Assumptions about the content of the input mapsource map:
#
# 1.  All features will have a __dbl_pnum, __dbl_lnum, or __dbl_anum
# tag.  These will be converted to m2db_pnum, m2db_lnum, and
# m2db_anum respectively.
#
# 2.  If a feature with a given m2db_pnum, m2db_lnum, or m2db_anum is
# present at more than one level, then the feature in the lower numbered
# level is a duplicate of the feature in the higher numbered level, but
# probably with more aggressive simplification of the geometry applied.
# Thus,the feature in the lower numbered level shall be discarded,
# except that it shall be used to update the m2db_minz value for the
# feature.  All attributes except m2db_geometry, m2db_level, m2db_minz,
# and m2db_maxz should be identical.
#
# 3.  If a feature with a given m2db_pnum, m2db_lnum, or m2db_anum is
# present in more than one tile at the same level, then all attributes
# except for m2db_geometry will be identical.  The final m2db_geometry
# of the feature shall be the union of all of the geometries from the
# individual tiles.
#
# 4.  A feature with a given m2db_pnum, m2db_lnum, or m2db_anum should
# have a single continuous zoom range.  Thus, if m2db_minz for a feature
# is not equal to the minzoom of the subfile where it is found, then
# that feature should not be found in tiles for lower level subfiles.
#
# Mapsforge maps created by forgemap with the dbl option will satisfy
# all of these assumptions.  Maps that do not satisfy these assumptions
# are only partially supported.

# See "docs/Specification-Binary-Map-File.md" at
# <https://www.github.com/mapsforge/mapsforge> for details of the
# mapsforge file format.


# If True, features from non-dbl maps will still be trimmed to tile
# boundaries, else they will not.
TRIM_NONUM = True


def _prepare_dbc(dbc):
    # Create metadata tables following FDO RFC 16
    # <https://trac.osgeo.org/fdo/wiki/FDORfc16>
    dbc.execute(
        """
        CREATE TABLE geometry_columns (
            f_table_name TEXT,
            f_geometry_column TEXT,
            geometry_type INTEGER,
            coord_dimension INTEGER,
            srid INTEGER,
            geometry_format TEXT);"""
    )
    dbc.execute(
        """
        CREATE TABLE spatial_ref_sys (
            srid INTEGER UNIQUE,
            auth_name TEXT,
            auth_srid INTEGER,
            srtext TEXT);"""
    )
    # populate spatial_ref_sys with EPSG:4326 as srid=1
    dbc.execute(
        """
        INSERT INTO spatial_ref_sys (srid, auth_name, auth_srid, srtext)
            VALUES (1, 'EPSG', 4326, ?);""",
        (
            "GEOGCS["
            '"WGS 84",'
            "DATUM["
            '"WGS_1984",'
            "SPHEROID["
            '"WGS 84",'
            "6378137,"
            "298.257223563,"
            'AUTHORITY["EPSG","7030"]],'
            'AUTHORITY["EPSG","6326"]],'
            'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
            'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
            'AXIS["Latitude",NORTH],'
            'AXIS["Longitude",EAST],'
            'AUTHORITY["EPSG","4326"]]',
        ),
    )

    # Create the data and custom metaddata tables
    dbc.execute("""CREATE TABLE README (desc TEXT, text TEXT);""")
    # Always create points, lines, and areas tables, even if some may
    # remain empty.  All geometry columns will be WKB (as opposed to
    # WKT).  Additional columns will be added as needed.  Columns are
    # arranged such that for a given m2db_pnum, m2db_lnum, or m2db_anum,
    # all columns after m2db_geometry should be identical for additional
    # appearances of this feature in other tiles.
    #
    # column names are chosen to reduce the likelihood of matching a tag
    # key name for which a new column should be added.
    dbc.execute(
        """
        CREATE TABLE points (
            m2db_pnum INTEGER PRIMARY KEY,
            m2db_level INTEGER,
            m2db_minz INTEGER,
            m2db_maxz INTEGER,
            m2db_geometry BLOB,
            m2db_layer_num INTEGER,
            m2db_tags TEXT);"""
    )
    dbc.execute(
        """
        INSERT INTO geometry_columns (
                f_table_name, f_geometry_column, geometry_type,
                coord_dimension, srid, geometry_format
            )
            VALUES ("points", "m2db_geometry", 4, 2, 1, "WKB");"""
    )
    dbc.execute(
        """
        CREATE TABLE lines (
            m2db_lnum INTEGER PRIMARY KEY,
            m2db_level INTEGER,
            m2db_minz INTEGER,
            m2db_maxz INTEGER,
            m2db_geometry BLOB,
            m2db_layer_num INTEGER,
            m2db_tags TEXT);"""
    )
    dbc.execute(
        """
        INSERT INTO geometry_columns (
                f_table_name, f_geometry_column, geometry_type,
                coord_dimension, srid, geometry_format
            )
            VALUES ("lines", "m2db_geometry", 5, 2, 1, "WKB");"""
    )
    dbc.execute(
        """
        CREATE TABLE areas (
            m2db_anum INTEGER PRIMARY KEY,
            m2db_level INTEGER,
            m2db_minz INTEGER,
            m2db_maxz INTEGER,
            m2db_geometry BLOB,
            m2db_layer_num INTEGER,
            m2db_tags TEXT);"""
    )
    dbc.execute(
        """
        INSERT INTO geometry_columns (
                f_table_name, f_geometry_column, geometry_type,
                coord_dimension, srid, geometry_format
            )
            VALUES ("areas", "m2db_geometry", 6, 2, 1, "WKB");"""
    )
    dbc.execute(
        """
        CREATE TABLE subfiles (
            subfile_num INTEGER PRIMARY KEY,
            level INTEGER,
            minzoom INTEGER,
            maxzoom INTEGER);"""
    )
    dbc.execute(
        """
        CREATE TABLE limits (
            key TEXT,
            value REAL);"""
    )


# Parse values of various types from the source file.
def _parse_u16(srcf):
    return unpack(">H", srcf.read(2))[0]


def _parse_s16(srcf):
    return unpack(">h", srcf.read(2))[0]


def _parse_u32(srcf):
    return unpack(">I", srcf.read(4))[0]


def _parse_s32(srcf):
    return unpack(">i", srcf.read(4))[0]


def _parse_u64(srcf):
    return unpack(">Q", srcf.read(8))[0]


def _parse_f(srcf):
    return unpack(">f", srcf.read(4))[0]


# Parse a variable length Unsigned Little Endian Base 128 (ULEB128)
# value of up to 8 bytes in length.
def _parse_vu128(srcf):
    output = 0
    for i in range(8):
        bval = ord(srcf.read(1))
        output += (bval & 0x7F) * 2 ** (7 * i)
        if not bval & 0x80:
            return output
    raise ValueError(f"unable to decode MF VBU-U at offset: {srcf.tell() - 8}")


# Parse a variable length Signed Little Endian Base 128 (SLEB128) value
# of up to 8 bytes in length.
def _parse_vs128(srcf):
    bval = ord(srcf.read(1))
    if not bval & 0x80:  # check to see if first byte is also last byte
        output = bval & 0x3F
        sign = -1 if bval & 0x40 else 1
        return sign * output
    output = bval & 0x7F  # first (where first byte is not last byte)

    for i in range(1, 8):
        bval = ord(srcf.read(1))
        if not bval & 0x80:  # last byte
            output += (bval & 0x3F) * 2 ** (7 * i)
            sign = -1 if (bval & 0x40) else 1
            return sign * output
        # continue to next byte
        output += (bval & 0x7F) * 2 ** (7 * i)
    raise ValueError(f"unable to decode MF VBU-S at offset: {srcf.tell() - 8}")


def _parse_lstr(srcf):
    # parse a variable length string
    return str(srcf.read(_parse_vu128(srcf)), encoding="utf-8")


# Conversions between EPSG:4326 lat/lon and tile coordinates.
def _xfromlon(tile_z, lon):
    return 2 ** (tile_z - 1) * (lon / 180.0 + 1.0)


def _yfromlat(tile_z, lat):
    return (2 ** (tile_z - 1)) * (
        2
        - (math.log(math.tan((0.25 + lat / 360.0) * math.pi)) / math.pi + 1.0)
    )


def _lonfromx(tile_z, tile_x):
    return (tile_x / (2**tile_z) * 2.0 - 1.0) * 180.0


def _latfromy(tile_z, tile_y):
    return (
        math.atan(
            math.exp(
                (((1 << tile_z) - tile_y) / (2 ** (tile_z - 1)) - 1.0)
                * math.pi
            )
        )
        / math.pi
        - 0.25
    ) * 360.0


def _parse_fileheader(srcf):
    srcf.seek(0)
    output = {}
    output["magic_byte"] = srcf.read(20)
    if output["magic_byte"] != b"mapsforge binary OSM":
        raise ValueError("Not a valid map file.  Wrong header bytes")
    srcf.read(24)  # skip past some unneeded data
    # bounding box decimal degrees
    output["minlat"] = round(_parse_s32(srcf) / 1000000.0, 6)
    output["minlon"] = round(_parse_s32(srcf) / 1000000.0, 6)
    output["maxlat"] = round(_parse_s32(srcf) / 1000000.0, 6)
    output["maxlon"] = round(_parse_s32(srcf) / 1000000.0, 6)
    _ = _parse_u16(srcf)  # tile width and height in px
    _ = _parse_lstr(srcf)  # projection
    flags = ord(srcf.read(1))
    output["debuginfo"] = bool(flags & 0x80)
    if flags & 0x40:  # startPosition
        output["startlat"] = round(_parse_s32(srcf) / 1000000.0, 6)
        output["startlon"] = round(_parse_s32(srcf) / 1000000.0, 6)
    if flags & 0x20:  # startZoom
        output["startzoom"] = ord(srcf.read(1))
    if flags & 0x10:  # langPref
        output["lang"] = _parse_lstr(srcf)
    if flags & 0x08:  # comment
        output["commstr"] = _parse_lstr(srcf)
    if flags & 0x04:  # createdBy
        output["createdby"] = _parse_lstr(srcf)
    # POI Tags
    ptags_count = _parse_u16(srcf)
    output["ptags"] = [_parse_lstr(srcf) for i in range(ptags_count)]
    # Way Tags
    wtags_count = _parse_u16(srcf)
    output["wtags"] = [_parse_lstr(srcf) for i in range(wtags_count)]
    # Zoom Intervals (subFiles)
    output["subfiles_count"] = ord(srcf.read(1))
    output["subfiles"] = []
    for i in range(output["subfiles_count"]):
        subfile = {}
        subfile["level"] = ord(srcf.read(1))
        subfile["minzoom"] = ord(srcf.read(1))
        subfile["maxzoom"] = ord(srcf.read(1))
        # absolute offset to start of subfile (bytes)
        subfile["offset"] = _parse_u64(srcf)
        subfile["length"] = _parse_u64(srcf)  # length of subfile (bytes)
        # calculate (integer) range of X and Y for Z=level
        subfile["minx"] = int(_xfromlon(subfile["level"], output["minlon"]))
        subfile["maxx"] = int(_xfromlon(subfile["level"], output["maxlon"]))
        subfile["miny"] = int(_yfromlat(subfile["level"], output["maxlat"]))
        subfile["maxy"] = int(_yfromlat(subfile["level"], output["minlat"]))
        subfile["y_count"] = subfile["maxy"] - subfile["miny"] + 1
        subfile["x_count"] = subfile["maxx"] - subfile["minx"] + 1
        subfile["tile_count"] = subfile["y_count"] * subfile["x_count"]
        output["subfiles"].append(subfile)
    output["levels"] = [subfile["level"] for subfile in output["subfiles"]]
    return output


# tile indicies are stored in a flat array for each subfile.
def _get_tilenum(tile_x, tile_y, subfile):
    if (
        tile_x < subfile["minx"]
        or tile_x > subfile["maxx"]
        or tile_y < subfile["miny"]
        or tile_y > subfile["maxy"]
    ):
        raise ValueError(
            f"Invalid tile coordinates: ({tile_x}, {tile_y},"
            f" {subfile['level']})"
        )
    return (tile_y - subfile["miny"]) * subfile["x_count"] + (
        tile_x - subfile["minx"]
    )


def _parse_tags(srcf, count, taglist):
    # Tags encoded in a mapsforge file may be fixed strings or may only
    # indicate the key and the type of the value.  In the latter case,
    # the value is encoded independently for each feature that uses such
    # a tag.  This is normally used where the value is unique (or nearly
    # so) to each feature.  Return a list of the fixed string tags, but
    # return the variable tags separately as a dict.  The fixed tags
    # will be stored in the database as a list, while a data column for
    # each variable tag will be added to the database table so that the
    # value can be explicitly set.
    tags = [taglist[_parse_vu128(srcf)] for i in range(count)]
    vtags = {}
    for tag in tags[:]:
        if tag.endswith("=%i"):
            if ":colour" in tag:
                vtags[tag[:-3]] = hex(_parse_s32(srcf))[2:]
            else:
                vtags[tag[:-3]] = str(_parse_s32(srcf))
        elif tag.endswith("=%f"):
            vtags[tag[:-3]] = str(_parse_f(srcf))
        elif tag.endswith("=%h"):
            vtags[tag[:-3]] = str(_parse_s16(srcf))
        elif tag.endswith("=%s"):
            vtags[tag[:-3]] = _parse_lstr(srcf)
        else:
            continue
        tags.remove(tag)
    return tags, vtags


def _parse_points(srcf, endoffset, origin_md, tagslist, debuginfo, poi_counts):
    output = []
    tile_z = 0
    for tile_z, count in enumerate(poi_counts):
        for _ in range(count):
            if debuginfo:
                if not srcf.read(32).startswith(b"***POIStart"):
                    raise ValueError(
                        f"POI Start (***POIStart) not found at offset:"
                        f" {hex(srcf.tell() - 32)}"
                    )
            point = {"m2db_minz": tile_z}
            point["lat"] = round(
                (origin_md["lat"] + _parse_vs128(srcf)) / 1000000.0, 6
            )
            point["lon"] = round(
                (origin_md["lon"] + _parse_vs128(srcf)) / 1000000.0, 6
            )
            layer_num = ord(srcf.read(1))
            tag_count = layer_num & 0x0F
            point["m2db_layer_num"] = layer_num // 16 - 5
            point["m2db_tags"], vtags = _parse_tags(srcf, tag_count, tagslist)
            if "__dbl_pnum" in vtags:
                point["m2db_pnum"] = int(vtags["__dbl_pnum"])
                del vtags["__dbl_pnum"]
            else:
                point["nonum"] = True
            flags = ord(srcf.read(1))
            # The following unnamed values are given keys that are
            # unlikely to conflict with the keys of named values in the
            # vtags dict or the output sqlite database.
            if flags & 0x80:  # POI name
                vtags["m2db_name"] = _parse_lstr(srcf)
            if flags & 0x40:  # house number
                vtags["m2db_house_num"] = _parse_lstr(srcf)
            if flags & 0x20:  # elevation in meters
                vtags["m2db_elevation"] = _parse_vu128(srcf)
            point["vtags"] = vtags
            output.append(point)
            if srcf.tell() > endoffset:
                raise ValueError(
                    f"Error reading points {hex(srcf.tell())} >"
                    f" {hex(endoffset)}."
                )
    if srcf.tell() != endoffset:
        raise ValueError(
            f"Error reading points {hex(srcf.tell())} !=" f" {hex(endoffset)}."
        )
    return output


def _parse_ways(srcf, endoffset, origin_md, tagslist, debuginfo, way_counts):
    output = []
    tile_z = 0
    for tile_z, count in enumerate(way_counts):
        for way_num in range(count):
            # For a non-dbl mapsforge file, if any coordinate block has
            # less than 4 nodes or is not closed, then rings will be
            # changed to False.  This will be used to determine whether
            # to put the feature into the lines table as a
            # MultiLineString or into the areas table as a MultiPolygon.
            # rings is not used when parsing ways from dbl mapsforge
            # files since these features are all uniquely identified as
            # line or area features.
            rings = True
            if debuginfo:
                if not srcf.read(32).startswith(b"---WayStartX"):
                    raise ValueError(
                        f"Way Start (---WayStartX) not found at offset:"
                        f" {hex(srcf.tell() - 32)}"
                    )
            way = {"m2db_minz": tile_z}
            # The difference between way_endoffset and the value of
            # srcf.tell() indicates the amount of data used to encode
            # this way.
            way_endoffset = _parse_vu128(srcf) + srcf.tell()
            if way_endoffset > endoffset:
                raise ValueError("Error reading ways")

            _parse_u16(srcf)  # subtile_bitmap (not needed here)
            layer_num = ord(srcf.read(1))
            tags_count = layer_num & 0x0F
            way["m2db_layer_num"] = layer_num // 16 - 5
            way["m2db_tags"], vtags = _parse_tags(srcf, tags_count, tagslist)
            if "__dbl_lnum" in vtags:
                way["m2db_lnum"] = int(vtags["__dbl_lnum"])
                del vtags["__dbl_lnum"]
            elif "__dbl_anum" in vtags:
                way["m2db_anum"] = int(vtags["__dbl_anum"])
                del vtags["__dbl_anum"]
            else:
                way["nonum"] = True

            flags = ord(srcf.read(1))
            # The following unnamed values are given keys that are
            # unlikely to conflict with the keys of named values in the
            # vtags dict or the output sqlite database.
            if flags & 0x80:  # Way name
                vtags["m2db_name"] = _parse_lstr(srcf)
            if flags & 0x40:  # house number
                vtags["m2db_house_num"] = _parse_lstr(srcf)
            if flags & 0x20:  # reference
                vtags["m2db_reference"] = _parse_lstr(srcf)
            if flags & 0x10:
                # label position as offsets from first way coordinates.
                label_coords = {
                    "dlat": round(_parse_vs128(srcf) / 1000000.0, 6),
                    "dlon": round(_parse_vs128(srcf) / 1000000.0, 6),
                }
            else:
                label_coords = None
            way["vtags"] = vtags
            if flags & 0x08:  # number of way data blocks.
                dblock_count = _parse_vu128(srcf)
            else:
                dblock_count = 1  # default is one block.
            # Uses double delta encoding rather than single delta.
            ddencode = bool(flags & 0x04)
            way["coordinates"] = []
            for _ in range(dblock_count):  # for each data block
                dblock = []
                cblock_count = _parse_vu128(srcf)
                for _ in range(cblock_count):
                    lat_prev_md = origin_md["lat"]
                    lon_prev_md = origin_md["lon"]
                    dlat_md = 0
                    dlon_md = 0
                    cblock = []
                    node_count = _parse_vu128(srcf)
                    # Note: with double delta encoding, the double delta
                    # is only implemented starting on the third point
                    for node_num in range(node_count):
                        if ddencode and node_num > 1:
                            dlat_md += _parse_vs128(srcf)
                            dlon_md += _parse_vs128(srcf)
                        else:
                            dlat_md = _parse_vs128(srcf)
                            dlon_md = _parse_vs128(srcf)
                        cblock.append(
                            [
                                round((lon_prev_md + dlon_md) / 1000000.0, 6),
                                round((lat_prev_md + dlat_md) / 1000000.0, 6),
                            ]
                        )
                        lon_prev_md += dlon_md
                        lat_prev_md += dlat_md
                    if len(cblock) > 0:
                        if "nonum" in way:
                            # Check whether cblock is a closed ring,
                            # allowing for a 1 microdegree coordinate
                            # error.
                            if (
                                len(cblock) < 4
                                or abs(cblock[0][0] - cblock[-1][0]) > 0.000001
                                or abs(cblock[0][1] - cblock[-1][1]) > 0.000001
                            ):
                                rings = False
                            elif (
                                cblock[0][0] != cblock[-1][0]
                                or cblock[0][1] != cblock[-1][1]
                            ):
                                # Force this nearly closed ring to be
                                # closed.
                                cblock[-1] = cblock[0][:]
                        dblock.append(cblock)
                if len(dblock) > 0:
                    way["coordinates"].append(dblock)
            if label_coords:
                vtags["m2db_label_lat"] = (
                    way["coordinates"][0][0][1] + label_coords["lat"]
                )
                vtags["m2db_label_lon"] = (
                    way["coordinates"][0][0][0] + label_coords["lon"]
                )
            if srcf.tell() != way_endoffset:
                raise ValueError(
                    f"Invalid way ({way_num}).  [{hex(srcf.tell())}"
                    f" != {hex(way_endoffset)}]"
                )

            if "nonum" in way:
                if rings:
                    way["nonum"] = "MULTIPOLYGON"
                else:
                    way["nonum"] = "MULTILINESTRING"
            if "m2db_anum" in way or (
                "nonum" in way and way["nonum"] == "MULTIPOLYGON"
            ):
                # The format of way["coordinates"] is similar to a
                # GeoJSON MultiPolygon.  Convert them instead to the
                # form required by the shapely MultiPolygon constructor.
                way["coordinates"] = [
                    (dblock[0], dblock[1:]) for dblock in way["coordinates"]
                ]
            else:
                # Convert coordinates to a form suitable to pass to the
                # shapely MultiLineString constructor.  This is also
                # similar to the form of coordinates of a GeoJSON
                # MultiLineString.
                # Because the Mapsforge specification is somewhat
                # ambiguous about the form of coordinates for
                # multilinestrings, this correctly converts coordinates
                # that may contain multiple coordinate blocks per data
                # block as well as multiple data blocks.
                way["coordinates"] = [
                    cblock
                    for dblock in way["coordinates"]
                    for cblock in dblock
                ]
            output.append(way)
    return output


def _parse_tile(tile_x, tile_y, tile_z, fileheader, srcf):
    if (
        tile_z not in fileheader["levels"]
        or tile_x < int(_xfromlon(tile_z, fileheader["minlon"]))
        or tile_x > int(_xfromlon(tile_z, fileheader["maxlon"]))
        or tile_y < int(_yfromlat(tile_z, fileheader["maxlat"]))
        or tile_y > int(_yfromlat(tile_z, fileheader["minlat"]))
    ):
        raise ValueError(
            f"Invalid tile coordinates: ({tile_x}, {tile_y}, {tile_z})"
        )
    subfile = fileheader["subfiles"][fileheader["levels"].index(tile_z)]
    tilenum = _get_tilenum(tile_x, tile_y, subfile)
    offset = subfile["offset"] + subfile["tile_offsets"][tilenum]

    if tilenum + 1 < len(subfile["tile_offsets"]):
        endoffset = subfile["offset"] + subfile["tile_offsets"][tilenum + 1]
    else:
        endoffset = subfile["offset"] + subfile["length"]

    if offset == endoffset:
        # This is an empty tile.
        return [], []
    srcf.seek(offset)

    if fileheader["debuginfo"]:
        if not srcf.read(32).startswith(b"###TileStart"):
            raise ValueError(
                f"Tile Header (###TileStartX,Y###) not found at offset:"
                f" {hex(offset)}"
            )

    poi_counts = [0] * subfile["minzoom"]
    way_counts = [0] * subfile["minzoom"]

    for _ in range(subfile["minzoom"], subfile["maxzoom"] + 1):
        poi_counts.append(_parse_vu128(srcf))
        way_counts.append(_parse_vu128(srcf))
    first_wayoffset = _parse_vu128(srcf) + srcf.tell()
    # EPSG:4326 (lat, lon) in microdegrees at the tile's northwest
    # corner where pixel coordinates are (0,0)
    origin_md = {
        "lat": int(round(_latfromy(tile_z, tile_y) * 1000000)),
        "lon": int(round(_lonfromx(tile_z, tile_x) * 1000000)),
    }
    points = _parse_points(
        srcf,
        first_wayoffset,
        origin_md,
        fileheader["ptags"],
        fileheader["debuginfo"],
        poi_counts,
    )
    ways = _parse_ways(
        srcf,
        endoffset,
        origin_md,
        fileheader["wtags"],
        fileheader["debuginfo"],
        way_counts,
    )
    return points, ways


def _progress(tilenum, tilecount, snum, scount):
    if (
        # A somewhat arbitrary choice for when to print a progress
        # indicator
        tilenum == 1
        or (tilecount < 20)
        or (tilecount < 200 and tilenum % 10 == 0)
        or (tilecount < 2000 and tilenum % 100 == 0)
        or (tilenum % 1000 == 0)
    ):
        print(
            f"Parsing tile {tilenum} of {tilecount} from"
            f" layer {snum + 1} of {scount}."
        )


def _update_table(dbc, subfile, table, clist, vlist, sgeom, feature):
    pkey = f"m2db_{table[0]}num"
    geom_index = clist.index("m2db_geometry")
    level_index = clist.index("m2db_level")
    minz_index = clist.index("m2db_minz")
    row = dbc.execute(
        f"""
        SELECT {",".join(clist)}
            FROM {table}
            WHERE {pkey} = ?;""",
        (feature[pkey],),
    ).fetchone()
    if not row:
        dbc.execute(
            f"""
            INSERT INTO {table} ({",".join(clist)})
                VALUES ({",".join(["?"]*len(vlist))});""",
            vlist,
        )
    elif row[level_index] == subfile["level"]:
        if row[:geom_index] != tuple(vlist[:geom_index]) or row[
            geom_index + 1 :
        ] != tuple(vlist[geom_index + 1 :]):
            dbc.commit()
            print(row[:geom_index] + row[geom_index + 1 :])
            print(vlist[:geom_index] + vlist[geom_index + 1 :])
            raise ValueError(
                f"Discrepancy in value for {pkey}={feature[pkey]} found"
            )
        merged_geom = sgeom.union(swkb.loads(row[geom_index]))
        merged_geom = merged_geom.normalize()  # useful for diagnostics
        dbc.execute(
            f"""
            UPDATE {table}
                SET m2db_geometry = ?
                WHERE {pkey} = ?;""",
            (merged_geom.wkb, feature[pkey]),
        )
    elif row[level_index] > subfile["level"]:
        if row[geom_index + 1 :] != tuple(vlist[geom_index + 1 :]):
            print(row[geom_index + 1 :])
            print(vlist[geom_index + 1 :])
            raise ValueError(
                f"Discrepancy in value for {pkey}={feature[pkey]} found"
            )
        if row[minz_index] == feature["m2db_minz"]:
            # m2db_minz has already been updated by another point
            # from this level
            return
        if row[minz_index] != subfile["maxzoom"] + 1:
            dbc.commit()
            print(f"{row[minz_index]} != {subfile['maxzoom'] + 1}")
            raise ValueError(
                f"Discontinuous zoom values for {pkey}={feature[pkey]}"
                f" found"
            )
        dbc.execute(
            f"""
            UPDATE {table}
                SET m2db_minz = ?
                WHERE {pkey} = ?;""",
            (feature["m2db_minz"], feature[pkey]),
        )
    else:
        dbc.commit()
        raise ValueError("Logic error in _update_table()")


def _tile_features(
    dbc,
    srcf,
    dbl,
    tile_x,
    tile_y,
    subfile,
    fileheader,
    clists,
    nonums,  # p, l, a
    vtagslist,
):
    maxz = subfile["maxzoom"]
    tile_box = sg.box(
        round(_lonfromx(subfile["level"], tile_x), 6),
        round(_latfromy(subfile["level"], tile_y), 6),
        round(_lonfromx(subfile["level"], tile_x + 1), 6),
        round(_latfromy(subfile["level"], tile_y + 1), 6),
    )
    points, ways = _parse_tile(
        tile_x, tile_y, subfile["level"], fileheader, srcf
    )
    for point in points:
        vlist = [None] * len(clists[0])
        vlist[clists[0].index("m2db_level")] = subfile["level"]
        vlist[clists[0].index("m2db_maxz")] = maxz
        if not dbl:
            vlist[clists[0].index("m2db_pnum")] = nonums[0]
            nonums[0] += 1
        sgeom = sg.MultiPoint([sg.Point(point["lon"], point["lat"])])
        if dbl or TRIM_NONUM:
            sgeom = tile_box.intersection(sgeom, grid_size=0.000001)
        if sgeom.is_empty:
            # This occurs when the point is in the buffer included
            # around the edges of a tile.
            continue
        sgeom = sgeom.normalize()  # useful for diagnostics
        vlist[clists[0].index("m2db_geometry")] = sgeom.wkb
        for key in point:
            if key in ["lat", "lon", "nonum"]:
                continue
            if key == "vtags":
                for vtag, value in point["vtags"].items():
                    if vtag not in vtagslist:
                        vtagslist.append(vtag)
                    if vtag not in clists[0]:
                        print("new points vtag", vtag)
                        # add the column to the table
                        dbc.execute(
                            f"""
                            ALTER TABLE points
                                ADD COLUMN {vtag};"""
                        )
                        clists[0].append(vtag)
                        vlist.append(value)
                    else:
                        vlist[clists[0].index(vtag)] = value
            elif key == "m2db_tags":
                vlist[clists[0].index("m2db_tags")] = json.dumps(
                    point[key], ensure_ascii=False
                )
            else:
                vlist[clists[0].index(key)] = point[key]
        if not dbl:
            dbc.execute(
                f"""
                INSERT INTO points ({",".join(clists[0])})
                    VALUES ({",".join(["?"]*len(vlist))});""",
                vlist,
            )
            continue
        _update_table(dbc, subfile, "points", clists[0], vlist, sgeom, point)

    for way in ways:
        if (
            "m2db_lnum" in way
            or "nonum" in way
            and way["nonum"] == "MULTILINESTRING"
        ):
            vlist = [None] * len(clists[1])
            vlist[clists[1].index("m2db_level")] = subfile["level"]
            vlist[clists[1].index("m2db_maxz")] = maxz
            if not dbl:
                vlist[clists[1].index("m2db_lnum")] = nonums[1]
                nonums[1] += 1
            # Using grid_size=0.000001 ensures that sgeom has all
            # coordinates rounded to the nearest microdegree, which was
            # already the case for the source geometry, but might not
            # otherwise still be the case after intersection()
            sgeom = sg.MultiLineString(way["coordinates"])
            if sgeom.is_empty:
                raise ValueError("empty?" + str(way))
            if (not tile_box.covers(sgeom)) and (dbl or TRIM_NONUM):
                # intersection() will modify any self-intersecting
                # linestrings within sgeom.
                # So, by not doing this for sgeom that do not need
                # trimming, it avoids many unnecessary changes to sgeom.
                sgeom = tile_box.intersection(sgeom, grid_size=0.000001)
            if sgeom.is_empty:
                # This occurs when the geometry is
                # entirely in the buffer included around
                # the edges of a tile.
                continue
            if sgeom.geom_type == "GeometryCollection":
                lslist = [
                    ls for ls in sgeom.geoms if ls.geom_type == "LineString"
                ]
                if not lslist:
                    dbc.commit()
                    print(sgeom.wkt)
                    raise ValueError(
                        "Unexpected geom_type from intersection for" " lines."
                    )
                sgeom = sg.MultiLineString(lslist)
            if sgeom.geom_type not in ["LineString", "MultiLineString"]:
                continue  # Discard intersection (Point or MultiPoint)
            sgeom = sgeom.normalize()  # useful for diagnostics
            vlist[clists[1].index("m2db_geometry")] = sgeom.wkb
            for key in way:
                if key in ["coordinates", "nonum"]:
                    continue
                if key == "vtags":
                    for vtag, value in way["vtags"].items():
                        if vtag not in vtagslist:
                            vtagslist.append(vtag)
                        if vtag not in clists[1]:
                            print("new lines vtag", vtag)
                            # add the column to the table
                            dbc.execute(
                                f"""
                                ALTER TABLE lines
                                    ADD COLUMN {vtag};"""
                            )
                            clists[1].append(vtag)
                            vlist.append(value)
                        else:
                            vlist[clists[1].index(vtag)] = value
                elif key == "m2db_tags":
                    vlist[clists[1].index("m2db_tags")] = json.dumps(
                        way[key], ensure_ascii=False
                    )
                else:
                    vlist[clists[1].index(key)] = way[key]
            if not dbl:
                dbc.execute(
                    f"""
                    INSERT INTO lines ({",".join(clists[1])})
                        VALUES ({",".join(["?"]*len(vlist))});""",
                    vlist,
                )
                continue
            _update_table(dbc, subfile, "lines", clists[1], vlist, sgeom, way)
        elif (
            "m2db_anum" in way
            or "nonum" in way
            and way["nonum"] == "MULTIPOLYGON"
        ):
            vlist = [None] * len(clists[2])
            vlist[clists[2].index("m2db_level")] = subfile["level"]
            vlist[clists[2].index("m2db_maxz")] = maxz
            if not dbl:
                vlist[clists[2].index("m2db_anum")] = nonums[2]
                nonums[2] += 1
            sgeom = sg.MultiPolygon(way["coordinates"])
            if not sgeom.is_valid:
                # Even though forgemap ensures that source geometries
                # are valid and of appropriate dimension, encoding them
                # with microdegree integer coordinates can make them
                # invalid.
                vgeom = sv.make_valid(sgeom)
                if vgeom.is_valid:
                    # Sometimes make_valid returns a geometry
                    # collection containing non-area elements.
                    if vgeom.geom_type == "GeometryCollection":
                        pglist = [
                            pg
                            for pg in vgeom.geoms
                            if pg.geom_type == "Polygon"
                        ]
                        if not pglist:
                            print("Invalid multipolygon found.  Discarding")
                            print(f"   {vgeom.wkt[:80]}...{vgeom.wkt[-80:]}")
                            continue
                        vgeom2 = sg.MultiPolygon(pglist)
                        if not vgeom2.is_valid:
                            print("Invalid multipolygon found.  Discarding")
                            print(f"   {vgeom.wkt[:80]}...{vgeom.wkt[-80:]}")
                            continue
                        vgeom = vgeom2
                    print("Invalid geometry found.  Modified to make valid:")
                    print(f"invalid: {sgeom.wkt[:80]}...{sgeom.wkt[-80:]}")
                    print(f"valid: {vgeom.wkt[:80]}...{vgeom.wkt[-80:]}")
                    sgeom = vgeom
                else:
                    print("Invalid geometry found.  Discarding")
                    print(f"   {sgeom.wkt[:80]}...{sgeom.wkt[-80:]}")
                    continue
            # Using grid_size=0.000001 ensures that sgeom has all
            # coordinates rounded to the nearest microdegree, which was
            # already the case for the source geometry, but might not
            # otherwise still be the case after intersection().  This
            # may make later use of linemerge() more robust.
            if (not tile_box.covers(sgeom)) and (dbl or TRIM_NONUM):
                sgeom = tile_box.intersection(sgeom, grid_size=0.000001)
            if sgeom.is_empty:
                # This occurs when the geometry is
                # entirely in the buffer included around
                # the edges of a tile.
                continue
            if sgeom.geom_type == "GeometryCollection":
                pglist = [
                    pg for pg in sgeom.geoms if pg.geom_type == "Polygon"
                ]
                if not pglist:
                    # Discard intersection (LineString, MultiLineString,
                    # Point, MultiPoint)
                    continue
                sgeom = sg.MultiPolygon(pglist)
            if sgeom.geom_type not in ["Polygon", "MultiPolygon"]:
                # Discard intersection (LineString, MultiLineString,
                # Point, MultiPoint)
                continue
            sgeom = sgeom.normalize()  # useful for diagnostics
            vlist[clists[2].index("m2db_geometry")] = sgeom.wkb
            for key in way:
                if key in ["coordinates", "nonum"]:
                    continue
                if key == "vtags":
                    for vtag, value in way["vtags"].items():
                        if vtag not in vtagslist:
                            vtagslist.append(vtag)
                        if vtag not in clists[2]:
                            print("new areas vtag", vtag)
                            # add the column to the table
                            dbc.execute(
                                f"""
                                ALTER TABLE areas
                                    ADD COLUMN {vtag};"""
                            )
                            clists[2].append(vtag)
                            vlist.append(value)
                        else:
                            vlist[clists[2].index(vtag)] = value
                elif key == "m2db_tags":
                    vlist[clists[2].index("m2db_tags")] = json.dumps(
                        way[key], ensure_ascii=False
                    )
                else:
                    vlist[clists[2].index(key)] = way[key]
            if not dbl:
                dbc.execute(
                    f"""
                    INSERT INTO areas ({",".join(clists[2])})
                        VALUES ({",".join(["?"]*len(vlist))});""",
                    vlist,
                )
                continue
            _update_table(dbc, subfile, "areas", clists[2], vlist, sgeom, way)
        else:
            dbc.commit()
            print(way)
            raise ValueError
        dbc.commit()


def _merge_lines(dbc):
    (max_lnum,) = dbc.execute("SELECT max(m2db_lnum) FROM lines;").fetchone()
    print("merging lines")
    # 2.1 microdegrees in lat and 2.1 in lon.
    # using 2.1 rather than 2 allows for rounding error
    tol = 0.0000021
    multi_count0 = 0
    multi_count1 = 0
    multi_count2 = 0
    loops = 0
    non_loops = []
    for m2db_lnum in range(max_lnum + 1):
        row = dbc.execute(
            """
            SELECT m2db_geometry
                FROM lines
                WHERE m2db_lnum = ?;""",
            (m2db_lnum,),
        ).fetchone()
        if row is None:
            continue
        (geom_wkb,) = row
        geom = swkb.loads(geom_wkb)
        if geom.geom_type == "MultiLineString" and len(geom.geoms) > 1:
            gcount = len(geom.geoms)
            multi_count0 += 1
            geom = linemerge(geom)
            if geom.geom_type == "MultiLineString" and len(geom.geoms) > 1:
                multi_count1 += 1
                # compare the start and end coordinates of each
                # linestring of the multilinestring with the end
                # and start coordinates of all preceeding
                # linestrings, adjusting the coordinates to exactly
                # match those that already match within tol.
                # Assume that direction of all linestrings in the
                # multilinestring is already correct.
                multi_coords = [ls.coords for ls in geom.geoms]
                for i in range(1, len(multi_coords)):
                    for j in range(i):
                        if (
                            abs(multi_coords[i][0][0] - multi_coords[j][-1][0])
                            <= tol
                            and abs(
                                multi_coords[i][0][1] - multi_coords[j][-1][1]
                            )
                            <= tol
                        ):
                            # replace multi_coords[i][0] with
                            # multi_coords[j][-1]
                            multi_coords[i] = (
                                multi_coords[j][-1:] + multi_coords[i][1:]
                            )
                        if (
                            abs(multi_coords[i][-1][0] - multi_coords[j][0][0])
                            <= tol
                            and abs(
                                multi_coords[i][-1][1] - multi_coords[j][0][1]
                            )
                            <= tol
                        ):
                            # replace multi_coords[i][-1] with
                            # multi_coords[j][0]
                            multi_coords[i] = (
                                multi_coords[i][:-1] + multi_coords[j][:1]
                            )
                # create new multilinestring from modified multi_coords,
                # and try linemerge again
                geom = linemerge(sg.MultiLineString(multi_coords))
                if geom.geom_type == "MultiLineString" and len(geom.geoms) > 1:
                    multi_count2 += 1
                    for line in geom.geoms:
                        if (
                            line.coords[0][0] == line.coords[-1][0]
                            and line.coords[0][1] == line.coords[-1][1]
                        ):
                            # This multilinestring contains at least one
                            # loop and at least two linestrings.  So, it
                            # cannot be converted to a single linestring
                            # by linemerge().
                            loops += 1
                            break
                    else:
                        non_loops.append(m2db_lnum)
            if geom.geom_type == "LineString" or len(geom.geoms) < gcount:
                geom = geom.normalize()  # useful for diagnostics
                dbc.execute(
                    """
                    UPDATE lines
                        SET m2db_geometry = ?
                        WHERE m2db_lnum = ?;""",
                    (geom.wkb, m2db_lnum),
                )
    dbc.commit()
    # Diagnostic
    print(
        "multi_counts",
        multi_count0,
        multi_count1,
        multi_count2,
        "loops =",
        loops,
    )
    # print("non_loops m2db_lnums", non_loops)


def _write_config(config_pathname, db_pathname, dblstr, fileheader, vtagslist):
    with open(config_pathname, "w", encoding="utf-8") as tomlf:
        tomlf.write(
            f"## This is a configuration file for forgemap.\n"
            f"## It was created by map2db.py to facilitate building a"
            f" new mapsforge map from the output database:"
            f' "{db_pathname}".\n'
            f"## See the forgemap documentation for a full explanation"
            f" of the content of this configuration file.\n\n"
        )
        tomlf.write(
            f"## name of the mapsforge file to produce.\n"
            f'output_pathname = "{db_pathname}.map"\n'
        )
        tomlf.write(f'dbl = """{dblstr}"""\n')
        tomlf.write(
            f'minlat = {fileheader["minlat"]}\n'
            f'maxlat = {fileheader["maxlat"]}\n'
            f'minlon = {fileheader["minlon"]}\n'
            f'maxlon = {fileheader["maxlon"]}\n'
        )
        if "startlat" in fileheader:
            tomlf.write(f'start_lat = {fileheader["startlat"]}\n')
        if "startlon" in fileheader:
            tomlf.write(f'start_lon = {fileheader["startlon"]}\n')
        if "startzoom" in fileheader:
            tomlf.write(f'start_zoom = {fileheader["startzoom"]}\n')
        if "commstr" in fileheader:
            tomlf.write(f'comment_str = """{fileheader["commstr"]}"""\n')
        if "createdby" in fileheader:
            tomlf.write(f'created_by = """{fileheader["createdby"]}"""\n')
        vtags_str = ""
        if vtagslist:
            vtags_str += "value_tags = [\n"
            tomlf.write(
                '## tags whose "key" is in variable_tags will have'
                ' their "value" encoded with each individual feature.'
                "  This is useful for widely used keys whose values are"
                " each used by only one (or a few) features.\n"
                "variable_tags = [\n"
            )
            for vtag in vtagslist:
                if vtag in [
                    "m2db_name",
                    "m2db_house_num",
                    "m2db_elevation",
                    "m2db_reference",
                    "m2db_label_lat",
                    "m2db_label_lon",
                ]:
                    continue
                vtags_str += f'    {{column_name = "{vtag}"}},\n'
                tomlf.write(f'    "{vtag}",\n')
            tomlf.write("]\n")
            vtags_str += "]\n"

        tomlf.write(
            "## [[map_data]] is a list of tables, each of which defines"
            " a subfile that will contain the vector map data for a"
            " specified zoom interval.  These must be in order from"
            " smallest to largest.  All values must be integers.\n"
            "## For each subfile: minzoom <= level <= maxzoom\n"
            "## map_data[i][maxzoom] < map_data[i+1][minzoom]\n"
            "## (Typically: map_data[i][maxzoom] + 1 ="
            " map_data[i+1][minzoom])\n"
        )
        for subfile in fileheader["subfiles"]:
            tomlf.write(
                f"[[map_data]]\n"
                f'minzoom = {subfile["minzoom"]}\n'
                f'level = {subfile["level"]}\n'
                f'maxzoom = {subfile["maxzoom"]}\n\n'
            )

        # Set nosimplify to the "level" of the last subfile since
        # database contents are suitable for at least this subfile
        # without simplification.
        tomlf.write(
            f"## [[sources]] is a list of tables, each of which defines"
            f" a datafile from which source data will be read, and how"
            f" that data will be interpreted.  More than one source may"
            f" use the same datafile, so that different subsets of data"
            f" from that datafile may be interpreted differently.  This"
            f" is done by using different layer or sql settings.\n"
            f"[[sources]]\n"
            f'pathname = "{db_pathname}"\n'
            f'layer = "points"\n'
            f'name_column = "m2db_name"\n'
            f'house_num_column = "m2db_house_num"\n'
            f'elevation_column = "m2db_elevation"\n'
            f'nosimplify = {fileheader["subfiles"][-1]["level"]}\n'
            f'minzoom_column = "m2db_minz"\n'
            f'maxzoom_column = "m2db_maxz"\n'
            f'json_tags = ["m2db_tags"]\n'
            f"{vtags_str}"
            f"\n"
            f"[[sources]]\n"
            f'pathname = "{db_pathname}"\n'
            f'layer = "lines"\n'
            f'name_column = "m2db_name"\n'
            f'house_num_column = "m2db_house_num"\n'
            f'reference_column = "m2db_reference"\n'
            f'label_lat_column = "m2db_label_lat"\n'
            f'label_lon_column = "m2db_label_lon"\n'
            f'nosimplify = {fileheader["subfiles"][-1]["level"]}\n'
            f'minzoom_column = "m2db_minz"\n'
            f'maxzoom_column = "m2db_maxz"\n'
            f'json_tags = ["m2db_tags"]\n'
            f"{vtags_str}"
            f"\n"
            f"[[sources]]\n"
            f'pathname = "{db_pathname}"\n'
            f'layer = "areas"\n'
            f'name_column = "m2db_name"\n'
            f'house_num_column = "m2db_house_num"\n'
            f'reference_column = "m2db_reference"\n'
            f'label_lat_column = "m2db_label_lat"\n'
            f'label_lon_column = "m2db_label_lon"\n'
            f'nosimplify = {fileheader["subfiles"][-1]["level"]}\n'
            f'minzoom_column = "m2db_minz"\n'
            f'maxzoom_column = "m2db_maxz"\n'
            f'json_tags = ["m2db_tags"]\n'
            f"{vtags_str}"
        )


def map2db(src_pathname, db_pathname):
    """
    Recreate the geospatial database used to build a mapsforge map.

    Parameters
    ----------
    src_pathname : str
        The pathname of the mapsforge map to read.
    db_pathname : str
        The pathname of the geospatial database to be created.
    """

    if os.path.exists(db_pathname):
        os.remove(db_pathname)
    with (
        sqlite3.connect(db_pathname) as dbc,
        open(src_pathname, "rb") as srcf,
    ):
        _prepare_dbc(dbc)
        # These columns are initially present in the points, lines, and
        # areas tables respectively.  Additional columns will be added
        # when/if they are needed.  m2db_pnum, m2db_lnum, and m2db_anum
        # are the primary keys to these tables.
        p_clist = [
            "m2db_pnum",
            "m2db_level",
            "m2db_minz",
            "m2db_maxz",
            "m2db_geometry",
            "m2db_layer_num",
            "m2db_tags",
        ]
        l_clist = [
            "m2db_lnum",
            "m2db_level",
            "m2db_minz",
            "m2db_maxz",
            "m2db_geometry",
            "m2db_layer_num",
            "m2db_tags",
        ]
        a_clist = [
            "m2db_anum",
            "m2db_level",
            "m2db_minz",
            "m2db_maxz",
            "m2db_geometry",
            "m2db_layer_num",
            "m2db_tags",
        ]

        fileheader = _parse_fileheader(srcf)
        if fileheader["wtags"][-1].startswith("_lbd_="):
            dbl = True
            dbldesc = (
                "Text that was encoded in the map file specifically so"
                " that it could be included here."
            )
            dblstr = fileheader["wtags"][-1][len("_lbd_=") :][-1::-1]
        else:
            dbl = False
            dbldesc = "Text generated by map2db.py"

            dblstr = (
                f"The map file from which this database was created,"
                f" '{src_pathname}', does not appear to have been"
                f" created by forgemap with the dbl option.  Thus, this"
                f" database includes the raw content from all tiles"
                f" rather than a properly reconstructed database with"
                f" various parts of each feature merged together and"
                f" more simplified versions of such features discarded"
                f" in favor of the less simplified versions.  The dbl"
                f" option of forgemap also explicitly encodes whether"
                f" the source geometry was a line or area feature."
                f"  Since this information is not available, way"
                f" features in which every coordinate block forms a"
                f" closed ring with at least four points are assumed to"
                f" be areas, while the remainder are assumed to be"
                f" lines.  This may misclassify some features as areas"
                f" that should have been lines.  Since not all closed"
                f" rings of at least four points are valid, this may"
                f" produce unexpected results."
                f"\nSince no database license statement was encoded in"
                f" the map file, the user is responsible for"
                f" determining what legal restrictions apply to use of"
                f" this database."
            )
        dbc.execute(
            """
            INSERT INTO README (desc, text)
                VALUES (?, ?);""",
            (dbldesc, dblstr),
        )
        dbc.execute(
            """
            INSERT INTO README (desc, text)
                VALUES (?, ?);""",
            (
                "Text generated by map2db.py",
                f"(This database was recreated from {src_pathname} by"
                f" map2db.py <https://www.github.com/pflarue/map2db>.)",
            ),
        )
        if "commstr" in fileheader and fileheader["commstr"]:
            dbc.execute(
                """
                INSERT INTO README (desc, text)
                    VALUES (?, ?);""",
                (
                    "Text from the 'comment' field in the file header"
                    " of the map file from which this database was"
                    " recreated. This often includes a copyright and/or"
                    " license statement for the map file (not for this"
                    " database!).",
                    fileheader["commstr"],
                ),
            )
        if "createdby" in fileheader and fileheader["createdby"]:
            dbc.execute(
                """
                INSERT INTO README (desc, text)
                    VALUES (?, ?);""",
                (
                    "Text from the 'created by' field in the file"
                    " header of the map file from which this database"
                    " was recreated. Per the specification for the"
                    " Mapsforge Binary Map File Format, this should be"
                    " the name of the application which created the"
                    " map file.",
                    fileheader["createdby"],
                ),
            )

        for key in ["minlat", "minlon", "maxlat", "maxlon"]:
            dbc.execute(
                """
                INSERT INTO limits (key, value)
                    VALUES (?, ?);""",
                (key, fileheader[key]),
            )
        dbc.commit()

        vtagslist = []  # This will be used for the config.toml file.
        nonums = [0, 0, 0]  # Used for non-dbl map files
        for subfile_num, subfile in enumerate(fileheader["subfiles"]):
            dbc.execute(
                """
                INSERT INTO subfiles
                    (subfile_num, level, minzoom, maxzoom)
                    VALUES (?, ?, ?, ?);""",
                (
                    subfile_num,
                    subfile["level"],
                    subfile["minzoom"],
                    subfile["maxzoom"],
                ),
            )
            srcf.seek(subfile["offset"])
            if fileheader["debuginfo"]:
                srcf.seek(subfile["offset"])
                if srcf.read(16) != b"+++IndexStart+++":
                    raise ValueError(
                        f"Tile Index Header (+++IndexStart+++) not"
                        f" found at offset: {hex(subfile['offset'])}"
                    )
            subfile["tile_offsets"] = []
            for _ in range(subfile["tile_count"]):
                # The msb of this byte undicates that the tile is
                # entirely covered with water.  This attribute will not
                # be encoded in the database.
                subfile["tile_offsets"].append(
                    (ord(srcf.read(1)) & 0x7F) * 2**32 + _parse_u32(srcf)
                )

        # Cycle through subfiles in reverse order so as to encounter
        # the less simplified geometry first.
        for snum, subfile in enumerate(reversed(fileheader["subfiles"])):
            tilecount = (subfile["maxy"] - subfile["miny"] + 1) * (
                subfile["maxx"] - subfile["minx"] + 1
            )
            print(f"level = {subfile['level']} with {tilecount} tiles.")
            for tile_y in range(subfile["miny"], subfile["maxy"] + 1):
                for tile_x in range(subfile["minx"], subfile["maxx"] + 1):
                    tilenum = _get_tilenum(tile_x, tile_y, subfile) + 1
                    _progress(
                        tilenum, tilecount, snum, len(fileheader["subfiles"])
                    )
                    # dbc.commit() at end of _tile_features()
                    _tile_features(
                        dbc,
                        srcf,
                        dbl,
                        tile_x,
                        tile_y,
                        subfile,
                        fileheader,
                        [p_clist, l_clist, a_clist],
                        nonums,
                        vtagslist,
                    )
        if dbl:
            _merge_lines(dbc)  # dbc.commit() at end of _merge_lines()
            # Also write a config.toml file for forgemap
            config_pathname = db_pathname + ".config.toml"
            _write_config(
                config_pathname, db_pathname, dblstr, fileheader, vtagslist
            )
            print(f"forgemap configuration file written to {config_pathname}")
        dbc.execute(
            """
            INSERT INTO README (desc, text)
                VALUES (?, ?);""",
            (
                "Text generated by map2db.py",
                "(database recreation completed.)",
            ),
        )
    dbc.commit()
    print(f"\ndatabase written to '{db_pathname}'")
    print("\nThe following is also in the README table of the database:")
    print(f'Description: """\n{dbldesc}\n"""')
    print(f'Text: """\n{dblstr}\n"""')


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("map pathname not provided")
    if not sys.argv[1].endswith(".map"):
        raise ValueError(
            f"{sys.argv[1]} does not appear to be a mapsforge map file"
        )
    src = sys.argv[1]
    if len(sys.argv) >= 3:
        dst = sys.argv[2]
    else:
        dst = src[:-4] + ".db"
    map2db(src, dst)
