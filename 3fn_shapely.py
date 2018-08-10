'''
2.shapely :(has geos_c.dll)
1.translate()	from shapely.affinity import translate	(affinity.py in shapely folder)
2.if speedups.available: speedups.enable()	shapely import speedups,	(_speedups.py)
3.rotate(bx, angle, origin=bx.centroid)		shapely.affinity import rotate		(affinity.py in shapely folder)	
4.Shapely : shapely.geometry import Point, LinearRing, LineString, box,  (geometry folder in shapely)
5.from shapely.geometry import Polygon, MultiPolygon (geometry folder)	Polygon(exCtr, inCtr)	
6.poly.is_valid:	polygon.py
8.poly.interiors	polygon.py
9.exCtr = tPolyC[tP[0]].exterior.coords[:]	polygon.py
'''
#1.translate()	from shapely.affinity import translate	(affinity.py in shapely folder)

def translate(geom, xoff=0.0, yoff=0.0, zoff=0.0):
    """Returns a translated geometry shifted by offsets along each dimension.

    The general 3D affine transformation matrix for translation is:

        / 1  0  0 xoff \ 
        | 0  1  0 yoff |
        | 0  0  1 zoff |
        \ 0  0  0   1  /
    """
    matrix = (1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0,
              xoff, yoff, zoff)
    return affine_transform(geom, matrix)

#from math import sin, cos, tan, pi

__all__ = ['affine_transform', 'rotate', 'scale', 'skew', 'translate']


def affine_transform(geom, matrix):
    """Returns a transformed geometry using an affine transformation matrix.

    The coefficient matrix is provided as a list or tuple with 6 or 12 items
    for 2D or 3D transformations, respectively.

    For 2D affine transformations, the 6 parameter matrix is::

        [a, b, d, e, xoff, yoff]

    which represents the augmented matrix::

                            / a  b xoff \ 
        [x' y' 1] = [x y 1] | d  e yoff |
                            \ 0  0   1  /

    or the equations for the transformed coordinates::

        x' = a * x + b * y + xoff
        y' = d * x + e * y + yoff

    For 3D affine transformations, the 12 parameter matrix is::

        [a, b, c, d, e, f, g, h, i, xoff, yoff, zoff]

    which represents the augmented matrix::

                                 / a  b  c xoff \ 
        [x' y' z' 1] = [x y z 1] | d  e  f yoff |
                                 | g  h  i zoff |
                                 \ 0  0  0   1  /

    or the equations for the transformed coordinates::

        x' = a * x + b * y + c * z + xoff
        y' = d * x + e * y + f * z + yoff
        z' = g * x + h * y + i * z + zoff
    """
    if geom.is_empty:
        return geom
    if len(matrix) == 6:
        ndim = 2
        a, b, d, e, xoff, yoff = matrix
        if geom.has_z:
            ndim = 3
            i = 1.0
            c = f = g = h = zoff = 0.0
            matrix = a, b, c, d, e, f, g, h, i, xoff, yoff, zoff
    elif len(matrix) == 12:
        ndim = 3
        a, b, c, d, e, f, g, h, i, xoff, yoff, zoff = matrix
        if not geom.has_z:
            ndim = 2
            matrix = a, b, d, e, xoff, yoff
    else:
        raise ValueError("'matrix' expects either 6 or 12 coefficients")
        
#--------------------------------------------------------------------------------------------------------
#2.if speedups.available: speedups.enable()	shapely import speedups,	(_speedups.py)

#import warnings

#from shapely.geometry import linestring, polygon
#from shapely import coords
#import shapely.affinity

try:
    from shapely.speedups import _speedups
    available = True
    import_error_msg = None
except ImportError:
    import sys
    available = False
    # TODO: This does not appear to do anything useful
    import_error_msg = sys.exc_info()[1]

from ..ftools import wraps
def method_wrapper(f):
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    return wraps(f)(wrapper)

__all__ = ['available', 'enable', 'disable']
_orig = {}

def enable():
    if not available:
        warnings.warn("shapely.speedups not available", RuntimeWarning)
        return
    
    if _orig:
        return
    
    _orig['CoordinateSequence.ctypes'] = coords.CoordinateSequence.ctypes
    coords.CoordinateSequence.ctypes = property(_speedups.coordseq_ctypes)
    
    _orig['CoordinateSequence.__iter__'] = coords.CoordinateSequence.__iter__
    coords.CoordinateSequence.__iter__ = method_wrapper(_speedups.coordseq_iter)

    _orig['geos_linestring_from_py'] = linestring.geos_linestring_from_py
    linestring.geos_linestring_from_py = _speedups.geos_linestring_from_py

    _orig['geos_linearring_from_py']  = polygon.geos_linearring_from_py
    polygon.geos_linearring_from_py = _speedups.geos_linearring_from_py
    
    _orig['affine_transform'] = shapely.affinity.affine_transform
    # copy docstring from original function
    def affine_transform(geom, matrix):
        return _speedups.affine_transform(geom, matrix)
    affine_transform.__doc__ = shapely.affinity.affine_transform.__doc__
    shapely.affinity.affine_transform = affine_transform

def disable():
    if not _orig:
        return

    coords.CoordinateSequence.ctypes = _orig['CoordinateSequence.ctypes']
    coords.CoordinateSequence.__iter__ = _orig['CoordinateSequence.__iter__']
    linestring.geos_linestring_from_py = _orig['geos_linestring_from_py']
    polygon.geos_linearring_from_py = _orig['geos_linearring_from_py']
    shapely.affinity.affine_transform = _orig['affine_transform']
    _orig.clear()
#--------------------------------------------------------------------------------------------------------------------------------
#3.rotate(bx, angle, origin=bx.centroid)		shapely.affinity import rotate		(affinity.py in shapely folder)

def rotate(geom, angle, origin='center', use_radians=False):
    """Returns a rotated geometry on a 2D plane.

    The angle of rotation can be specified in either degrees (default) or
    radians by setting ``use_radians=True``. Positive angles are
    counter-clockwise and negative are clockwise rotations.

    The point of origin can be a keyword 'center' for the bounding box
    center (default), 'centroid' for the geometry's centroid, a Point object
    or a coordinate tuple (x0, y0).

    The affine transformation matrix for 2D rotation is:

      / cos(r) -sin(r) xoff \ 
      | sin(r)  cos(r) yoff |
      \   0       0      1  /

    where the offsets are calculated from the origin Point(x0, y0):

        xoff = x0 - x0 * cos(r) + y0 * sin(r)
        yoff = y0 - x0 * sin(r) - y0 * cos(r)
    """
    if not use_radians:  # convert from degrees
        angle *= pi/180.0
    cosp = cos(angle)
    sinp = sin(angle)
    if abs(cosp) < 2.5e-16:
        cosp = 0.0
    if abs(sinp) < 2.5e-16:
        sinp = 0.0
    x0, y0 = interpret_origin(geom, origin, 2)

    matrix = (cosp, -sinp, 0.0,
              sinp,  cosp, 0.0,
              0.0,    0.0, 1.0,
              x0 - x0 * cosp + y0 * sinp, y0 - x0 * sinp - y0 * cosp, 0.0)
    return affine_transform(geom, matrix)
    
    
def interpret_origin(geom, origin, ndim):
    """Returns interpreted coordinate tuple for origin parameter.

    This is a helper function for other transform functions.

    The point of origin can be a keyword 'center' for the 2D bounding box
    center, 'centroid' for the geometry's 2D centroid, a Point object or a
    coordinate tuple (x0, y0, z0).
    """
    # get coordinate tuple from 'origin' from keyword or Point type
    if origin == 'center':
        # bounding box center
        minx, miny, maxx, maxy = geom.bounds
        origin = ((maxx + minx)/2.0, (maxy + miny)/2.0)
    elif origin == 'centroid':
        origin = geom.centroid.coords[0]
    elif isinstance(origin, str):
        raise ValueError("'origin' keyword %r is not recognized" % origin)
    elif hasattr(origin, 'type') and origin.type == 'Point':
        origin = origin.coords[0]

    # origin should now be tuple-like
    if len(origin) not in (2, 3):
        raise ValueError("Expected number of items in 'origin' to be "
                         "either 2 or 3")
    if ndim == 2:
        return origin[0:2]
    else:  # 3D coordinate
        if len(origin) == 2:
            return origin + (0.0,)
        else:
            return origin
            
#-----------------------------------------------------------------------------------------------------------------------------------
#Shapely : shapely.geometry import Point, LinearRing, LineString, box,  (geometry folder in shapely)

class LineString(BaseGeometry):
    """
    A one-dimensional figure comprising one or more line segments

    A LineString has non-zero length and zero area. It may approximate a curve
    and need not be straight. Unlike a LinearRing, a LineString is not closed.
    """

    def __init__(self, coordinates=None):
        """
        Parameters
        ----------
        coordinates : sequence
            A sequence of (x, y [,z]) numeric coordinate pairs or triples or
            an object that provides the numpy array interface, including
            another instance of LineString.

        Example
        -------
        Create a line with two segments

          >>> a = LineString([[0, 0], [1, 0], [1, 1]])
          >>> a.length
          2.0
        """
        BaseGeometry.__init__(self)
        if coordinates is not None:
            self._set_coords(coordinates)

    @property
    def __geo_interface__(self):
        return {
            'type': 'LineString',
            'coordinates': tuple(self.coords)
            }

    def svg(self, scale_factor=1., stroke_color=None):
        """Returns SVG polyline element for the LineString geometry.

        Parameters
        ==========
        scale_factor : float
            Multiplication factor for the SVG stroke-width.  Default is 1.
        stroke_color : str, optional
            Hex string for stroke color. Default is to use "#66cc99" if
            geometry is valid, and "#ff3333" if invalid.
        """
        if self.is_empty:
            return '<g />'
        if stroke_color is None:
            stroke_color = "#66cc99" if self.is_valid else "#ff3333"
        pnt_format = " ".join(["{0},{1}".format(*c) for c in self.coords])
        return (
            '<polyline fill="none" stroke="{2}" stroke-width="{1}" '
            'points="{0}" opacity="0.8" />'
            ).format(pnt_format, 2. * scale_factor, stroke_color)

    @property
    def ctypes(self):
        if not self._ctypes_data:
            self._ctypes_data = self.coords.ctypes
        return self._ctypes_data

    def array_interface(self):
        """Provide the Numpy array protocol."""
        return self.coords.array_interface()

    __array_interface__ = property(array_interface)

    # Coordinate access
    def _set_coords(self, coordinates):
        self.empty()
        self._geom, self._ndim = geos_linestring_from_py(coordinates)

    coords = property(BaseGeometry._get_coords, _set_coords)

    @property
    def xy(self):
        """Separate arrays of X and Y coordinate values

        Example:

          >>> x, y = LineString(((0, 0), (1, 1))).xy
          >>> list(x)
          [0.0, 1.0]
          >>> list(y)
          [0.0, 1.0]
        """
        return self.coords.xy

    def parallel_offset(
            self, distance, side='right',
            resolution=16, join_style=JOIN_STYLE.round, mitre_limit=5.0):

        """Returns a LineString or MultiLineString geometry at a distance from
        the object on its right or its left side.

        The side parameter may be 'left' or 'right' (default is 'right'). The
        resolution of the buffer around each vertex of the object increases by
        increasing the resolution keyword parameter or third positional
        parameter. Vertices of right hand offset lines will be ordered in
        reverse.

        The join style is for outside corners between line segments. Accepted
        values are JOIN_STYLE.round (1), JOIN_STYLE.mitre (2), and
        JOIN_STYLE.bevel (3).

        The mitre ratio limit is used for very sharp corners. It is the ratio
        of the distance from the corner to the end of the mitred offset corner.
        When two line segments meet at a sharp angle, a miter join will extend
        far beyond the original geometry. To prevent unreasonable geometry, the
        mitre limit allows controlling the maximum length of the join corner.
        Corners with a ratio which exceed the limit will be beveled.
        """
        if mitre_limit == 0.0:
            raise ValueError(
                'Cannot compute offset from zero-length line segment')
        try:
            return geom_factory(self.impl['parallel_offset'](
                self, distance, resolution, join_style, mitre_limit, side))
        except OSError:
            raise TopologicalError()
            
#-------------------------------------------------------------------------------------------
#LinearRing(in polygon.py)

class LinearRing(LineString):
    """
    A closed one-dimensional feature comprising one or more line segments

    A LinearRing that crosses itself or touches itself at a single point is
    invalid and operations on it may fail.
    """

    def __init__(self, coordinates=None):
        """
        Parameters
        ----------
        coordinates : sequence
            A sequence of (x, y [,z]) numeric coordinate pairs or triples

        Rings are implicitly closed. There is no need to specific a final
        coordinate pair identical to the first.

        Example
        -------
        Construct a square ring.

          >>> ring = LinearRing( ((0, 0), (0, 1), (1 ,1 ), (1 , 0)) )
          >>> ring.is_closed
          True
          >>> ring.length
          4.0
        """
        BaseGeometry.__init__(self)
        if coordinates is not None:
            self._set_coords(coordinates)

    @property
    def __geo_interface__(self):
        return {
            'type': 'LinearRing',
            'coordinates': tuple(self.coords)
            }

    # Coordinate access

    _get_coords = BaseGeometry._get_coords

    def _set_coords(self, coordinates):
        self.empty()
        self._geom, self._ndim = geos_linearring_from_py(coordinates)

    coords = property(_get_coords, _set_coords)

    @property
    def is_ccw(self):
        """True is the ring is oriented counter clock-wise"""
        return bool(self.impl['is_ccw'](self))

    @property
    def is_simple(self):
        """True if the geometry is simple, meaning that any self-intersections
        are only at boundary points, else False"""
        return LineString(self).is_simple


            
            
            




