import numpy as np

from osgeo import ogr, gdal, osr
import fiona

from shapely.geometry import shape, Point, Polygon, MultiLineString, MultiPoint, MultiPolygon, LineString




def load_shapefile(filename, parameters = []):
  
    c = fiona.open(filename)

    if c[0]['geometry']['type'] == 'Polygon':
        shp = MultiPolygon([shape(pol['geometry']) for pol in c])

    elif c[0]['geometry']['type'] == 'LineString':
        shp = MultiLineString([shape(pol['geometry']) for pol in c])

    elif c[0]['geometry']['type'] == 'Point':
        shp = MultiPoint([shape(pol['geometry']) for pol in c])

    else:
        shp = [shape(pol['geometry']) for pol in c]
        
        
    
    if parameters is 'all':
        parameters = c[0]['properties'].keys()
        
    if type(parameters) is not list:
    
        parameters = list(parameters)


    shp_params = {}

    for param in parameters:
        shp_params[param] = [line['properties'][param] for line in c]

    c = None

    return shp, shp_params
    
    
    
def read_tiff_as_array(filename, get_info = True, normalize = False):
    
    src_ds = gdal.Open(filename)

    try:
        srcband = src_ds.GetRasterBand(1)
    except RuntimeError, e:
        # for example, try GetRasterBand(10)
        print 'Band ( %i ) not found' % band_num
        print e
        sys.exit(1)

    raster = srcband.ReadAsArray()
    
    if normalize:
        raster = (raster - np.min(raster)) / (np.max(raster) - np.min(raster))
    
    if get_info:
    
        ulx, xres, xskew, uly, yskew, yres  = src_ds.GetGeoTransform()
        lrx = ulx + (src_ds.RasterXSize * xres)
        lry = uly + (src_ds.RasterYSize * yres)
        
        srcband = src_ds = None
        
        return raster, (ulx, lrx, lry, uly, np.abs(xres), np.abs(yres))
    
    else:
    
        srcband = src_ds = None
        
        return raster
    
    
def create_shapefile_from_shapely_multi(features, filename,
                                        fields = {}, field_type = {},
                                        buffer_width = 0, spatial_ref = 32645):
    '''
    Creates a shapefile from a
    Shapely MultiPolygon, MultiLineString, or MultiPoint
    '''


    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(filename)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(spatial_ref)

    layer = ds.CreateLayer('', srs, ogr.wkbPolygon)

    for f in fields.keys():
        fieldDefn = ogr.FieldDefn(f, field_type[f])
        layer.CreateField(fieldDefn)

    defn = layer.GetLayerDefn()


    for i in range(len(features)):

        poly = features[i].buffer(buffer_width)

        # Create a new feature (attribute and geometry)
        feat = ogr.Feature(defn)

        for f in fields.keys():
            feat.SetField(f, fields[f][i])

        # Make a geometry from Shapely object
        geom = ogr.CreateGeometryFromWkb(poly.wkb)
        feat.SetGeometry(geom)

        layer.CreateFeature(feat)
        feat = geom = None  # destroy these


    # Save and close everything
    ds = layer = feat = geom = None




def outline_to_mask(line, x, y):
    """
    Create mask from outline contour

    Parameters
    ----------
    line: array-like (N, 2)
    x, y: 1-D grid coordinates (input for meshgrid)

    Returns
    -------
    mask : 2-D boolean array (True inside)

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> poly = Point(0,0).buffer(1)
    >>> x = np.linspace(-5,5,100)
    >>> y = np.linspace(-5,5,100)
    >>> mask = outline_to_mask(poly.boundary, x, y)
    
    Modified from https://gist.github.com/perrette/a78f99b76aed54b6babf3597e0b331f8
    
    """
    import matplotlib.path as mplp
    mpath = mplp.Path(line)
    X, Y = np.meshgrid(x, y)
    points = np.array((X.flatten(), Y.flatten())).T
    mask = mpath.contains_points(points).reshape(X.shape)
    return mask


def _grid_bbox(x, y):
    dx = dy = 0
    return x[0]-dx/2, x[-1]+dx/2, y[0]-dy/2, y[-1]+dy/2

def _bbox_to_rect(bbox):
    l, r, b, t = bbox
    return Polygon([(l, b), (r, b), (r, t), (l, t)])

def shp_mask(shp, x, y, m=None):
    """
    Use recursive sub-division of space and shapely contains method to create a raster mask on a regular grid.

    Parameters
    ----------
    shp : shapely's Polygon (or whatever with a "contains" method and intersects method)
    x, y : 1-D numpy arrays defining a regular grid
    m : mask to fill, optional (will be created otherwise)

    Returns
    -------
    m : boolean 2-D array, True inside shape.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> poly = Point(0,0).buffer(1)
    >>> x = np.linspace(-5,5,100)
    >>> y = np.linspace(-5,5,100)
    >>> mask = shp_mask(poly, x, y)
    """
    rect = _bbox_to_rect(_grid_bbox(x, y))
    
    if m is None:
        m = np.zeros((y.size, x.size), dtype=bool)
               
    if not shp.intersects(rect):
        m[:] = False
    
    elif shp.contains(rect):
        m[:] = True
    
    else:
        k, l = m.shape
        
        if k == 1 and l == 1:
            m[:] = shp.contains(Point(x[0], y[0]))
            
        elif k == 1:
            m[:, :l//2] = shp_mask(shp, x[:l//2], y, m[:, :l//2])
            m[:, l//2:] = shp_mask(shp, x[l//2:], y, m[:, l//2:])
            
        elif l == 1:
            m[:k//2] = shp_mask(shp, x, y[:k//2], m[:k//2])
            m[k//2:] = shp_mask(shp, x, y[k//2:], m[k//2:])
        
        else:
            m[:k//2, :l//2] = shp_mask(shp, x[:l//2], y[:k//2], m[:k//2, :l//2])
            m[:k//2, l//2:] = shp_mask(shp, x[l//2:], y[:k//2], m[:k//2, l//2:])
            m[k//2:, :l//2] = shp_mask(shp, x[:l//2], y[k//2:], m[k//2:, :l//2])
            m[k//2:, l//2:] = shp_mask(shp, x[l//2:], y[k//2:], m[k//2:, l//2:])
        
    return m


def Polygon_axes(polygon):
    '''
    Calculates major and minor axes, and major axis azimuth, for a Shapely polygon
    '''

    rect = polygon.minimum_rotated_rectangle.exterior

    len1 = LineString(rect.coords[0:2]).length
    len2 = LineString(rect.coords[1:3]).length

    if len1 > len2:
        direction = azimuth(Point(rect.coords[0]), Point(rect.coords[1]))
    else:
        direction = azimuth(Point(rect.coords[1]), Point(rect.coords[2]))
        
    if direction > 180:
        direction = direction - 180

    axminor, axmajor = np.sort([len1, len2])

    return axminor, axmajor, direction



def azimuth(point1, point2):
    '''azimuth between 2 shapely points (interval 0 - 360), from vertical'''
    
    angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
    az = np.degrees(angle) if angle > 0 else np.degrees(angle) + 360
    
    return az