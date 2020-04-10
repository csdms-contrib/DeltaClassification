import numpy as np
import cPickle as pickle
from scipy.ndimage import morphology

from osgeo import ogr, gdal, osr
import fiona

from shapely.geometry import shape, Point, Polygon, MultiLineString, MultiPoint, MultiPolygon, LineString
from shapely import affinity

from utilities import *


# Set calculate_params to True to calculate all parameters (slow)
# Set to False to load previously saved parameters

calculate_params = False

network_filepath = '_input/network.shp'
island_filepath = '_input/islands.shp'
patch_filepath = '_input/patches.shp'




''' LOAD INPUT DATA '''

print 'Loading input data'

# load island shapefile
islands, _ = load_shapefile(island_filepath)

# load patch shapefile
# A patch is the area between centerlines of bounding channels for each island
patches, _ = load_shapefile(patch_filepath)

# load network shapefile
network_lines, params = load_shapefile(network_filepath, parameters = ['Width'])
network_widths = params['Width']




''' FIND BOUNDING AND INTERIOR CHANNELS

Identify network lines that surround islands (bounding channels)
or that drain islands (interior channels)
'''

if calculate_params:

    print 'Finding bounding and interior channels'

    # get midpoints of all network lines
    # add buffers so they touch both neighbors
    midpts = [l.interpolate(0.5, normalized=True).buffer(5) for l in network_lines]

    bounding_channels = []
    interior_channels = []

    # compare the location of midpts to island patches
    for polygon in patches:

        # if midpoint intersects patch, the line is a bounding channel
        touch = [i for i,l in enumerate(midpts) if polygon.exterior.intersects(l)]
        bounding_channels.append(touch)

        # if midpoint is completely within patch,
        # then the line is an interior channel
        touch = [i for i,l in enumerate(midpts) if polygon.contains(l)]
        interior_channels.append(touch)

    pickle.dump(bounding_channels,
                open( '_input_processed/bounding_channels.p', "wb" ) )
    pickle.dump(interior_channels,
                open( '_input_processed/interior_channels.p', "wb" ) )
    
else:
    
    bounding_channels = pickle.load(
                        open( '_input_processed/bounding_channels.p', "rb" ))
    interior_channels = pickle.load(
                        open( '_input_processed/interior_channels.p', "rb" ))



''' CALCULATE BASIC GEOMETRIC METRICS '''

print 'Calculating geometries'

interior_lengths = [sum([network_lines[j].length for j in interior_channels[i]]) if len(interior_channels[i])>0 else 0 for i in range(len(islands))] 

perimeter = np.array([i.boundary.length for i in islands])
wetted_perimeter = perimeter + 2 * np.array(interior_lengths)   
perimeter_convex_hull = np.array([i.convex_hull.exterior.length for i in islands])

area = np.array([i.area for i in islands])
area_convex_hull = np.array([i.convex_hull.area for i in islands])

a = np.array(map(Polygon_axes, islands))
major_axis = a[:,1]
minor_axis = a[:,0]
aspect_ratio = major_axis / minor_axis

circularity = 4 * np.pi * area / perimeter**2
solidity = area / area_convex_hull
concavity = area_convex_hull - area
convexity = perimeter_convex_hull / perimeter
dry_shape_factor = perimeter / np.sqrt(area)
wet_shape_factor = wetted_perimeter / np.sqrt(area)

polygon_metrics = {'Area': area,
                'Perimeter': perimeter,
                'WetPerim': wetted_perimeter,
                'CH_Area': area_convex_hull,
                'CH_Perim': perimeter_convex_hull,
                'AspectR': aspect_ratio,
                'Circular': circularity,
                'Solidity': solidity,
                'Concavity': concavity,
                'Convexity': convexity,
                'DryShapeF': dry_shape_factor,
                'WetShapeF': wet_shape_factor}

pickle.dump(polygon_metrics,
            open( '_metrics/metrics__base_metrics.p', "wb" ) ) 





''' CALCULATE MAXIMUM DISTANCE FROM ANY WATER BODY '''

if calculate_params:

    print 'Calculating maximum distance from water bodies'

    maximum_edge_distance = np.zeros((len(islands),))
    cellsize = 30

    for n,i in enumerate(islands):

        print n

        minx, miny, maxx, maxy = i.bounds

        minx = np.floor(minx) - 1 * cellsize
        maxx = np.ceil(maxx) + 1 * cellsize
        miny = np.floor(miny) - 1 * cellsize
        maxy = np.ceil(maxy) + 1 * cellsize

        x = np.arange(minx, maxx , cellsize)
        y = np.arange(miny, maxy , cellsize)

        mask = outline_to_mask(i.exterior, x, y)
        distmap = morphology.distance_transform_edt(mask)

        maximum_edge_distance[n] = distmap.max() * cellsize

        mask = dist = None

    pickle.dump(maximum_edge_distance,
                open( '_metrics/metrics__edge_distance.p', "wb" ) )

else:
    
    maximum_edge_distance = pickle.load(
                            open( '_metrics/metrics__edge_distance.p', "rb" ))   

    
polygon_metrics['EdgeDist'] = maximum_edge_distance
    
    
    
    
    

''' CALCULATE BOUNDING CHANNEL WIDTH STATISTICS '''

if calculate_params:

    print 'Calculating statistics of bounding channels'

    network_min_widths = np.zeros((len(islands),))
    network_avg_widths = np.zeros((len(islands),))
    network_max_widths = np.zeros((len(islands),))

    for n in range(len(islands)):

        i = islands[n]

        # network lines off coast (for closing patches) have width 9999 - ignore
        channels = [network_lines[b] for b in bounding_channels[n] if network_widths[b] != 9999]
        widths = [network_widths[b] for b in bounding_channels[n] if network_widths[b] != 9999]
        lengths = [c.length for c in channels]


        tot_length = sum(lengths)
        network_avg_widths[n] = sum([widths[b] * lengths[b] for b in range(len(widths))]) / tot_length
        network_max_widths[n] = max(widths)
        network_min_widths[n] = min(widths)

    pickle.dump([network_min_widths, network_avg_widths, network_max_widths],
                open( '_metrics/metrics__bounding_channel_widths.p', "wb" ) )

else:
        
    network_min_widths, network_avg_widths, network_max_widths = pickle.load(
                    open( '_metrics/metrics__bounding_channel_widths.p', "rb" ))
    
polygon_metrics['Max_Width'] = network_max_widths
polygon_metrics['Min_Width'] = network_min_widths
polygon_metrics['Avg_Width'] = network_avg_widths




''' FIND NUMBER OF OUTFLOW CHANNELS

Count the number of outflow channels by counting the number of interior
channels that cross the boundary of the island.
We use patches instead of islands to avoid false positives where wide
channels are preserved in the outline of islands.
'''

if calculate_params:

    print 'Count number of outflow channels'

    num_outflow = np.zeros((len(islands),), dtype = 'int')

    for n in range(len(islands)):

        lines = [i for i in interior_channels[n]]

        # interior channels that touch the boundary
        outflow = []

        for l in lines:
            if network_lines[l].intersects(patches[n].exterior):
                outflow.append(l)

        num_outflow[n] = len(outflow)

    pickle.dump(num_outflow,
                open( '_metrics/metrics__outflow_channels.p', "wb" ) )

else:
    
    num_outflow = pickle.load(
                        open( '_metrics/metrics__outflow_channels.p', "rb" ))
    
polygon_metrics['NumOutflow'] = num_outflow






''' MEASURE THE FRACTAL DIMENSION OF BOUNDING CHANNELS '''


def fractal_dimension(Z, threshold=0.9):

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
                               np.arange(0, Z.shape[1], k), axis=1)

        # Count non-empty (0) and non-full boxes (k*k)
        return len(np.where((S > 0) & (S < k*k))[0])

    p = min(Z.shape)
    n = 2**np.floor(np.log(p)/np.log(2))
    n = int(np.log(n)/np.log(2))
    sizes = 2**np.arange(n, 1, -1)

    # Box counting with decreasing size
    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    # Fit the successive log(sizes) with log (counts)
    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    return -coeffs[0]



if calculate_params:

    print 'Calculating fractal dimension'

    fractal_dimensions = np.zeros((len(islands),))

    for j in range(len(islands)):

        print j

        outline = patches[j].exterior
        _,_,angle = Polygon_axes(outline)
        outline = affinity.rotate(outline, angle, origin='centroid')

        minx, miny, maxx, maxy = outline.bounds

        cellsize = 5

        if (maxy - miny > 10000) or (maxx - minx > 10000):
            cellsize = 30

        if (maxy - miny > 50000) or (maxx - minx > 50000):
            cellsize = 60

        minx = np.floor(minx) - 1 * cellsize
        maxx = np.ceil(maxx) + 1 * cellsize
        miny = np.floor(miny) - 1 * cellsize
        maxy = np.ceil(maxy) + 1 * cellsize

        x = np.arange(minx, maxx , cellsize)
        y = np.arange(miny, maxy , cellsize)

        mask = outline_to_mask(outline, x, y)
        fractal_dimensions[j] = fractal_dimension(mask)

        mask = None

    pickle.dump(fractal_dimensions,
                open( '_metrics/metrics__fractal_dimensions.p', "wb" ) )    
    
else:
    
    fractal_dimensions = pickle.load(
                        open( '_metrics/metrics__fractal_dimensions.p', "rb" ))
    
    
polygon_metrics['FractalD'] = fractal_dimensions







''' PICKLE PARAMETERS '''

pickle.dump(polygon_metrics,
            open( '_metrics/metrics__all_metrics.p', "wb" ) )





''' SAVE PARAMETERS TO SHAPEFILE '''

field_type = {}

for k in polygon_metrics.keys():
    
    if polygon_metrics[k][0].dtype == 'float':
        field_type[k] = ogr.OFTReal
        
    if polygon_metrics[k][0].dtype == 'int':
        field_type[k] = ogr.OFTInteger

        
create_shapefile_from_shapely_multi(islands,
                            '_output/islands_properties.shp',
                            fields = polygon_metrics,
                            field_type = field_type)

