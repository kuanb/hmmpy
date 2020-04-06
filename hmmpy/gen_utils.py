from .trellis import *

import osmnx as ox
from shapely.geometry import Polygon, box

from pyproj import CRS
from pyproj import Transformer
from shapely.geometry import LineString
from shapely.ops import transform

def make_graph_and_trace():
	# north, south, east, west = 37.79, 37.78, -122.41, -122.43
	p = Polygon([[-122.27465629577637,37.7983230783235],[-122.26096630096434,37.7983230783235],[-122.26096630096434,37.80761398306056],[-122.27465629577637,37.80761398306056],[-122.27465629577637,37.7983230783235]])
	west, south, east, north = p.bounds
	G1 = ox.graph_from_bbox(north, south, east, west)

	# make meter projected
	G1 = ox.project_graph(G1)

	print('Converting to this projection:', G1.graph['crs'])
	project_to_meter = Transformer.from_crs(
	    CRS.from_epsg(4326), G1.graph['crs'], always_xy=True
	)

	# create trace path and reproject
	trace = LineString([[-122.2701072692871,37.798865645037765],[-122.26969957351685,37.79971339755073],[-122.26935625076295,37.800222044388384],[-122.26892709732056,37.80079850657009],[-122.26845502853392,37.80135800967933],[-122.26856231689453,37.801646236899785],[-122.26939916610718,37.80208705282608],[-122.27008581161499,37.80235832285767],[-122.27113723754883,37.8026635004523]])
	trace = transform(project_to_meter.transform, LineString(trace))

	return G1, trace