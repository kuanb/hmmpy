import copy

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString

DEFAULT_TRACE_NODE_RADIUS = 100  # meters


def graph_edges_to_gdf(G):
    rows = []
    for node_fr, node_to, data in G.edges(data=True):
        if 'geometry' in data:
            line = data['geometry']
        else:
            # if it doesn't have a geometry attribute, the edge is a straight
            # line from node to node
            x1 = G.nodes[node_fr]['x']
            y1 = G.nodes[node_fr]['y']
            x2 = G.nodes[node_to]['x']
            y2 = G.nodes[node_to]['y']
            line = LineString([(x1, y1), (x2, y2)])

        rows.append({
            'from': node_fr,
            'to': node_to,
            'geometry': line,
        })

    gdf = gpd.GeoDataFrame(rows)
    gdf = gdf.drop_duplicates(['from', 'to'])
    return gdf.set_index(['from', 'to'])


def graph_nodes_to_gdf(G):
    rows = []
    for node_id, data in G.nodes(data=True):
        rows.append({
            'id': node_id,
            'geometry': Point(data['x'], data['y']),
        })

    gdf = gpd.GeoDataFrame(rows)
    gdf = gdf.drop_duplicates(['id'])
    return gdf.set_index('id')


def emission_probability(gps_point, road_edge):
    point_distance = road_edge.project(gps_point)
    interpolated_point = road_edge.interpolate(point_distance)
    dist = interpolated_point.distance(gps_point)
    return dist ** 2


def generate_phase_candidates(
        graph_gdf: gpd.GeoDataFrame,
        trace_linestring: LineString,
        search_radius: int=DEFAULT_TRACE_NODE_RADIUS
    ):
    """Create an array of arrays, each containing the candidate edges related to a given GPS trace point.
    
    Parameters:
        graph_gdf: GeoDataFrame containing geometries of edges of graph
        trace_linestring: shape object of GPS path
        search_radius: distance in meters to search within for candidate edges from a GPS point
    
    Returns:
        phase_candidates: an array containing arrays of candidate edgees for each GPS trace point
    """
    # initialize results candidates list
    phase_candidates = []

    # extract Rtree spatial index for downstream reference
    spatial_index = graph_gdf.sindex

    # iterate through each point in trace linestring
    xys = list(zip(*trace_linestring.xy))
    for i, trace_xy in enumerate(xys):
        p = Point(*trace_xy)
        area = p.buffer(DEFAULT_TRACE_NODE_RADIUS)

        # utilize spatial index for faster intersections
        possible_matches_index = list(spatial_index.intersection(area.bounds))
        possible_matches = graph_gdf.iloc[possible_matches_index]
        precise_matches = possible_matches[possible_matches.intersects(area)]
        
        emission_probability = []
        for i, row in precise_matches.iterrows():
            line_geom = row.geometry
            node_fr_id, node_to_id = row.name
            emission_probability.append(
                emission_probability(p, line_geom))

        edge_contents = []
        for edge_id, shape, emission_prob in zip(
            precise_matches.index.tolist(),
            precise_matches.geometry.tolist(),
            emission_probability):
            edge_contents.append({
                'edge_id': edge_id,
                'shape': shape,
                'emission_probability': emission_prob,
            })
        
            
        # add to queue; one set of candidates per trace point
        phase_candidates.append({
            'stage': i,
            'point': p,
            'is_last_point': i == (len(xys) - 1),
            'edges': edge_contents,
        })
    
    return phase_candidates


def make_trellis_node_id(stage, node_fr, node_to):
    return 'stage={}&fr={}&to={}'.format(stage, node_fr, node_to)

def read_trellis_node_id(node_id):
    a, b, c = node_id.split('&')
    stage = a.replace('stage=', '')
    node_fr = b.replace('fr=', '')
    node_to = c.replace('to=', '')
    return (stage, node_fr, node_to)


def generate_trellis_graph(phase_candidates: list):
    # initalize trellis
    trellis = nx.DiGraph()

    for phase_candidate in phase_candidates:
        stage = phase_candidate['stage']
        gps_point = phase_candidate['point']
        for node_data in phase_candidate['edges']:
            node_id = make_trellis_node_id(stage, *node_data['edge_id'])
            
            trellis.add_node(
                node_id,
                graph_node_fr=node_data['edge_id'][0],
                graph_node_to=node_data['edge_id'][1],
                gps_point=gps_point,
                edge_shape=node_data['shape'],
                emission_probability=node_data['emission_probability'])
        
    for candidates_fr, candidates_to in zip(phase_candidates[:-1], phase_candidates[1:]):
        stage_fr = candidates_fr['stage']
        stage_to = candidates_to['stage']

        # add edges for all "from-to" pairs between each phase to the next
        for nodes_data_cf in candidates_fr['edges']:
            nodes_cf = nodes_data_cf['edge_id']
            node_from_id = make_trellis_node_id(stage_fr, *nodes_cf)

            for nodes_data_to in candidates_to['edges']:
                nodes_ct = nodes_data_to['edge_id']
                node_to_id = make_trellis_node_id(stage_to, *nodes_ct)
                
                trellis.add_edge(node_from_id, node_to_id)
    
    return trellis


def transition_probability(
        road_network_graph,
        source_id,
        dest_id,
        gps_point_fr,
        gps_point_to,
        beta=DEFAULT_BETA):

    network_distance = nx.shortest_path_length(
        road_network_graph,
        source=source_id,
        target=dest_id,
        weight='length')
    dist_between_points = gps_point_fr.distance(gps_point_to)

    return network_distance**2 / dist_between_points**2


def path_probability(
        road_network_graph,
        trellis,
        trellis_node_fr,
        trellis_node_to,
        is_last_point):

    node_fr = trellis.nodes[trellis_node_fr]
    node_to = trellis.nodes[trellis_node_to]

    start_edge = node_fr['edge_shape']
    end_edge = node_to['edge_shape']
    
    emission_weight = node_fr['emission_probability']
    
    if is_last_point:
        emission_weight += node_to['emission_probability']
        
    transition_weight = transition_probability(
        road_network_graph,
        node_fr['graph_node_fr'],
        node_to['graph_node_to'],
        node_fr['gps_point'],
        node_to['gps_point'])

    # negative weight so log-likelihoods are non-negative in graph search
    return -transition_weight + -emission_weight


