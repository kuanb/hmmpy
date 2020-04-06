import networkx as nx

from .trellis import *
from .gen_utils import *

G1, trace = make_graph_and_trace()


# def match(graph: nx.DiGraph, trace: LineString):

road_network_graph = G1.copy()

# take networkx graph and convert to reference geodataframe
edges_gdf = graph_edges_to_gdf(road_network_graph)
nodes_gdf = graph_nodes_to_gdf(road_network_graph)

phase_candidates = generate_phase_candidates(edges_gdf, trace)
trellis = generate_trellis_graph(phase_candidates)

start_node_ids = []
for phase_candidate in phase_candidates:
    stage = phase_candidate['stage']
    gps_point = phase_candidate['point']
    for node_data in phase_candidate['edges']:
        start_node_ids.append(make_trellis_node_id(stage, *node_data['edge_id']))

paths = []
for start_node_id in start_node_ids:
    try:
        nx.single_source_dijkstra(
            trellis,
            start_node_id,
            weight=lambda start, end, _: path_probability(road_network_graph, trellis, start, end, False))
    except nx.NetworkXNoPath:
        paths.append([])

paths