"""Combinatorial type of a polyhedron given by face cycles: the vertex-edge-face incidence graph, compared by isomorphism."""
import networkx as nx
def incidence(face_cycles):
    G = nx.Graph(); E = {}
    for f, cyc in enumerate(face_cycles):
        G.add_node(('f', f), kind='f')
        for i in range(len(cyc)):
            a, b = cyc[i], cyc[(i + 1) % len(cyc)]
            e = E.setdefault(frozenset((a, b)), len(E))
            G.add_node(('v', a), kind='v'); G.add_node(('v', b), kind='v'); G.add_node(('e', e), kind='e')
            G.add_edge(('v', a), ('e', e)); G.add_edge(('v', b), ('e', e))            # an edge is incident with BOTH endpoints
            G.add_edge(('e', e), ('f', f)); G.add_edge(('v', a), ('f', f))
    return G
def same_type(c1, c2):
    return nx.is_isomorphic(incidence(c1), incidence(c2), node_match=lambda x, y: x['kind'] == y['kind'])
def doubled_graph(face_cycles):
    """faces as nodes, an edge for every pair of faces sharing two edges"""
    n = len(face_cycles); owner = {}
    for f, cyc in enumerate(face_cycles):
        for i in range(len(cyc)): owner.setdefault(frozenset((cyc[i], cyc[(i + 1) % len(cyc)])), []).append(f)
    cnt = {}
    for e, fs in owner.items(): cnt[frozenset(fs)] = cnt.get(frozenset(fs), 0) + 1
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(tuple(k) for k, v in cnt.items() if v == 2); return G
