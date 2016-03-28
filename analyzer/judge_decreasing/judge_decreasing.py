import networkx as nx
import json
import sys

g = nx.read_weighted_edgelist( sys.argv[1] )
assortativity = nx.degree_assortativity_coefficient(g)

# calculate degree histo
degrees = g.degree().values()
degree_set = sorted( set(degrees) )
hist = { x:list(degrees).count(x) for x in degree_set }
ratio_1 = hist[1] / max( hist.values() )
max_at = max( hist, key = lambda i : hist[i] )

output = { "assortativity": assortativity, "p1_ratio": ratio_1, "pk_max_at": max_at}
print( json.dumps( output ) )

