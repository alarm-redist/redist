#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, unicode_literals,\
    absolute_import, division

import sys
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_node_cut

def minmax(e0, e1):
#    return (min(e[0], e[1]), max(e[0], e[1]))
    return (min(e0, e1), max(e0, e1))

def get_farthest_two_vertices(G):

    ecc = nx.eccentricity(G)

    max_dist = nx.diameter(G, ecc)

    s = [key for key in ecc.keys() if ecc[key] == max_dist][0]

    sp = nx.single_source_shortest_path_length(G, s)

    e = [key for key in sp.keys() if sp[key] == max_dist][0]

    return (s, e)

def choose_next(G, vertex_order_dict, order_v_list):
    scores = [0] * len(order_v_list)
    for i in range(len(order_v_list)):
        vs = list(G.neighbors(order_v_list[i]))
        count = 0
        for v in vs:
            if vertex_order_dict[v] >= 1:
                scores[i] += 1
    max_score = max(scores)
    n_scores = [9999999] * len(order_v_list)
    for i in range(len(order_v_list)):
        if scores[i] == max_score:
            n_scores[i] = 0
            vs = list(G.neighbors(order_v_list[i]))
            for v in vs:
                if vertex_order_dict[v] >= 1:
                    n_scores[i] += vertex_order_dict[v]
    min_n_score = min(n_scores)
    return order_v_list[n_scores.index(min_n_score)]

def order_by_NDS(G, vertex_order_dict, order_vertex_set, result_edge_list):
    order_v_list = list(order_vertex_set)

    v_count = 1
    for v in vertex_order_dict:
        if vertex_order_dict[v] >= 1:
            v_count += 1

    while len(order_v_list) > 0:
        v = choose_next(G, vertex_order_dict, order_v_list)
        order_v_list.remove(v)
        vertex_order_dict[v] = v_count
        v_count += 1
        vs = list(G.neighbors(v))
        for w in sorted(vs, key=(lambda x: vertex_order_dict[x])):
            for i in range(G.number_of_edges(v, w)): # coping with parallel edges
                if vertex_order_dict[w] >= 1:
                    result_edge_list.append((v, w))

def split_graph(G, s, t, left_vertex_set, right_vertex_set,
                cut_set,
                vertex_order_dict, result_edge_list):
    H = G.copy()

    #print "left =", left_vertex_set
    #print "right =", right_vertex_set
    #print "cut_set =", cut_set

    #if len(right_vertex_set) >= 10:
    #    exit(1)
    
    for v in left_vertex_set:
        if v != s and v != t:
            H = nx.contracted_nodes(H, s, v)
    for v in right_vertex_set:
        if v != s and v != t:
            H = nx.contracted_nodes(H, t, v)

    if len(H.nodes()) <= 10 or H.has_edge(s, t):
        order_vertex_set = set(G.nodes())
        order_vertex_set -= left_vertex_set
        order_vertex_set -= right_vertex_set
        order_vertex_set |= cut_set

        #print "order_by_NDS start"
        #print "order_vertex_set =", order_vertex_set
        order_by_NDS(G, vertex_order_dict, order_vertex_set,
                     result_edge_list)
        #print "vertex_order_dict =", vertex_order_dict
        #print "result_edge_list =", result_edge_list
        return

    cut = minimum_st_node_cut(H, s, t)
    #print "cut =", cut
    Hc = H.copy()
    Hc.remove_nodes_from(cut)

    cc_list = list(nx.connected_components(Hc))
    ccs = None
    for cc in cc_list:
        if s in cc:
            ccs = cc
            cc_list.remove(cc)
            break
    #print "cc_list =", cc_list
    #print "ccs =", ccs
    new_right_vertex_set = set()
    new_right_vertex_set |= right_vertex_set
    new_right_vertex_set |= cut
    for cc in cc_list:
        new_right_vertex_set |= cc

    #print "new_right_vertex_set =", new_right_vertex_set
    split_graph(G, s, t, left_vertex_set, new_right_vertex_set, cut,
                vertex_order_dict, result_edge_list)

    new_left_vertex_set = set()
    new_left_vertex_set |= left_vertex_set
    new_left_vertex_set |= ccs
    new_left_vertex_set |= cut

    #print "new_left_vertex_set =", new_left_vertex_set
    split_graph(G, s, t, new_left_vertex_set, right_vertex_set, cut_set,
                vertex_order_dict, result_edge_list)

# return (G, deg1_edges, deg2_dict, deg2_cycle, is_cycle, is_tree)
# is_cycle: whether G is a cycle
# is_tree: whether G is a tree
def remove_deg12(G):

    deg1_edges = []
    deg2_dict = {}
    deg2_cycle = {}

    found = True
    while found:
        found = False
        for n in G.nodes():
            if G.degree(n) == 1:
                ns = list(G.neighbors(n))
                G.remove_node(n)
                deg1_edges.append((ns[0], n))
                found = True
                break

    if len(G.edges()) == 0: # case for trees
        e = deg1_edges[-1]
        G.add_edge(*e)
        deg1_edges.remove(e)
        return (G, deg1_edges, deg2_dict, deg2_cycle, False, True)

    #print "d", deg1_edges
            
    found = True
    cycle = False
    while found:
        found = False
        for n in G.nodes():
            if G.degree(n) == 2:
                walk = [n]
                ns = list(G.neighbors(n))
                # print(n, ns)
                c = ns[0]
                walk.insert(0, c)
                while G.degree(c) == 2:
                    if n == c:
                        return (G, deg1_edges, deg2_dict, deg2_cycle, True, False)
                    nc = list(G.neighbors(c))
                    if nc[0] != walk[1]:
                        c = nc[0]
                    else:
                        c = nc[1]
                    walk.insert(0, c)
                c = ns[1]
                walk.append(c)
                while G.degree(c) == 2:
                    nc = list(G.neighbors(c))
                    if nc[0] != walk[-2]:
                        c = nc[0]
                    else:
                        c = nc[1]
                    walk.append(c)
                #print walk

                G.remove_nodes_from(walk[1:-1])
                if walk[0] == walk[-1]: # cycle
                    deg2_cycle.setdefault(walk[0], []).append(walk)
                    #deg2_cycle[walk[0]] = walk
                else:
                    has_e = G.has_edge(walk[0], walk[-1])
                    G.add_edge(walk[0], walk[-1])
                    deg2_dict.setdefault(minmax(walk[0], walk[-1]), []).append((walk, has_e))
                #print G.edges()
                found = True
                break

    return (G, deg1_edges, deg2_dict, deg2_cycle, False, False)

def get_cycle(edge_list):
    result_edge_list = []
    v = edge_list[0][0]
    prev_v = -1
    while len(result_edge_list) < len(edge_list):
        vs = [x for x in edge_list if x[0] == v or x[1] == v]
        next_v = (set(vs[0]).union(set(vs[1])) - set([prev_v]) - set([v])).pop()
        result_edge_list.append(minmax(v, next_v))
        prev_v = v
        v = next_v
    return result_edge_list

def check_connected_order(edge_list):
    touched_vertices = set()
    for i in range(len(edge_list)):
        if i > 0 and (edge_list[i][0] not in touched_vertices) and (edge_list[i][1] not in touched_vertices):
            #sys.stderr.write(str(i) + " " + str(edge_list[i][0]) + " " + str(edge_list[i][1]))
            return False
        touched_vertices.add(edge_list[i][0])
        touched_vertices.add(edge_list[i][1])
    return True

def recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list):

    #print deg1_edges
    #print deg2_dict
    #print deg2_cycle
    found = True
    while found:
        #print "res ", result_edge_list
        found = False
        for i in range(len(result_edge_list)):
            e = minmax(result_edge_list[i][0], result_edge_list[i][1])
            if e in deg2_dict:
                found = True
                #print "recover ", e
                for x in deg2_dict[e]:
                    if e in result_edge_list:
                        result_edge_list.remove(e)
                    else:
                        result_edge_list.remove((e[1], e[0]))
                    touched_vertices = [y for l in result_edge_list[0:i] for y in l] # flatten
                    #print touched_vertices
                    if x[0][0] in touched_vertices:
                        for j in reversed(range(len(x[0]) - 1)):
                            result_edge_list.insert(i, (x[0][j], x[0][j + 1]))
                    else: #elif x[0][-1] in touched_vertices:
                        for j in range(len(x[0]) - 1):
                            result_edge_list.insert(i, (x[0][j], x[0][j + 1]))
                    #else:
                    #    print x
                    #    raise Exception("graph disconnected?")
                deg2_dict.pop(e)
                break # for i
        if not found and len(deg2_cycle) > 0:
            for v in deg2_cycle:
                for i in range(len(result_edge_list)):
                    if v in result_edge_list[i]:
                        #print "recover ", v
                        found = True
                        for cy in deg2_cycle[v]:
                            for j in range(len(cy) - 1):
                                c0 = cy[j]
                                c1 = cy[j + 1]
                                #print "insert", c0, c1
                                result_edge_list.insert(i + 1, minmax(c0, c1))
                        deg2_cycle.pop(v)
                        break # for i
                if found:
                    break # for v

    #print result_edge_list

    found = True
    while len(deg1_edges) > 0:
        found = False
        for e in deg1_edges:
            for i in range(len(result_edge_list)):
                if e[0] in result_edge_list[i]:
                    result_edge_list.insert(i + 1, e)
                    found = True
                    deg1_edges.remove(e)
                    break
            if found:
                break
        if not found:
            print("not found error!")
            exit(1)

def get_order_by_cut(edge_list):

    G = nx.MultiGraph()
    G.add_edges_from(edge_list)

    (H, deg1_edges, deg2_dict, deg2_cycle, is_cycle, is_tree) = remove_deg12(G)

    if is_cycle:
        #print "edges ", G.edges()
        result_edge_list = get_cycle(list(G.edges()))
        recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
    elif is_tree:
        result_edge_list = list(G.edges())
        recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
    else:
        G = H
        s, t = get_farthest_two_vertices(G) #print "s, t = ", s, t

        result_edge_list = []
        vertex_order_dict = {}
        for v in G.nodes():
            vertex_order_dict[v] = -1

        vertex_order_dict[s] = 1
        left = set([s])
        right = set([t])
        cut_set = right.copy()
        split_graph(G, s, t, left, right, cut_set, vertex_order_dict, result_edge_list)
        #order_vertex_list = G.nodes()
        #order_vertex_list.remove(left_vertex_set[0])
        #order_by_NDS(G, vertex_order_dict, order_vertex_list, result_edge_list

        #return split_graph(G, [s], edge_list, [t])

        recover_deg12(deg1_edges, deg2_dict, deg2_cycle, result_edge_list)
        
    return result_edge_list

def get_order_by_cut_with_check(edge_list):

    result_edge_list = get_order_by_cut(edge_list)

    # check start
    G1 = nx.MultiGraph()
    G2 = nx.MultiGraph()

    #for e in edge_list:
    #    G1.add_edge(*minmax(e[0], e[1]))
    #for e in result_edge_list:
    #    G2.add_edge(*minmax(e[0], e[1]))
    
    G1.add_edges_from(edge_list)
    G2.add_edges_from(result_edge_list)

    if not nx.is_isomorphic(G1, G2):
        sys.stderr.write("not isomorphic!")
        #print "G1 =", edge_list
        #print "G2 =", result_edge_list
        print("G1 =", G1.edges())
        print("G2 =", G2.edges())

        s1 = []
        s2 = []
        for e in edge_list:
            s1.append(minmax(e[0], e[1]))
        for e in result_edge_list:
            s2.append(minmax(e[0], e[1]))
        ss1 = set(s1)
        ss2 = set(s2)

        print("ss1 - ss2 =", (ss1 - ss2))
        print("ss2 - ss1 =", (ss2 - ss1))

        exit(1)

    #for e in result_edge_list:
    #    c = minmax(e[0], e[1])
    #    print c[0], c[1]
    if not check_connected_order(result_edge_list):
        sys.stderr.write("not connected order")
        exit(1)
    # check end

    return result_edge_list


if __name__ == '__main__':

    edge_list = []

    for line in sys.stdin:
        ar = line.strip().split()
        if len(ar) < 2:
            sys.stderr.write("Each edge must have two vertices! " + line.strip() +
                " does not.")
            exit(1)
        edge_list.append((int(ar[0]), int(ar[1])))

    if len(edge_list) == 0:
        print("The input graph is empty.")
        exit(1)

    new_edge_list = get_order_by_cut(edge_list)

    G1 = nx.MultiGraph()
    G2 = nx.MultiGraph()

    G1.add_edges_from(edge_list)
    G2.add_edges_from(new_edge_list)

    if not nx.is_isomorphic(G1, G2):
        sys.stderr.write("not isomorphic!")
        #exit(1)

    for e in new_edge_list:
        c = minmax(e[0], e[1])
        print(c[0], c[1])
    if not check_connected_order(new_edge_list):
        sys.stderr.write("not connected order")
