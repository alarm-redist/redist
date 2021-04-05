/*
 * enumpart
 * Copyright (c) 2018 -- 2019 Jun Kawahara
 */

#include <climits>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <random>
#include <array>
#include <algorithm>
#include <functional>

#include <tdzdd/DdSpecOp.hpp>
#include <tdzdd/DdStructure.hpp>

#include "GraphPartitionSpec.hpp"
#include "GraphPartitionSpecL31.hpp"
#include "ComponentWeightSpec.hpp"
#include "ComponentRatioSpec.hpp"
#include "Graph.hpp"
#include "testGP.hpp"

#include "SAPPOROBDD/ZBDD.h"
#include <tdzdd/spec/SapporoZdd.hpp>
#include <tdzdd/eval/ToZBDD.hpp>

#include "SBDD_helper.h"

using namespace sbddh;

#include "BigInteger.hpp"
#include "RandomSample.hpp"

using namespace tdzdd;

std::string options[][2] = { //
        {"nola", "Do not use lookahead"}, //
        {"a", "Read <graph_file> as an adjacency list"}, //
        {"k <n>", "<n> components"}, //
        {"lower <n>", "<n> lower bound"}, //
        {"upper <n>", "<n> upper bound"}, //
        {"ratio <f>", "<f> ratio"}, //
        {"noloop", "no loop"}, //
        {"noreduce", "not reduce result ZDD"}, //
        {"L31", "Optimization when frontier size is less than 31"}, //
        {"M32", "Use the algorithm for frontier size at least 32"}, //
        {"test", "Test this program"}, //
        {"count", "Report the number of solutions"}, //
        {"graph", "Dump input graph to STDOUT in DOT format"}, //
        {"allsols", "Output all solutions in the edge list format"}, //
        {"solutions <n>", "Output the first <n> solutions in the edge list format"}, //
        {"comp", "Make output format component number"}, //
        {"drawsol <n>", "Output the <n>-th solution to STDOUT in DOT format"}, //
        {"sample <n>", "Output <n> solutions uniformly and randomly"}, //
        {"sample-old <n>", "Output <n> solutions uniformly and randomly (old version)"}, //
        {"zdd", "Dump result ZDD to STDOUT in DOT format"}, //
        {"export", "Dump result ZDD to STDOUT"}, //
        {"readfile", "Read ZDD from a file instead of running the algorithm"}};

std::map<std::string,bool> opt;
std::map<std::string,int> optNum;
std::map<std::string,double> optDouble;
std::map<std::string,std::string> optStr;

void usage(char const* cmd) {
    std::cerr << "usage: " << cmd
              << " [ <option>... ] [ <graph_file> [ <vertex_group_file> ]]\n";
    std::cerr << "options\n";
    for (unsigned i = 0; i < sizeof(options) / sizeof(options[0]); ++i) {
        std::cerr << "  -" << options[i][0];
        for (unsigned j = options[i][0].length(); j < 10; ++j) {
            std::cerr << " ";
        }
        std::cerr << ": " << options[i][1] << "\n";
    }
}

class EdgeDecorator {
    int const n;
    std::set<int> const& levels;

public:
    EdgeDecorator(int n, std::set<int> const& levels) :
            n(n), levels(levels) {
    }

    std::string operator()(Graph::EdgeNumber a) const {
        return levels.count(n - a) ?
                "[style=bold]" : "[style=dotted,color=gray]";
    }
};


BigInteger computeMap(const DdStructure<2>& dd, NodeId node, std::map<NodeId, BigInteger>& solutionMap)
{
    if (solutionMap.find(node) != solutionMap.end()) { // found
        return solutionMap[node];
    } else { // not found
        solutionMap[node] = computeMap(dd, dd.child(node, 0), solutionMap)
            + computeMap(dd, dd.child(node, 1), solutionMap);
        return solutionMap[node];
    }
}

void Replace(std::vector<int>& vec, int n, int src, int dest)
{
    for (int i = 1; i <= n; ++i) {
        if (vec[i] == src) {
            vec[i] = dest;
        }
    }
}

int translateVertex(int v, const Graph& g)
{
    std::string s = g.vertexName(v);
    std::istringstream iss(s);
    int vn;
    iss >> vn;
    return vn;
}

/*
 * Check whether partition meets weight bounds
 * Author: Cory McCartan
 */
bool check_print(const std::set<int>& s, const Graph& g, int k,
                 const std::vector<unsigned int> weights, int lower, int upper) {
    // convert to list of assignments (copied from OutputDistrict)
    std::vector<int> vec;
    for (std::set<int>::const_reverse_iterator itor = s.rbegin();
         itor != s.rend(); ++itor) {
        vec.push_back(g.edgeSize() - *itor); // 0-origin edge list
    }

    const int n = g.vertexSize();;
    std::vector<int> verarray(n + 1);

    for (int i = 1; i <= n; ++i) {
        verarray[i] = i + n;
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        const Graph::EdgeInfo& edge = g.edgeInfo(vec[i]);
        int c1 = verarray[translateVertex(edge.v1, g)];
        int c2 = verarray[translateVertex(edge.v2, g)];
        if (c1 > c2) {
            std::swap(c1, c2);
        }
        Replace(verarray, n, c2, c1);
    }

    // Aggregate by component
    std::vector<unsigned int> aggr(k, 0);
    int count = 0;
    for (int i = 1; i <= n; ++i) {
        if (verarray[i] > n) {
            Replace(verarray, n, verarray[i], count);
            ++count;
        }
        aggr[verarray[i]] += weights[i-1];
        // early exit
        if (aggr[verarray[i]] > upper) return false;
    }

    // check bounds
    for (int i = 0; i < k; i++) {
        if (aggr[i] < lower || aggr[i] > upper) return false;
    }

    // if reached here we're OK
    for (int i = 1; i <= n; ++i) {
        std::cout << verarray[i];
        if (i < n) {
            std::cout << " ";
        }
    }
    std::cout << std::endl;

    return true;
}


void OutputDistrict(const std::vector<int>& vec, const Graph& graph)
{
    const int n = graph.vertexSize();;
    std::vector<int> verarray(n + 1);

    for (int i = 1; i <= n; ++i) {
        verarray[i] = i + n;
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        const Graph::EdgeInfo& edge = graph.edgeInfo(vec[i]);
        int c1 = verarray[translateVertex(edge.v1, graph)];
        int c2 = verarray[translateVertex(edge.v2, graph)];
        if (c1 > c2) {
            std::swap(c1, c2);
        }
        Replace(verarray, n, c2, c1);
    }

    int count = 0;
    for (int i = 1; i <= n; ++i) {
        if (verarray[i] > n) {
            Replace(verarray, n, verarray[i], count);
            ++count;
        }
        std::cout << verarray[i];
        if (i < n) {
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}

// comp_format: true -> component number format,
//              false -> edge number format
void PrintEdgeSet(const std::set<int>& s, const Graph& g, bool comp_format)
{
    if (comp_format) {
        std::vector<int> vec;
        for (std::set<int>::const_reverse_iterator itor = s.rbegin();
             itor != s.rend(); ++itor) {
            vec.push_back(g.edgeSize() - *itor); // 0-origin edge list
        }
        OutputDistrict(vec, g);
    } else {
        for (std::set<int>::const_reverse_iterator itor = s.rbegin();
             itor != s.rend(); ++itor) {
            std::cout << (g.edgeSize() - *itor + 1) << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[]) {

    srand(time(NULL));

    for (unsigned i = 0; i < sizeof(options) / sizeof(options[0]); ++i) {
        opt[options[i][0]] = false;
    }

    std::string graphFileName;
    std::string zddFileName;
    std::string weightFileName;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string s = argv[i];
            if (s[0] == '-') {
                s = s.substr(1);

                if (opt.count(s)) {
                    opt[s] = true;
                }
                else if (i + 1 < argc && opt.count(s + " <n>")) {
                    opt[s] = true;
                    optNum[s] = std::atoi(argv[++i]);
                }
                else if (i + 1 < argc && opt.count(s + " <f>")) {
                    opt[s] = true;
                    optDouble[s] = std::atof(argv[++i]);
                }
                else if (i + 1 < argc && opt.count(s + " " + argv[i + 1])) {
                    opt[s] = true;
                    optStr[s] = argv[++i];
                }
                else {
                    throw std::exception();
                }
            }
            else if (graphFileName.empty()) {
                graphFileName = s;
            }
            else if (zddFileName.empty()) {
                zddFileName = s;
            }
            else if (weightFileName.empty()) {
                weightFileName = s;
            }
            else {
                throw std::exception();
            }
        }
    }
    catch (std::exception& e) {
        usage(argv[0]);
        return 1;
    }

    MessageHandler::showMessages();
    MessageHandler mh;
    mh.begin("started");

    if (opt["test"]) {
        testGP();
        testGPL31();
        std::cout << "test passed!" << std::endl;
        return 0;
    }

    Graph g;
    try {
        if (!graphFileName.empty()) {
            if (opt["a"]) {
                g.readAdjacencyList(graphFileName);
            }
            else {
                g.readEdges(graphFileName);
            }
        }
        else {
            g.addEdge("v1", "v2");
            g.addEdge("v1", "v3");
            g.addEdge("v1", "v4");
            g.addEdge("v2", "v4");
            g.addEdge("v2", "v5");
            g.addEdge("v3", "v4");
            g.addEdge("v3", "v6");
            g.addEdge("v4", "v5");
            g.addEdge("v4", "v6");
            g.addEdge("v4", "v7");
            g.addEdge("v5", "v7");
            g.addEdge("v6", "v7");
            g.update();
        }

        int const m = g.vertexSize();
        int const n = g.edgeSize();

        mh << "#vertex = " << m << ", #edge = " << n << "\n";

        if (g.edgeSize() == 0)
            throw std::runtime_error("ERROR: The graph is empty!!!");

        if (opt["graph"]) {
            g.dump(std::cout);
            return 0;
        }
        int k = (opt["k"] ? optNum["k"] : -1);

        // If we use "-comp", the vertex numbers should be 1,2,...,|V|.
        // We check it.
        if (opt["comp"]) {
            for (int v = 1; v <= m; ++v) {
                int vn = translateVertex(v, g);
                if (!(1 <= vn && vn <= m)) {
                    std::cerr << "Vertex numbers should be in {1,...,n}." << std::endl;
                    exit(1);
                }
            }
        }

        std::vector<unsigned int> weight_list;
        if (!opt["readfile"] && weightFileName.empty() && !zddFileName.empty())
            weightFileName = zddFileName;

        if (weightFileName.empty()) {
            for (int jj = 0; jj < m; ++jj) {
                weight_list.push_back(1);
            }
        } else {
            std::ifstream ifs(weightFileName.c_str());
            if (!ifs) {
                std::cerr << weightFileName << " not found!" << std::endl;
            }
            for (int jj = 0; jj < m; ++jj) {
                int c;
                ifs >> c;
                weight_list.push_back(c);
            }
        }
        // negative => not binding
        int lower = -2;
        int upper = -1;

        DdStructure<2> dd;
        ZBDD dd_s;
        bool dd_s_initialized = false;

        if (opt["readfile"]) { // read ZDD from file
            if (!dd_s_initialized) {
                BDD_Init(1024, 1024 * 1024 * 1024);
            }
            FILE* fp = fopen(zddFileName.c_str(), "r");
            if (fp == NULL) {
                std::cerr << "File " << zddFileName << " cannot be opened." << std::endl;
                exit(1);
            }
            dd_s = ZBDD_Import(fp);
            dd_s_initialized = true;
            fclose(fp);
            SapporoZdd szdd(dd_s);

            dd = DdStructure<2>(szdd);
        } else { // run frontier-based search

            if (opt["L31"]) {
                GraphPartitionSpecL31 gpspec(g, k, opt["noloop"], !opt["nola"]);
                dd = DdStructure<2>(gpspec);
            } else if (opt["M32"]) {
                GraphPartitionSpec gpspec(g, k, opt["noloop"], !opt["nola"], false);
                dd = DdStructure<2>(gpspec);
            } else {
                if (g.maxFrontierSize() < 31) {
                    GraphPartitionSpecL31 gpspec(g, k, opt["noloop"], !opt["nola"]);
                    dd = DdStructure<2>(gpspec);
                } else {
                    GraphPartitionSpec gpspec(g, k, opt["noloop"], !opt["nola"], false);
                    dd = DdStructure<2>(gpspec);
                }
            }

            if (!opt["noreduce"]) {
                dd.zddReduce();
            }

            if (opt["ratio"]) {
                int sum = 0;
                for (size_t i = 0; i < weight_list.size(); ++i) {
                    sum += weight_list[i];
                }
                double ratio = optDouble["ratio"];
                lower = static_cast<int>(floor(static_cast<double>(sum) /
                                        (ratio * (k - 1) + 1)));
                upper = static_cast<int>(ceil(ratio * static_cast<double>(sum) /
                                             (ratio + (k - 1))));
                std::cerr << "lower = " << lower << ", upper = " << upper << std::endl;
                //ComponentRatioSpec crspec(g, weight_list,
                //                          lower, upper,
                //                          ratio,
                //                          opt["noloop"], !opt["nola"]);

                //dd = DdStructure<2>(crspec);
                //dd.zddReduce();
                //dd.zddSubset(gpspec);

                //dd = DdStructure<2>(zddIntersection(gpspec, crspec));
                //dd.zddSubset(crspec);
                //dd.zddReduce();
            } else if (opt["lower"] || opt["upper"]) {
                lower = (opt["lower"] ? optNum["lower"] : 0);
                upper = (opt["upper"] ? optNum["upper"] : INT_MAX);
                //ComponentWeightSpec cwspec(g, weight_list,
                //                           lower, upper,
                //                           opt["noloop"], !opt["nola"]);
                //dd.zddSubset(cwspec);
                //dd.zddReduce();
            }
        }

        mh << "\n#node = " << dd.size() << ", #solution = "
                << std::setprecision(10)
                << dd.evaluate(ZddCardinality<double>())
                << "\n";

        mh << "lower = " << lower << ", upper = " << upper << "\n";

        if (opt["count"]) {
            MessageHandler mh;
            mh.begin("counting solutions") << " ...";
            mh << "\n#solution = " << dd.evaluate(ZddCardinality<>());
            mh.end();
        }

        if (opt["zdd"]) dd.dumpDot(std::cout, "ZDD");
        if (opt["export"]) dd.dumpSapporo(std::cout);

        bool check_flag = lower > 0 && upper > lower && opt["comp"];
        if (opt["solutions"] || opt["allsols"]) {
            int count = optNum["solutions"];
            if (opt["allsols"]) {
                count = 1;
            }

            for (DdStructure<2>::const_iterator t = dd.begin(); t != dd.end();
                    ++t) {
                // check bounds if necessary
                if (check_flag) {
                    bool ok = check_print(*t, g, k, weight_list, lower, upper);
                    if (!ok) continue;
                } else {
                    PrintEdgeSet(*t, g, opt["comp"]);
                }

                //EdgeDecorator edges(n, *t);
                //g.dump(std::cout, edges);
                if (opt["solutions"]) {
                    --count;
                }
                if (count == 0) break;
            }
        }
        if (opt["sample"]) {
            if (!dd_s_initialized) {
                BDD_Init(1024, 2199023255552ll); // max usage 2 TiB
                dd_s = dd.evaluate(ToZBDD());
                dd_s_initialized = true;
            }
            int const sampleNum = optNum["sample"];

            ZBDD_CountMap cmap;
            BigInteger bi = ZBDD_CountSolutions(dd_s, &cmap);

            BigIntegerRandom random;
            for (int i = 0; i < sampleNum; ++i) {
                std::set<bddvar> s = ZBDD_SampleRandomly(dd_s, cmap, random);
                // change unsigned -> signed
                std::set<int> r;
                std::set<bddvar>::iterator itor = s.begin();
                for ( ; itor != s.end(); ++itor) {
                    r.insert(*itor);
                }

                // check bounds if necessary
                if (check_flag) {
                    bool ok = check_print(r, g, k, weight_list, lower, upper);
                    if (!ok) i--;
                } else {
                    PrintEdgeSet(r, g, opt["comp"]);
                }
            }

        } else if (opt["sample-old"]) {
            int const sampleNum = optNum["sample"];

            std::map<NodeId, BigInteger> solutionMap;
            solutionMap[NodeId(0)] = BigInteger(0);
            solutionMap[NodeId(1)] = BigInteger(1);

            computeMap(dd, dd.root(), solutionMap);

            std::random_device rd;
            std::array<int, std::mt19937::state_size> seed_data;
            std::generate_n(seed_data.data(), seed_data.size(), std::ref(rd));
            std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
            std::mt19937 engine(seq);

            std::uniform_real_distribution<double> dis(0.0, 1.0);
            for (int i = 0; i < sampleNum; ++i) {
                std::set<int> edgeSet;
                NodeId node = dd.root();
                while (node != NodeId(0) && node != NodeId(1)) {
                    NodeId child0 = dd.child(node, 0);
                    NodeId child1 = dd.child(node, 1);
                    double p = dis(engine);
                    double th = (double)solutionMap[child0] /
                        (double)(solutionMap[child0] + solutionMap[child1]);
                    if (p < th) { // 0-child
                        node = child0;
                    } else { // 1-child
                        //edgeSet.push_back(n - node.row() + 1);
                        edgeSet.insert(node.row());
                        node = child1;
                    }
                }
                assert(node != NodeId(0));
                PrintEdgeSet(edgeSet, g, opt["comp"]);
                //std::cout << std::endl;
            }
        }
        if (opt["drawsol"]) {
            int const n = g.edgeSize();
            int num = optNum["drawsol"];

            for (DdStructure<2>::const_iterator t = dd.begin(); t != dd.end();
                    ++t) {
                --num;
                if (num > 0) { // output the <num>-th solution
                    continue;
                }
                EdgeDecorator edges(n, *t);
                g.dump(std::cout, edges);
                break;
            }
        }
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    mh.end("finished");
    return 0;
}
