#define test(b) testfunc((int64_t)(b), __LINE__)
#define test_eq(v1, v2) testfunc_eq((int64_t)(v1), (int64_t)(v2), __LINE__)

void testfunc(int64_t b, int error_line)
{
    if (b == 0) {
        fprintf(stderr, "not expected value at line %d\n", error_line);
        exit(1);
    }
}

void testfunc_eq(int64_t v1, int64_t v2, int error_line)
{
    if (v1 != v2) {
        fprintf(stderr, "%lld != %lld at line %d\n", v1, v2, error_line);
        exit(1);
    }
}

Graph vectorToGraph(int* edge_vector, int n)
{
    Graph g;
    char buff1[255];
    char buff2[255];
    for (size_t i = 0; i < n; i += 2) {
        sprintf(buff1, "%d", edge_vector[i]);
        sprintf(buff2, "%d", edge_vector[i + 1]);
        g.addEdge(buff1, buff2);
    }
    g.update();
    return g;
}

Graph makeSmallGraph()
{
    int ga[] = {1, 2, 1, 4, 2, 3, 2, 5, 3, 6, 4, 5, 5, 6};
    return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
}

Graph makeNotConnectedGraph()
{
    int ga[] = {1, 2, 1, 3, 2, 3, 3, 4, 4, 7, 5, 6};
    return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
}

Graph makeCliqueGraph(int n)
{
    Graph g;
    char buff1[255];
    char buff2[255];
    for (int i = 1; i <= n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            sprintf(buff1, "%d", i);
            sprintf(buff2, "%d", j);
            g.addEdge(buff1, buff2);
        }
    }
    g.update();
    return g;
}

Graph makePathGraph(int n)
{
    Graph g;
    char buff1[255];
    char buff2[255];
    for (int i = 1; i < n; ++i) {
        sprintf(buff1, "%d", i);
        sprintf(buff2, "%d", i + 1);
        g.addEdge(buff1, buff2);
    }
    g.update();
    return g;
}

Graph makeGridGraph(int n)
{
    if (n == 2) {
        int ga[] = {1, 2, 1, 3, 2, 4, 3, 4};
        return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
    } else if (n == 3) {
        int ga[] = {1, 2, 1, 4, 2, 3, 2, 5, 4, 5, 4, 7, 3, 6, 5, 6, 5, 8, 7, 8, 6, 9, 8, 9};
        return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
    } else if (n == 4) {
        int ga[] = {1, 2, 1, 5, 2, 3, 2, 6, 5, 6, 5, 9, 3, 4, 3, 7, 6, 7, 6, 10, 9, 10, 9,
                    13, 4, 8, 7, 8, 7, 11, 10, 11, 10, 14, 13, 14, 8, 12, 11, 12, 11, 15,
                    14, 15, 12, 16, 15, 16};
        return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
    } else if (n == 5) {
        int ga[] = {1, 2, 1, 6, 2, 3, 2, 7, 6, 7, 6, 11, 3, 8, 7, 8, 3, 4, 7, 12, 11, 12,
                    11, 16, 4, 9, 8, 9, 8, 13, 12, 13, 4, 5, 12, 17, 16, 17, 16, 21, 5, 10,
                    9, 10, 9, 14, 13, 14, 13, 18, 17, 18, 17, 22, 21, 22, 10, 15, 14, 15,
                    14, 19, 18, 19, 18, 23, 22, 23, 15, 20, 19, 20, 19, 24, 23, 24, 20,
                    25, 24, 25};
        return vectorToGraph(ga, sizeof(ga) / sizeof(ga[0]));
    } else {
        std::cerr << "not supported" << std::endl;
        exit(1);
        return Graph();
    }
}

void testGPL31() {
    {
        int answer[] = {1, 15, 29, 21, 7, 1, 0};
        for (int k = 1; k <= 7; ++k) {
            Graph G = makeSmallGraph();
            GraphPartitionSpecL31 gpspec(G, k, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {0, 1, 6, 13, 13, 6, 1, 0};
        for (int k = 1; k <= 8; ++k) {
            Graph G = makeNotConnectedGraph();
            GraphPartitionSpecL31 gpspec(G, k, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {1, 15, 25, 10, 1, 0};
        for (int k = 1; k <= 6; ++k) {
            Graph G = makeCliqueGraph(5);
            GraphPartitionSpecL31 gpspec(G, k, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {1, 5, 10, 10, 5, 1, 0};
        for (int k = 1; k <= 7; ++k) {
            Graph G = makePathGraph(6);
            GraphPartitionSpecL31 gpspec(G, k, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    int answer2[] = {6, 53, 627, 16213};
    int answer3[] = {4, 258, 10830, 709351};
    for (int k = 2; k <= 3; ++k) {
        for (int n = 2; n <= 5; ++n) {
            Graph G = makeGridGraph(n);
            GraphPartitionSpecL31 gpspec(G, k, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, (k == 2 ? answer2[n - 2] : answer3[n - 2]));
        }
    }
}

void testGP() {
    {
        int answer[] = {1, 15, 29, 21, 7, 1, 0};
        for (int k = 1; k <= 7; ++k) {
            Graph G = makeSmallGraph();
            GraphPartitionSpec gpspec(G, k, false, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {0, 1, 6, 13, 13, 6, 1, 0};
        for (int k = 1; k <= 7; ++k) {
            Graph G = makeNotConnectedGraph();
            GraphPartitionSpec gpspec(G, k, false, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {1, 15, 25, 10, 1, 0};
        for (int k = 1; k <= 6; ++k) {
            Graph G = makeCliqueGraph(5);
            GraphPartitionSpec gpspec(G, k, false, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    {
        int answer[] = {1, 5, 10, 10, 5, 1, 0};
        for (int k = 1; k <= 7; ++k) {
            Graph G = makePathGraph(6);
            GraphPartitionSpec gpspec(G, k, false, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, answer[k - 1]);
        }
    }

    int answer2[] = {6, 53, 627, 16213};
    int answer3[] = {4, 258, 10830, 709351};
    for (int k = 2; k <= 3; ++k) {
        for (int n = 2; n <= 5; ++n) {
            Graph G = makeGridGraph(n);
            GraphPartitionSpec gpspec(G, k, false, false, false);
            tdzdd::DdStructure<2> dd = tdzdd::DdStructure<2>(gpspec);
            int64_t sol = dd.evaluate(tdzdd::ZddCardinality<int64_t>());
            test_eq(sol, (k == 2 ? answer2[n - 2] : answer3[n - 2]));
        }
    }
}
