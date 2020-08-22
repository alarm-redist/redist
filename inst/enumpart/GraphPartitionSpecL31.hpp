/*
 * TdZdd: a Top-down/Breadth-first Decision Diagram Manipulation Framework
 * by Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2014 ERATO MINATO Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#pragma once

#include <cassert>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <tdzdd/DdSpec.hpp>

#include "Graph.hpp"

struct GraphPartitionSpecL31Count {
    int16_t comp; ///< uncolored edge component counter.

    GraphPartitionSpecL31Count()
            : comp(0) {
    }

    GraphPartitionSpecL31Count(int16_t uncoloredEdgeComponents)
            : comp(uncoloredEdgeComponents) {
    }

    size_t hash() const {
        return comp;
    }

    bool operator==(GraphPartitionSpecL31Count const& o) const {
        return comp == o.comp;
    }

    friend std::ostream& operator<<(std::ostream& os,
            GraphPartitionSpecL31Count const& o) {
        return os << o.comp;
    }
};

class GraphPartitionSpecL31Mate {
public:
    typedef int32_t Offset;
    typedef uint32_t uOffset;

private:
    Offset hoc; ///< offset to head or FPS.
    Offset nxt; ///< offset to next connected vertex.

public:
    GraphPartitionSpecL31Mate(Offset hoc = 0)
            : hoc(hoc), nxt(0) {
    }

    bool operator==(GraphPartitionSpecL31Mate const& o) const {
        return this == &o;
    }

    bool operator!=(GraphPartitionSpecL31Mate const& o) const {
        return this != &o;
    }

    void clear() {
        hoc = 0;
        nxt = 0;
    }

    bool isHead() const {
        return hoc >= 0;
    }

    bool isTail() const {
        return nxt == 0;
    }

    bool isIsolated() const {
        return isHead() && isTail();
    }

    GraphPartitionSpecL31Mate& head() {
        return isHead() ? *this : *(this + hoc);
    }

    GraphPartitionSpecL31Mate const& head() const {
        return isHead() ? *this : *(this + hoc);
    }

    GraphPartitionSpecL31Mate& next() {
        return *(this + nxt);
    }

    GraphPartitionSpecL31Mate const& next() const {
        return *(this + nxt);
    }

    void addToFPS(int offset) {
        hoc |= (1 << offset);
    }

    void eraseFromFPS(int offset) {
        uOffset v = static_cast<uOffset>(hoc) & ~(1u << offset);
        hoc = static_cast<Offset>(v);
    }

    bool isInFPS(int offset) const {
        return ((hoc >> offset) & 1) != 0;
    }

    void shiftHoc(int d) {
        hoc >>= d;
    }

    void print() const {
        std::cerr << "(" << hoc << "," << nxt << ")";
    }

    void print(std::ostream& ost) const {
        ost << "(" << hoc << "," << nxt << ")";
    }

    void print(FILE* fp) const {
        fprintf(fp, "(%d,%d)", hoc, nxt);
    }

    void mergeLists(GraphPartitionSpecL31Mate& o1, GraphPartitionSpecL31Mate& o2,
                    GraphPartitionSpecL31Mate* mate) {
        GraphPartitionSpecL31Mate* p1 = &o1.head();
        GraphPartitionSpecL31Mate* p2 = &o2.head();
        if (p1 == p2) return;
        if (p1 > p2) std::swap(p1, p2);

        // update fps start

        assert(p2->hoc >= 0);
        p1->hoc |= p2->hoc;

        for (GraphPartitionSpecL31Mate* q = this; q < p1; ++q) {
            if (q->hoc > 0 && q->isInFPS(p2 - mate)) {
                q->eraseFromFPS(p2 - mate);
                q->addToFPS(p1 - mate);
            }
        }

        for (GraphPartitionSpecL31Mate* q = p1 + 1; q < p2; ++q) {
            if (q->hoc > 0 && q->isInFPS(p2 - mate)) {
                q->eraseFromFPS(p2 - mate);
                p1->addToFPS(q - this);
            }
        }

        // update fps end

        for (GraphPartitionSpecL31Mate* q = p2;; q += q->nxt) {
            q->hoc = p1 - q;
            if (q->nxt == 0) break;
        }

        GraphPartitionSpecL31Mate* p = p1;
        GraphPartitionSpecL31Mate* q = p2;

        while (true) {
            assert(p != q);
            GraphPartitionSpecL31Mate* pp = p + p->nxt;
            assert(p <= pp && pp != q);

            while (p < pp && pp < q) {
                p = pp;
                pp += pp->nxt;
                assert(p <= pp && pp != q);
            }

            assert(p == pp || q < pp);
            p->nxt = q - p;
            if (p == pp) break;
            p = q, q = pp;
        }
    }

    void replaceHeadWith(GraphPartitionSpecL31Mate& newHead, GraphPartitionSpecL31Mate* mate) const {
        GraphPartitionSpecL31Mate const* p = &head();
        GraphPartitionSpecL31Mate* q = &newHead;

        Offset v = p->hoc;
        assert(v >= 0);

        for (int i = 0; i < q - mate; ++i) {
            if (((v >> i) & 1) != 0) {
                v = static_cast<Offset>(static_cast<uOffset>(v) & ~(1u << i));
                (mate + i)->addToFPS(q - mate);
            }
        }
        q->hoc = v;

        while (q->nxt > 0) {
            q += q->nxt;
            q->hoc = &newHead - q;
        }
    }

    void removeFromList(GraphPartitionSpecL31Mate const& o) {
        if (o.nxt == 0) {
            for (GraphPartitionSpecL31Mate* p = this; p <= &o; ++p) {
                if (p + p->nxt == &o) p->nxt = 0;
            }
        }
        else {
            for (GraphPartitionSpecL31Mate* p = this; p <= &o; ++p) {
                if (p + p->nxt == &o) p->nxt += o.nxt;
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os,
            GraphPartitionSpecL31Mate const& o) {
        return os << "<" << o.hoc << "," << o.nxt << ">";
    }
};

class GraphPartitionSpecL31: public tdzdd::HybridDdSpec<GraphPartitionSpecL31,
        GraphPartitionSpecL31Count,GraphPartitionSpecL31Mate,2> {
    typedef GraphPartitionSpecL31Count Count;
    typedef GraphPartitionSpecL31Mate Mate;

    Graph const& graph;
    int const m;
    int const n;
    int const mateSize;
    std::vector<Mate> initialMate;
    int numCOMP;
    bool const noLoop;
    bool const lookahead;

    int takable(Count& c, Mate const* mate, Graph::EdgeInfo const& e) const {
        Mate const& w1 = mate[e.v1 - e.v0];
        Mate const& w2 = mate[e.v2 - e.v0];

        // don't connect again
        if (noLoop && w1.head() == w2.head()) return false;
        //if (w1.head() == w2.head()) return false;

        if (w1.head() != w2.head()) { // v1 and v2 are in distinct components
            if (&w1.head() < &w2.head()) {
                if (w1.head().isInFPS(&w2.head() - mate)) {
                    return false;
                }
            } else {
                if (w2.head().isInFPS(&w1.head() - mate)) {
                    return false;
                }
            }
        }

        if (e.v1final && e.v2final) {
            if (w1.isIsolated() && w2.isIsolated()) { // new component leaves immediately
                if (c.comp == 0) return false;
                if (c.comp > 0) --c.comp;
            }
            else if (w1.isHead() && w2 == w1.next() && w2.isTail()) { // existing component leaves
                if (c.comp == 0) return false;
                if (c.comp > 0) --c.comp;
            }
        }

        if (e.finalEdge && c.comp > 0) return false;
        return true;
    }

    bool doTake(Count& count, Mate* mate, Graph::EdgeInfo const& e) const {
        Count c = count;

        if (!takable(c, mate, e)) return false;

        count = c;
        mate[0].mergeLists(mate[e.v1 - e.v0], mate[e.v2 - e.v0], mate);
        return true;
    }

    bool doNotTake(Count& count, Mate* mate, Graph::EdgeInfo const& e) const {
        Count c = count;
        Mate& w1 = mate[e.v1 - e.v0];
        Mate& w2 = mate[e.v2 - e.v0];

        if (w1.head() == w2.head()) {
            return false;
        }
        if (&w1.head() < &w2.head()) {
            w1.head().addToFPS(&w2.head() - mate);
        } else {
            w2.head().addToFPS(&w1.head() - mate);
        }

        if (e.v1final && w1.isIsolated()) {
            if (c.comp >= 0) {
                if (c.comp == 0) return false;
                --c.comp;
            }
        }

        if (e.v2final && w2.isIsolated()) {
            if (c.comp >= 0) {
                if (c.comp == 0) return false;
                --c.comp;
            }
        }

        if (e.v1final && e.v2final && w1.isHead() && w2 == w1.next()
                && w2.isTail()) { // existing component leaves) {
            if (c.comp == 0) return false;
            if (c.comp > 0) --c.comp;
        }

        if (e.finalEdge && c.comp > 0) return false;
        count = c;
        return true;
    }

    void update(Mate* mate, Graph::EdgeInfo const& e,
            Graph::EdgeInfo const& ee) const {
        int const d = ee.v0 - e.v0;
        assert(d >= 0);
        Mate* p1 = &mate[e.v1 - e.v0];
        Mate* p2 = &mate[e.v2 - e.v0];
        Mate* pd = p1 + d;

        for (Mate* q = p1; q < pd; ++q) {
            Mate* qq = &q->next();
            if (qq >= pd) {
                q->replaceHeadWith(*qq, mate);
            }
        }

        // update fps
        if (d > 0) {
            for (Mate* q = mate + d; q < mate + mateSize; ++q) {
                if (q->isHead() > 0) {
                    q->shiftHoc(d);
                }
            }
        }

        if (e.v2final) {
            mate[0].removeFromList(*p2);
            p2->clear();
        }

        if (e.v1final) {
            mate[0].removeFromList(*p1);
            p1->clear();
        }

        if (d > 0) {
            std::memmove(p1, pd, (mateSize - d) * sizeof(*mate));
            for (int i = mateSize - d; i < mateSize; ++i) {
                p1[i] = initialMate[ee.v0 + i];
            }
        }
    }

public:
    GraphPartitionSpecL31(Graph const& graph, int numCOMP = -1,
            bool noLoop = false, bool lookahead = true)
            : graph(graph), m(graph.vertexSize()), n(graph.edgeSize()),
              mateSize(graph.maxFrontierSize()), initialMate(1 + m + mateSize),
              numCOMP(numCOMP), noLoop(noLoop), lookahead(lookahead) {
        this->setArraySize(mateSize);

        int s = sizeof(GraphPartitionSpecL31Mate::Offset) * 8 - 1;
        if (mateSize > s) {
            std::cerr << "mateSize = " << mateSize << ". It must not exceed "
                      << s << "." << std::endl;
            exit(1);
        }

        for (int v = 1; v <= m; ++v) {
            initialMate[v] = Mate(0);
        }
    }

    int getRoot(Count& count, Mate* mate) const {
        int const v0 = graph.edgeInfo(0).v0;

        count = Count(numCOMP);

        for (int i = 0; i < mateSize; ++i) {
            mate[i] = initialMate[v0 + i];
        }

        return n;
    }

    int getChild(Count& count, Mate* mate, int level, int take) const {
        assert(1 <= level && level <= n);
        int i = n - level;
        Graph::EdgeInfo const* e = &graph.edgeInfo(i);

        if (take) {
            if (!doTake(count, mate, *e)) return 0;
        }
        else {
            if (!doNotTake(count, mate, *e)) return 0;
        }

        if (++i == n) return -1;

        Graph::EdgeInfo const* ee = &graph.edgeInfo(i);
        update(mate, *e, *ee);

        while (lookahead) {
            e = ee;

            Count c = count;
            if (takable(c, mate, *e)) break;
            if (!doNotTake(count, mate, *e)) return 0;

            if (++i == n) return -1;

            ee = &graph.edgeInfo(i);
            update(mate, *e, *ee);
        }

        assert(i < n);
        return n - i;
    }

    size_t hashCode(Count const& count) const {
        return count.hash();
    }
};
