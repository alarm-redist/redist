/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Link-Cut Tree Implementation
 *
 * Based on Sleator & Tarjan's Link-Cut Trees.
 * Translated from CycleWalk.jl's splaytrees.jl and linkcuttrees.jl
 ********************************************************/

#include "cw_lct.h"
#include <stdexcept>
#include <algorithm>

// ============================================================
// LCTNode Implementation
// ============================================================

LCTNode::LCTNode(int v)
    : vertex(v), parent(nullptr), path_parent(nullptr), reversed(false) {
    children[0] = nullptr;
    children[1] = nullptr;
}

int LCTNode::child_index() const {
    if (parent == nullptr) return -1;
    return (parent->children[1] == this) ? 1 : 0;
}

bool LCTNode::is_splay_root() const {
    return parent == nullptr;
}

// ============================================================
// LinkCutTree Implementation
// ============================================================

LinkCutTree::LinkCutTree(int n) {
    nodes.reserve(n);
    for (int i = 0; i < n; i++) {
        nodes.emplace_back(i);
    }
}

LCTNode* LinkCutTree::node(int v) {
    return &nodes[v];
}

int LinkCutTree::size() const {
    return static_cast<int>(nodes.size());
}

// ============================================================
// Splay Tree Operations
// ============================================================

void LinkCutTree::set_child(LCTNode* p, int idx, LCTNode* c) {
    if (p != nullptr) {
        p->children[idx] = c;
    }
    if (c != nullptr) {
        c->parent = p;
    }
}

void LinkCutTree::push_reversed(LCTNode* n) {
    if (n == nullptr) return;

    // Build stack of ancestors (iterative instead of recursive to prevent stack overflow)
    std::vector<LCTNode*> stack;
    LCTNode* cur = n;
    int guard = 0;
    while (cur != nullptr) {
        if (++guard > 10000) {
            throw std::runtime_error("LCT push_reversed: parent cycle detected");
        }
        stack.push_back(cur);
        cur = cur->parent;
    }

    // Process from top (root) down to n
    for (int i = (int)stack.size() - 1; i >= 0; i--) {
        LCTNode* node = stack[i];
        if (node->reversed) {
            std::swap(node->children[0], node->children[1]);
            if (node->children[0] != nullptr) {
                node->children[0]->reversed = !node->children[0]->reversed;
            }
            if (node->children[1] != nullptr) {
                node->children[1]->reversed = !node->children[1]->reversed;
            }
            node->reversed = false;
        }
    }
}

void LinkCutTree::rotate_up(LCTNode* n) {
    int i = n->child_index();
    LCTNode* p = n->parent;
    LCTNode* g = p->parent;

    // Update n's parent
    n->parent = g;
    if (g != nullptr) {
        int j = p->child_index();
        g->children[j] = n;
    } else {
        // n becomes new splay root, inherit path parent from p
        n->path_parent = p->path_parent;
        p->path_parent = nullptr;
        if (n->path_parent != nullptr) {
            n->path_parent->path_children.erase(p);
            n->path_parent->path_children.insert(n);
        }
    }

    // Move n's appropriate child to p
    set_child(p, i, n->children[1 - i]);

    // Make p a child of n
    set_child(n, 1 - i, p);
}

void LinkCutTree::splay(LCTNode* n) {
    if (n == nullptr) return;

    // First push all reversal flags down to n
    push_reversed(n);

    int guard = 0;
    while (n->parent != nullptr) {
        if (++guard > 10000) {
            throw std::runtime_error("LCT splay: infinite loop detected");
        }
        LCTNode* p = n->parent;

        if (p->parent == nullptr) {
            // Zig: p is the splay root
            rotate_up(n);
        } else if (n->child_index() == p->child_index()) {
            // Zig-zig: n and p are same-side children
            rotate_up(p);
            rotate_up(n);
        } else {
            // Zig-zag: n and p are opposite-side children
            rotate_up(n);
            rotate_up(n);
        }
    }
}

void LinkCutTree::replace_right_subtree(LCTNode* n, LCTNode* r) {
    LCTNode* old_right = n->children[1];

    // Move old right child to be a path child
    if (old_right != nullptr) {
        old_right->path_parent = n;
        old_right->parent = nullptr;
        n->path_children.insert(old_right);
    }

    // Set new right child
    n->children[1] = r;
    if (r != nullptr) {
        r->parent = n;
        r->path_parent = nullptr;
        n->path_children.erase(r);
    }
}

// ============================================================
// Core LCT Operations
// ============================================================

void LinkCutTree::expose(LCTNode* n) {
    if (n == nullptr) return;

    splay(n);
    replace_right_subtree(n);

    int guard = 0;
    while (n->path_parent != nullptr) {
        if (++guard > 10000) {
            throw std::runtime_error("LCT expose: infinite loop detected");
        }
        LCTNode* p = n->path_parent;
        splay(p);
        replace_right_subtree(p, n);
        splay(n);
    }
}

void LinkCutTree::expose(int v) {
    expose(node(v));
}

void LinkCutTree::link(LCTNode* u, LCTNode* v) {
    expose(u);

    // u must be the root of its represented tree
    // After expose, u has no left child iff it's the root
    if (u->children[0] != nullptr) {
        throw std::invalid_argument("link: u must be the root of its represented tree");
    }

    expose(v);

    // u and v must be in different trees
    if (u->parent != nullptr || u->path_parent != nullptr) {
        throw std::invalid_argument("link: cannot link nodes in the same tree");
    }

    u->path_parent = v;
    v->path_children.insert(u);
}

void LinkCutTree::link(int u, int v) {
    link(node(u), node(v));
}

void LinkCutTree::cut(LCTNode* u) {
    expose(u);

    // u must not be the root (must have a left child after expose)
    if (u->children[0] == nullptr) {
        throw std::invalid_argument("cut: cannot cut the root of the represented tree");
    }

    LCTNode* left = u->children[0];
    left->parent = nullptr;
    u->children[0] = nullptr;
}

void LinkCutTree::cut(int u) {
    cut(node(u));
}

void LinkCutTree::evert(LCTNode* u) {
    expose(u);
    u->reversed = true;
}

void LinkCutTree::evert(int u) {
    evert(node(u));
}

LCTNode* LinkCutTree::find_root(LCTNode* u) {
    expose(u);

    // The root is the leftmost node in the splay tree
    LCTNode* r = u;
    push_reversed(r);
    int guard = 0;
    while (r->children[0] != nullptr) {
        if (++guard > 10000) {
            throw std::runtime_error("LCT find_root: infinite loop detected");
        }
        r = r->children[0];
        push_reversed(r);
    }

    // Splay the root to amortize cost
    splay(r);
    return r;
}

int LinkCutTree::find_root(int u) {
    return find_root(node(u))->vertex;
}

bool LinkCutTree::same_tree(int u, int v) {
    return find_root(u) == find_root(v);
}

void LinkCutTree::traverse_inorder(LCTNode* n, std::vector<int>& result, bool reversed) {
    if (n == nullptr) return;

    bool current_reversed = reversed ^ n->reversed;
    int left_idx = current_reversed ? 1 : 0;
    int right_idx = 1 - left_idx;

    traverse_inorder(n->children[left_idx], result, current_reversed);
    result.push_back(n->vertex);
    traverse_inorder(n->children[right_idx], result, current_reversed);
}

std::vector<int> LinkCutTree::find_path(LCTNode* u) {
    expose(u);

    // Traverse the splay tree in-order to get the path from root to u
    std::vector<int> path;
    LCTNode* root = u;
    int guard = 0;
    while (root->parent != nullptr) {
        if (++guard > 10000) {
            throw std::runtime_error("LCT find_path: infinite loop detected");
        }
        root = root->parent;
    }
    traverse_inorder(root, path, false);
    return path;
}

std::vector<int> LinkCutTree::find_path(int u) {
    return find_path(node(u));
}

std::vector<int> LinkCutTree::find_path_to_root(int u) {
    std::vector<int> path = find_path(u);
    std::reverse(path.begin(), path.end());
    return path;
}
