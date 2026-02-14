/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Link-Cut Tree Data Structure
 *
 * Implementation of Sleator & Tarjan's Link-Cut Trees for
 * dynamic tree operations in O(log n) amortized time.
 * Reference: https://dl.acm.org/doi/10.1145/3828.3835
 ********************************************************/

#ifndef CW_LCT_H
#define CW_LCT_H

#include <vector>
#include <set>

/*
 * LCTNode: A node in the Link-Cut Tree.
 *
 * Each node represents a vertex in the represented forest.
 * The LCT uses splay trees as auxiliary trees to represent
 * preferred paths.
 *
 * Key concepts:
 * - parent: Parent in the splay tree (auxiliary tree)
 * - path_parent: Parent in the represented tree when this node
 *   is the root of its splay tree but not the root of the
 *   represented tree
 * - children[0/1]: Left/right children in splay tree
 * - reversed: Lazy reversal flag for evert operation
 */
struct LCTNode {
    int vertex;                     // Vertex ID (0-indexed)
    LCTNode* parent;                // Splay tree parent
    LCTNode* path_parent;           // Path parent (for non-preferred edges)
    LCTNode* children[2];           // Left (0) and right (1) children
    bool reversed;                  // Lazy reversal flag
    std::set<LCTNode*> path_children; // Set of path children (for traversal)

    LCTNode(int v);

    // Get child index (0 for left, 1 for right) in parent's children
    int child_index() const;

    // Check if this node is the root of its splay tree
    bool is_splay_root() const;
};

/*
 * LinkCutTree: Container for LCT nodes.
 *
 * Provides an interface to the LCT operations using vertex IDs.
 * Nodes are stored in a vector and accessed by vertex ID.
 */
class LinkCutTree {
public:
    std::vector<LCTNode> nodes;

    // Constructor: Create n isolated nodes
    explicit LinkCutTree(int n);

    // Access node by vertex ID (0-indexed)
    LCTNode* node(int v);

    // Number of nodes
    int size() const;

    // === Core LCT Operations ===

    /*
     * expose(v): Make the path from v to the root a preferred path.
     * After this operation, v is the root of its splay tree and
     * the rightmost node on the preferred path.
     */
    void expose(int v);
    void expose(LCTNode* n);

    /*
     * link(u, v): Make u a child of v in the represented tree.
     * Precondition: u must be the root of its represented tree.
     */
    void link(int u, int v);
    void link(LCTNode* u, LCTNode* v);

    /*
     * cut(u): Cut u from its parent in the represented tree.
     * Precondition: u must not be the root of its represented tree.
     */
    void cut(int u);
    void cut(LCTNode* u);

    /*
     * evert(u): Make u the root of its represented tree.
     * This reverses the path from u to the current root.
     */
    void evert(int u);
    void evert(LCTNode* u);

    /*
     * find_root(u): Find the root of the represented tree containing u.
     * Returns the vertex ID of the root.
     */
    int find_root(int u);
    LCTNode* find_root(LCTNode* u);

    /*
     * same_tree(u, v): Check if u and v are in the same represented tree.
     */
    bool same_tree(int u, int v);

    /*
     * find_path(u): Get the path from the root to u in the represented tree.
     * Returns a vector of vertex IDs in order from root to u.
     */
    std::vector<int> find_path(int u);
    std::vector<int> find_path(LCTNode* u);

    /*
     * find_path_to_root(u): Get the path from u to the root.
     * Returns a vector of vertex IDs in order from u to root.
     */
    std::vector<int> find_path_to_root(int u);

    // === Splay Tree Operations (internal) ===

private:
    // Push lazy reversal down to children
    void push_reversed(LCTNode* n);

    // Rotate n up in the splay tree
    void rotate_up(LCTNode* n);

    // Splay n to the root of its auxiliary tree
    void splay(LCTNode* n);

    // Replace right subtree of n, moving old right to path child
    void replace_right_subtree(LCTNode* n, LCTNode* r = nullptr);

    // Set child at index (0=left, 1=right)
    void set_child(LCTNode* parent, int idx, LCTNode* child);

    // Helper to traverse splay tree in order
    void traverse_inorder(LCTNode* n, std::vector<int>& result, bool reversed);
};

#endif // CW_LCT_H
