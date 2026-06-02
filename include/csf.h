#ifndef CSF_H
#define CSF_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ── Graph primitives ─────────────────────────────────────────────── */

typedef struct {
    int u, v;            /* edge endpoints */
    double weight;       /* edge weight (conductance) */
} CSFEdge;

typedef struct {
    int n;               /* number of vertices */
    int m;               /* number of edges */
    CSFEdge *edges;      /* edge list */
    /* adjacency: adj[v] = { edge indices } */
    int **adj;
    int  *adj_count;
    int  *adj_cap;
} CSFGraph;

/* Build a graph from an edge list. Takes ownership of edges array. */
CSFGraph csf_graph_create(int n, int m, CSFEdge *edges);
CSFGraph csf_graph_copy(const CSFGraph *g);
void     csf_graph_free(CSFGraph *g);

/* Preset topologies */
CSFGraph csf_graph_path(int n);
CSFGraph csf_graph_cycle(int n);
CSFGraph csf_graph_tree(int n, const CSFEdge *tree_edges); /* n-1 edges */
CSFGraph csf_graph_complete(int n);

/* ── Dense matrix utilities ───────────────────────────────────────── */

typedef struct {
    int    n;
    double *data;       /* row-major, n*n */
} CSFMatrix;

CSFMatrix csf_mat_alloc(int n);
void      csf_mat_free(CSFMatrix *m);
void      csf_mat_zero(CSFMatrix *m);
void      csf_mat_identity(CSFMatrix *m);
double   *csf_mat_at(CSFMatrix *m, int i, int j);

/* Multiply: C = A*B */
void csf_mat_mul(const CSFMatrix *A, const CSFMatrix *B, CSFMatrix *C);
/* Add: C = alpha*A + beta*B */
void csf_mat_add(double alpha, const CSFMatrix *A,
                 double beta,  const CSFMatrix *B, CSFMatrix *C);

/* ── Conservation sheaf ───────────────────────────────────────────── */

/*
 * A conservation sheaf assigns a stalk (vector space of conserved
 * quantities) to each vertex and a restriction map (flux operator)
 * to each edge. The sheaf Laplacian encodes conservation: flux in
 * equals flux out at every vertex.
 */

typedef struct {
    int stalk_dim;            /* dimension of each stalk */
    CSFGraph graph;

    /* restriction maps: one per edge, size stalk_dim x stalk_dim */
    CSFMatrix *restrictions;  /* length = graph.m */

    /* coboundary adjoints (transpose of restriction), precomputed */
    CSFMatrix *restrictions_T;

    /* sheaf Laplacian: n*stalk_dim x n*stalk_dim */
    CSFMatrix laplacian;
    bool      laplacian_valid;

    /* Whether the sheaf is static (restriction maps fixed).
     * Only static sheaves satisfy the non-decreasing spectral gap theorem.
     * Evolving sheaves (restriction maps change with flow) are EXPERIMENTAL. */
    bool      is_static;
} ConservationSheaf;

/* Construct a conservation sheaf with identity restriction maps */
ConservationSheaf csf_sheaf_create(int stalk_dim, CSFGraph graph);
/* Construct with custom restriction maps (copies them) */
ConservationSheaf csf_sheaf_create_with_restrictions(
    int stalk_dim, CSFGraph graph, const CSFMatrix *restrictions);
void csf_sheaf_free(ConservationSheaf *s);

/* Build/rebuild the sheaf Laplacian */
void csf_sheaf_build_laplacian(ConservationSheaf *s);

/* Verify conservation: check flux balance at a vertex */
bool csf_verify_conservation_at_vertex(const ConservationSheaf *s,
                                       const double *quantities, /* n*stalk_dim */
                                       int vertex);
/* Check global conservation: total quantity preserved */
bool csf_global_conservation_check(const ConservationSheaf *s,
                                   const double *quantities);

/* ── Flow dynamics ────────────────────────────────────────────────── */

typedef struct {
    int stalk_dim;
    int n_vertices;
    double *quantities;   /* n_vertices * stalk_dim */
    double dt;            /* time step */
    double time;          /* accumulated time */
} FlowState;

/* Create a flow state (copies initial quantities) */
FlowState csf_flow_create(int stalk_dim, int n_vertices,
                          const double *initial, double dt);
FlowState csf_flow_copy(const FlowState *fs);
void csf_flow_free(FlowState *fs);

/* One step of sheaf heat equation: dq/dt = -L q */
void csf_flow_step(FlowState *fs, const ConservationSheaf *s);

/* Iterate until steady state or max_steps.
 * Returns true if converged. */
bool csf_flow_converge(FlowState *fs, const ConservationSheaf *s,
                       int max_steps, double tol);

/* Shannon entropy of flow state (per-vertex average, normalized to [0,1]) */
double csf_flow_entropy(const FlowState *fs);

/* Spectral gap of conservation sheaf Laplacian (smallest non-zero eigenvalue) */
double csf_flow_spectral_gap(const ConservationSheaf *s);

/* ── Optimal transport on conservation sheaves ────────────────────── */

/* Wasserstein-like cost between two flow states (Euclidean) */
double csf_transport_cost(const FlowState *a, const FlowState *b);

/*
 * Compute optimal transport plan between two flow configurations.
 * For scalar quantities on vertices, this reduces to a matching problem.
 * Returns the transport matrix (n x n, caller must free).
 */
double *csf_transport_plan(const FlowState *src, const FlowState *dst);

/* Barycenter: find flow state minimizing total transport cost to targets */
FlowState csf_barycenter_flow(const FlowState **targets, int k);

/* ── Theorem verification ─────────────────────────────────────────── */

typedef struct {
    int    n_steps;
    double *time_points;    /* length n_steps */
    double *spectral_gaps;  /* length n_steps */
    double *entropies;      /* length n_steps */
    bool    theorem_holds;  /* non-decreasing spectral gap */
    int     violation_step; /* first violation, -1 if none */
} SpectralEvolution;

/* Track spectral gap evolution during flow */
SpectralEvolution csf_track_spectral_gap(const ConservationSheaf *s,
                                         const FlowState *initial,
                                         int n_steps);

/* Track spectral gap for evolving sheaf (EXPERIMENTAL — not a proven theorem) */
SpectralEvolution csf_track_spectral_gap_evolving(
    const ConservationSheaf *s_template,
    const FlowState *initial,
    int n_steps);

/* Verify the theorem: spectral gap is non-decreasing during flow */
bool csf_verify_non_decreasing_gap(const SpectralEvolution *ev);

/* Compute critical time: when gap stops increasing significantly */
double csf_compute_critical_time(const SpectralEvolution *ev);

void csf_spectral_evolution_free(SpectralEvolution *ev);

/* ── Eigenvalue computation (power iteration) ─────────────────────── */

/* Compute smallest eigenvalue of M (for positive semi-definite matrices) */
double csf_eigenvalue_min(const CSFMatrix *M, int iterations, double tol);

/* Compute k smallest eigenvalues via deflation */
double *csf_eigenvalues_k_smallest(const CSFMatrix *M, int k,
                                   int iterations, double tol);

/* Matrix-vector multiply: y = M*x */
void csf_mat_vec(const CSFMatrix *M, const double *x, double *y);

#ifdef __cplusplus
}
#endif

#endif /* CSF_H */
