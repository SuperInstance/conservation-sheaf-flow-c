#include "csf.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── Graph primitives ─────────────────────────────────────────────── */

static void graph_add_adj(CSFGraph *g, int v, int eidx) {
    if (g->adj_count[v] >= g->adj_cap[v]) {
        g->adj_cap[v] = g->adj_cap[v] ? g->adj_cap[v] * 2 : 4;
        g->adj[v] = realloc(g->adj[v], g->adj_cap[v] * sizeof(int));
    }
    g->adj[v][g->adj_count[v]++] = eidx;
}

CSFGraph csf_graph_create(int n, int m, CSFEdge *edges) {
    CSFGraph g;
    g.n = n;
    g.m = m;
    g.edges = edges;
    g.adj = calloc(n, sizeof(int *));
    g.adj_count = calloc(n, sizeof(int));
    g.adj_cap = calloc(n, sizeof(int));
    for (int e = 0; e < m; e++) {
        graph_add_adj(&g, edges[e].u, e);
        graph_add_adj(&g, edges[e].v, e);
    }
    return g;
}

void csf_graph_free(CSFGraph *g) {
    free(g->edges);
    for (int i = 0; i < g->n; i++) free(g->adj[i]);
    free(g->adj);
    free(g->adj_count);
    free(g->adj_cap);
    memset(g, 0, sizeof(*g));
}

CSFGraph csf_graph_copy(const CSFGraph *g) {
    CSFGraph c;
    c.n = g->n;
    c.m = g->m;
    c.edges = malloc(g->m * sizeof(CSFEdge));
    memcpy(c.edges, g->edges, g->m * sizeof(CSFEdge));
    c.adj = calloc(g->n, sizeof(int *));
    c.adj_count = calloc(g->n, sizeof(int));
    c.adj_cap = calloc(g->n, sizeof(int));
    for (int v = 0; v < g->n; v++) {
        c.adj_count[v] = g->adj_count[v];
        c.adj_cap[v] = g->adj_count[v];
        c.adj[v] = malloc(g->adj_count[v] * sizeof(int));
        memcpy(c.adj[v], g->adj[v], g->adj_count[v] * sizeof(int));
    }
    return c;
}

CSFGraph csf_graph_path(int n) {
    CSFEdge *e = malloc((n > 1 ? n - 1 : 0) * sizeof(CSFEdge));
    for (int i = 0; i < n - 1; i++) {
        e[i].u = i; e[i].v = i + 1; e[i].weight = 1.0;
    }
    return csf_graph_create(n, n > 1 ? n - 1 : 0, e);
}

CSFGraph csf_graph_cycle(int n) {
    if (n < 2) return csf_graph_path(n);
    CSFEdge *e = malloc(n * sizeof(CSFEdge));
    for (int i = 0; i < n - 1; i++) {
        e[i].u = i; e[i].v = i + 1; e[i].weight = 1.0;
    }
    e[n - 1].u = n - 1; e[n - 1].v = 0; e[n - 1].weight = 1.0;
    return csf_graph_create(n, n, e);
}

CSFGraph csf_graph_tree(int n, const CSFEdge *tree_edges) {
    CSFEdge *e = malloc((n - 1) * sizeof(CSFEdge));
    memcpy(e, tree_edges, (n - 1) * sizeof(CSFEdge));
    return csf_graph_create(n, n - 1, e);
}

CSFGraph csf_graph_complete(int n) {
    int m = n * (n - 1) / 2;
    CSFEdge *e = malloc(m * sizeof(CSFEdge));
    int idx = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            e[idx].u = i; e[idx].v = j; e[idx].weight = 1.0;
            idx++;
        }
    return csf_graph_create(n, m, e);
}

/* ── Dense matrix utilities ───────────────────────────────────────── */

CSFMatrix csf_mat_alloc(int n) {
    CSFMatrix m;
    m.n = n;
    m.data = calloc((size_t)n * n, sizeof(double));
    return m;
}

void csf_mat_free(CSFMatrix *m) {
    free(m->data);
    m->data = NULL;
    m->n = 0;
}

void csf_mat_zero(CSFMatrix *m) {
    memset(m->data, 0, (size_t)m->n * m->n * sizeof(double));
}

void csf_mat_identity(CSFMatrix *m) {
    csf_mat_zero(m);
    for (int i = 0; i < m->n; i++) m->data[i * m->n + i] = 1.0;
}

double *csf_mat_at(CSFMatrix *m, int i, int j) {
    return &m->data[i * m->n + j];
}

void csf_mat_mul(const CSFMatrix *A, const CSFMatrix *B, CSFMatrix *C) {
    int n = A->n;
    csf_mat_zero(C);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) {
            double a = A->data[i * n + k];
            for (int j = 0; j < n; j++)
                C->data[i * n + j] += a * B->data[k * n + j];
        }
}

void csf_mat_add(double alpha, const CSFMatrix *A,
                 double beta,  const CSFMatrix *B, CSFMatrix *C) {
    int n = A->n;
    for (int i = 0; i < n * n; i++)
        C->data[i] = alpha * A->data[i] + beta * B->data[i];
}

void csf_mat_vec(const CSFMatrix *M, const double *x, double *y) {
    int n = M->n;
    for (int i = 0; i < n; i++) {
        y[i] = 0;
        for (int j = 0; j < n; j++)
            y[i] += M->data[i * n + j] * x[j];
    }
}

/* ── Conservation sheaf ───────────────────────────────────────────── */

ConservationSheaf csf_sheaf_create(int stalk_dim, CSFGraph graph) {
    ConservationSheaf s;
    s.stalk_dim = stalk_dim;
    /* Deep-copy the graph so the sheaf owns its own memory.
     * This prevents double-free when both the sheaf and the caller's
     * graph are freed independently. */
    s.graph = csf_graph_copy(&graph);
    s.is_static = true;
    s.restrictions = malloc(graph.m * sizeof(CSFMatrix));
    s.restrictions_T = malloc(graph.m * sizeof(CSFMatrix));

    /* Default: identity restriction maps (standard diffusion sheaf) */
    for (int e = 0; e < graph.m; e++) {
        s.restrictions[e] = csf_mat_alloc(stalk_dim);
        csf_mat_identity(&s.restrictions[e]);
        s.restrictions_T[e] = csf_mat_alloc(stalk_dim);
        csf_mat_identity(&s.restrictions_T[e]);
    }

    int N = graph.n * stalk_dim;
    s.laplacian = csf_mat_alloc(N);
    s.laplacian_valid = false;
    csf_sheaf_build_laplacian(&s);
    return s;
}

ConservationSheaf csf_sheaf_create_with_restrictions(
    int stalk_dim, CSFGraph graph, const CSFMatrix *restrictions) {
    ConservationSheaf s;
    s.stalk_dim = stalk_dim;
    /* Deep-copy the graph so the sheaf owns its own memory. */
    s.graph = csf_graph_copy(&graph);
    s.is_static = true;
    s.restrictions = malloc(graph.m * sizeof(CSFMatrix));
    s.restrictions_T = malloc(graph.m * sizeof(CSFMatrix));

    for (int e = 0; e < graph.m; e++) {
        s.restrictions[e] = csf_mat_alloc(stalk_dim);
        memcpy(s.restrictions[e].data, restrictions[e].data,
               (size_t)stalk_dim * stalk_dim * sizeof(double));
        s.restrictions_T[e] = csf_mat_alloc(stalk_dim);
        /* Transpose */
        for (int i = 0; i < stalk_dim; i++)
            for (int j = 0; j < stalk_dim; j++)
                *csf_mat_at(&s.restrictions_T[e], i, j) =
                    restrictions[e].data[j * stalk_dim + i];
    }

    int N = graph.n * stalk_dim;
    s.laplacian = csf_mat_alloc(N);
    s.laplacian_valid = false;
    csf_sheaf_build_laplacian(&s);
    return s;
}

void csf_sheaf_free(ConservationSheaf *s) {
    for (int e = 0; e < s->graph.m; e++) {
        csf_mat_free(&s->restrictions[e]);
        csf_mat_free(&s->restrictions_T[e]);
    }
    free(s->restrictions);
    free(s->restrictions_T);
    csf_mat_free(&s->laplacian);
    csf_graph_free(&s->graph);
}

void csf_sheaf_build_laplacian(ConservationSheaf *s) {
    int d = s->stalk_dim;
    int N = s->graph.n * d;
    csf_mat_zero(&s->laplacian);

    /*
     * Sheaf Laplacian: L = B^T W B
     * where B is the coboundary (restriction map on oriented edges)
     * and W is the edge weight diagonal.
     *
     * L_{vv} = sum_{e ~ v} w_e * R_e^T R_e   (diagonal block)
     * L_{vu} = -w_e * R_e^T R_e               (off-diagonal, for edge v-u)
     *
     * For conservation sheaf with identity restrictions:
     * This reduces to the standard graph Laplacian, which is correct.
     */

    for (int e = 0; e < s->graph.m; e++) {
        int u = s->graph.edges[e].u;
        int v = s->graph.edges[e].v;
        double w = s->graph.edges[e].weight;

        /* Compute R^T R for this edge */
        CSFMatrix RtR = csf_mat_alloc(d);
        CSFMatrix Rt = s->restrictions_T[e];
        CSFMatrix R = s->restrictions[e];
        /* RtR = R^T * R */
        csf_mat_mul(&Rt, &R, &RtR);

        /* Add w * RtR to diagonal blocks (u,u) and (v,v) */
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++) {
                double val = w * RtR.data[i * d + j];
                s->laplacian.data[(u * d + i) * N + (u * d + j)] += val;
                s->laplacian.data[(v * d + i) * N + (v * d + j)] += val;
                /* Subtract from off-diagonal blocks */
                s->laplacian.data[(u * d + i) * N + (v * d + j)] -= val;
                s->laplacian.data[(v * d + i) * N + (u * d + j)] -= val;
            }
        csf_mat_free(&RtR);
    }
    s->laplacian_valid = true;
}

bool csf_verify_conservation_at_vertex(const ConservationSheaf *s,
                                       const double *quantities, int vertex) {
    int d = s->stalk_dim;
    int n = s->graph.n;

    /* For conservation: sum of Laplacian applied to quantities at this vertex
       should be ~0. More precisely, flux in = flux out. */
    double flux = 0;
    for (int i = 0; i < d; i++) {
        double lap_val = 0;
        for (int j = 0; j < n * d; j++) {
            lap_val += s->laplacian.data[(vertex * d + i) * (n * d) + j] * quantities[j];
        }
        flux += lap_val;
    }
    return fabs(flux) < 1e-9;
}

bool csf_global_conservation_check(const ConservationSheaf *s,
                                   const double *quantities) {
    int d = s->stalk_dim;
    int n = s->graph.n;

    /* Total conserved quantity per stalk component should be constant.
       Check that the all-ones vector is in the kernel of the Laplacian,
       which means the sum of quantities is preserved under flow. */
    for (int c = 0; c < d; c++) {
        double total = 0;
        for (int v = 0; v < n; v++)
            total += quantities[v * d + c];
        if (fabs(total) < 1e-12) continue; /* skip zero components */

        /* Check L * 1_c = 0 (where 1_c is 1 on component c at each vertex) */
        for (int v = 0; v < n; v++) {
            double lap_val = 0;
            for (int j = 0; j < n * d; j++)
                lap_val += s->laplacian.data[(v * d + c) * (n * d) + j] *
                           (j % d == c ? 1.0 : 0.0);
            if (fabs(lap_val) > 1e-9) return false;
        }
    }
    return true;
}

/* ── Eigenvalue computation ───────────────────────────────────────── */

/* Solve M * x = b via Gaussian elimination (for small dense matrices).
   Overwrites b with solution x. M is destroyed. */
static void csf_solve_system(int n, double *M, double *b) {
    /* Forward elimination with partial pivoting */
    for (int col = 0; col < n; col++) {
        /* Find pivot */
        int pivot = col;
        double max_val = fabs(M[col * n + col]);
        for (int row = col + 1; row < n; row++) {
            double v = fabs(M[row * n + col]);
            if (v > max_val) { max_val = v; pivot = row; }
        }
        /* Swap rows */
        if (pivot != col) {
            for (int j = 0; j < n; j++) {
                double tmp = M[col * n + j];
                M[col * n + j] = M[pivot * n + j];
                M[pivot * n + j] = tmp;
            }
            double tmp = b[col]; b[col] = b[pivot]; b[pivot] = tmp;
        }
        if (fabs(M[col * n + col]) < 1e-15) {
            /* Singular — zero out to skip */
            M[col * n + col] = 1e-15;
        }
        /* Eliminate below */
        for (int row = col + 1; row < n; row++) {
            double factor = M[row * n + col] / M[col * n + col];
            for (int j = col + 1; j < n; j++)
                M[row * n + j] -= factor * M[col * n + j];
            M[row * n + col] = 0;
            b[row] -= factor * b[col];
        }
    }
    /* Back substitution */
    for (int row = n - 1; row >= 0; row--) {
        for (int j = row + 1; j < n; j++)
            b[row] -= M[row * n + j] * b[j];
        b[row] /= M[row * n + row];
    }
}

double csf_eigenvalue_min(const CSFMatrix *M, int iterations, double tol) {
    /* For PSD matrices: find smallest eigenvalue via inverse iteration.
       Inverse iteration on A converges to the eigenvector for the eigenvalue
       closest to 0, which for a PSD matrix is λ_min = 0.
       To find λ_1 (smallest non-zero), we project out the null space first. */
    int n = M->n;
    double *v = malloc(n * sizeof(double));
    double *w = malloc(n * sizeof(double));
    double *Acopy = malloc((size_t)n * n * sizeof(double));

    /* Initialize with non-constant vector */
    for (int i = 0; i < n; i++) v[i] = (i % 2 == 0) ? 1.0 : -1.0;

    /* Project out constant (kernel of graph Laplacian) */
    double mean = 0;
    for (int i = 0; i < n; i++) mean += v[i];
    mean /= n;
    for (int i = 0; i < n; i++) v[i] -= mean;

    /* Normalize */
    double norm = 0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm > 1e-15)
        for (int i = 0; i < n; i++) v[i] /= norm;

    /* Inverse iteration: solve (A + shift*I) * w = v, then v = w/||w||
       The shift ensures the matrix is non-singular and finds eigenvalue near 0.
       Since we projected out the kernel, this converges to λ_1. */
    double shift = 1e-6;
    for (int it = 0; it < iterations; it++) {
        /* Copy shifted matrix */
        for (int i = 0; i < n * n; i++) Acopy[i] = M->data[i];
        for (int i = 0; i < n; i++) Acopy[i * n + i] += shift;

        /* Solve (A + shift*I) * w = v */
        memcpy(w, v, n * sizeof(double));
        csf_solve_system(n, Acopy, w);

        /* Project out constant again */
        mean = 0;
        for (int i = 0; i < n; i++) mean += w[i];
        mean /= n;
        for (int i = 0; i < n; i++) w[i] -= mean;

        norm = 0;
        for (int i = 0; i < n; i++) norm += w[i] * w[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = w[i] / norm;
    }

    /* Rayleigh quotient: v^T A v = eigenvalue */
    csf_mat_vec(M, v, w);
    double rq = 0;
    for (int i = 0; i < n; i++) rq += v[i] * w[i];

    free(v); free(w); free(Acopy);
    return rq;
}

double *csf_eigenvalues_k_smallest(const CSFMatrix *M, int k,
                                   int iterations, double tol) {
    int n = M->n;
    double *eigs = malloc(k * sizeof(double));
    double *v = malloc(n * sizeof(double));
    double *w = malloc(n * sizeof(double));
    double *Acopy = malloc((size_t)n * n * sizeof(double));

    /* lambda_0 = 0 for connected graph Laplacian */
    eigs[0] = 0.0;

    if (k <= 1) { free(v); free(w); free(Acopy); return eigs; }

    /* Find λ_1 via inverse iteration with kernel projected out */
    for (int i = 0; i < n; i++) v[i] = (i % 2 == 0) ? 1.0 : -1.0;

    /* Project out constant */
    double mean = 0;
    for (int i = 0; i < n; i++) mean += v[i];
    mean /= n;
    for (int i = 0; i < n; i++) v[i] -= mean;
    double norm = 0;
    for (int i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm > 1e-15)
        for (int i = 0; i < n; i++) v[i] /= norm;

    double shift = 1e-6;
    for (int it = 0; it < iterations; it++) {
        for (int i = 0; i < n * n; i++) Acopy[i] = M->data[i];
        for (int i = 0; i < n; i++) Acopy[i * n + i] += shift;
        memcpy(w, v, n * sizeof(double));
        csf_solve_system(n, Acopy, w);
        mean = 0;
        for (int i = 0; i < n; i++) mean += w[i];
        mean /= n;
        for (int i = 0; i < n; i++) w[i] -= mean;
        norm = 0;
        for (int i = 0; i < n; i++) norm += w[i] * w[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (int i = 0; i < n; i++) v[i] = w[i] / norm;
    }

    csf_mat_vec(M, v, w);
    double rq = 0;
    for (int i = 0; i < n; i++) rq += v[i] * w[i];
    eigs[1] = rq;

    /* For k > 2, use deflation with inverse iteration */
    /* Store found eigenvectors for deflation */
    double *evecs = NULL;
    if (k > 2) {
        evecs = malloc((size_t)(k - 1) * n * sizeof(double));
        memcpy(evecs, v, n * sizeof(double));  /* first non-trivial eigenvector */
    }

    for (int eig = 2; eig < k; eig++) {
        /* Initialize with diverse vector */
        for (int i = 0; i < n; i++) v[i] = sin((double)(i * (eig + 2)));

        /* Project out all previously found eigenvectors and the constant */
        for (int def = 0; def < eig - 1; def++) {
            double dot = 0;
            for (int i = 0; i < n; i++) dot += v[i] * evecs[def * n + i];
            for (int i = 0; i < n; i++) v[i] -= dot * evecs[def * n + i];
        }
        /* Project out constant */
        mean = 0;
        for (int i = 0; i < n; i++) mean += v[i];
        mean /= n;
        for (int i = 0; i < n; i++) v[i] -= mean;

        norm = 0;
        for (int i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm > 1e-15)
            for (int i = 0; i < n; i++) v[i] /= norm;

        /* Shifted inverse iteration to find next eigenvalue */
        double target_shift = eigs[eig - 1] + 1e-6;
        for (int it = 0; it < iterations; it++) {
            for (int i = 0; i < n * n; i++) Acopy[i] = M->data[i];
            for (int i = 0; i < n; i++) Acopy[i * n + i] -= target_shift;
            memcpy(w, v, n * sizeof(double));
            csf_solve_system(n, Acopy, w);

            /* Project out all previous eigenvectors + constant */
            for (int def = 0; def < eig - 1; def++) {
                double dot = 0;
                for (int i = 0; i < n; i++) dot += w[i] * evecs[def * n + i];
                for (int i = 0; i < n; i++) w[i] -= dot * evecs[def * n + i];
            }
            mean = 0;
            for (int i = 0; i < n; i++) mean += w[i];
            mean /= n;
            for (int i = 0; i < n; i++) w[i] -= mean;

            norm = 0;
            for (int i = 0; i < n; i++) norm += w[i] * w[i];
            norm = sqrt(norm);
            if (norm < 1e-15) break;
            for (int i = 0; i < n; i++) v[i] = w[i] / norm;
        }

        csf_mat_vec(M, v, w);
        rq = 0;
        for (int i = 0; i < n; i++) rq += v[i] * w[i];
        eigs[eig] = rq;
        memcpy(evecs + (eig - 1) * n, v, n * sizeof(double));
    }

    free(v); free(w); free(Acopy); free(evecs);
    return eigs;
}

/* ── Flow state ───────────────────────────────────────────────────── */

FlowState csf_flow_create(int stalk_dim, int n_vertices,
                          const double *initial, double dt) {
    FlowState fs;
    fs.stalk_dim = stalk_dim;
    fs.n_vertices = n_vertices;
    fs.quantities = malloc((size_t)n_vertices * stalk_dim * sizeof(double));
    memcpy(fs.quantities, initial, (size_t)n_vertices * stalk_dim * sizeof(double));
    fs.dt = dt;
    fs.time = 0;
    return fs;
}

void csf_flow_free(FlowState *fs) {
    free(fs->quantities);
    fs->quantities = NULL;
}

FlowState csf_flow_copy(const FlowState *fs) {
    FlowState c;
    c.stalk_dim = fs->stalk_dim;
    c.n_vertices = fs->n_vertices;
    c.dt = fs->dt;
    c.time = fs->time;
    size_t sz = (size_t)fs->n_vertices * fs->stalk_dim * sizeof(double);
    c.quantities = malloc(sz);
    memcpy(c.quantities, fs->quantities, sz);
    return c;
}

double csf_flow_entropy(const FlowState *fs) {
    int n = fs->n_vertices;
    int d = fs->stalk_dim;

    /* Compute entropy per component, then average */
    double total_entropy = 0;
    int valid_components = 0;

    for (int c = 0; c < d; c++) {
        double sum = 0;
        for (int v = 0; v < n; v++)
            sum += fabs(fs->quantities[v * d + c]);

        if (sum < 1e-15) continue;
        valid_components++;

        double entropy = 0;
        for (int v = 0; v < n; v++) {
            double p = fabs(fs->quantities[v * d + c]) / sum;
            if (p > 1e-15)
                entropy -= p * log(p);
        }
        /* Normalize: max entropy = log(n) */
        if (n > 1)
            entropy /= log(n);
        total_entropy += entropy;
    }

    return valid_components > 0 ? total_entropy / valid_components : 0;
}

/* ── Spectral gap ─────────────────────────────────────────────────── */

double csf_flow_spectral_gap(const ConservationSheaf *s) {
    double *eigs = csf_eigenvalues_k_smallest(&s->laplacian, 2, 200, 1e-10);
    double gap = eigs[1];
    free(eigs);
    return gap;
}
