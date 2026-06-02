#include "csf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        printf("  FAIL: %s (line %d)\n", msg, __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define ASSERT_FEQ(a, b, tol, msg) ASSERT(fabs((a) - (b)) < (tol), msg)
#define TASSERT(cond, msg) do { if (!(cond)) { printf("  FAIL: %s\n", msg); tests_failed++; } else { tests_passed++; } } while(0)

/* ══════════════════════════════════════════════════════════════════════
 * Helper: create random initial quantities
 * ══════════════════════════════════════════════════════════════════════ */
static double *random_quantities(int n, int d) {
    double *q = malloc((size_t)n * d * sizeof(double));
    for (int i = 0; i < n * d; i++)
        q[i] = (double)(i % 7 + 1) / 7.0;  /* deterministic "random" */
    return q;
}

/* ══════════════════════════════════════════════════════════════════════
 * 1. Graph construction tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_graph_path(void) {
    printf("test_graph_path...\n");
    CSFGraph g = csf_graph_path(5);
    ASSERT(g.n == 5, "path n=5");
    ASSERT(g.m == 4, "path m=4");
    csf_graph_free(&g);
}

static void test_graph_cycle(void) {
    printf("test_graph_cycle...\n");
    CSFGraph g = csf_graph_cycle(4);
    ASSERT(g.n == 4, "cycle n=4");
    ASSERT(g.m == 4, "cycle m=4");
    csf_graph_free(&g);
}

static void test_graph_complete(void) {
    printf("test_graph_complete...\n");
    CSFGraph g = csf_graph_complete(4);
    ASSERT(g.n == 4, "complete n=4");
    ASSERT(g.m == 6, "complete m=6");
    csf_graph_free(&g);
}

static void test_graph_tree(void) {
    printf("test_graph_tree...\n");
    CSFEdge edges[] = {{0,1,1},{0,2,1},{1,3,1},{1,4,1}};
    CSFGraph g = csf_graph_tree(5, edges);
    ASSERT(g.n == 5, "tree n=5");
    ASSERT(g.m == 4, "tree m=4");
    csf_graph_free(&g);
}

static void test_graph_single_vertex(void) {
    printf("test_graph_single_vertex...\n");
    CSFGraph g = csf_graph_path(1);
    ASSERT(g.n == 1, "single n=1");
    ASSERT(g.m == 0, "single m=0");
    csf_graph_free(&g);
}

/* ══════════════════════════════════════════════════════════════════════
 * 2. Matrix tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_mat_identity(void) {
    printf("test_mat_identity...\n");
    CSFMatrix m = csf_mat_alloc(3);
    csf_mat_identity(&m);
    ASSERT_FEQ(m.data[0], 1.0, 1e-12, "I[0][0]=1");
    ASSERT_FEQ(m.data[4], 1.0, 1e-12, "I[1][1]=1");
    ASSERT_FEQ(m.data[8], 1.0, 1e-12, "I[2][2]=1");
    ASSERT_FEQ(m.data[1], 0.0, 1e-12, "I[0][1]=0");
    csf_mat_free(&m);
}

static void test_mat_multiply(void) {
    printf("test_mat_multiply...\n");
    CSFMatrix a = csf_mat_alloc(2);
    a.data[0]=1; a.data[1]=2; a.data[2]=3; a.data[3]=4;
    CSFMatrix b = csf_mat_alloc(2);
    b.data[0]=2; b.data[1]=0; b.data[2]=1; b.data[3]=3;
    CSFMatrix c = csf_mat_alloc(2);
    csf_mat_mul(&a, &b, &c);
    ASSERT_FEQ(c.data[0], 4.0, 1e-12, "AB[0][0]=4");
    ASSERT_FEQ(c.data[1], 6.0, 1e-12, "AB[0][1]=6");
    ASSERT_FEQ(c.data[2], 10.0, 1e-12, "AB[1][0]=10");
    ASSERT_FEQ(c.data[3], 12.0, 1e-12, "AB[1][1]=12");
    csf_mat_free(&a); csf_mat_free(&b); csf_mat_free(&c);
}

/* ══════════════════════════════════════════════════════════════════════
 * 3. Conservation sheaf tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_sheaf_create_free(void) {
    printf("test_sheaf_create_free...\n");
    CSFGraph g = csf_graph_path(3);
    ConservationSheaf s = csf_sheaf_create(2, g);
    ASSERT(s.stalk_dim == 2, "stalk_dim=2");
    ASSERT(s.laplacian_valid, "laplacian valid");
    csf_sheaf_free(&s);
}

static void test_sheaf_laplacian_structure(void) {
    printf("test_sheaf_laplacian_structure...\n");
    /* Path graph with stalk_dim=1 should give standard graph Laplacian */
    CSFGraph g = csf_graph_path(3);
    ConservationSheaf s = csf_sheaf_create(1, g);
    CSFMatrix *L = &s.laplacian;
    /* L = [[1,-1,0],[-1,2,-1],[0,-1,1]] */
    ASSERT_FEQ(*csf_mat_at(L,0,0), 1.0, 1e-12, "L[0][0]=1");
    ASSERT_FEQ(*csf_mat_at(L,0,1), -1.0, 1e-12, "L[0][1]=-1");
    ASSERT_FEQ(*csf_mat_at(L,1,1), 2.0, 1e-12, "L[1][1]=2");
    ASSERT_FEQ(*csf_mat_at(L,2,2), 1.0, 1e-12, "L[2][2]=1");
    csf_sheaf_free(&s);
}

static void test_sheaf_global_conservation(void) {
    printf("test_sheaf_global_conservation...\n");
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {1.0, 2.0, 3.0, 4.0};
    ASSERT(csf_global_conservation_check(&s, q), "global conservation holds");
    csf_sheaf_free(&s);
}

static void test_sheaf_vertex_conservation(void) {
    printf("test_sheaf_vertex_conservation...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    /* Uniform distribution: Laplacian gives 0 everywhere */
    double q[] = {1.0, 1.0, 1.0, 1.0};
    ASSERT(csf_verify_conservation_at_vertex(&s, q, 0), "vertex 0 conserved");
    ASSERT(csf_verify_conservation_at_vertex(&s, q, 2), "vertex 2 conserved");
    csf_sheaf_free(&s);
}

static void test_sheaf_kernel_constant(void) {
    printf("test_sheaf_kernel_constant...\n");
    /* Constant vector should be in kernel of Laplacian */
    CSFGraph g = csf_graph_cycle(5);
    ConservationSheaf s = csf_sheaf_create(1, g);
    int N = 5;
    double ones[5] = {1,1,1,1,1};
    double result[5];
    csf_mat_vec(&s.laplacian, ones, result);
    double norm = 0;
    for (int i = 0; i < N; i++) norm += result[i] * result[i];
    ASSERT(sqrt(norm) < 1e-10, "constant in kernel");
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * 4. Flow dynamics tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_flow_create_free(void) {
    printf("test_flow_create_free...\n");
    double q[] = {1.0, 2.0, 3.0};
    FlowState fs = csf_flow_create(1, 3, q, 0.01);
    ASSERT(fs.n_vertices == 3, "n_vertices=3");
    ASSERT(fs.stalk_dim == 1, "stalk_dim=1");
    ASSERT_FEQ(fs.quantities[1], 2.0, 1e-12, "q[1]=2");
    csf_flow_free(&fs);
}

static void test_flow_step_diffusion(void) {
    printf("test_flow_step_diffusion...\n");
    CSFGraph g = csf_graph_path(3);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {0.0, 3.0, 0.0};
    FlowState fs = csf_flow_create(1, 3, q, 0.01);

    /* After one step, middle should decrease, ends should increase */
    double mid_before = fs.quantities[1];
    csf_flow_step(&fs, &s);
    ASSERT(fs.quantities[1] < mid_before, "middle decreased");
    ASSERT(fs.quantities[0] > 0, "left increased");
    ASSERT(fs.quantities[2] > 0, "right increased");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

static void test_flow_converges(void) {
    printf("test_flow_converges...\n");
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {0.0, 0.0, 4.0, 0.0};
    FlowState fs = csf_flow_create(1, 4, q, 0.05);

    bool converged = csf_flow_converge(&fs, &s, 5000, 1e-8);
    ASSERT(converged, "flow converged");

    /* Should converge to uniform */
    double mean = 1.0;
    ASSERT_FEQ(fs.quantities[0], mean, 0.01, "near uniform [0]");
    ASSERT_FEQ(fs.quantities[3], mean, 0.01, "near uniform [3]");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

static void test_flow_entropy_increases(void) {
    printf("test_flow_entropy_increases...\n");
    CSFGraph g = csf_graph_path(5);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {5.0, 0.0, 0.0, 0.0, 0.0};
    FlowState fs = csf_flow_create(1, 5, q, 0.01);

    double e0 = csf_flow_entropy(&fs);
    for (int i = 0; i < 100; i++) csf_flow_step(&fs, &s);
    double e1 = csf_flow_entropy(&fs);

    ASSERT(e1 > e0 - 1e-10, "entropy increased after 100 steps");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

static void test_flow_total_preserved(void) {
    printf("test_flow_total_preserved...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {1.0, 2.0, 3.0, 4.0};
    FlowState fs = csf_flow_create(1, 4, q, 0.05);

    double total_before = 0;
    for (int i = 0; i < 4; i++) total_before += fs.quantities[i];

    for (int i = 0; i < 50; i++) csf_flow_step(&fs, &s);

    double total_after = 0;
    for (int i = 0; i < 4; i++) total_after += fs.quantities[i];

    ASSERT_FEQ(total_before, total_after, 0.01, "total quantity preserved");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * 5. Spectral gap tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_spectral_gap_path(void) {
    printf("test_spectral_gap_path...\n");
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double gap = csf_flow_spectral_gap(&s);
    /* Path(4) spectral gap = 2(1 - cos(pi/4)) = 2 - sqrt(2) ≈ 0.5858 */
    ASSERT(gap > 0, "spectral gap positive");
    ASSERT_FEQ(gap, 0.5858, 0.02, "path(4) spectral gap ≈ 0.586");
    printf("    path(4) spectral gap = %.6f (expected ≈ 0.5858)\n", gap);
    csf_sheaf_free(&s);
}

static void test_spectral_gap_cycle(void) {
    printf("test_spectral_gap_cycle...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double gap = csf_flow_spectral_gap(&s);
    /* Cycle(4) spectral gap = 2(1 - cos(pi/2)) = 2 */
    ASSERT(gap > 0, "spectral gap positive");
    ASSERT_FEQ(gap, 2.0, 0.05, "cycle(4) spectral gap ≈ 2.0");
    printf("    cycle(4) spectral gap = %.6f (expected ≈ 2.0)\n", gap);
    csf_sheaf_free(&s);
}

static void test_spectral_gap_complete(void) {
    printf("test_spectral_gap_complete...\n");
    CSFGraph g = csf_graph_complete(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double gap = csf_flow_spectral_gap(&s);
    /* K4 spectral gap = 4 */
    ASSERT(gap > 0, "spectral gap positive");
    ASSERT_FEQ(gap, 4.0, 0.05, "K4 spectral gap ≈ 4.0");
    printf("    K4 spectral gap = %.6f (expected ≈ 4.0)\n", gap);
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * 6. Transport tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_transport_cost_identity(void) {
    printf("test_transport_cost_identity...\n");
    double q[] = {1.0, 2.0, 3.0};
    FlowState a = csf_flow_create(1, 3, q, 0.01);
    FlowState b = csf_flow_create(1, 3, q, 0.01);
    double cost = csf_transport_cost(&a, &b);
    ASSERT_FEQ(cost, 0.0, 1e-12, "zero cost for identical states");
    csf_flow_free(&a); csf_flow_free(&b);
}

static void test_transport_cost_symmetric(void) {
    printf("test_transport_cost_symmetric...\n");
    double q1[] = {1.0, 0.0, 0.0};
    double q2[] = {0.0, 0.0, 1.0};
    FlowState a = csf_flow_create(1, 3, q1, 0.01);
    FlowState b = csf_flow_create(1, 3, q2, 0.01);
    double c1 = csf_transport_cost(&a, &b);
    double c2 = csf_transport_cost(&b, &a);
    ASSERT_FEQ(c1, c2, 1e-12, "transport cost symmetric");
    csf_flow_free(&a); csf_flow_free(&b);
}

static void test_transport_plan(void) {
    printf("test_transport_plan...\n");
    double q1[] = {2.0, 0.0};
    double q2[] = {0.0, 2.0};
    FlowState a = csf_flow_create(1, 2, q1, 0.01);
    FlowState b = csf_flow_create(1, 2, q2, 0.01);
    double *plan = csf_transport_plan(&a, &b);
    /* Should transfer from vertex 0 to vertex 1 */
    ASSERT(plan[0 * 2 + 1] > 0, "transport from 0 to 1");
    csf_flow_free(&a); csf_flow_free(&b);
    free(plan);
}

static void test_barycenter(void) {
    printf("test_barycenter...\n");
    double q1[] = {1.0, 0.0};
    double q2[] = {0.0, 1.0};
    FlowState a = csf_flow_create(1, 2, q1, 0.01);
    FlowState b = csf_flow_create(1, 2, q2, 0.01);
    const FlowState *targets[] = {&a, &b};
    FlowState bc = csf_barycenter_flow(targets, 2);
    ASSERT_FEQ(bc.quantities[0], 0.5, 1e-12, "barycenter[0]=0.5");
    ASSERT_FEQ(bc.quantities[1], 0.5, 1e-12, "barycenter[1]=0.5");
    csf_flow_free(&a); csf_flow_free(&b); csf_flow_free(&bc);
}

/* ══════════════════════════════════════════════════════════════════════
 * 7. THE THEOREM — Non-decreasing spectral gap
 * ══════════════════════════════════════════════════════════════════════ */
static void test_theorem_path(void) {
    printf("test_theorem_path (THE THEOREM)...\n");
    CSFGraph g = csf_graph_path(5);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {5.0, 0.0, 0.0, 0.0, 0.0};
    FlowState fs = csf_flow_create(1, 5, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 50);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "THEOREM HOLDS on path(5)");
    printf("    spectral gap = %.6f (constant = non-decreasing ✓)\n",
           ev.spectral_gaps[0]);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_theorem_cycle(void) {
    printf("test_theorem_cycle (THE THEOREM)...\n");
    CSFGraph g = csf_graph_cycle(6);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {3.0, 0.0, 0.0, 1.0, 0.0, 2.0};
    FlowState fs = csf_flow_create(1, 6, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 50);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "THEOREM HOLDS on cycle(6)");
    printf("    spectral gap = %.6f ✓\n", ev.spectral_gaps[0]);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_theorem_tree(void) {
    printf("test_theorem_tree (THE THEOREM)...\n");
    CSFEdge edges[] = {{0,1,1},{0,2,1},{0,3,1},{1,4,1},{1,5,1}};
    CSFGraph g = csf_graph_tree(6, edges);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {0.0, 2.0, 0.0, 0.0, 1.0, 3.0};
    FlowState fs = csf_flow_create(1, 6, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 50);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "THEOREM HOLDS on tree");
    printf("    spectral gap = %.6f ✓\n", ev.spectral_gaps[0]);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_theorem_complete(void) {
    printf("test_theorem_complete (THE THEOREM)...\n");
    CSFGraph g = csf_graph_complete(5);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {4.0, 0.0, 0.0, 0.0, 1.0};
    FlowState fs = csf_flow_create(1, 5, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 50);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "THEOREM HOLDS on K5");
    printf("    spectral gap = %.6f ✓\n", ev.spectral_gaps[0]);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_theorem_multidim(void) {
    printf("test_theorem_multidim (THE THEOREM, stalk_dim=2)...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(2, g);
    double q[] = {1.0, 0.5, 0.0, 0.0, 2.0, 0.3, 0.0, 0.0};
    FlowState fs = csf_flow_create(2, 4, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 50);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "THEOREM HOLDS multidim cycle");
    printf("    spectral gap = %.6f ✓\n", ev.spectral_gaps[0]);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * 8. Critical time tests
 * ══════════════════════════════════════════════════════════════════════ */
static void test_critical_time(void) {
    printf("test_critical_time...\n");
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {0.0, 0.0, 4.0, 0.0};
    FlowState fs = csf_flow_create(1, 4, q, 0.05);

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 200);
    double tc = csf_compute_critical_time(&ev);
    ASSERT(tc >= 0, "critical time non-negative");
    printf("    critical time = %.4f\n", tc);

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * 9. Edge cases
 * ══════════════════════════════════════════════════════════════════════ */
static void test_single_vertex(void) {
    printf("test_single_vertex...\n");
    CSFGraph g = csf_graph_path(1);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {5.0};
    FlowState fs = csf_flow_create(1, 1, q, 0.01);

    /* No edges, Laplacian is 1x1 zero matrix */
    ASSERT_FEQ(s.laplacian.data[0], 0.0, 1e-12, "L is zero for single vertex");

    csf_flow_step(&fs, &s);
    ASSERT_FEQ(fs.quantities[0], 5.0, 1e-12, "quantity unchanged (no neighbors)");

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 10);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "theorem holds trivially");

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_two_vertices(void) {
    printf("test_two_vertices...\n");
    CSFGraph g = csf_graph_path(2);
    ConservationSheaf s = csf_sheaf_create(1, g);

    double q[] = {1.0, 0.0};
    FlowState fs = csf_flow_create(1, 2, q, 0.01);

    for (int i = 0; i < 200; i++) csf_flow_step(&fs, &s);

    /* Should converge to (0.5, 0.5) */
    ASSERT_FEQ(fs.quantities[0], 0.5, 0.05, "converged to uniform [0]");
    ASSERT_FEQ(fs.quantities[1], 0.5, 0.05, "converged to uniform [1]");

    SpectralEvolution ev = csf_track_spectral_gap(&s, &fs, 10);
    ASSERT(csf_verify_non_decreasing_gap(&ev), "theorem holds on 2 vertices");

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

static void test_uniform_initial(void) {
    printf("test_uniform_initial...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {1.0, 1.0, 1.0, 1.0};
    FlowState fs = csf_flow_create(1, 4, q, 0.01);

    /* Uniform is already steady state */
    for (int i = 0; i < 10; i++) csf_flow_step(&fs, &s);
    for (int i = 0; i < 4; i++)
        ASSERT_FEQ(fs.quantities[i], 1.0, 0.01, "uniform stays uniform");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

static void test_laplacian_psd(void) {
    printf("test_laplacian_psd...\n");
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);

    /* Check positive semi-definite: x^T L x >= 0 for random x */
    double x[] = {1.0, -2.0, 3.0, -1.0};
    double Lx[4];
    csf_mat_vec(&s.laplacian, x, Lx);
    double xtLx = 0;
    for (int i = 0; i < 4; i++) xtLx += x[i] * Lx[i];
    ASSERT(xtLx >= -1e-10, "Laplacian is PSD");

    csf_sheaf_free(&s);
}

static void test_zero_vector_steady(void) {
    printf("test_zero_vector_steady...\n");
    CSFGraph g = csf_graph_path(3);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {0.0, 0.0, 0.0};
    FlowState fs = csf_flow_create(1, 3, q, 0.01);

    for (int i = 0; i < 10; i++) csf_flow_step(&fs, &s);
    for (int i = 0; i < 3; i++)
        ASSERT_FEQ(fs.quantities[i], 0.0, 1e-12, "zero stays zero");

    csf_flow_free(&fs);
    csf_sheaf_free(&s);
}

static void test_flow_spectral_gap_positive(void) {
    printf("test_flow_spectral_gap_positive...\n");
    /* Spectral gap should be positive for connected graphs */
    CSFGraph g = csf_graph_cycle(5);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double gap = csf_flow_spectral_gap(&s);
    ASSERT(gap > 1e-6, "spectral gap positive for connected graph");
    printf("    gap = %.6f\n", gap);
    csf_sheaf_free(&s);
}

static void test_entropy_bounds(void) {
    printf("test_entropy_bounds...\n");
    /* Uniform should have entropy ~1.0, concentrated ~0 */
    double uniform[] = {1.0, 1.0, 1.0, 1.0};
    FlowState fu = csf_flow_create(1, 4, uniform, 0.01);
    double eu = csf_flow_entropy(&fu);
    ASSERT(eu > 0.99, "uniform entropy near 1");
    printf("    uniform entropy = %.6f\n", eu);

    double conc[] = {4.0, 0.0, 0.0, 0.0};
    FlowState fc = csf_flow_create(1, 4, conc, 0.01);
    double ec = csf_flow_entropy(&fc);
    ASSERT(ec < 0.5, "concentrated entropy < 0.5");
    printf("    concentrated entropy = %.6f\n", ec);

    ASSERT(eu > ec, "uniform entropy > concentrated");

    csf_flow_free(&fu);
    csf_flow_free(&fc);
}

/* ══════════════════════════════════════════════════════════════════════
 * 10. Double-free safety test
 * ══════════════════════════════════════════════════════════════════════ */
static void test_no_double_free(void) {
    printf("test_no_double_free...\n");
    /* Sheaf takes a deep copy of the graph, so freeing both should be safe */
    CSFGraph g = csf_graph_path(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    /* Free the original graph — sheaf has its own copy */
    csf_graph_free(&g);
    /* Verify sheaf still has valid data */
    ASSERT(s.graph.n == 4, "sheaf graph still valid after original freed");
    ASSERT(s.laplacian_valid, "sheaf laplacian still valid");
    /* Now free the sheaf — should not double-free */
    csf_sheaf_free(&s);
    ASSERT(1, "no crash from double-free");
}

/* ══════════════════════════════════════════════════════════════════════
 * 11. Evolving sheaf — EXPERIMENTAL (not a proven theorem)
 * ══════════════════════════════════════════════════════════════════════ */
static void test_evolving_sheaf_experimental(void) {
    printf("test_evolving_sheaf_experimental (OPEN QUESTION)...\n");
    CSFGraph g = csf_graph_cycle(4);
    ConservationSheaf s = csf_sheaf_create(1, g);
    double q[] = {3.0, 0.0, 1.0, 0.0};
    FlowState fs = csf_flow_create(1, 4, q, 0.01);

    SpectralEvolution ev = csf_track_spectral_gap_evolving(&s, &fs, 50);
    printf("    evolving sheaf: gap[0]=%.3f, gap[last]=%.3f\n",
           ev.spectral_gaps[0], ev.spectral_gaps[ev.n_steps - 1]);
    printf("    theorem_holds=%s (EXPERIMENTAL — not a proven theorem)\n",
           ev.theorem_holds ? "true" : "false");
    /* We do NOT assert theorem_holds here — it's an open question */
    ASSERT(ev.n_steps == 50, "evolving sheaf tracked all steps");

    csf_flow_free(&fs);
    csf_spectral_evolution_free(&ev);
    csf_sheaf_free(&s);
}

/* ══════════════════════════════════════════════════════════════════════
 * MAIN
 * ══════════════════════════════════════════════════════════════════════ */
int main(void) {
    printf("══════════════════════════════════════════════════════════\n");
    printf("  Conservation Sheaf Flow — Test Suite\n");
    printf("══════════════════════════════════════════════════════════\n\n");

    printf("── Graph construction ──\n");
    test_graph_path();
    test_graph_cycle();
    test_graph_complete();
    test_graph_tree();
    test_graph_single_vertex();

    printf("\n── Matrix utilities ──\n");
    test_mat_identity();
    test_mat_multiply();

    printf("\n── Conservation sheaf ──\n");
    test_sheaf_create_free();
    test_sheaf_laplacian_structure();
    test_sheaf_global_conservation();
    test_sheaf_vertex_conservation();
    test_sheaf_kernel_constant();

    printf("\n── Flow dynamics ──\n");
    test_flow_create_free();
    test_flow_step_diffusion();
    test_flow_converges();
    test_flow_entropy_increases();
    test_flow_total_preserved();

    printf("\n── Spectral gap ──\n");
    test_spectral_gap_path();
    test_spectral_gap_cycle();
    test_spectral_gap_complete();
    test_flow_spectral_gap_positive();

    printf("\n── Optimal transport ──\n");
    test_transport_cost_identity();
    test_transport_cost_symmetric();
    test_transport_plan();
    test_barycenter();

    printf("\n════════════════════════════════════════════════════════\n");
    printf("  THE THEOREM — Non-decreasing spectral gap\n");
    printf("════════════════════════════════════════════════════════\n\n");

    printf("── Theorem verification on multiple topologies ──\n");
    test_theorem_path();
    test_theorem_cycle();
    test_theorem_tree();
    test_theorem_complete();
    test_theorem_multidim();

    printf("\n── Critical time ──\n");
    test_critical_time();

    printf("\n── Edge cases ──\n");
    test_single_vertex();
    test_two_vertices();
    test_uniform_initial();
    test_laplacian_psd();
    test_zero_vector_steady();
    test_entropy_bounds();

    printf("\n── Double-free safety ──\n");
    test_no_double_free();

    printf("\n══════════════════════════════════════════════════════════\n");
    printf("  Evolving Sheaf — EXPERIMENTAL (Open Question)\n");
    printf("══════════════════════════════════════════════════════════\n\n");
    test_evolving_sheaf_experimental();

    printf("\n══════════════════════════════════════════════════════════\n");
    printf("  Results: %d passed, %d failed, %d total\n",
           tests_passed, tests_failed, tests_passed + tests_failed);
    printf("══════════════════════════════════════════════════════════\n");

    return tests_failed > 0 ? 1 : 0;
}
