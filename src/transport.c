#include "csf.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── Optimal transport on conservation sheaves ────────────────────── */

double csf_transport_cost(const FlowState *a, const FlowState *b) {
    int N = a->n_vertices * a->stalk_dim;
    double cost = 0;
    for (int i = 0; i < N; i++) {
        double d = a->quantities[i] - b->quantities[i];
        cost += d * d;
    }
    return sqrt(cost);
}

/*
 * Transport plan for scalar (d=1) quantities:
 * Sinkhorn-style or simple greedy matching.
 * For d=1, optimal transport on a graph reduces to matching
 * surpluses to deficits.
 *
 * Returns an n x n matrix (row = source, col = destination).
 */
double *csf_transport_plan(const FlowState *src, const FlowState *dst) {
    int n = src->n_vertices;
    int d = src->stalk_dim;
    double *plan = calloc((size_t)n * n, sizeof(double));

    if (d != 1) {
        /* For multi-dimensional, use identity transport (each to each) */
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int c = 0; c < d; c++) s += src->quantities[i * d + c];
            double t = 0;
            for (int c = 0; c < d; c++) t += dst->quantities[i * d + c];
            plan[i * n + i] = fabs(s - t);
        }
        return plan;
    }

    /* Greedy matching of surpluses to deficits */
    double *surplus = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        surplus[i] = src->quantities[i] - dst->quantities[i];

    /* Simple north-west corner style: match positive to negative */
    for (int i = 0; i < n; i++) {
        if (surplus[i] <= 0) continue;
        double remaining = surplus[i];
        for (int j = 0; j < n && remaining > 1e-15; j++) {
            if (surplus[j] >= 0) continue;
            double transfer = fmin(remaining, -surplus[j]);
            plan[i * n + j] = transfer;
            remaining -= transfer;
            surplus[j] += transfer;
        }
    }

    free(surplus);
    return plan;
}

/*
 * Barycenter: for each vertex/component, take the mean of all target states.
 * This minimizes sum of squared distances (which equals Wasserstein for
 * the Euclidean ground metric).
 */
FlowState csf_barycenter_flow(const FlowState **targets, int k) {
    if (k == 0) {
        FlowState fs = {0};
        return fs;
    }

    int n = targets[0]->n_vertices;
    int d = targets[0]->stalk_dim;
    int N = n * d;

    double *avg = calloc(N, sizeof(double));
    for (int t = 0; t < k; t++)
        for (int i = 0; i < N; i++)
            avg[i] += targets[t]->quantities[i];
    for (int i = 0; i < N; i++)
        avg[i] /= k;

    FlowState result = csf_flow_create(d, n, avg, targets[0]->dt);
    free(avg);
    return result;
}
