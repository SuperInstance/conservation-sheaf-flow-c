#include "csf.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ── Flow dynamics ────────────────────────────────────────────────── */

void csf_flow_step(FlowState *fs, const ConservationSheaf *s) {
    int N = fs->n_vertices * fs->stalk_dim;

    /* q_{t+1} = q_t - dt * L * q_t  (forward Euler on heat equation) */
    double *Lq = malloc(N * sizeof(double));
    csf_mat_vec(&s->laplacian, fs->quantities, Lq);

    for (int i = 0; i < N; i++)
        fs->quantities[i] -= fs->dt * Lq[i];

    fs->time += fs->dt;
    free(Lq);
}

bool csf_flow_converge(FlowState *fs, const ConservationSheaf *s,
                       int max_steps, double tol) {
    int N = fs->n_vertices * fs->stalk_dim;
    double *prev = malloc(N * sizeof(double));

    for (int step = 0; step < max_steps; step++) {
        memcpy(prev, fs->quantities, N * sizeof(double));
        csf_flow_step(fs, s);

        /* Check convergence: max change */
        double max_change = 0;
        for (int i = 0; i < N; i++) {
            double delta = fabs(fs->quantities[i] - prev[i]);
            if (delta > max_change) max_change = delta;
        }
        if (max_change < tol) {
            free(prev);
            return true;
        }
    }
    free(prev);
    return false;
}
