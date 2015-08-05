
#include <stdlib.h>

struct cell {
    unsigned int i,j,k;
};

typedef struct {
    unsigned int num_cells;
    struct cell *cells;
    // index is an index into the cells array for when storing the cell indexes
    unsigned int index;
} plane;

// Compute the order of the sweep for the first octant
plane *compute_sweep_order(int nx, int ny, int nz)
{
    unsigned int nplanes = nx + ny + nz - 2;
    plane *planes = (plane *)malloc(sizeof(plane)*nplanes);
    for (unsigned int i = 0; i < nplanes; i++)
    {
        planes[i].num_cells = 0;
    }

    // Cells on each plane have equal co-ordinate sum
    for (unsigned int k = 0; k < nz; k++)
    {
        for (unsigned int j = 0; j < ny; j++)
        {
            for (unsigned int i = 0; i < nx; i++)
            {
                unsigned int n = i + j + k;
                planes[n].num_cells++;
            }
        }
    }

    // Allocate the memory for each plane
    for (unsigned int i = 0; i < nplanes; i++)
    {
        planes[i].cells = (struct cell *)malloc(sizeof(struct cell)*planes[i].num_cells);
        planes[i].index = 0;
    }

    // Store the cell indexes in the plane array
    for (unsigned int k = 0; k < nz; k++)
    {
        for (unsigned int j = 0; j < ny; j++)
        {
            for (unsigned int i = 0; i < nx; i++)
            {
                unsigned int n = i + j + k;
                unsigned int idx = planes[n].index;
                planes[n].cells[idx].i = i;
                planes[n].cells[idx].j = j;
                planes[n].cells[idx].k = k;
                planes[n].index += 1;
            }
        }
    }

    return planes;
}

#define qtot(m,i,j,k,g) qtot[(m)+(cmom*(i))+(cmom*nx*(j))+(cmom*nx*ny*(k))+(cmom*nx*ny*nz*(g))]
#define ec(a,l,o) ec[(a)+(nang*(l))+(nang*cmom*(o))]

#define flux_i(a,g,j,k) flux_i[(a)+(nang*(g))+(nang*ng*(j))+(nang*ng*ny*(k))]
// #define flux_i(idx,j,k) flux_i[(idx)+(nang*ng*(j))+(nang*ng*ny*(k))]

#define flux_j(a,g,i,k) flux_j[(a)+(nang*(g))+(nang*ng*(i))+(nang*ng*nx*(k))]
// #define flux_j(idx,i,k) flux_j[(idx)+(nang*ng*(i))+(nang*ng*nx*(k))]

#define flux_k(a,g,i,j) flux_k[(a)+(nang*(g))+(nang*ng*(i))+(nang*ng*nx*(j))]
// #define flux_k(idx,i,j) flux_k[(idx)+(nang*ng*(i))+(nang*ng*nx*(j))]

#define mu(a) mu[(a)]
#define hj(a) hj[(a)]
#define hk(a) hk[(a)]
#define vdelt(g) vdelt[(g)]

#define flux_in(a,g,i,j,k,o) flux_in[(a)+(nang*(g))+(nang*ng*(i))+(nang*ng*nx*(j))+(nang*ng*nx*ny*(k))+(nang*ng*nx*ny*nz*(o))]
// #define flux_in(idx,i,j,k,o) flux_in[(idx)+(nang*ng*(i))+(nang*ng*nx*(j))+(nang*ng*nx*ny*(k))+(nang*ng*nx*ny*nz*(o))]

#define flux_out(a,g,i,j,k,o) flux_out[(a)+(nang*(g))+(nang*ng*(i))+(nang*ng*nx*(j))+(nang*ng*nx*ny*(k))+(nang*ng*nx*ny*nz*(o))]
// #define flux_out(idx,i,j,k,o) flux_out[(idx)+(nang*ng*(i))+(nang*ng*nx*(j))+(nang*ng*nx*ny*(k))+(nang*ng*nx*ny*nz*(o))]

#define dinv(a,i,j,k,g) dinv[(a)+(nang*(i))+(nang*nx*(j))+(nang*nx*ny*(k))+(nang*nx*ny*nz*(g))]

#define t_xs(i,j,k,g) t_xs[(i)+(nx*(j))+(nx*ny*(k))+(nx*ny*nz*(g))]

void omp_sweep_c_(int *p_, int *ng_, int *nang_, int *nx_, int *ny_, int *nz_, int *istep_, int *jstep_, int *kstep_, int *o_, int *noct_, int *cmom_,
                         const double * restrict qtot, const double * restrict ec, const double * restrict mu, const double * restrict hi, const double * restrict hj, const double * restrict hk, const double * restrict vdelt, const double * restrict dinv,
                         const double * restrict t_xs,
                         double * restrict flux_i, double * restrict flux_j, double * restrict flux_k, double * restrict flux_in, double * restrict flux_out
    )
{

    const int nx = *nx_;
    const int ny = *ny_;
    const int nz = *nz_;
    const int ng = *ng_;
    const int nang = *nang_;
    const int noct = *noct_;
    const int cmom = *cmom_;
    const int istep = *istep_;
    const int jstep = *jstep_;
    const int kstep = *kstep_;

    const int p = *p_ - 1;
    const int o = *o_ - 1;

    plane *sweep_order = compute_sweep_order(nx, ny, nz);

    #pragma novector
    #pragma omp parallel for
    for (int c = 0; c < sweep_order[p].num_cells; c++)
    {
        const int i = (istep > 0) ? sweep_order[p].cells[c].i : nx - sweep_order[p].cells[c].i - 1;
        const int j = (jstep > 0) ? sweep_order[p].cells[c].j : ny - sweep_order[p].cells[c].j - 1;
        const int k = (kstep > 0) ? sweep_order[p].cells[c].k : nz - sweep_order[p].cells[c].k - 1;

        #pragma ivdep
        for (int g = 0; g < ng; g++)
            //#pragma ivdep
            #pragma omp simd
            for (int a = 0; a < nang; a++)
            {
                double source = qtot(0,i,j,k,g);
                #pragma novector
                for (int l = 1; l < cmom; l++)
                    source += ec(a,l,o) * qtot(l,i,j,k,g);

                double psi = source + (flux_i(a,g,j,k) * mu(a) * (*hi)) + (flux_j(a,g,i,k) * hj(a)) + (flux_k(a,g,i,j) * hk(a));

                if (vdelt(g) != 0.0)
                    psi += vdelt(g) * flux_in(a,g,i,j,k,o);

                psi *= dinv(a,i,j,k,g);

                double tmp_flux_i = 2.0 * psi - flux_i(a,g,j,k);
                double tmp_flux_j = 2.0 * psi - flux_j(a,g,i,k);
                double tmp_flux_k = 2.0 * psi - flux_k(a,g,i,j);

                if (vdelt(g) != 0.0)
                    psi = 2.0*psi - flux_in(a,g,i,j,k,o);

                // FIXUP

                int num_to_fix = 4;

                double zeros[4];// = {1.0, 1.0, 1.0, 1.0};

                #pragma novector
                for (int fix = 0; fix < 4; fix++)
                    zeros[fix] = 1.0;

                #pragma unroll (4)
                #pragma novector
                for (int fix = 0; fix < 4; fix++)
                {
                    if (tmp_flux_i < 0.0)
                        zeros[0] = 0.0;
                    else
                        zeros[0] = 1.0;

                    if (tmp_flux_j < 0.0)
                        zeros[1] = 0.0;
                    else
                        zeros[1] = 1.0;

                    if (tmp_flux_k < 0.0)
                        zeros[2] = 0.0;
                    else
                        zeros[2] = 1.0;

                    if (psi < 0.0)
                        zeros[3] = 0.0;
                    else
                        zeros[3] = 1.0;

                    if (num_to_fix == zeros[0] + zeros[1] + zeros[2] + zeros[3])
                        continue;

                    num_to_fix = zeros[0] + zeros[1] + zeros[2] + zeros[3];

                    psi =
                        flux_i(a,g,j,k)*mu(a)*(*hi)*(1.0+zeros[0]) +
                        flux_j(a,g,j,k)*hj(a)*(1.0+zeros[1]) +
                        flux_k(a,g,i,j)*hk(a)*(1.0+zeros[2]);

                    if (vdelt(g) != 0.0)
                        psi += vdelt(g) * flux_in(a,g,i,j,k,o) *  (1.0+zeros[3]);

                    psi = 0.5*psi + source;
                    double recalc_denom = t_xs(i,j,k,g);
                    recalc_denom += mu(a) * (*hi) * zeros[0];
                    recalc_denom += hj(a) * zeros[1];
                    recalc_denom += hk(a) * zeros[2];
                    recalc_denom += vdelt(g) * zeros[3];

                    if (recalc_denom > 1.0E-12)
                        psi /= recalc_denom;
                    else
                        psi = 0.0;

                    tmp_flux_i = 2.0 * psi - flux_i(a,g,j,k);
                    tmp_flux_j = 2.0 * psi - flux_j(a,g,i,k);
                    tmp_flux_k = 2.0 * psi - flux_k(a,g,i,j);

                    if (vdelt(g) != 0.0)
                        psi = 2.0*psi - flux_in(a,g,i,j,k,o);

                }

                tmp_flux_i *= zeros[0];
                tmp_flux_j *= zeros[1];
                tmp_flux_k *= zeros[2];
                psi *= zeros[3];

                flux_i(a,g,j,k) = tmp_flux_i;
                flux_j(a,g,i,k) = tmp_flux_j;
                flux_k(a,g,i,j) = tmp_flux_k;
                flux_out(a,g,i,j,k,o) = psi;
            }
    }

}
