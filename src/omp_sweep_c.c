
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

#define qtot(m,i,j,k,g) qtot[(m)+((*cmom)*(i))+((*cmom)*(*nx)*(j))+((*cmom)*(*nx)*(*ny)*(k))+((*cmom)*(*nx)*(*ny)*(*nz)*(g))]
#define ec(a,l,o) ec[(a)+((*nang)*(l))+((*nang)*(*cmom)*(o))]
#define flux_i(a,g,j,k) flux_i[(a)+((*nang)*(g))+((*nang)*(*ng)*(j))+((*nang)*(*ng)*(*ny)*(k))]
#define flux_j(a,g,i,k) flux_j[(a)+((*nang)*(g))+((*nang)*(*ng)*(i))+((*nang)*(*ng)*(*nx)*(k))]
#define flux_k(a,g,i,j) flux_k[(a)+((*nang)*(g))+((*nang)*(*ng)*(i))+((*nang)*(*ng)*(*nx)*(j))]
#define mu(a) mu[(a)]
#define hj(a) hj[(a)]
#define hk(a) hk[(a)]
#define vdelt(g) vdelt[(g)]
#define flux_in(a,g,i,j,k,o) flux_in[(a)+((*nang)*(g))+((*nang)*(*ng)*(i))+((*nang)*(*ng)*(*nx)*(j))+((*nang)*(*ng)*(*nx)*(*ny)*(k))+((*nang)*(*ng)*(*nx)*(*ny)*(*nz)*(o))]
#define flux_out(a,g,i,j,k,o) flux_out[(a)+((*nang)*(g))+((*nang)*(*ng)*(i))+((*nang)*(*ng)*(*nx)*(j))+((*nang)*(*ng)*(*nx)*(*ny)*(k))+((*nang)*(*ng)*(*nx)*(*ny)*(*nz)*(o))]
#define dinv(a,i,j,k,g) dinv[(a)+((*nang)*(i))+((*nang)*(*nx)*(j))+((*nang)*(*nx)*(*ny)*(k))+((*nang)*(*nx)*(*ny)*(*nz)*(g))]

void omp_sweep_c_(const int *p, const int *ng, const int *nang, const int *nx, const int *ny, const int *nz, const int *istep, const int *jstep, const int *kstep, const int *o, const int *noct, const int *cmom,
                         const double *qtot, const double *ec, const double *mu, const double *hi, const double *hj, const double *hk, const double *vdelt, const double *dinv,
                         double * restrict flux_i, double * restrict flux_j, double * restrict flux_k, double * restrict flux_in, double * restrict flux_out
    )
{

    plane *sweep_order = compute_sweep_order(*nx, *ny, *nz);
    for (int c = 0; c < sweep_order[*p-1].num_cells; c++)
    {
        int i, j, k;
        if (*istep > 0)
            i = sweep_order[*p-1].cells[c].i;
        else
            i = (*nx) - sweep_order[*p-1].cells[c].i - 1;
        if (*jstep > 0)
            j = sweep_order[*p-1].cells[c].j;
        else
            j = (*ny) - sweep_order[*p-1].cells[c].j - 1;
        if (*kstep > 0)
            k = sweep_order[*p-1].cells[c].k;
        else
            k = (*nz) - sweep_order[*p-1].cells[c].k - 1;

        for (int g = 0; g < (*ng); g++)
            for (int a = 0; a < (*nang); a++)
            {
                double source = qtot(0,i,j,k,g);
                #pragma novector
                for (int l = 1; l < *cmom; l++)
                    source += ec(a,l,(*o)-1) * qtot(l,i,j,k,g);

                double psi = source + (flux_i(a,g,j,k) * mu(a) * (*hi)) + (flux_j(a,g,i,k) * hj(a)) + (flux_k(a,g,i,j) * hk(a));

                if (vdelt(g) != 0.0)
                    psi += vdelt(g) * flux_in(a,g,i,j,k,(*o)-1);

                psi *= dinv(a,i,j,k,g);

                flux_i(a,g,j,k) = 2.0 * psi - flux_i(a,g,j,k);
                flux_j(a,g,i,k) = 2.0 * psi - flux_j(a,g,i,k);
                flux_k(a,g,i,j) = 2.0 * psi - flux_k(a,g,i,j);

                if (vdelt(g) != 0.0)
                    psi = 2.0*psi - flux_in(a,g,i,j,k,(*o)-1);

                flux_out(a,g,i,j,k,(*o)-1) = psi;
            }
    }

}
