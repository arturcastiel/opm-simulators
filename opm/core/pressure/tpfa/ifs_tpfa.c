#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/linalg/sparse_sys.h>

#include <opm/core/newwells.h>
#include <opm/core/pressure/flow_bc.h>
#include <opm/core/pressure/tpfa/ifs_tpfa.h>


struct ifs_tpfa_impl {
    double *fgrav;              /* Accumulated grav contrib/face */

    /* Linear storage */
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
impl_deallocate(struct ifs_tpfa_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        free(pimpl->ddata);
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
static struct ifs_tpfa_impl *
impl_allocate(struct UnstructuredGrid *G,
              struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_impl *new;

    size_t nnu;
    size_t ddata_sz;

    nnu = G->number_of_cells;
    if (W != NULL) {
        nnu += W->number_of_wells;
    }

    ddata_sz  = 2 * nnu;                 /* b, x */
    ddata_sz += 1 * G->number_of_faces;  /* fgrav */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->ddata = malloc(ddata_sz * sizeof *new->ddata);

        if (new->ddata == NULL) {
            impl_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct CSRMatrix *
ifs_tpfa_construct_matrix(struct UnstructuredGrid *G,
                          struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2, w, i, nc, nnu;
    size_t nnz;

    struct CSRMatrix *A;

    nc = nnu = G->number_of_cells;
    if (W != NULL) {
        nnu += W->number_of_wells;
    }

    A = csrmatrix_new_count_nnz(nnu);

    if (A != NULL) {
        /* Self connections */
        for (c1 = 0; c1 < nnu; c1++) {
            A->ia[ c1 + 1 ] = 1;
        }

        /* Other connections */
        for (f = 0; f < G->number_of_faces; f++) {
            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            if ((c1 >= 0) && (c2 >= 0)) {
                A->ia[ c1 + 1 ] += 1;
                A->ia[ c2 + 1 ] += 1;
            }
        }

        if (W != NULL) {
            /* Well <-> cell connections */
            for (w = i = 0; w < W->number_of_wells; w++) {
                for (; i < W->well_connpos[w + 1]; i++) {
                    c1 = W->well_cells[i];

                    A->ia[ 0  + c1 + 1 ] += 1; /* c -> w */
                    A->ia[ nc + w  + 1 ] += 1; /* w -> c */
                }
            }
        }

        nnz = csrmatrix_new_elms_pushback(A);
        if (nnz == 0) {
            csrmatrix_delete(A);
            A = NULL;
        }
    }

    if (A != NULL) {
        /* Fill self connections */
        for (i = 0; i < nnu; i++) {
            A->ja[ A->ia[ i + 1 ] ++ ] = i;
        }

        /* Fill other connections */
        for (f = 0; f < G->number_of_faces; f++) {
            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            if ((c1 >= 0) && (c2 >= 0)) {
                A->ja[ A->ia[ c1 + 1 ] ++ ] = c2;
                A->ja[ A->ia[ c2 + 1 ] ++ ] = c1;
            }
        }

        if (W != NULL) {
            /* Fill well <-> cell connections */
            for (w = i = 0; w < W->number_of_wells; w++) {
                for (; i < W->well_connpos[w + 1]; i++) {
                    c1 = W->well_cells[i];

                    A->ja[ A->ia[ 0  + c1 + 1 ] ++ ] = nc + w;
                    A->ja[ A->ia[ nc + w  + 1 ] ++ ] = c1    ;
                }
            }
        }

        assert ((size_t) A->ia[ nnu ] == nnz);

        /* Guarantee sorted rows */
        csrmatrix_sortrows(A);
    }

    return A;
}


/* ---------------------------------------------------------------------- */
/* fgrav = accumarray(cf(j), grav(j).*sgn(j), [nf, 1]) */
/* ---------------------------------------------------------------------- */
static void
compute_grav_term(struct UnstructuredGrid *G, const double *gpress,
                  double *fgrav)
/* ---------------------------------------------------------------------- */
{
    int    c, i, f, c1, c2;
    double s;

    vector_zero(G->number_of_faces, fgrav);

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f  = G->cell_faces[i];

            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            s  = 2.0*(c1 == c) - 1.0;

            if ((c1 >= 0) && (c2 >= 0)) {
                fgrav[f] += s * gpress[i];
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
assemble_bc_contrib(struct UnstructuredGrid             *G    ,
                    const struct FlowBoundaryConditions *bc   ,
                    const double                        *trans,
                    struct ifs_tpfa_data                *h    )
/* ---------------------------------------------------------------------- */
{
    int is_neumann, is_outflow;
    int f, c1, c2;

    size_t i, j, ix;
    double s, t;

    is_neumann = 1;

    for (i = 0; i < bc->nbc; i++) {
        if (bc->type[ i ] == BC_PRESSURE) {
            is_neumann = 0;

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[i + 1]; j++) {
                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0)); /* BCs on ext. faces only */

                is_outflow = c1 >= 0;

                t  = trans[ f ];
                s  = 2.0*is_outflow - 1.0;
                c1 = is_outflow ? c1 : c2;
                ix = csrmatrix_elm_index(c1, c1, h->A);

                h->A->sa[ ix ] += t;
                h->b    [ c1 ] += t * bc->value[ i ];
                h->b    [ c1 ] -= s * t * h->pimpl->fgrav[ f ];
            }
        }

        else if (bc->type[ i ] == BC_FLUX_TOTVOL) {
            /* We currently support individual flux faces only. */
            assert (bc->cond_pos[i + 1] - bc->cond_pos[i] == 1);

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[i + 1]; j++) {
                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0)); /* BCs on ext. faces only */

                c1 = (c1 >= 0) ? c1 : c2;

                /* Interpret BC as flow *INTO* cell */
                h->b[ c1 ] += bc->value[ i ];
            }
        }

        /* Other types currently not handled */
    }

    return is_neumann;
}


/* ---------------------------------------------------------------------- */
static void
boundary_fluxes(struct UnstructuredGrid             *G     ,
                const struct FlowBoundaryConditions *bc    ,
                const double                        *trans ,
                const double                        *cpress,
                const struct ifs_tpfa_data          *h     ,
                double                              *fflux )
/* ---------------------------------------------------------------------- */
{
    int    f, c1, c2;
    size_t i, j;
    double s, dh;

    for (i = 0; i < bc->nbc; i++) {
        if (bc->type[ i ] == BC_PRESSURE) {
            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[ i + 1 ]; j++) {

                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0));

                if (c1 < 0) {   /* Environment -> c2 */
                    dh = bc->value[ i ] - cpress[c2];
                }
                else {          /* c1 -> environment */
                    dh = cpress[c1] - bc->value[ i ];
                }

                fflux[f] = trans[f] * (dh + h->pimpl->fgrav[f]);
            }
        }

        else if (bc->type[ i ] == BC_FLUX_TOTVOL) {
            assert (bc->cond_pos[i+1] - bc->cond_pos[i] == 1);

            for (j = bc->cond_pos[ i ]; j < bc->cond_pos[ i + 1 ]; j++) {
                    
                f  = bc->face[ j ];
                c1 = G->face_cells[2*f + 0];
                c2 = G->face_cells[2*f + 1];

                assert ((c1 < 0) ^ (c2 < 0));

                /* BC flux is positive into reservoir. */
                s = 2.0*(c1 < 0) - 1.0;

                fflux[f] = s * bc->value[ i ];
            }
        }
    }
}


/* ======================================================================
 * Public interface below separator.
 * ====================================================================== */

/* ---------------------------------------------------------------------- */
struct ifs_tpfa_data *
ifs_tpfa_construct(struct UnstructuredGrid *G,
                   struct Wells            *W)
/* ---------------------------------------------------------------------- */
{
    struct ifs_tpfa_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->pimpl = impl_allocate(G, W);
        new->A     = ifs_tpfa_construct_matrix(G, W);

        if ((new->pimpl == NULL) || (new->A == NULL)) {
            ifs_tpfa_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        new->b = new->pimpl->ddata;
        new->x = new->b             + new->A->m;

        new->pimpl->fgrav = new->x  + new->A->m;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_assemble(struct UnstructuredGrid      *G     ,
                  const struct ifs_tpfa_forces *F     ,
                  const double                 *trans ,
                  const double                 *gpress,
                  struct ifs_tpfa_data         *h     )
/* ---------------------------------------------------------------------- */
{
    int c1, c2, c, i, f, j1, j2;

    int is_neumann;

    double s;

    csrmatrix_zero(         h->A);
    vector_zero   (h->A->m, h->b);

    compute_grav_term(G, gpress, h->pimpl->fgrav);

    for (c = i = 0; c < G->number_of_cells; c++) {
        j1 = csrmatrix_elm_index(c, c, h->A);

        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            c1 = G->face_cells[2*f + 0];
            c2 = G->face_cells[2*f + 1];

            s  = 2.0*(c1 == c) - 1.0;
            c2 = (c1 == c) ? c2 : c1;

            h->b[c] -= trans[f] * (s * h->pimpl->fgrav[f]);

            if (c2 >= 0) {
                j2 = csrmatrix_elm_index(c, c2, h->A);

                h->A->sa[j1] += trans[f];
                h->A->sa[j2] -= trans[f];
            }
        }
    }

    is_neumann = 1;
    if ((F != NULL) && (F->bc != NULL)) {
        is_neumann = assemble_bc_contrib(G, F->bc, trans, h);
    }

    if ((F != NULL) && (F->src != NULL)) {
        for (c = 0; c < G->number_of_cells; c++) {
            h->b[c] += F->src[c];
        }
    }

    if (is_neumann) {
        /* Remove zero eigenvalue associated to constant pressure */
        h->A->sa[0] *= 2;
    }
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_press_flux(struct UnstructuredGrid      *G    ,
                    const struct ifs_tpfa_forces *F    ,
                    const double                 *trans,
                    struct ifs_tpfa_data         *h    ,
                    struct ifs_tpfa_solution     *soln )
/* ---------------------------------------------------------------------- */
{
    int    c1, c2, f;
    double dh;

    double *cpress, *fflux;

    assert (soln             != NULL);
    assert (soln->cell_press != NULL);
    assert (soln->face_flux  != NULL);

    cpress = soln->cell_press;
    fflux  = soln->face_flux ;

    /* Assign cell pressure directly from solution vector */
    memcpy(cpress, h->x, G->number_of_cells * sizeof *cpress);

    for (f = 0; f < G->number_of_faces; f++) {
        c1 = G->face_cells[2*f + 0];
        c2 = G->face_cells[2*f + 1];

        if ((c1 >= 0) && (c2 >= 0)) {
            dh = cpress[c1] - cpress[c2] + h->pimpl->fgrav[f];
            fflux[f] = trans[f] * dh;
        } else {
            fflux[f] = 0.0;
        }
    }

    if (F != NULL) {
        if (F->bc != NULL) {
            boundary_fluxes(G, F->bc, trans, cpress, h, fflux);
        }
    }
}


/* ---------------------------------------------------------------------- */
void
ifs_tpfa_destroy(struct ifs_tpfa_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        csrmatrix_delete(h->A);
        impl_deallocate (h->pimpl);
    }

    free(h);
}
