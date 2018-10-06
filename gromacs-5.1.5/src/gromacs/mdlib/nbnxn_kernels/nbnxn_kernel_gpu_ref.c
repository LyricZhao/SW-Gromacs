# include "gmxpre.h"
# include "nbnxn_kernel_gpu_ref.h"
# include "config.h"

# include <math.h>
# include <assert.h>

# include "gromacs/legacyheaders/force.h"
# include "gromacs/legacyheaders/typedefs.h"
# include "gromacs/legacyheaders/types/simple.h"
# include "gromacs/math/utilities.h"
# include "gromacs/math/vec.h"
# include "gromacs/mdlib/nb_verlet.h"
# include "gromacs/mdlib/nbnxn_consts.h"
# include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
# include "gromacs/pbcutil/ishift.h"

# include <athread.h>

# define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
# define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

# define LYRIC_DEBUG

extern void slave_sw_computing_core();

static inline void addTrans(void* &dest, void* src, int siz) {
  memcpy(dest, src, siz);
  dest += siz;
}

void
nbnxn_kernel_gpu_ref(const nbnxn_pairlist_t     *nbl,
                     const nbnxn_atomdata_t     *nbat,
                     const interaction_const_t  *iconst,
                     rvec                       *shift_vec,
                     int                         force_flags,
                     int                         clearF,
                     real  *                     f,
                     real  *                     fshift,
                     real  *                     Vc,
                     real  *                     Vvdw) {

    /* Variables Init */
    if (clearF == enbvClearFYes) clear_f(nbat, 0, f);

    /* Global */
    gmx_bool bEner = (force_flags & GMX_FORCE_ENERGY);
    real facel = iconst->epsfac;
    const real *x = nbat->x;
    const nbnxn_excl_t *excl[2];

    /* Attributes */
    const nbnxn_sci_t *nbln;
    int n, ish3, cj4_ind0, cj4_ind1, sci;
    real shX, shY, shZ;

    /* Temp Vars */
    int im, ci, ic, ia, jm, im, cj;
    int is, ifs, js, jfs;
    real iq;

    void **startPoint = (void **) malloc(nbl->nsci * sizeof(void *));
    void *transferData = malloc(56 * 1024);
    void *transferDataGlobal = malloc();
    void *tempPtr = transferData;
    real *vctot = (real *) malloc(sizeof(real) * nbl->nsci);
    real *Vvdwtot = (real *) malloc(sizeof(real) * nbl->nsci);
    real *shiftvec = shift_vec[0];

    for (n = 0; n < nbl->nsci; n++) {

        startPoint[n] = tempPtr;
        nbln = &nbl->sci[n];
        ish3             = 3*nbln->shift;
        shX              = shiftvec[ish3];
        shY              = shiftvec[ish3+1];
        shZ              = shiftvec[ish3+2];
        cj4_ind0         = nbln->cj4_ind_start;
        cj4_ind1         = nbln->cj4_ind_end;
        sci              = nbln->sci;
        vctot[n]         = 0;
        Vvdwtot[n]       = 0;
        addTrans(tempPtr, (void *) nbln, sizeof(nbnxn_sci_t));
        addTrans(tempPtr, (void *) &shiftvec[ish3], sizeof(real) * 3);

        /* MPE Works */
        if (nbln->shift == CENTRAL && nbl->cj4[cj4_ind0].cj[0] == sci*NCL_PER_SUPERCL) {
            for (im = 0; im < NCL_PER_SUPERCL; im++) {
                ci = sci*NCL_PER_SUPERCL + im;
                for (ic = 0; ic < CL_SIZE; ic++) {
                    ia = ci*CL_SIZE + ic;
                    iq = x[ia*nbat->xstride+3];
                    vctot[n] += iq*iq;
                }
            }
            vctot[n] *= -facel*iconst->ewaldcoeff_q*M_1_SQRTPI;
        }

        /* Data Package */
        for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {

            excl[0] = &nbl->excl[nbl->cj4[cj4_ind].imei[0].excl_ind];
            excl[1] = &nbl->excl[nbl->cj4[cj4_ind].imei[1].excl_ind];

            addTrans(tempPtr, (void *) excl[0]->pair, 32 * sizeof(unsigned int));
            addTrans(tempPtr, (void *) excl[1]->pair, 32 * sizeof(unsigned int));
            addTrans(tempPtr, (void *) &nbl->cj4[cj4_ind].imei[0].imask, sizeof(unsigned int));

            /* kernel starts here */
            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) { // NBNXN_GPU_JGROUP_SIZE = 4

                cj = nbl->cj4[cj4_ind].cj[jm];
                addTrans(tempPtr, (void *) &cj, sizeof(int));

                for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
                    if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

                        ci = sci * NCL_PER_SUPERCL + im;

                        /* ic: 0 to 8 */
                        addTrans(tempPtr, (void *) &x[(ci * CL_SIZE) * nbat->xstride], 8 * CL_SIZE * nbat->xstride * sizeof(real));
                        addTrans(tempPtr, (void *) &type[ci * CL_SIZE], 8 * sizeof(int));
                        addTrans(tempPtr, (void *) &f[(ci * CL_SIZE) * nbat->fstride], 8 * CL_SIZE * nbat->fstride * sizeof(real));

                        /* jc: 0 to 8 */
                        addTrans(tempPtr, (void *) &x[(cj * CL_SIZE) * nbat->xstride], 8 * CL_SIZE * nbat->xstride * sizeof(real));
                        addTrans(tempPtr, (void *) &type[cj * CL_SIZE], 8 * sizeof(int));
                        addTrans(tempPtr, (void *) &f[(cj * CL_SIZE)] * nbat->fstride, 8 * CL_SIZE * nbat->fstride * sizeof(real));
                    }
                }
            }
        }
    }

    /* Calling Kernel */
    long paras[3] = {(long) nbl->nsci, (long) transferData, (long) transferDataGlobal};
    athread_spawn(sw_computing_core, &paras);
    athread_join();
    while(1);

# ifndef LYRIC_DEBUG
    gmx_bool bEner = (force_flags & GMX_FORCE_ENERGY);
    /* Reduce */
    if(bEner) {
      int n;
      for (n = 0; n < nbl->nsci; n++) {
        Vc[0]         = Vc[0]   + vctot[n];
        Vvdw[0]       = Vvdw[0] + Vvdwtot[n];
      }
    }
# endif
    free(transferData);
}
