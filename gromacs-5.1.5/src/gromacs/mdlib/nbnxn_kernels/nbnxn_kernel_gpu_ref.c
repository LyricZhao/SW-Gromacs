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

void *tempPtrG;
int ttsize = 0;

static inline void setP(void *dest, void *pa) {
  memcpy(dest, pa, 4);
}

static inline void addTrans(void* useless, void* src, int siz) {
  ttsize += siz;
  assert(ttsize < 8 * 1024 * 1024);
  memcpy(tempPtrG, src, siz);
  tempPtrG += siz;
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
    const real *x = nbat->x;
    const nbnxn_excl_t *excl[2];

    /* Attributes */
    const nbnxn_sci_t *nbln;
    int n, ish3, cj4_ind, cj4_ind0, cj4_ind1, sci;
    real shX, shY, shZ;
    real rcut2 = iconst->rcoulomb*iconst->rcoulomb;
    real rvdw2 = iconst->rvdw*iconst->rvdw;
    int *type = nbat->type;

    /* Temp Vars */
    int im, ci, ic, ia, jm, cj;
    int is, ifs, js, jfs;
    real iq, facel = iconst->epsfac;

    /* Loop Buffer Data */
    void *transferData = malloc(8 * 1024 * 1024); tempPtrG = transferData; /* 8MB */

    /* Global Data */
    int vdwparam_size = nbat->ntype * nbat->ntype * 2 * sizeof(real);
    real *vdwparam = (real *) malloc(vdwparam_size);
    memcpy((void *) vdwparam, (void *)nbat->nbfp, vdwparam_size);
    int Ftab_size = iconst->tabq_size * sizeof(real);
    real *Ftab = (real *) malloc(Ftab_size);
    memcpy((void *) Ftab, (void *) iconst->tabq_coul_F, Ftab_size);

    real *vctot = (real *) malloc(sizeof(real) * nbl->nsci);
    real *Vvdwtot = (real *) malloc(sizeof(real) * nbl->nsci);
    real *shiftvec = shift_vec[0];

    for (n = 0; n < nbl->nsci; n++) {
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
        addTrans(tempPtrG, (void *) nbln, sizeof(nbnxn_sci_t));
        addTrans(tempPtrG, (void *) &shiftvec[ish3], sizeof(real) * 3);

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

            addTrans(tempPtrG, (void *) excl[0]->pair, 32 * sizeof(unsigned int));
            addTrans(tempPtrG, (void *) excl[1]->pair, 32 * sizeof(unsigned int));
            addTrans(tempPtrG, (void *) &(nbl->cj4[cj4_ind].imei[0].imask), sizeof(unsigned int));

            /* kernel starts here */
            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) { // NBNXN_GPU_JGROUP_SIZE = 4

                cj = nbl->cj4[cj4_ind].cj[jm];
                addTrans(tempPtrG, (void *) &cj, sizeof(int));

                continue;
                for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
                    if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

                        ci = sci * NCL_PER_SUPERCL + im;

                        /* ic: 0 to 8 */
                        addTrans(tempPtrG, (void *) &x[(ci * CL_SIZE) * nbat->xstride], (7 * CL_SIZE * nbat->xstride + 4) * sizeof(real));
                        addTrans(tempPtrG, (void *) &type[ci * CL_SIZE], 8 * sizeof(int));
                        addTrans(tempPtrG, (void *) &f[(ci * CL_SIZE) * nbat->fstride], (7 * CL_SIZE * nbat->fstride + 3) * sizeof(real));

                        /* jc: 0 to 8 */
                        addTrans(tempPtrG, (void *) &x[(cj * CL_SIZE) * nbat->xstride], (7 * CL_SIZE * nbat->xstride + 4) * sizeof(real));
                        addTrans(tempPtrG, (void *) &type[cj * CL_SIZE], 8 * sizeof(int));
                        addTrans(tempPtrG, (void *) &f[(cj * CL_SIZE) * nbat->fstride], (7 * CL_SIZE * nbat->fstride + 3 )* sizeof(real));
                    }
                }
            }
        }
    }

    long long_nsci = nbl->nsci;
    long long_nxstride = nbat->xstride;
    long long_nfstride = nbat->fstride;
    long long_ntype = nbat->ntype;

    /* Calling Kernel */
    char paras[18 * 8];
    setP((void *) &paras[ 0 * 8], (void *) &(long_nsci));
    setP((void *) &paras[ 1 * 8], (void *) &(transferData));
    setP((void *) &paras[ 2 * 8], (void *) &(long_nxstride));
    setP((void *) &paras[ 3 * 8], (void *) &(long_nfstride));
    setP((void *) &paras[ 4 * 8], (void *) &(long_ntype));
    setP((void *) &paras[ 5 * 8], (void *) &(rcut2));
    setP((void *) &paras[ 6 * 8], (void *) &(rvdw2));
    setP((void *) &paras[ 7 * 8], (void *) &(vctot));
    setP((void *) &paras[ 8 * 8], (void *) &(Vvdwtot));
    setP((void *) &paras[ 9 * 8], (void *) &(vdwparam));
    setP((void *) &paras[10 * 8], (void *) &(iconst->tabq_coul_F));
    setP((void *) &paras[11 * 8], (void *) &(iconst->k_rf));
    setP((void *) &paras[12 * 8], (void *) &(iconst->c_rf));
    setP((void *) &paras[13 * 8], (void *) &(iconst->tabq_scale));
    setP((void *) &paras[14 * 8], (void *) &(iconst->ewaldcoeff_q));
    setP((void *) &paras[15 * 8], (void *) &(iconst->sh_ewald));
    setP((void *) &paras[16 * 8], (void *) &(iconst->sh_invrc6));
    setP((void *) &paras[17 * 8], (void *) &(iconst->epsfac));


    athread_spawn(sw_computing_core, &paras);
    athread_join();
    while(1);

/*
# ifndef LYRIC_DEBUG
    gmx_bool bEner = (force_flags & GMX_FORCE_ENERGY);
    if(bEner) {
      int n;
      for (n = 0; n < nbl->nsci; n++) {
        Vc[0]         = Vc[0]   + vctot[n];
        Vvdw[0]       = Vvdw[0] + Vvdwtot[n];
      }
    }
# endif
*/
    free(transferData);
    free(vdwparam);
    free(Ftab);
    free(vctot);
    free(Vvdwtot);
}
