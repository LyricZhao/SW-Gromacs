# include "gmxpre.h"
# include "nbnxn_kernel_gpu_ref.h"
# include "config.h"

# include <math.h>

# include "gromacs/legacyheaders/force.h"
# include "gromacs/legacyheaders/typedefs.h"
# include "gromacs/legacyheaders/types/simple.h"
# include "gromacs/math/utilities.h"
# include "gromacs/math/vec.h"
# include "gromacs/mdlib/nb_verlet.h"
# include "gromacs/mdlib/nbnxn_consts.h"
# include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
# include "gromacs/pbcutil/ishift.h"

# include "athread.h"
# include "spinlock.h"

# define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
# define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

# define vc_max 128
# define thread_max 64

extern SLAVE_FUN(sw_computing_core)();

tMPI_Spinlock_t bEner_lock[vc_max];
tMPI_Spinlock_t thread_lock[thread_max];
tMPI_Spinlock_t para_pass_lock;
int isEnd, my_looper;

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
    const nbnxn_sci_t  *nbln;
    const real         *x;
    gmx_bool            bEner;
    gmx_bool            bEwald;
    const real         *Ftab = NULL;
    real                rcut2, rvdw2, rlist2;
    int                 ntype;
    real                facel;
    int                 n;
    int                 ish3;
    int                 sci;
    int                 cj4_ind0, cj4_ind1, cj4_ind;
    int                 ci, cj;
    int                 ic, ia, im, jm;
    real                shX, shY, shZ;
    real                iq;
    real                vctot[vc_max];
    real                Vvdwtot[vc_max];
    const real         *shiftvec;
    real               *vdwparam;
    int                *shift;
    int                *type;
    const nbnxn_excl_t *excl[2];

    if (clearF == enbvClearFYes) clear_f(nbat, 0, f);

    bEner = (force_flags & GMX_FORCE_ENERGY);

    bEwald = EEL_FULL(iconst->eeltype);
    if (bEwald) Ftab = iconst->tabq_coul_F;

    rcut2               = iconst->rcoulomb*iconst->rcoulomb;
    rvdw2               = iconst->rvdw*iconst->rvdw;
    rlist2              = nbl->rlist*nbl->rlist;
    type                = nbat->type;
    facel               = iconst->epsfac;
    shiftvec            = shift_vec[0];
    vdwparam            = nbat->nbfp;
    ntype               = nbat->ntype;

    x = nbat->x;

    long paras_package[26];
    long paras_array[4] = {(long) thread_lock, (long) &para_pass_lock, (long) &isEnd, (long) paras_package};
    isEnd = 0;
    tMPI_Spinlock_init(&para_pass_lock);
    for(my_looper = 0; my_looper < thread_max; ++ my_looper) {
      tMPI_Spinlock_init(&thread_lock[my_looper]);
    }
    athread_spawn64(SLAVE_FUN(sw_computing_core), &paras_array);
    assert(nbl -> nsci < vc_max);

    for (n = 0; n < nbl->nsci; n++) { // nsci changes (1 to 64)

        tMPI_Spinlock_init(&bEner_lock[n]);

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

        if (nbln->shift == CENTRAL && nbl->cj4[cj4_ind0].cj[0] == sci*NCL_PER_SUPERCL) {
            for (im = 0; im < NCL_PER_SUPERCL; im++) {
                ci = sci*NCL_PER_SUPERCL + im;
                for (ic = 0; ic < CL_SIZE; ic++) {
                    ia = ci*CL_SIZE + ic;
                    iq = x[ia*nbat->xstride+3];
                    vctot[n] += iq*iq;
                }
            }
            if (!bEwald) vctot[n] *= -facel*0.5*iconst->c_rf;
            else vctot[n] *= -facel*iconst->ewaldcoeff_q*M_1_SQRTPI;
        }

        for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {

            excl[0] = &nbl->excl[nbl->cj4[cj4_ind].imei[0].excl_ind];
            excl[1] = &nbl->excl[nbl->cj4[cj4_ind].imei[1].excl_ind];

            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) { // NBNXN_GPU_JGROUP_SIZE = 4

                cj = nbl->cj4[cj4_ind].cj[jm];

                /* load task in the thread pool */
                while(true) {
                  for(my_looper = 0; my_looper < thread_max; ++ my_looper) {
                    if(tMPI_Spinlock_islocked(&thread_lock[my_looper]) == 0) {
                      /* found, pass and break */
                      tMPI_Spinlock_lock(&thread_lock[my_looper]);

                      paras_package[0] = (long) nbl -> cj4[cj4_ind].imei[0].imask;
                      paras_package[1] = (long) jm; paras_package[2] = (long) sci; paras_package[3] = shX; paras_package[4] = shY; paras_package[5] = shZ;
                      paras_package[6] = (long) ntype; paras_package[7] = (long) cj; paras_package[8] = (long) (nbat -> xstride); paras_package[9] = (long) (nbat -> fstride);
                      paras_package[10] = (long) facel; paras_package[11] = (long) rcut2; paras_package[12] = (long) rvdw2;
                      paras_package[13] = (long) type;
                      paras_package[14] = (long) excl;
                      paras_package[15] = (long) bEwald; paras_package[16] = (long) bEner;
                      paras_package[17] = (long) iconst;
                      paras_package[18] = (long) Ftab; paras_package[19] = (long) x;
                      paras_package[20] = (long) vdwparam; paras_package[21] = (long) f; paras_package[22] = (long) fshift;
                      paras_package[23] = (long) &vctot[n]; paras_package[24] = (long) &Vvdwtot[n];
                      paras_package[25] = (long) &bEner_lock[n];
                      paras_package[26] = (long) nbln -> shift;
                      tMPI_Spinlock_lock(&para_pass_lock);

                      while(tMPI_Spinlock_islocked(&para_pass_lock));
                      goto wd_done;
                    }
                  }
                }
                wd_done: /* nothing, next task */
            }
        }
    }

    isEnd = 1;
    athread_join64();

    if(bEner) {
      for (n = 0; n < nbl->nsci; n++) {
        Vc[0]         = Vc[0]   + vctot[n];
        Vvdw[0]       = Vvdw[0] + Vvdwtot[n];
      }
    }
}
