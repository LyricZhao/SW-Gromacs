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

# define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
# define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

extern SLAVE_FUN(sw_computing_core)();

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
    int                 ic, jc, ia, ja, is, ifs, js, jfs, im, jm;
    int                 n0;
    real                shX, shY, shZ;
    real                fscal, tx, ty, tz;
    real                rinvsq;
    real                iq;
    real                qq, vcoul = 0, krsq, vctot;
    int                 nti;
    int                 tj;
    real                rt, r, eps;
    real                rinvsix;
    real                Vvdwtot;
    real                Vvdw_rep, Vvdw_disp;
    real                ix, iy, iz, fix, fiy, fiz;
    real                jx, jy, jz;
    real                dx, dy, dz, rsq, rinv;
    int                 int_bit;
    real                fexcl;
    real                c6, c12, cexp1, cexp2, br;
    const real       *  shiftvec;
    real       *        vdwparam;
    int       *         shift;
    int       *         type;
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

    for (n = 0; n < nbl->nsci; n++) { // nsci changes (1 to 64)

        nbln = &nbl->sci[n];
        ish3             = 3*nbln->shift;
        shX              = shiftvec[ish3];
        shY              = shiftvec[ish3+1];
        shZ              = shiftvec[ish3+2];
        cj4_ind0         = nbln->cj4_ind_start;
        cj4_ind1         = nbln->cj4_ind_end;
        sci              = nbln->sci;
        vctot            = 0;
        Vvdwtot          = 0;

        if (nbln->shift == CENTRAL && nbl->cj4[cj4_ind0].cj[0] == sci*NCL_PER_SUPERCL) {
            for (im = 0; im < NCL_PER_SUPERCL; im++) {
                ci = sci*NCL_PER_SUPERCL + im;
                for (ic = 0; ic < CL_SIZE; ic++) {
                    ia = ci*CL_SIZE + ic;
                    iq = x[ia*nbat->xstride+3];
                    vctot += iq*iq;
                }
            }
            if (!bEwald) vctot *= -facel*0.5*iconst->c_rf;
            else vctot *= -facel*iconst->ewaldcoeff_q*M_1_SQRTPI;
        }

        for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {

            excl[0] = &nbl->excl[nbl->cj4[cj4_ind].imei[0].excl_ind];
            excl[1] = &nbl->excl[nbl->cj4[cj4_ind].imei[1].excl_ind];

            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) { // NBNXN_GPU_JGROUP_SIZE = 4

                cj = nbl->cj4[cj4_ind].cj[jm];

                /* using sw thread pool here */
                
                /* sw slave function call */
            }
        }

        if (bEner) {
            Vc[0]         = Vc[0]   + vctot;
            Vvdw[0]       = Vvdw[0] + Vvdwtot;
        }
    }
}
