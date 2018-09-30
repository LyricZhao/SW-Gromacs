/* SW Slave Computing Kernel */

# include "gmxpre.h"
# include "nbnxn_kernel_gpu_ref.h"
# include "config.h"

# include <math.h>

# include "gromacs/legacyheaders/force.h"
# include "gromacs/legacyheaders/typedefs.h"
# include "gromacs/legacyheaders/types/simple.h"
# include "gromacs/mdlib/nb_verlet.h"
# include "gromacs/mdlib/nbnxn_consts.h"
# include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
# include "gromacs/pbcutil/ishift.h"

# include "slave.h"
# include "unifunc.h"

# define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
# define BLOCK_HIH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
# define BLOCK_SIZE(id, p, n) (BLOCK_HIH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
# define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

# define CL_SIZE 8
# define NCL_PER_SUPERCL 8
# define THREAD_TOT 64

/* Variables */
__thread_kernel("compute") int nsci, threadID, threadST, threadED, force_flags;
__thread_kernel("compute") const nbnxn_pairlist_t *nbl;
__thread_kernel("compute") const nbnxn_atomdata_t *nbat;
__thread_kernel("compute") const interaction_const_t *iconst;
__thread_kernel("compute") rvec *shift_vec;
__thread_kernel("compute") real *f, *fshift, *vctot, *Vvdwtot;


static inline void setParas(long *paras) {
  nbl = (const nbnxn_pairlist_t *) paras[0]; nbat = (const nbnxn_atomdata_t *) paras[1]; iconst = (const interaction_const_t *) paras[2];
  shift_vec = (rvec *) paras[3]; force_flags = (int) paras[4]; f = (real *) paras[5];
  fshift = (real *) paras[6]; vctot = (real *) paras[7]; Vvdwtot = (real *) paras[8];

  nsci = nbl -> nsci;
  threadID = athread_get_id(-1);
  threadST = BLOCK_LOW(threadID, THREAD_TOT, nsci);
  threadED = BLOCK_HIH(threadID, THREAD_TOT, nsci);
}

void sw_computing_core(long *paras) {
  setParas(paras);

  /* Definition */
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
  real                qq, vcoul = 0, krsq;
  int                 nti;
  int                 tj;
  real                rt, r, eps;
  real                rinvsix;
  real                Vvdw_rep, Vvdw_disp;
  real                ix, iy, iz, fix, fiy, fiz;
  real                jx, jy, jz;
  real                dx, dy, dz, rsq, rinv;
  int                 int_bit;
  real                fexcl;
  real                c6, c12, cexp1, cexp2, br;
  const real         *shiftvec;
  real               *vdwparam;
  int                *shift;
  int                *type;
  const nbnxn_excl_t *excl[2];

  /* Variables */
  bEner = (force_flags & GMX_FORCE_ENERGY);
  bEwald = EEL_FULL(iconst->eeltype);
  if (bEwald) Ftab    = iconst->tabq_coul_F;
  rcut2               = iconst->rcoulomb*iconst->rcoulomb;
  rvdw2               = iconst->rvdw*iconst->rvdw;
  rlist2              = nbl->rlist*nbl->rlist;
  type                = nbat->type;
  facel               = iconst->epsfac;
  shiftvec            = shift_vec[0];
  vdwparam            = nbat->nbfp;
  ntype               = nbat->ntype;
  x = nbat->x;

  /* Kernel */
  for(n = threadST; n <= threadED; ++ n) {

    nbln             = &nbl->sci[n];
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

            for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
                if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

                    ci               = sci * NCL_PER_SUPERCL + im;

                    for (ic = 0; ic < CL_SIZE; ic++) { // CL_SIZE = 8

                        ia               = ci*CL_SIZE + ic;
                        is               = ia*nbat->xstride;
                        ifs              = ia*nbat->fstride;
                        ix               = shX + x[is+0];
                        iy               = shY + x[is+1];
                        iz               = shZ + x[is+2];
                        iq               = facel*x[is+3];
                        nti              = ntype*2*type[ia];
                        fix              = 0;
                        fiy              = 0;
                        fiz              = 0;

                        for (jc = 0; jc < CL_SIZE; jc++) { // CL_SIZE = 8
                            /* core computing */

                            ja               = cj*CL_SIZE + jc;
                            if (nbln->shift == CENTRAL && ci == cj && ja <= ia) continue;
                            int_bit = ((excl[jc>>2]->pair[(jc & 3)*CL_SIZE+ic] >> (jm*NCL_PER_SUPERCL+im)) & 1);

                            js               = ja*nbat->xstride;
                            jfs              = ja*nbat->fstride;
                            jx               = x[js+0];
                            jy               = x[js+1];
                            jz               = x[js+2];
                            dx               = ix - jx;
                            dy               = iy - jy;
                            dz               = iz - jz;
                            rsq              = dx*dx + dy*dy + dz*dz;
                            if (rsq >= rcut2) continue;

                            // avoid NaN for excluded pairs at r=0
                            rsq             += (1.0 - int_bit)*NBNXN_AVOID_SING_R2_INC;

                            rinv             = gmx_software_invsqrt(rsq);
                            rinvsq           = rinv * rinv;
                            fscal            = 0;

                            qq               = iq*x[js+3];
                            if (!bEwald) {
                                // Reaction-field
                                krsq  = iconst->k_rf*rsq;
                                fscal = qq*(int_bit*rinv - 2*krsq)*rinvsq;
                                if (bEner) vcoul = qq*(int_bit*rinv + krsq - iconst->c_rf);
                            }
                            else {
                                r     = rsq*rinv;
                                rt    = r*iconst->tabq_scale;
                                n0    = rt;
                                eps   = rt - n0;
                                fexcl = (1 - eps)*Ftab[n0] + eps*Ftab[n0+1];
                                fscal = qq*(int_bit*rinvsq - fexcl)*rinv;
                                if (bEner) vcoul = qq*((int_bit - gmx_erff(iconst->ewaldcoeff_q*r))*rinv - int_bit*iconst->sh_ewald);
                            }
                            if (rsq < rvdw2) {
                                tj        = nti + 2*type[ja];

                                // Vanilla Lennard-Jones cutoff
                                c6        = vdwparam[tj];
                                c12       = vdwparam[tj+1];

                                rinvsix   = int_bit*rinvsq*rinvsq*rinvsq;
                                Vvdw_disp = c6*rinvsix;
                                Vvdw_rep  = c12*rinvsix*rinvsix;
                                fscal    += (Vvdw_rep - Vvdw_disp)*rinvsq;

                                if (bEner) {
                                    vctot[n]   += vcoul;
                                    Vvdwtot[n] +=
                                        (Vvdw_rep - int_bit*c12*iconst->sh_invrc6*iconst->sh_invrc6)/12 -
                                        (Vvdw_disp - int_bit*c6*iconst->sh_invrc6)/6;
                                }
                            }

                            tx        = fscal*dx;
                            ty        = fscal*dy;
                            tz        = fscal*dz;
                            fix       = fix + tx;
                            fiy       = fiy + ty;
                            fiz       = fiz + tz;
                            f[jfs+0] -= tx;
                            f[jfs+1] -= ty;
                            f[jfs+2] -= tz;
                        }

                        f[ifs+0]        += fix;
                        f[ifs+1]        += fiy;
                        f[ifs+2]        += fiz;
                        fshift[ish3]     = fshift[ish3]   + fix;
                        fshift[ish3+1]   = fshift[ish3+1] + fiy;
                        fshift[ish3+2]   = fshift[ish3+2] + fiz;
                    }
                }
            }
        } /* JGroup_loop End */
    } /* cj4_loop End */
    /* Reduce vctot/vvdwtot in the main core */
  } /* End Main Loop */
}
