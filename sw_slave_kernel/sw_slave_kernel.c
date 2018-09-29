/* SW Slave Computing Kernel */

# include <math.h>

# include "gromacs/legacyheaders/force.h"
# include "gromacs/legacyheaders/typedefs.h"
# include "gromacs/legacyheaders/types/simple.h"
# include "gromacs/math/utilities.h"
# include "gromacs/math/vec.h"
# include "gromacs/mdlib/nbnxn_consts.h"
# include "gromacs/pbcutil/ishift.h"

# include "slave.h"
# include "spinlock.h"

# define CL_SIZE 8
# define NCL_PER_SUPERCL 8

/* Thread */
__thread_kernel("compute") int threadID;
__thread_kernel("compute") tMPI_Spinlock_t *para_pass_lock, *thread_lock;
__thread_kernel("compute") int* isEnd;

/* Parameters */
__thread_kernel("compute") long *paras_package;
__thread_kernel("compute") unsigned int imask;
__thread_kernel("compute") int jm, sci, shX, shY, shZ, ntype, cj, xstride, fstride;
__thread_kernel("compute") real facel, rcut2, rvdw2;
__thread_kernel("compute") int *type;
__thread_kernel("compute") const unsigned int **excl;
__thread_kernel("compute") gmx_bool bEwald, bEner;
__thread_kernel("compute") const interaction_const_t *iconst;
__thread_kernel("compute") const real *Ftab, *x;
__thread_kernel("compute") real *vdwparam, *f, *fshift, *vctot, *Vvdwtot;
__thread_kernel("compute") tMPI_Spinlock_t *bEner_lock;
__thread_kernel("compute") int nbln_shift;

static inline void setParas(long *paras) {
  threadID = athread_get_id(-1);
  thread_lock = &(((tMPI_Spinlock_t *) paras[0])[threadID]);
  para_pass_lock = (tMPI_Spinlock_t *) paras[1];
  isEnd = (int *) paras[2]; paras_package = (long *) paras[3];
}

static inline void transferParas() {
  imask = (unsigned int) paras_package[0];
  jm = (int) paras_package[1]; sci = (int) paras_package[2]; shX = (int) paras_package[3]; shY = (int) paras_package[4]; shZ = (int) paras_package[5];
  ntype = (int) paras_package[6]; cj = (int) paras_package[7]; xstride = (int) paras_package[8]; fstride = (int) paras_package[9];
  facel = (real) paras_package[10]; rcut2 = (real) paras_package[11]; rvdw2 = (real) paras_package[12];
  type = (int *) paras_package[13];
  excl = (const unsigned int **) paras_package[14];
  bEwald = (gmx_bool) paras_package[15]; bEner = (gmx_bool) paras_package[16];
  iconst = (const interaction_const_t *) paras_package[17];
  Ftab = (const real *) paras_package[18]; x = (const real *) paras_package[19];
  vdwparam = (real *) paras_package[20]; f = (real *) paras_package[21]; fshift = (real *) paras_package[22];
  vctot = (real *) paras_package[23]; Vvdwtot = (real *) paras_package[24];
  bEner_lock = (tMPI_Spinlock_t *) paras_package[25];
  nbln_shift = (int) paras_package[26];
}

static void real_computing_core() {
  /* Stack Variables */
  int ci, im, ic, ia, is, ifs, ix, iy, iz, nti, jc, ja, int_bit, js, jfs, n0, tj, ish3;
  real iq, fix, fiy, fiz, jx, jy, jz, dx, dy, dz, rsq, rinv, rinvsq, fscal, qq;
  real krsq, vcoul = 0, r, rt, eps, fexcl, c6, c12, rinvsix, Vvdw_disp, Vvdw_rep, tx, ty, tz;
  for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
      if ((imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

          ci               = sci * NCL_PER_SUPERCL + im;

          for (ic = 0; ic < CL_SIZE; ic++) { // CL_SIZE = 8

              ia               = ci*CL_SIZE + ic;
              is               = ia*xstride;
              ifs              = ia*fstride;
              ix               = shX + x[is+0];
              iy               = shY + x[is+1];
              iz               = shZ + x[is+2];
              iq               = facel*x[is+3];
              nti              = ntype*2*type[ia];
              fix              = 0;
              fiy              = 0;
              fiz              = 0;

              for (jc = 0; jc < CL_SIZE; jc++) { // CL_SIZE = 8

                  ja               = cj*CL_SIZE + jc;
                  if (nbln_shift == CENTRAL && ci == cj && ja <= ia) continue;
                  int_bit = ((excl[jc>>2][(jc & 3)*CL_SIZE+ic] >> (jm*NCL_PER_SUPERCL+im)) & 1);

                  js               = ja*xstride;
                  jfs              = ja*fstride;
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

                  rinv             = gmx_invsqrt(rsq);
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
                      if (bEner) vcoul = qq*((int_bit - gmx_erf(iconst->ewaldcoeff_q*r))*rinv - int_bit*iconst->sh_ewald);
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
                          tMPI_Spinlock_lock(bEner_lock);
                          (*vctot)   += vcoul;
                          (*Vvdwtot) +=
                              (Vvdw_rep - int_bit*c12*iconst->sh_invrc6*iconst->sh_invrc6)/12 -
                              (Vvdw_disp - int_bit*c6*iconst->sh_invrc6)/6;
                          tMPI_Spinlock_unlock(bEner_lock);
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
}

void sw_computing_core(long *paras) {
  setParas(paras);
  while(1) {
    while(1) {
      if(*isEnd) return; /* Finish */
      if(tMPI_Spinlock_islocked(thread_lock) == 1) {

        /* Pass Parameters */
        while(tMPI_Spinlock_islocked(para_pass_lock) == 0);
        tMPI_Spinlock_unlock(para_pass_lock);
        transferParas();

        /* Computing Kernel */
        real_computing_core();

        /* Done */
        break;
      }
    }
  }
}
