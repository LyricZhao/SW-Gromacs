/* SW Slave Computing Kernel */

# include "gmxpre.h"
# include "nbnxn_kernel_gpu_ref.h"
# include "config.h"

# include <math.h>
# include <stdio.h>
# include <assert.h>

# include "gromacs/legacyheaders/force.h"
# include "gromacs/legacyheaders/typedefs.h"
# include "gromacs/legacyheaders/types/simple.h"
# include "gromacs/mdlib/nb_verlet.h"
# include "gromacs/mdlib/nbnxn_consts.h"
# include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
# include "gromacs/pbcutil/ishift.h"

# include <slave.h>
# include "unifunc.h"

# define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
# define BLOCK_HIH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
# define BLOCK_SIZE(id, p, n) (BLOCK_HIH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
# define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

# define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
# define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

# define THREAD_TOT 64

__thread_kernel("compute") real vdwparam[2] __attribute__((aligned(32)));
__thread_kernel("compute") real vctot[384], Vvdwtot[384] __attribute__((aligned(32)));
__thread_kernel("compute") real Ftab[2] __attribute__((aligned(32)));
__thread_kernel("compute") real xii[256], fii[256] __attribute__((aligned(32)));
__thread_kernel("compute") real xjj[256], fjj[256] __attribute__((aligned(32)));
__thread_kernel("compute") int typeii[8], typejj[8] __attribute__((aligned(32)));
__thread_kernel("compute") void *tempPtr __attribute__((aligned(32)));

static void cpe_printf(const char *fmt, ...){
  volatile long vprintf_addr = (long)vprintf;
  int (*vprintf_ptr)(const char *, va_list) = (void*)vprintf_addr;
  va_list vlist;
  va_start(vlist, fmt);
  vprintf_ptr(fmt, vlist);
  va_end(vlist);
}

static inline void getTrans(void *useless, void *dest_cpe, int siz) {
  volatile int getReply = 0;
  athread_get(PE_MODE, tempPtr, dest_cpe, siz, (void*) &getReply, 0, 0, 0);
  while(getReply != 1);
  tempPtr += siz;
}

static inline void getAll(void *src_mpe, void *dest_cpe, int siz) {
  volatile int getReply = 0;
  athread_get(PE_MODE, src_mpe, dest_cpe, siz, (void*) &getReply, 0, 0, 0);
  while(getReply != 1);
  tempPtr += siz;
}

static inline void setP(void *src, void *dest) {
  memcpy(dest, src, 8);
}

void sw_computing_core(char *paras) {
  int nsci = ((long *)paras)[0];
  int threadID = athread_get_id(-1);
  int threadST = BLOCK_LOW(threadID, THREAD_TOT, nsci);
  int threadED = BLOCK_HIH(threadID, THREAD_TOT, nsci) + 1;
  if(threadST >= threadED) return;

  /* SetParas */
  void **startPoint_addr;
  int *LoopBufferSize_addr;
  int nbat_xstride, nbat_fstride, ntype;
  long long_nbat_xstride, long_nbat_fstride, long_ntype;
  real rcut2, rvdw2;
  real *vctot_addr, *Vvdwtot_addr, *vdwparam_addr, *Ftab_addr;
  int vdwparam_size, Ftab_size;
  long long_vdwparam_size, long_Ftab_size;
  real k_rf, c_rf, tabq_scale, ewaldcoeff_q, sh_ewald, sh_invrc6, facel;

  setP((void *) &paras[ 1 * 8], (void *) &(tempPtr));
  setP((void *) &paras[ 2 * 8], (void *) &(long_nbat_xstride));
  setP((void *) &paras[ 3 * 8], (void *) &(long_nbat_fstride));
  setP((void *) &paras[ 4 * 8], (void *) &(long_ntype));
  setP((void *) &paras[ 5 * 8], (void *) &(rcut2));
  setP((void *) &paras[ 6 * 8], (void *) &(rvdw2));
  setP((void *) &paras[ 7 * 8], (void *) &(vctot_addr));
  setP((void *) &paras[ 8 * 8], (void *) &(Vvdwtot_addr));
  setP((void *) &paras[ 9 * 8], (void *) &(vdwparam_addr));
  setP((void *) &paras[10 * 8], (void *) &(Ftab_addr));
  setP((void *) &paras[11 * 8], (void *) &(k_rf));
  setP((void *) &paras[12 * 8], (void *) &(c_rf));
  setP((void *) &paras[13 * 8], (void *) &(tabq_scale));
  setP((void *) &paras[14 * 8], (void *) &(ewaldcoeff_q));
  setP((void *) &paras[15 * 8], (void *) &(sh_ewald));
  setP((void *) &paras[16 * 8], (void *) &(sh_invrc6));
  setP((void *) &paras[17 * 8], (void *) &(facel));
  nbat_xstride = long_nbat_xstride; nbat_fstride = long_nbat_fstride; ntype = long_ntype;

  /* Temp Vars */
  int n, sci, cj, im, ic, ia, ish3, jm, is, ifs, ja, jfs, jc, int_bit, n0, nti, js, ci, tj;
  real shX, shY, shZ, ix, iy, iz, fix, fiy, fiz, jx, jy, jz, dx, dy, dz, rsq, rinv, eps, rt, r, qq, c6, c12, cexp1, cexp2, rinvsix, Vvdw_rep, Vvdw_disp, iq, rinvsq, fscal, fexcl, tx, ty, tz;

  /* Attributes */
  real vcoul = 0, fshift[3];
  int cj4_ind, cj4_ind0, cj4_ind1;
  nbnxn_sci_t nbln;
  unsigned int excl_pair[2][32], imask;

  for(n = threadST; n < threadED; ++ n) {

    getTrans(tempPtr, &nbln, sizeof(nbnxn_sci_t));
    cj4_ind0         = nbln.cj4_ind_start;
    cj4_ind1         = nbln.cj4_ind_end;
    sci              = nbln.sci;
    ish3             = nbln.shift * 3;
    getTrans(tempPtr, (void *) &shX, sizeof(real));
    getTrans(tempPtr, (void *) &shY, sizeof(real));
    getTrans(tempPtr, (void *) &shZ, sizeof(real));
    fshift[0] = fshift[1] = fshift[2] = 0;


    for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {

      getTrans(tempPtr, (void *) excl_pair[0], 32 * sizeof(unsigned int));
      getTrans(tempPtr, (void *) excl_pair[1], 32 * sizeof(unsigned int));
      getTrans(tempPtr, (void *) &imask, sizeof(unsigned int));

      for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) {

        getTrans(tempPtr, (void *) &cj, sizeof(int));

        for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
          if ((imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

            /* 8KB Here */
            getTrans(tempPtr, (void *) xii, 8 * CL_SIZE * nbat_xstride * sizeof(real));
            getTrans(tempPtr, (void *) typeii, 8 * sizeof(int));
            getTrans(tempPtr, (void *) fii, 8 * CL_SIZE * nbat_fstride * sizeof(real));

            getTrans(tempPtr, (void *) xjj, 8 * CL_SIZE * nbat_xstride * sizeof(real));
            getTrans(tempPtr, (void *) typejj, 8 * sizeof(int));
            getTrans(tempPtr, (void *) fjj, 8 * CL_SIZE * nbat_fstride * sizeof(real));

            for (ic = 0; ic < CL_SIZE; ic++) { // CL_SIZE = 8
              ia               = ic;
              is               = ia*nbat_xstride;
              ifs              = ia*nbat_fstride;
              ix               = shX + xii[is+0];
              iy               = shY + xii[is+1];
              iz               = shZ + xii[is+2];
              iq               = facel*xii[is+3];
              nti              = ntype*2*typeii[ia];
              fix              = 0;
              fiy              = 0;
              fiz              = 0;

              for (jc = 0; jc < CL_SIZE; jc++) { // CL_SIZE = 8

                  ja               = cj*CL_SIZE + jc;
                  if (nbln.shift == CENTRAL && ci == cj && ja + cj * CL_SIZE <= ia + ci * CL_SIZE) continue;
                  int_bit = ((excl_pair[jc>>2][(jc & 3)*CL_SIZE+ic] >> (jm*NCL_PER_SUPERCL+im)) & 1);

                  js               = ja*nbat_xstride;
                  jfs              = ja*nbat_fstride;
                  jx               = xjj[js+0];
                  jy               = xjj[js+1];
                  jz               = xjj[js+2];
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

                  qq               = iq*xjj[js+3];
                  r     = rsq*rinv;
                  rt    = r*tabq_scale;
                  n0    = rt;
                  eps   = rt - n0;
                  getAll((void *) &Ftab_addr[n0], (void *) Ftab, sizeof(real) * 2);
                  fexcl = (1 - eps)*Ftab[0] + eps*Ftab[1];
                  fscal = qq*(int_bit*rinvsq - fexcl)*rinv;
                  // if (bEner) vcoul = qq*((int_bit - gmx_erf(ewaldcoeff_q*r))*rinv - int_bit*sh_ewald); !!!
                  if (rsq < rvdw2) {
                      tj        = nti + 2*typejj[ja];

                      // Vanilla Lennard-Jones cutoff
                      getAll((void *) &vdwparam_addr[tj], (void *) vdwparam, sizeof(real) * 2);
                      c6        = vdwparam[0];
                      c12       = vdwparam[1];

                      rinvsix   = int_bit*rinvsq*rinvsq*rinvsq;
                      Vvdw_disp = c6*rinvsix;
                      Vvdw_rep  = c12*rinvsix*rinvsix;
                      fscal    += (Vvdw_rep - Vvdw_disp)*rinvsq;

                      /*
                      if (bEner) {
                          vctot[n]   += vcoul;
                          Vvdwtot[n] +=
                              (Vvdw_rep - int_bit*c12*sh_invrc6*sh_invrc6)/12 -
                              (Vvdw_disp - int_bit*c6*sh_invrc6)/6;
                      }
                      */
                  }

                  tx        = fscal*dx;
                  ty        = fscal*dy;
                  tz        = fscal*dz;
                  fix       = fix + tx;
                  fiy       = fiy + ty;
                  fiz       = fiz + tz;
                  fjj[jfs+0] -= tx;
                  fjj[jfs+1] -= ty;
                  fjj[jfs+2] -= tz;
              } /* jc */

              fii[ifs+0]        += fix;
              fii[ifs+1]        += fiy;
              fii[ifs+2]        += fiz;
              fshift[0]   = fshift[0] + fix;
              fshift[1]   = fshift[1] + fiy;
              fshift[2]   = fshift[2] + fiz;
            } /* ic */

            /* put f back */
          } /* if */
        } /* im_Loop */
      } /* jm_Loop */
    } /* cj4_Loop */
    /* put fshift back */
  }
}
