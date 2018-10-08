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

# define buffer_size (256 * 1024 * 1024) // 256MB

void *tempPtrG;
int first_time_run = 0;

static inline void setP(void *dest, void *pa) {
  memcpy(dest, pa, 8);
}

static inline void addTrans(void* useless, void* src, int siz) {
  memcpy(tempPtrG, src, siz);
  tempPtrG += siz;
}

static inline void fakTrans(void* useless, void* src, int siz) {
  tempPtrG += siz;
}

char *transferData;
real *vctotCopy, *VvdwtotCopy, *vdwparamCopy, *FtabCopy, *fshift_sum;
void **startPoint;

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
    real *x = nbat->x;
    nbnxn_excl_t *excl[2];

    /* Attributes */
    const nbnxn_sci_t *nbln;
    int n, ish3, cj4_ind, cj4_ind0, cj4_ind1, sci;
    real shX, shY, shZ;
    real rcut2 = iconst->rcoulomb*iconst->rcoulomb;
    real rvdw2 = iconst->rvdw*iconst->rvdw;
    int *type = nbat->type;
    gmx_bool bEner;
    gmx_bool bEwald;
    bEner = (force_flags & GMX_FORCE_ENERGY);

    /* Array */
    if(first_time_run == 0) {
      transferData = (char *) malloc(buffer_size);
      vctotCopy = (real *) malloc(sizeof(real) * 512);
      VvdwtotCopy = (real *) malloc(sizeof(real) * 512);
      startPoint = (void *) malloc(sizeof(void *) * 512);
      vdwparamCopy = (real *) malloc(4096 * sizeof(real));
      FtabCopy = (real *) malloc(4096 * sizeof(real));
      fshift_sum = (real *) malloc(2048 * sizeof(real));
      first_time_run = 1;
    }

    /* Global Data */
    int vdwparam_size = nbat->ntype * nbat->ntype * 2 * sizeof(real);
    memcpy((void *) vdwparamCopy, (void *) nbat->nbfp, vdwparam_size);
    int Ftab_size = iconst->tabq_size * sizeof(real);
    memcpy((void *) FtabCopy, (void *) iconst->tabq_coul_F, Ftab_size);

    /* Temp Vars */
    int im, ci, ic, ia, jm, cj;
    int is, ifs, js, jfs, jc, ja;
    real iq, facel = iconst->epsfac;

    /* Loop Buffer Data */
    tempPtrG = transferData;

    real *shiftvec = shift_vec[0];

    // printf("mpe_start_addr = %p\n", transferData);
    for (n = 0; n < nbl->nsci; n++) {
        nbln = &nbl->sci[n];
        ish3             = 3*nbln->shift;
        shX              = shiftvec[ish3];
        shY              = shiftvec[ish3+1];
        shZ              = shiftvec[ish3+2];
        cj4_ind0         = nbln->cj4_ind_start;
        cj4_ind1         = nbln->cj4_ind_end;
        sci              = nbln->sci;
        vctotCopy[n]         = 0;
        VvdwtotCopy[n]       = 0;
        startPoint[n] = tempPtrG;
        addTrans(tempPtrG, (void *) &ish3, sizeof(int));
        addTrans(tempPtrG, (void *) &cj4_ind0, sizeof(int));
        addTrans(tempPtrG, (void *) &cj4_ind1, sizeof(int));
        addTrans(tempPtrG, (void *) &sci, sizeof(int));
        addTrans(tempPtrG, (void *) &shiftvec[ish3], sizeof(real) * 3);
        // printf("?? mpe: n = %d %d %d %.2lf %.2lf %.2lf\n", n, cj4_ind0, cj4_ind1, shX, shY, shZ);
        // int *ppp = (int *) startPoint[n];
        // printf("^^ %p %d %d %d %d\n", ppp, ppp[0], ppp[1], ppp[2], ppp[3]);

        /* MPE Works */
        if (nbln->shift == CENTRAL && nbl->cj4[cj4_ind0].cj[0] == sci*NCL_PER_SUPERCL) {
            for (im = 0; im < NCL_PER_SUPERCL; im++) {
                ci = sci*NCL_PER_SUPERCL + im;
                for (ic = 0; ic < CL_SIZE; ic++) {
                    ia = ci*CL_SIZE + ic;
                    iq = x[ia*nbat->xstride+3];
                    vctotCopy[n] += iq*iq;
                }
            }
            vctotCopy[n] *= -facel*iconst->ewaldcoeff_q*M_1_SQRTPI;
        }

        /* Data Package */
        for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {

            excl[0] = &nbl->excl[nbl->cj4[cj4_ind].imei[0].excl_ind];
            excl[1] = &nbl->excl[nbl->cj4[cj4_ind].imei[1].excl_ind];

            addTrans(tempPtrG, (void *) &(excl[0]->pair[0]), 32 * sizeof(unsigned int));
            addTrans(tempPtrG, (void *) &(excl[1]->pair[0]), 32 * sizeof(unsigned int));
            addTrans(tempPtrG, (void *) &(nbl->cj4[cj4_ind].imei[0].imask), sizeof(unsigned int));

            // printf("imask mpe: %d %u\n", cj4_ind, nbl->cj4[cj4_ind].imei[0].imask);

            /* kernel starts here */
            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) { // NBNXN_GPU_JGROUP_SIZE = 4

                cj = nbl->cj4[cj4_ind].cj[jm];
                // printf("mpe cj = %d\n", cj);
                addTrans(tempPtrG, (void *) &cj, sizeof(int));

                for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
                    if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

                        ci = sci * NCL_PER_SUPERCL + im;

                        /* ic: 0 to 8 */
                        addTrans(tempPtrG, (void *) &x[(ci * CL_SIZE) * nbat->xstride], (7 * nbat->xstride + 4) * sizeof(real));
                        addTrans(tempPtrG, (void *) &type[ci * CL_SIZE], 8 * sizeof(int));
                        fakTrans(tempPtrG, (void *) &f[(ci * CL_SIZE) * nbat->fstride], (7 * nbat->fstride + 3) * sizeof(real));

                        /* jc: 0 to 8 */
                        addTrans(tempPtrG, (void *) &x[(cj * CL_SIZE) * nbat->xstride], (7 * nbat->xstride + 4) * sizeof(real));
                        addTrans(tempPtrG, (void *) &type[cj * CL_SIZE], 8 * sizeof(int));
                        fakTrans(tempPtrG, (void *) &f[(cj * CL_SIZE) * nbat->fstride], (7 * nbat->fstride + 3 )* sizeof(real));
                    }
                }
            }
        }
    }

    long long_nsci = nbl->nsci;
    long long_nxstride = nbat->xstride;
    long long_nfstride = nbat->fstride;
    long long_ntype = nbat->ntype;
    long long_bEner = bEner;

    void **bsp = startPoint;

    /* Calling Kernel */
    char paras[20 * 8];
    // printf("immmp: %d %d %d %d %lf %lf\n", nbl->nsci, nbat->xstride, nbat->fstride, nbat->ntype, rcut2, rvdw2);
    setP((void *) &paras[ 0 * 8], (void *) &(long_nsci));
    setP((void *) &paras[ 1 * 8], (void *) &(bsp));
    setP((void *) &paras[ 2 * 8], (void *) &(long_nxstride));
    setP((void *) &paras[ 3 * 8], (void *) &(long_nfstride));
    setP((void *) &paras[ 4 * 8], (void *) &(long_ntype));
    setP((void *) &paras[ 5 * 8], (void *) &(rcut2));
    setP((void *) &paras[ 6 * 8], (void *) &(rvdw2));
    setP((void *) &paras[ 7 * 8], (void *) &(vctotCopy));
    setP((void *) &paras[ 8 * 8], (void *) &(VvdwtotCopy));
    setP((void *) &paras[ 9 * 8], (void *) &(vdwparamCopy));
    setP((void *) &paras[10 * 8], (void *) &(FtabCopy));
    setP((void *) &paras[11 * 8], (void *) &(iconst->k_rf));
    setP((void *) &paras[12 * 8], (void *) &(iconst->c_rf));
    setP((void *) &paras[13 * 8], (void *) &(iconst->tabq_scale));
    setP((void *) &paras[14 * 8], (void *) &(iconst->ewaldcoeff_q));
    setP((void *) &paras[15 * 8], (void *) &(iconst->sh_ewald));
    setP((void *) &paras[16 * 8], (void *) &(iconst->sh_invrc6));
    setP((void *) &paras[17 * 8], (void *) &(iconst->epsfac));
    setP((void *) &paras[18 * 8], (void *) &(long_bEner));
    setP((void *) &paras[19 * 8], (void *) &(fshift_sum));

    printf("mpe_addr: %p\n", f);
    athread_spawn(sw_computing_core, &paras);
    athread_join();

    if(bEner) {

      for (n = 0; n < nbl->nsci; n++) {
        nbln = &nbl->sci[n];
        ish3             = 3*nbln->shift;
        fshift[ish3 + 0] += fshift_sum[n * 3 + 0];
        fshift[ish3 + 1] += fshift_sum[n * 3 + 1];
        fshift[ish3 + 2] += fshift_sum[n * 3 + 2];
        Vc[0]         = Vc[0]   + vctotCopy[n];
        Vvdw[0]       = Vvdw[0] + VvdwtotCopy[n];
      }
    }

    tempPtrG = transferData;
    for (n = 0; n < nbl->nsci; ++ n) {
      nbln = &nbl->sci[n];
      cj4_ind0         = nbln->cj4_ind_start;
      cj4_ind1         = nbln->cj4_ind_end;
      sci              = nbln->sci;
      tempPtrG += 4 * sizeof(int) + 3 * sizeof(real);
      for (cj4_ind = cj4_ind0; (cj4_ind < cj4_ind1); cj4_ind++) {
          tempPtrG += 65 * sizeof(unsigned int);
          for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++) {
            cj = nbl->cj4[cj4_ind].cj[jm];
            tempPtrG += sizeof(int);

            for (im = 0; im < NCL_PER_SUPERCL; im++) { // NCL_PER_SUPERCL = 2 * 2 * 2 = 8
                if ((nbl->cj4[cj4_ind].imei[0].imask >> (jm*NCL_PER_SUPERCL+im)) & 1) {

                    ci = sci * NCL_PER_SUPERCL + im;

                    /* putback here */
                    tempPtrG += (7 * nbat->xstride + 4) * sizeof(real) + 8 * sizeof(int);
                    real *iip = tempPtrG;
                    for(ic = 0; ic < 8; ++ ic) {
                      ia = ci * CL_SIZE + ic;
                      ifs = ia*nbat->fstride;
                      f[ifs+0] += iip[ifs-ic*nbat->fstride];
                      f[ifs+1] += iip[ifs-ic*nbat->fstride];
                      f[ifs+2] += iip[ifs-ic*nbat->fstride];
                    }
                    tempPtrG += (7 * nbat->fstride + 3) * sizeof(real);

                    /* putback here */
                    tempPtrG += (7 * nbat->xstride + 4) * sizeof(real) + 8 * sizeof(int);
                    real *jjp = tempPtrG;
                    for(jc = 0; jc < 8; ++ jc) {
                      ja = cj * CL_SIZE + jc;
                      jfs = ja*nbat->fstride;
                      f[jfs+0] += iip[jfs-jc*nbat->fstride];
                      f[jfs+1] += iip[jfs-jc*nbat->fstride];
                      f[jfs+2] += iip[jfs-jc*nbat->fstride];
                    }
                    tempPtrG += (7 * nbat->fstride + 3) * sizeof(real);
                }
            }

          }
      }
    }
    while(1);
}
