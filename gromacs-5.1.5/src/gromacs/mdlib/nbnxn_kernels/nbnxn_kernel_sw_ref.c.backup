/* SunWay TaihuLight - Gromacs Kernel */

#include "gmxpre.h"

#include "nbnxn_kernel_sw_ref.h"

#include "config.h"

#include <assert.h>
#include <math.h>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Typedefs for declaring lookup tables of kernel functions.
 */

typedef void (*p_nbk_func_noener)(const nbnxn_pairlist_t     *nbl,
                                  const nbnxn_atomdata_t     *nbat,
                                  const interaction_const_t  *ic,
                                  rvec                       *shift_vec,
                                  real                       *f,
                                  real                       *fshift);

typedef void (*p_nbk_func_ener)(const nbnxn_pairlist_t     *nbl,
                                const nbnxn_atomdata_t     *nbat,
                                const interaction_const_t  *ic,
                                rvec                       *shift_vec,
                                real                       *f,
                                real                       *fshift,
                                real                       *Vvdw,
                                real                       *Vc);

/* Analytical reaction-field kernels */
#define CALC_COUL_RF
#define LJ_CUT
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
#undef CALC_COUL_RF


/* Tabulated exclusion interaction electrostatics kernels */
#define CALC_COUL_TAB
#define LJ_CUT
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
/* Twin-range cut-off kernels */
#define VDW_CUTOFF_CHECK
#define LJ_CUT
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_CUT
#define LJ_FORCE_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_FORCE_SWITCH
#define LJ_POT_SWITCH
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_POT_SWITCH
#define LJ_EWALD
#define LJ_CUT
#define LJ_EWALD_COMB_GEOM
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_GEOM
#define LJ_EWALD_COMB_LB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_includes.h"
#undef LJ_EWALD_COMB_LB
#undef LJ_CUT
#undef LJ_EWALD
#undef VDW_CUTOFF_CHECK
#undef CALC_COUL_TAB


enum {
    coultRF, coultTAB, coultTAB_TWIN, coultNR
};

enum {
    vdwtCUT, vdwtFSWITCH, vdwtPSWITCH, vdwtEWALDGEOM, vdwtEWALDLB, vdwtNR
};

p_nbk_func_noener p_nbk_c_noener[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ElecRF_VdwLJ_F_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_F_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_F_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_ref           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_F_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_F_ref        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_F_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_F_ref }
};

p_nbk_func_ener p_nbk_c_ener[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ElecRF_VdwLJ_VF_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_VF_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_VF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_ref            },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VF_ref         },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VF_ref  }
};

p_nbk_func_ener p_nbk_c_energrp[coultNR][vdwtNR] =
{
    { nbnxn_kernel_ElecRF_VdwLJ_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJFsw_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJPsw_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VgrpF_ref,           nbnxn_kernel_ElecRF_VdwLJEwCombLB_VgrpF_ref           },
    { nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJFsw_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJPsw_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombGeom_VgrpF_ref,        nbnxn_kernel_ElecQSTab_VdwLJEwCombLB_VgrpF_ref        },
    { nbnxn_kernel_ElecQSTabTwinCut_VdwLJ_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJFsw_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJPsw_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombGeom_VgrpF_ref, nbnxn_kernel_ElecQSTabTwinCut_VdwLJEwCombLB_VgrpF_ref }
};

#include "athread.h"
extern SLAVE_FUN(sw_slave_kernel)();

void
nbnxn_kernel_sw_ref(const nbnxn_pairlist_set_t *nbl_list,
                 const nbnxn_atomdata_t     *nbat,
                 const interaction_const_t  *ic,
                 rvec                       *shift_vec,
                 int                         force_flags,
                 int                         clearF,
                 real                       *fshift,
                 real                       *Vc,
                 real                       *Vvdw)
{
    int                nnbl;
    nbnxn_pairlist_t **nbl;
    int                coult;
    int                vdwt;
    int                nb;
    int                nthreads gmx_unused;

    nnbl = nbl_list->nnbl;
    nbl  = nbl_list->nbl;

    if (EEL_RF(ic->eeltype) || ic->eeltype == eelCUT)
    {
        coult = coultRF;
    }
    else
    {
        if (ic->rcoulomb == ic->rvdw)
        {
            coult = coultTAB;
        }
        else
        {
            coult = coultTAB_TWIN;
        }
    }

    if (ic->vdwtype == evdwCUT)
    {
        switch (ic->vdw_modifier)
        {
            case eintmodPOTSHIFT:
            case eintmodNONE:
                vdwt = vdwtCUT;
                break;
            case eintmodFORCESWITCH:
                vdwt = vdwtFSWITCH;
                break;
            case eintmodPOTSWITCH:
                vdwt = vdwtPSWITCH;
                break;
            default:
                gmx_incons("Unsupported VdW modifier");
                break;
        }
    }
    else if (ic->vdwtype == evdwPME)
    {
        if (ic->ljpme_comb_rule == ljcrGEOM)
        {
            assert(nbat->comb_rule == ljcrGEOM);
            vdwt = vdwtEWALDGEOM;
        }
        else
        {
            assert(nbat->comb_rule == ljcrLB);
            vdwt = vdwtEWALDLB;
        }
    }
    else
    {
        gmx_incons("Unsupported vdwtype in nbnxn reference kernel");
    }

    /* here use athread instead of openmp */
    /* original: for (nb = 0; nb < nnbl; nb++) */
    /*
      Paras Summary:
        nnbl, nbat, clearF, force_flags, fshift, coult, vdwt, ic, nbl, shift_vec

      Function Summary:
        clear_f, clear_fshift, p_nbk_c_noener[][], p_nbk_c_ener[][], p_nbk_c_energrp[][]
     */
    nthreads = 64; /* fixed */
    athread_enter64();
    long paras[10] = {(long) nnbl, (long) nbat, (long) clearF, (long) force_flags, (long) fshift, (long) coult, (long) vdwt, (long) ic, (long) nbl, (long) shift_vec};
	  athread_spawn64(SLAVE_FUN(sw_slave_kernel), &paras);
	  athread_join64();

    if (force_flags & GMX_FORCE_ENERGY)
    {
        reduce_energies_over_lists(nbat, nnbl, Vvdw, Vc);
    }
}
