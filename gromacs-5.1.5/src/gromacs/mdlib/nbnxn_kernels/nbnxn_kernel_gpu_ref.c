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

# define vc_max 512

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

    /* Defined Variables */
    real vctot[vc_max];
    real Vvdwtot[vc_max];
    assert(vc_max >= nbl -> nsci);

    /* Variables Init */
    if (clearF == enbvClearFYes) clear_f(nbat, 0, f);

    /* Calling Kernel */
    long paras[9] = {(long) nbl, (long) nbat, (long) iconst, (long) shift_vec, (long) force_flags, (long) f, (long) fshift, (long) vctot, (long) Vvdwtot};
    athread_spawn64(SLAVE_FUN(sw_computing_core), paras);
    athread_join64();

    /* Reduce */
    if(bEner) {
      int n;
      for (n = 0; n < nbl->nsci; n++) {
        Vc[0]         = Vc[0]   + vctot[n];
        Vvdw[0]       = Vvdw[0] + Vvdwtot[n];
      }
    }
}
