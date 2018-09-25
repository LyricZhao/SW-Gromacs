/* SunWay TaihuLight - Gromacs Kernel */

#ifndef _nbnxn_kernel_sw_ref_h
#define _nbnxn_kernel_sw_ref_h

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Wrapper call for the non-bonded n vs n reference kernels */
void
nbnxn_kernel_sw_ref(const nbnxn_pairlist_set_t *nbl_list,
                 const nbnxn_atomdata_t     *nbat,
                 const interaction_const_t  *ic,
                 rvec                       *shift_vec,
                 int                         force_flags,
                 int                         clearF,
                 real                       *fshift,
                 real                       *Vc,
                 real                       *Vvdw);

#ifdef __cplusplus
}
#endif

#endif
