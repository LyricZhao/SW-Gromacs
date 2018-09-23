/*
  Paras Summary:
    nnbl(int)
    nbat(const nbnxn_atomdata_t *)
    clearF(int)
    force_flags(int)
    fshift(real *)
    coult(int)
    vdwt(int)
    ic(const interaction_const_t *)
    nbl(nbnxn_pairlist_t **)
    shift_vec(rvec *)

  Function Summary:
    clear_f, clear_fshift, p_nbk_c_noener[][], p_nbk_c_ener[][], p_nbk_c_energrp[][]
 */

 __thread_kernel("compute") int nnbl, clearF, force_flags, coult, vdwt;
 __thread_kernel("compute") const nbnxn_atomdata_t *nbat;
 __thread_kernel("compute") real *fshift;
 __thread_kernel("compute") const interaction_const_t *ic;
 __thread_kernel("compute") nbnxn_pairlist_t **nbl;
 __thread_kernel("compute") rvec *shift_vec;

void set_sw_paras(long *paras) {
  nnbl = (int) paras[0]; nbat = (const nbnxn_atomdata_t *) paras[1]; clearF = (int) paras[2];
  force_flags = (int) paras[3]; fshift = (real *) paras[4]; coult = (int) paras[5];
  vdwt = (int) paras[6]; ic = (const interaction_const_t *) paras[7]; nbl = (nbnxn_pairlist_t **) paras[8];
  shift_vec = (rvec *) paras[9];
}

void sw_slave_kernel(long *paras) {
  set_sw_paras(paras);
  int threadID = athread_get_id(-1);
  int threadST = threadID * (nnbl / 64);
  int threadED = (threadID + 1) * (nnbl / 64);
  if(threadID == 63) threadED = nnbl;

  int nb;
  for(nb = threadST; nb < threadED; ++ nb) {
    nbnxn_atomdata_output_t *out;
    real                    *fshift_p;

    out = &nbat->out[nb]; // nbat

    if (clearF == enbvClearFYes) // clearF
    {
        clear_f(nbat, nb, out->f); // clear_f
    }

    if ((force_flags & GMX_FORCE_VIRIAL) && nnbl == 1) // force_flags
    {
        fshift_p = fshift; // fshift
    }
    else
    {
        fshift_p = out->fshift;

        if (clearF == enbvClearFYes)
        {
            clear_fshift(fshift_p); // clear_fshift
        }
    }

    if (!(force_flags & GMX_FORCE_ENERGY))
    {
        /* Don't calculate energies */
        p_nbk_c_noener[coult][vdwt](nbl[nb], nbat,
                                    ic,
                                    shift_vec,
                                    out->f,
                                    fshift_p); // all the functions
                                    // coult, vdwt, ic, nbl, shift_vec
    }
    else if (out->nV == 1)
    {
        /* No energy groups */
        out->Vvdw[0] = 0;
        out->Vc[0]   = 0;

        p_nbk_c_ener[coult][vdwt](nbl[nb], nbat,
                                  ic,
                                  shift_vec,
                                  out->f,
                                  fshift_p,
                                  out->Vvdw,
                                  out->Vc); // all the functions
    }
    else
    {
        /* Calculate energy group contributions */
        int i;

        for (i = 0; i < out->nV; i++)
        {
            out->Vvdw[i] = 0;
        }
        for (i = 0; i < out->nV; i++)
        {
            out->Vc[i] = 0;
        }

        p_nbk_c_energrp[coult][vdwt](nbl[nb], nbat,
                                     ic,
                                     shift_vec,
                                     out->f,
                                     fshift_p,
                                     out->Vvdw,
                                     out->Vc);
    }
  }
}
