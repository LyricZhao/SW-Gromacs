void sw_computing_core(long *paras) {
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
                          vctot   += vcoul;
                          Vvdwtot +=
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
}
