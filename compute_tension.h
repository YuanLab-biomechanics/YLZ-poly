/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hongyan Yuan (SUSTech) hongyan6@outlook.com
                        Xin Zhang (SUSTech) 12432359@mail.sustech.edu.cn
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(tension,ComputeTension);
// clang-format on
#else

#ifndef LMP_COMPUTE_TENSION_H
#define LMP_COMPUTE_TENSION_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTension : public Compute {
 public:
  ComputeTension(class LAMMPS *, int, char **);
  ~ComputeTension() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  void reset_extra_compute_fix(const char *) override;

 protected:
  double boltz,  inv_natomecof,  N_atom;
  int nvirial;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];    // ordering: xx,yy,zz,xy,xz,yz
  int pairhybridflag;
  class Pair *pairhybrid;
  int keflag, pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int fixflag, kspaceflag;

  void virial_compute(int, int);

 private:
  char *pstyle;
  int nsub;
};

}    // namespace LAMMPS_NS

#endif
#endif
