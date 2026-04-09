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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ylz/poly,PairYLZPoly);
// clang-format on

#else

#ifndef LMP_PAIR_YLZ_POLY_H
#define LMP_PAIR_YLZ_POLY_H

#include "pair.h"

namespace LAMMPS_NS {

class PairYLZPoly : public Pair {
 public:
  PairYLZPoly(LAMMPS *lmp);
  ~PairYLZPoly() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void *extract(const char *, int &) override;
 
 
 protected:
  double cut_global;
  double **epsilon, **sigma, **cut, **zeta, **mu, **beta, **lambda;

  class AtomVecEllipsoid *avec;

  void allocate();
  double ylz_poly_analytic(const int i, const int j, double a1[3][3], double a2[3][3], double *r12,
                      const double rsq, double *fforce, double *ttor, double *rtor);
  
};
}    // namespace LAMMPS_NS
#endif
#endif
