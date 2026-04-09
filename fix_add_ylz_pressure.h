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

#ifdef FIX_CLASS
// clang-format off
FixStyle(add/ylz/pressure,FixAddYLZPressure);
// clang-format on
#else

#ifndef LMP_FIX_ADD_YLZ_PRESSURE_H
#define LMP_FIX_ADD_YLZ_PRESSURE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddYLZPressure : public Fix {
 public:
  FixAddYLZPressure(class LAMMPS *, int, char **);
  ~FixAddYLZPressure() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_vector(int) override;
  double memory_usage() override;

  enum { CONSTANT, EQUAL };

 protected:
  double pressure_value;
  char *pstr;
  int pstyle;
  int pvar;

  class AtomVecEllipsoid *avec;
  int *ellipsoid;

  char *idregion;
  class Region *region;

  double foriginal[3], foriginal_all[3];
  int force_flag;
  int ilevel_respa;
};

} // namespace LAMMPS_NS

#endif
#endif
