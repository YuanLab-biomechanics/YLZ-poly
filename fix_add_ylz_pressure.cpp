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

#include "fix_add_ylz_pressure.h"
#include "atom.h"
#include "atom_masks.h"
#include "atom_vec_ellipsoid.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_extra.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"
#include "utils.h"
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAddYLZPressure::FixAddYLZPressure(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), avec(nullptr), ellipsoid(nullptr), pstr(nullptr), 
    idregion(nullptr), region(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix add/ylz/pressure", error);

  dynamic_group_allow = 1;
  scalar_flag = 0;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  energy_global_flag = 0;
  virial_global_flag = 0;
  respa_level_support = 1;
  ilevel_respa = 0;

  // 解析压力参数
  if (utils::strmatch(arg[3], "^v_")) {
    pstr = utils::strdup(arg[3] + 2);
    pstyle = EQUAL;
  } else {
    pressure_value = utils::numeric(FLERR, arg[3], false, lmp);
    pstyle = CONSTANT;
  }

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix add/ylz/pressure region", error);
      delete[] idregion;
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown fix add/ylz/pressure keyword: {}", arg[iarg]);
    }
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixAddYLZPressure::~FixAddYLZPressure()
{
  delete[] pstr;
  delete[] idregion;
}

/* ---------------------------------------------------------------------- */

int FixAddYLZPressure::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::init()
{
  // 确认原子类型为 ellipsoid
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR, "Fix add/ylz/pressure requires atom style ellipsoid");
  ellipsoid = atom->ellipsoid;

  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0) error->all(FLERR, "Variable {} for fix add/ylz/pressure does not exist", pstr);
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR, "Variable {} for fix add/ylz/pressure must be equal-style", pstr);
  }

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix add/ylz/pressure does not exist", idregion);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels - 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::post_force(int /* vflag */)
{
  double pressure = 0.0;
  if (pstyle == EQUAL) {
    modify->clearstep_compute();
    pressure = input->variable->compute_equal(pvar);
    modify->addstep_compute(update->ntimestep + 1);
  } else {
    pressure = pressure_value;
  }
  
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;
  
  if (region) region->prematch();
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (ellipsoid[i] < 0) continue;
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
      
      double *q = bonus[ellipsoid[i]].quat;
      double nx = 1.0 - 2.0 * (q[2]*q[2] + q[3]*q[3]);
      double ny = 2.0 * (q[1]*q[2] + q[0]*q[3]);
      double nz = 2.0 * (q[1]*q[3] - q[0]*q[2]);
      
      foriginal[0] += f[i][0];
      foriginal[1] += f[i][1];
      foriginal[2] += f[i][2];
      
      f[i][0] -= pressure * nx;
      f[i][1] -= pressure * ny;
      f[i][2] -= pressure * nz;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddYLZPressure::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixAddYLZPressure::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  if (n >= 0 && n < 3) return foriginal_all[n];
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixAddYLZPressure::memory_usage()
{
  return 0.0;
}