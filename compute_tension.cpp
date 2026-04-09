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

#include "compute_tension.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "dihedral.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "update.h"

#include <cctype>
#include <cstring>
using namespace LAMMPS_NS;

ComputeTension::ComputeTension(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), vptr(nullptr), id_temp(nullptr), pstyle(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR,"compute tension", error);
  if (igroup) error->all(FLERR,"Compute tension must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  tensionflag = 1;
  timeflag = 1;
  // store temperature ID used by tension computation
  // ensure it is valid for temperature computation

  if (strcmp(arg[3],"NULL") == 0) {
    id_temp = nullptr;
  } else {
    id_temp = utils::strdup(arg[3]);
    auto icompute = modify->get_compute_by_id(id_temp);
    if (!icompute)
      error->all(FLERR,"Could not find compute tension temperature ID {}", id_temp);
    if (!icompute->tempflag)
      error->all(FLERR,"Compute tension temperature ID {} does not compute temperature", id_temp);
  }

  // process optional args

  pairhybridflag = 0;
  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"ke") == 0) keflag = 1;

      else if (strcmp(arg[iarg],"pair/hybrid") == 0) {
        if (lmp->suffix)
          pstyle = utils::strdup(fmt::format("{}/{}",arg[++iarg],lmp->suffix));
        else
          pstyle = utils::strdup(arg[++iarg]);

        nsub = 0;

        if (narg > iarg) {
          if (isdigit(arg[iarg][0])) {
            nsub = utils::inumeric(FLERR,arg[iarg],false,lmp);
            ++iarg;
            if (nsub <= 0)
              error->all(FLERR,"Illegal compute tension command");
          }
        }

        // check if pair style with and without suffix exists

        pairhybrid = (Pair *) force->pair_match(pstyle,1,nsub);
        if (!pairhybrid && lmp->suffix) {
          pstyle[strlen(pstyle) - strlen(lmp->suffix) - 1] = '\0';
          pairhybrid = (Pair *) force->pair_match(pstyle,1,nsub);
        }

        if (!pairhybrid)
          error->all(FLERR,"Unrecognized pair style in compute tension command");

        pairhybridflag = 1;
      }
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
        kspaceflag = fixflag = 1;
      } else error->all(FLERR,"Illegal compute tension command");
      iarg++;
    }
  }

  // error check

  if (keflag && id_temp == nullptr)
    error->all(FLERR,"Compute tension requires temperature ID to include kinetic energy");

  vector = new double[size_vector];
  nvirial = 0;
  vptr = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeTension::~ComputeTension()
{
  delete[] id_temp;
  delete[] vector;
  delete[] vptr;
  delete[] pstyle;
}


void ComputeTension::init()
{
  boltz = force->boltz;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (keflag) {
    temperature = modify->get_compute_by_id(id_temp);
    if (!temperature)
      error->all(FLERR,"Could not find compute tension temperature ID {}", id_temp);
  }

  // recheck if pair style with and without suffix exists

  if (pairhybridflag) {
    pairhybrid = (Pair *) force->pair_match(pstyle,1,nsub);
    if (!pairhybrid && lmp->suffix) {
      strcat(pstyle,"/");
      strcat(pstyle,lmp->suffix);
      pairhybrid = (Pair *) force->pair_match(pstyle,1,nsub);
    }

    if (!pairhybrid)
      error->all(FLERR,"Unrecognized pair style in compute tension command");
  }

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = nullptr;

  if (pairhybridflag && force->pair) nvirial++;
  if (pairflag && force->pair) nvirial++;
  if (atom->molecular != Atom::ATOMIC) {
    if (bondflag && force->bond) nvirial++;
    if (angleflag && force->angle) nvirial++;
    if (dihedralflag && force->dihedral) nvirial++;
    if (improperflag && force->improper) nvirial++;
  }
  if (fixflag)
    for (auto &ifix : modify->get_fix_list())
      if (ifix->thermo_virial) nvirial++;

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (pairhybridflag && force->pair) {
      auto ph = dynamic_cast<PairHybrid *>(force->pair);
      ph->no_virial_fdotr_compute = 1;
      vptr[nvirial++] = pairhybrid->virial;
    }
    if (pairflag && force->pair) vptr[nvirial++] = force->pair->virial;
    if (bondflag && force->bond) vptr[nvirial++] = force->bond->virial;
    if (angleflag && force->angle) vptr[nvirial++] = force->angle->virial;
    if (dihedralflag && force->dihedral)
      vptr[nvirial++] = force->dihedral->virial;
    if (improperflag && force->improper)
      vptr[nvirial++] = force->improper->virial;
    if (fixflag)
    for (auto &ifix : modify->get_fix_list())
      if (ifix->virial_global_flag && ifix->thermo_virial)
          vptr[nvirial++] = ifix->virial;
  }

  // flag Kspace contribution separately, since not summed across procs

  if (kspaceflag && force->kspace) kspace_virial = force->kspace->virial;
  else kspace_virial = nullptr;
}

/* ----------------------------------------------------------------------
   compute total tension, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputeTension::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      temperature->compute_scalar();
  }
  // Added by ZX ******************************************************
  N_atom = (double) atom->natoms;
  inv_natomecof = 0.5 / N_atom;
  // Added by ZX ******************************************************
  virial_compute(3,3);
  if (keflag)
    scalar = (temperature->dof * boltz * temperature->scalar +
                virial[0] + virial[1] + virial[2]) / 3.0  * inv_natomecof ;
                
  else
    scalar = (virial[0] + virial[1] + virial[2]) / 3.0  * inv_natomecof ;
  
  return scalar;
}
/* ----------------------------------------------------------------------
   compute tension tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */
void ComputeTension::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_tension_flag)
    error->all(FLERR,"Must use 'kspace_modify tension/scalar no' for "
               "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double *ke_tensor;
  if (keflag) { 
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
    ke_tensor = temperature->vector;
  }

  // Added by ZX ******************************************************
  N_atom = (double) atom->natoms;
  inv_natomecof = 0.5 / N_atom;
  // **************************************************************************
  virial_compute(6,3);
  if (keflag) {
    for (int i = 0; i < 6; i++)
      vector[i] = ( ke_tensor[i] + virial[i])  * inv_natomecof ;
  } else
    for (int i = 0; i < 6; i++)
      vector[i] = virial[i] * inv_natomecof ;
  
}


void ComputeTension::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

  // KSpace virial contribution is already summed across procs

  if (kspace_virial)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction, only if pair contributions are included

  if (force->pair && pairflag && force->pair->tail_flag)
    for (i = 0; i < ndiag; i++) virial[i] += force->pair->ptail * inv_natomecof;
}

/* ---------------------------------------------------------------------- */

void ComputeTension::reset_extra_compute_fix(const char *id_new)
{
  delete[] id_temp;
  id_temp = utils::strdup(id_new);
}
