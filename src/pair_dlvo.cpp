/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_dlvo.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define PI 3.141592653589793238462643383279

// Note: we assume that all DVO particles have the same Hamaker, Radius, Psi, and SpringK.
// The code doesn't check that you obey this.  Generalizing this will be left as a future project.

/* ---------------------------------------------------------------------- */

PairDLVO::PairDLVO(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairDLVO::~PairDLVO()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(hamaker);
    memory->destroy(radius);
    memory->destroy(psisqrd);
    memory->destroy(debyeinv);
    memory->destroy(springk);
    memory->destroy(offset_coul);
    memory->destroy(offset_vdwl);
  }
}

/* ---------------------------------------------------------------------- */

void PairDLVO::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,forcelj,factor_lj;
  double r,rinv,h,hinv,h2inv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double epsr = force->dielectric;
  double epso = 1.0/(4*PI*force->qqr2e);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r = sqrt(rsq);
	rinv = 1.0/r;
 	h = r - 2*radius[itype][jtype]; 
	hinv  = 1.0/h;
	h2inv = hinv*hinv;
	  
        fpair = -epsr*epso*radius[itype][jtype]*psisqrd[itype][jtype]*debyeinv[itype][jtype]/(exp(-debyeinv[itype][jtype]*h)+1)*0.5
	        + radius[itype][jtype]*hamaker[itype][jtype]*0.08333333333333333333333333333*h2inv;
	if(h<0)
	  fpair += -springk[itype][jtype]*h;

	fpair *= factor_lj*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = epsr*(-epso*radius[itype][jtype]*psisqrd[itype][jtype]*0.5*log(1+exp(-debyeinv[itype][jtype]*h)) 
			+ offset_coul[itype][jtype]) 
	    + radius[itype][jtype]*hamaker[itype][jtype]*0.083333333333333*hinv 
	    + offset_vdwl[itype][jtype];
	  if(h<0)
	    evdwl += 0.5*springk[itype][jtype]*h*h;

          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairDLVO::compute_inner()
{
}

/* ---------------------------------------------------------------------- */

void PairDLVO::compute_middle()
{
}

/* ---------------------------------------------------------------------- */

void PairDLVO::compute_outer(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDLVO::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(hamaker,n+1,n+1,"pair:hamaker");
  memory->create(radius,n+1,n+1,"pair:radius");
  memory->create(psisqrd,n+1,n+1,"pair:psisqrd");
  memory->create(debyeinv,n+1,n+1,"pair:debyeinv");
  memory->create(springk,n+1,n+1,"pair:springk");
  memory->create(offset_coul,n+1,n+1,"pair:offset_coul");
  memory->create(offset_vdwl,n+1,n+1,"pair:offset_vdwl");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDLVO::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDLVO::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 87)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double hamaker_one = force->numeric(FLERR,arg[2]);
  double radius_one = force->numeric(FLERR,arg[3]);
  double psi_one = force->numeric(FLERR,arg[4]);
  double debyeinv_one = force->numeric(FLERR,arg[5]);
  double springk_one = force->numeric(FLERR,arg[6]);

  double cut_one = cut_global;
  if (narg == 8) cut_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      hamaker[i][j] = hamaker_one;
      radius[i][j] = radius_one;
      psisqrd[i][j] = psi_one*psi_one;
      debyeinv[i][j] = debyeinv_one;
      springk[i][j] = springk_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDLVO::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this,instance_me);
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this,instance_me);

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairDLVO::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDLVO::init_one(int i, int j)
{
  double r,h,hinv;
  double epso = 1.0/(4*PI*force->qqr2e);

  /*  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  */
  if (offset_flag) {
    h = cut[i][j];
    hinv  = 1.0/h;
    offset_coul[i][j] = epso*radius[i][j]*psisqrd[i][j]*0.5*log(1+exp(-debyeinv[i][j]*h));
    offset_vdwl[i][j] = -radius[i][j]*hamaker[i][j]*0.083333333333333*hinv;
  } else {
    offset_coul[i][j] = 0.0;
    offset_vdwl[i][j] = 0.0;
  }
  /*
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];
  */
  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

  }

  return cut[i][j] + 2*radius[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDLVO::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&hamaker[i][j],sizeof(double),1,fp); 
	fwrite(&radius[i][j],sizeof(double),1,fp);
	fwrite(&psisqrd[i][j],sizeof(double),1,fp);
	fwrite(&debyeinv[i][j],sizeof(double),1,fp);
	fwrite(&springk[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDLVO::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&hamaker[i][j],sizeof(double),1,fp);
          fread(&radius[i][j],sizeof(double),1,fp);
	  fread(&psisqrd[i][j],sizeof(double),1,fp);
	  fread(&debyeinv[i][j],sizeof(double),1,fp);
	  fread(&springk[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&hamaker[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&radius[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&psisqrd[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&debyeinv[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&springk[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDLVO::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDLVO::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDLVO::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,hamaker[i][i],radius[i][i],psisqrd[i][i],debyeinv[i][i],springk[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDLVO::write_data_all(FILE *fp)
{

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,hamaker[i][j],radius[i][j],psisqrd[i][j],debyeinv[i][j],springk[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDLVO::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double rinv,r,h,hinv,h2inv,fpair,philj;
  double epsr = force->dielectric;
  double epso = 1.0/(4*PI*force->qqr2e);

  r = sqrt(rsq);
  rinv = 1.0/r;
  h = r - 2*radius[itype][jtype]; 
  hinv  = 1.0/h;
  h2inv = hinv*hinv;
  
  fpair = -epsr*epso*radius[itype][jtype]*psisqrd[itype][jtype]*debyeinv[itype][jtype]/(exp(-debyeinv[itype][jtype]*h)+1)*0.5
    + radius[itype][jtype]*hamaker[itype][jtype]*0.08333333333333333333333333333*h2inv;
  if(h<0)
    fpair += -springk[itype][jtype]*h;
  
  fforce *= factor_lj*rinv;
  
  philj = epsr*(-epso*radius[itype][jtype]*psisqrd[itype][jtype]*0.5*log(1+exp(-debyeinv[itype][jtype]*h))
		+ offset_coul[itype][jtype])
    + radius[itype][jtype]*hamaker[itype][jtype]*0.083333333333333*hinv
    + offset_vdwl[itype][jtype];
  if(h<0)
    philj += 0.5*springk[itype][jtype]*h*h;
  
  return philj*factor_lj;
}

/* ---------------------------------------------------------------------- */

void *PairDLVO::extract(const char *str, int &dim)
{
  dim = 5;
  if (strcmp(str,"hamaker") == 0) return (void *) hamaker;
  if (strcmp(str,"radius") == 0) return (void *) radius;
  if (strcmp(str,"psisqrd") == 0) return (void *) psisqrd;
  if (strcmp(str,"debyeinv") == 0) return (void *) debyeinv;
  if (strcmp(str,"springk") == 0) return (void *) springk;
  return NULL;
}

