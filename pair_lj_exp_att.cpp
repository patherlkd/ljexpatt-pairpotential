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
   Contributing author: Luke Kristopher Davis (UCL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lj_exp_att.h"
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

/* ---------------------------------------------------------------------- */

PairLEA::PairLEA(LAMMPS *lmp) : Pair(lmp)
{
  extracted = 0;
  lj_switch = 1;
  respa_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLEA::~PairLEA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut_exp2);
    
    memory->destroy(cut);
    memory->destroy(cut_exp);
    memory->destroy(epsilon);
    memory->destroy(epsilon_exp);
    memory->destroy(sigma);
    memory->destroy(lambda);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(exp_sl);
    memory->destroy(exp_rcl);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLEA::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double r,rsq,r2inv,r6inv,forceexp,forcelj,factor_lj;
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

      if(rsq >= cut_exp2[itype][jtype] &&  rsq < cutsq[itype][jtype])
	{
	  r=sqrt(rsq);
	  //printf("LAMBDA[%d][%d]= %f\n",itype,jtype,lambda[itype][jtype]);
	  forceexp = (exp_sl[itype][jtype]*exp_rcl[itype][jtype] 
		    - exp_sl[itype][jtype]*exp(-r/lambda[itype][jtype]))*epsilon_exp[itype][jtype]/lambda[itype][jtype];

	  delx = delx/fabs(delx);
	  dely = dely/fabs(dely);
	  delz = delz/fabs(delz);
	  fpair=forceexp*factor_lj;
	  
	  f[i][0] += delx*fpair;
	  f[i][1] += dely*fpair;
	  f[i][2] += delz*fpair;
	  if(newton_pair || j < nlocal)
	  {
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;
	  }

	if(eflag)
	  {
	    evdwl = epsilon_exp[itype][jtype]*(exp_sl[itype][jtype]*exp_rcl[itype][jtype] -(exp_sl[itype][jtype]*exp(-r/lambda[itype][jtype])) - exp_sl[itype][jtype]*exp_rcl[itype][jtype]*(r-cut[itype][jtype])/lambda[itype][jtype]);
	  }
	if(evflag)
	  ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
	}
      else if (rsq < cut_exp2[itype][jtype] && lj_switch==1) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
      else;
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLEA::compute_inner()
{

}

/* ---------------------------------------------------------------------- */

void PairLEA::compute_middle()
{

}

/* ---------------------------------------------------------------------- */

void PairLEA::compute_outer(int eflag, int vflag)
{
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLEA::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut_exp2,n+1,n+1,"pair:cut_exp2");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_exp,n+1,n+1,"pair:cut_exp");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(epsilon_exp,n+1,n+1,"pair:epsilon_exp");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lambda,n+1,n+1,"pair:lambda");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(exp_sl,n+1,n+1,"pair:exp_sl");
  memory->create(exp_rcl,n+1,n+1,"pair:exp_rcl");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLEA::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  cut_lj= force->numeric(FLERR,arg[1]);
  min_global = force->numeric(FLERR,arg[2]);
  lj_switch = (unsigned int)force-> numeric(FLERR,arg[3]);

  if(lj_switch)
    printf("0==] In lj/exp/att  LJ repulsion switched on.\n");
  else
    printf("0==] In lj/exp/att LJ repulstion switched off.\n"); 

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) 
	  {
	    cut_exp[i][j] = cut_lj;
	    cut[i][j] = cut_global;
	  }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLEA::coeff(int narg, char **arg)
{
  if (narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double check = force->numeric(FLERR,arg[3]);
  double lambda_one = force -> numeric(FLERR,arg[4]);
  double cut_one;
  double cut_two;

if (narg>5)
     cut_one= force->numeric(FLERR,arg[5]);
  else
     cut_one = cut_lj;

  if(narg > 6)
    cut_two = force->numeric(FLERR,arg[6]);
  else
    cut_two = cut_global;


  if(cut_one > cut_two)
    error->all(FLERR, "Inner LJ cut off cannot be > to exp_att cut off");

  double epsilon_two,sigma_one = cut_one/pow(2,1.0/6.0);
  int count = 0;
  double pot=0,esl,er1l,er2l;
  esl=exp(cut_one/lambda_one);
  er1l=exp(-cut_one/lambda_one);
  er2l=exp(-cut_two/lambda_one);

  pot=esl*er2l -(esl*er1l) - esl*er2l*(cut_one-cut_two)/lambda_one;

  if(pot==0.0)
    epsilon_two=0.0;
  else
    epsilon_two=-min_global/pot;

  if(check>0.0)
    {
      min_global=check;
      printf("==========    lj/exp/att      =============");
      printf("Energy minimum: %f\n",min_global);
      
      printf("lj/exp/att | coeff: epsilon_one %f \n",epsilon_one);   
      printf("lj/exp/att | coeff: epsilon_two %f \n",epsilon_two);   
      printf("lj/exp/att | coeff: lambda_one %f \n",lambda_one);   
      printf("lj/exp/att | coeff: cut_one %f \n",cut_one);   
      printf("lj/exp/att | coeff: cut_two %f \n",cut_two);   
      printf("lj/exp/att | coeff: sigma_one %f \n\n",sigma_one);   
      printf("==================================");
    }
  else;

  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      epsilon_exp[i][j] = epsilon_two;
      sigma[i][j] = sigma_one;
      cut_exp[i][j] = cut_one;
      cut[i][j] = cut_two;
      lambda[i][j] = lambda_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLEA::init_style()
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
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLEA::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    
  }

  cut[j][i]=cut[i][j];
  
  lambda[j][i]=lambda[i][j];
  cut_exp[i][j] = mix_distance(cut_exp[i][i],cut_exp[j][j]);
  exp_rcl[i][j]=exp(-cut[i][j]/lambda[i][j]);
  exp_sl[i][j]=exp(cut_exp[i][j]/lambda[i][j]);

  cut_exp[j][i]=cut_exp[i][j];
  cut_exp2[i][j]=cut_exp[i][j]*cut_exp[i][j];
  cut_exp2[j][i]=cut_exp2[i][j];
  epsilon_exp[j][i]=epsilon_exp[i][j];

  exp_rcl[j][i]=exp_rcl[i][j];
  exp_sl[j][i]=exp_sl[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);


  offset[i][j] =  -epsilon_exp[i][j] + min_global;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

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

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLEA::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLEA::read_restart(FILE *fp)
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
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLEA::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLEA::read_restart_settings(FILE *fp)
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

void PairLEA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLEA::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLEA::single(int i, int j, int itype, int jtype, double rsq,double factor_coul, double factor_lj,double &fforce)
{
  double r2inv,r6inv,forceexp, forcelj,philj;
  if(rsq<cut_exp2[itype][jtype] && lj_switch==1)
    {
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;
      forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
      fforce = factor_lj*forcelj*r2inv;
      
      philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype])-offset[itype][jtype];
  return factor_lj*philj;
    }
  else if(rsq >= cut_exp2[itype][jtype] && rsq < cutsq[itype][jtype])
    {
      double r;
      r = sqrt(rsq);
      forceexp = (exp_sl[itype][jtype]*exp_rcl[itype][jtype] - 
		exp_sl[itype][jtype]*exp(-r/lambda[itype][jtype]))*epsilon_exp[itype][jtype]/lambda[itype][jtype];
      fforce = factor_lj*forceexp;
      philj = epsilon_exp[itype][jtype]*(exp_sl[itype][jtype]*exp_rcl[itype][jtype] 
       -(exp_sl[itype][jtype]*exp(-r/lambda[itype][jtype])) 
       - exp_sl[itype][jtype]*exp_rcl[itype][jtype]*(r-cut[itype][jtype])/lambda[itype][jtype]);
      return philj;
    }
}

/* ---------------------------------------------------------------------- */

void *PairLEA::extract(const char *str, int &dim)
{
  dim = 2;
 
  if (strcmp(str,"cut") == 0) return (void *) cut;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  if (strcmp(str,"epsilon_att") == 0) return (void *) epsilon_exp;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  
  return NULL;
}
