/*
 *  ---------------------------------------------------------------------------
 *  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
 *  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
 *  with point-matched, patched, overset, or embedded connectivities. Both
 *  multigrid and mesh sequencing are available in time-accurate or
 *  steady-state modes.
 *
 *  Copyright 2001 United States Government as represented by the Administrator
 *  of the National Aeronautics and Space Administration. All Rights Reserved.
 * 
 *  The CFL3D platform is licensed under the Apache License, Version 2.0 
 *  (the "License"); you may not use this file except in compliance with the 
 *  License. You may obtain a copy of the License at 
 *  http://www.apache.org/licenses/LICENSE-2.0. 
 * 
 *  Unless required by applicable law or agreed to in writing, software 
 *  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
 *  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
 *  License for the specific language governing permissions and limitations 
 *  under the License.
 *  ---------------------------------------------------------------------------
 */

/*  $Id$ */

#include <math.h>


void j0_(double *x, double *ans)
{
  *ans = j0(*x);
}

void j1_(double *x, double *ans)
{
  *ans = j1(*x);
}

void jn_(int *n, double *x, double *ans)
{
  *ans = jn(*n,*x);
}

void y0_(double *x, double *ans)
{
  *ans = y0(*x);
}

void y1_(double *x, double *ans)
{
  *ans = y1(*x);
}

void yn_(int *n, double *x, double *ans)
{
  *ans = yn(*n,*x);
}

void j0sp_(float *x, float *ans)
{
  *ans = j0((double)*x);
}

void j1sp_(float *x, float *ans)
{
  *ans = j1(*x);
}

void jnsp_(int *n, float *x, float *ans)
{
  *ans = jn(*n,*x);
}

void y0sp_(float *x, float *ans)
{
  *ans = y0(*x);
}

void y1sp_(float *x, float *ans)
{
  *ans = y1(*x);
}

void ynsp_(int *n, float *x, float *ans)
{
  *ans = yn(*n,*x);
}
