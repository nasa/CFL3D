
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
