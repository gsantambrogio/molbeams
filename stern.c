#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "stern.h"

static const double halflambda = 2.81; //halflambda/28 amu in (um/us)^2
static const double mu = 49.3; //half mu/28 amu in (um/us)^2/(V/um)

int jac (double t, const double y[], double *dfdy,double dfdt[], void *params){
  //temp:
  (void)(t); /* avoid unused parameter warning */
  //  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 6, 6);
  gsl_matrix * m = &dfdy_mat.matrix;
  //The Jacobian in this system is the identity
  gsl_matrix_set_identity (m);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;
  dfdt[4] = 0.0;
  dfdt[5] = 0.0;

  return GSL_SUCCESS;
}

int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */
  //  struct mol cc = *(struct mol *)params;
  double El, dEndEl;

//Force_x is - D[Energy,x] = dEnergy / dElectric * dElectric/dx
  El= quello che mi manda Giacomo (y[0],y[1],y[2]);
  //dEndEl(x,y,z)= El(x,y,z)*mu*mu / hypot(halflambda,(El(x,y,z)*mu));
  
  dEndEl = El*mu*mu / hypot(halflambda,(El*mu)); 
  
//y[] = x,y,z,vx,vy,vz
//f[] = vx,vy,vz,ax,ay,az
  f[0] = y[3];
  f[1] = y[4];
  f[2] = y[5];
  f[3] = -dEndEl * dEldx(y[0],y[1],y[2]); 
  f[4] = -dEndEl * dEldy(y[0],y[1],y[2]); 
  f[5] = -dEndEl * dEldz(y[0],y[1],y[2]); 

  
  return GSL_SUCCESS;
}
