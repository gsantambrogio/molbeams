#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "dip_dynamics.h"
//#include "params.h"
#include "molecules.h"
#include "lasers.h"
#include "electrodes.h"

#define a3PiHalfLambda 2.80747 //0.5 h 394 MHz / 28 amu  in um^2/us^2

int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */ //Forse serve per un laser impulsato
  struct mol cc = *(struct mol *)params; //non so se mi serve. Forse si', per conoscere eta' delle molecole
  double relx, rely, relz, w, ElDC, dEndEl;
  

    relx = y[3]-IPG.x0;
    rely = y[4]-IPG.y0;
    relz = y[5]-IPG.z0;

    w = IPG.w0*sqrt(1.0 +  (((y[4]-IPG.y0)*(y[4]-IPG.y0)*IPG.lambda0*IPG.lambda0)  /   (M_PI*M_PI*IPG.w0*IPG.w0*IPG.w0*IPG.w0))   );    
    dEndEl= gsPol*IPG.El0*(IPG.w0/w)*exp(-  (  relx*relx+relz*relz  ) / (w*w) ); 
    
    f[0] = dEndEl*(-2.)*IPG.El0*IPG.w0*(relx)*exp(-(relx*relx+relz*relz)/(w*w))/(w*w*w); //acceleration_x
    f[1] = dEndEl*2.*IPG.El0*(relx*relx+relz*relz)*rely*IPG.lambda0*IPG.lambda0*exp(-(relx*relx+relz*relz)/(w*w))/(M_PI*M_PI*IPG.w0*w*w*w*w*w)-IPG.El0*rely*IPG.lambda0*IPG.lambda0*exp(-(relx*relx+relz*relz)/(w*w))/(M_PI*M_PI*IPG.w0*w*w*w); //acceleration_y
    f[2] = dEndEl*(-2.)*IPG.El0*IPG.w0*relz*exp(-(relx*relx+relz*relz)/(w*w))/(w*w*w); //acceleration_z
    f[3] = y[0]; // v_x
    f[4] = y[1]; // v_y
    f[5] = y[2]; // v_z
    printf("ok\n");

    return GSL_SUCCESS;
}

  
//  fields(&t, &cc);
  /*
    f[0]=a_x;  dEn(El)/dx = dEn/dEl dEl/dx 
    f[1]=a_y;  etc
    f[2]=a_z;
    f[3]=v_x;
    f[4]=v_y;
    f[5]=v_z; */
