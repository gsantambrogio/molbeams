// gcc dipole01.c -o dipole01 dip_dynamics.c -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lpthread -lconfig

#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <pthread.h>
#include <math.h>
#include <libconfig.h>
#include "dip_dynamics.h"
#include "molecules.h"
#include "lasers.h"

#ifdef NAN
/* NAN is supported */
#endif
#define NUM_THREADS 1

#define hbar 0.063507799 //units: amu um^2/us
#define c 2.9979246e+08  //units um/us = m/s

pthread_mutex_t readmutex, writemutex;

//struct laser IPG;
//struct electrode StEl;

double t1;
double allconst, s, Gamma, delta, omegal, omegaz, w0, pos[3], k[3];
int axis=1;
void *trajcalc (void *);

int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */ //Forse serve per un laser impulsato
  //  struct mol cc = *(struct mol *)params; //non so se mi serve. Forse si', per conoscere eta' delle molecole
  switch(axis){
  case 1:
      if ((y[4]*y[4]+y[5]*y[5])<(w0*w0)){
	f[0] = k[0]*allconst / (1.0 + s + ( ( 2.0 * (delta + omegal * y[0]*k[0]/c))/ Gamma) * ( ( 2.0 * (delta + omegal * y[0]*k[0]/c))/ Gamma)   );
        f[1] = k[1]*allconst / (1.0 + s + ( ( 2.0 * (delta + omegal * y[1]*k[1]/c))/ Gamma) * ( ( 2.0 * (delta + omegal * y[1]*k[1]/c))/ Gamma)   );  //a_y
	f[2] = k[2]*allconst / (1.0 + s + ( ( 2.0 * (delta + omegal * y[2]*k[2]/c))/ Gamma) * ( ( 2.0 * (delta + omegal * y[2]*k[2]/c))/ Gamma)   );  //a_z
      }
      else{
	f[0] = 0.0;  //a_x
	f[1] = 0.0;  //a_y
	f[2] = 0.0;  //a_z
      }
      break;
    case 2:
      if ((y[3]*y[3]+y[5]*y[5])<(w0*w0)){
	f[0] = 0.0;  //a_x
        f[1] = allconst / (1.0 + s + ( ( 2.0 * (delta + omegal * y[0]/c))/ Gamma) * ( ( 2.0 * (delta + omegal * y[0]/c))/ Gamma)   );
	f[2] = 0.0;  //a_z
      }
      else{
	f[0] = 0.0;  //a_x
	f[1] = 0.0;  //a_y
	f[2] = 0.0;  //a_z
      }
      break;
    case 3:
      if ((y[4]*y[4]+y[3]*y[3])<(w0*w0)){
	f[0] = 0.0;  //a_x
        f[1] = 0.0;  //a_y
	f[2] = allconst / (1.0 + s + ( ( 2.0 * (delta + omegal * y[0]/c))/ Gamma) * ( ( 2.0 * (delta + omegal * y[0]/c))/ Gamma)   );
      }
      else{
	f[0] = 0.0;  //a_x
	f[1] = 0.0;  //a_y
	f[2] = 0.0;  //a_z
      }
      break;
  }

  f[3] = y[0]; // v_x
  f[4] = y[1]; // v_y
  f[5] = y[2]; // v_z
  //    printf("ok\n");

    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
  (void)(t);
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 6, 6);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, - k[0]*allconst *8.0 * omegal *(delta + y[0]*k[0] * omegal/c)/( c * Gamma * Gamma * (1.0 + s + 4.0*(delta + y[0]*k[0]*omegal/c)*(delta + y[0]*k[0]*omegal/c)/(Gamma*Gamma))*(1.0 + s + 4.0*(delta + y[0]*k[0]*omegal/c)*(delta + y[0]*k[0]*omegal/c)/(Gamma*Gamma)) )  );
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, - k[1]*allconst *8.0 * omegal *(delta + y[1]*k[1] * omegal/c)/( c * Gamma * Gamma * (1.0 + s + 4.0*(delta + y[1]*k[1]*omegal/c)*(delta + y[1]*k[1]*omegal/c)/(Gamma*Gamma))*(1.0 + s + 4.0*(delta + y[1]*k[1]*omegal/c)*(delta + y[1]*k[1]*omegal/c)/(Gamma*Gamma)) )  );
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 0.0);
  gsl_matrix_set (m, 1, 5, 0.0);
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, - k[2]*allconst *8.0 * omegal *(delta + y[2]*k[2] * omegal/c)/( c * Gamma * Gamma * (1.0 + s + 4.0*(delta + y[2]*k[2]*omegal/c)*(delta + y[2]*k[2]*omegal/c)/(Gamma*Gamma))*(1.0 + s + 4.0*(delta + y[2]*k[2]*omegal/c)*(delta + y[2]*k[2]*omegal/c)/(Gamma*Gamma)) )  );
  gsl_matrix_set (m, 2, 3, 0.0);
  gsl_matrix_set (m, 2, 4, 0.0);
  gsl_matrix_set (m, 2, 5, 0.0);
  gsl_matrix_set (m, 3, 0, 0.0);
  gsl_matrix_set (m, 3, 1, 0.0);
  gsl_matrix_set (m, 3, 2, 0.0);
  gsl_matrix_set (m, 3, 3, 0.0);
  gsl_matrix_set (m, 3, 4, 0.0);
  gsl_matrix_set (m, 3, 5, 0.0);
  gsl_matrix_set (m, 4, 0, 0.0);
  gsl_matrix_set (m, 4, 1, 0.0);
  gsl_matrix_set (m, 4, 2, 0.0);
  gsl_matrix_set (m, 4, 3, 0.0);
  gsl_matrix_set (m, 4, 4, 0.0);
  gsl_matrix_set (m, 4, 5, 0.0);
  gsl_matrix_set (m, 5, 0, 0.0);
  gsl_matrix_set (m, 5, 1, 0.0);
  gsl_matrix_set (m, 5, 2, 0.0);
  gsl_matrix_set (m, 5, 3, 0.0);
  gsl_matrix_set (m, 5, 4, 0.0);
  gsl_matrix_set (m, 5, 5, 0.0);  
  dfdt[0]= 0.0;
  dfdt[1]= 0.0;
  dfdt[2]= 0.0;
  dfdt[3]= 0.0;
  dfdt[4]= 0.0;
  dfdt[5]= 0.0;
}

int main (int argc, char *argv[]){
  if (argc != 2) {
    fprintf(stderr, "Usage: transcooling <evolution t in us>  \n");
    return 1;
  }

  t1 = atof(argv[1]);

  config_t cfg;
  
  config_init(&cfg);
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }

  Gamma = config_setting_get_float(config_lookup(&cfg, "CaF.Gamma"));
  omegaz = config_setting_get_float(config_lookup(&cfg, "CaF.omegaz"));
  double mass =  config_setting_get_float(config_lookup(&cfg, "CaF.mass"));
  s = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.s"));
  omegal = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.omegal"));
  delta = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.delta"));
  w0 =config_setting_get_float(config_lookup(&cfg, "lasers.nm531.w0"));
  k[0] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.propagation"),0);
  k[1] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.propagation"),1);
  k[2] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.propagation"),2);
  pos[0] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.position"),0);
  pos[1] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.position"),1);
  pos[2] =config_setting_get_float_elem(config_lookup(&cfg, "lasers.nm531.position"),2);
  double hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.hbarkappa"));
  allconst = hbarkappa * 0.5 * Gamma *  s / mass;
  //  printf("%f\t%f\t%f",k[0],k[1],k[2]);
  //printf("%f\t%f\t%f",pos[0],pos[1],pos[2]);

  pthread_t thread[NUM_THREADS];
  pthread_attr_t attr;
  void *joinstatus;
  int n;
  
  pthread_attr_init (&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&readmutex, NULL);
  pthread_mutex_init(&writemutex, NULL);
  for (n = 0; n < NUM_THREADS; ++n) {
    pthread_create (&thread[n], &attr, trajcalc, NULL);
  }
  pthread_attr_destroy(&attr);
  for (n = 0; n < NUM_THREADS; ++n) {
    pthread_join (thread[n], &joinstatus);
  }
  pthread_mutex_destroy(&readmutex);
  pthread_mutex_destroy(&writemutex);
  
  return 0;
}

void *trajcalc (void *arg){
  struct mol cc;
  char buf[0x1000];

  gsl_odeiv2_system sys = {func, jac, 6, &cc};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
  //gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);

  double y[6];

  for (;;){
    pthread_mutex_lock (&readmutex);
    if (fgets(buf, sizeof(buf), stdin) == NULL){
      pthread_mutex_unlock (&readmutex);
      break;
    }
    pthread_mutex_unlock (&readmutex);	 
    if (buf[0] == '#') {
      fputs(buf,stdout); //puts() would produce another \n at the end of the buffer
      continue;
    }

    sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf",
	   &cc.molname, &cc.tinit, &cc.pos[0], &cc.pos[1], &cc.pos[2], &cc.vel[0], &cc.vel[1], &cc.vel[2]);
    y[0] = cc.vel[0];
    y[1] = cc.vel[1];
    y[2] = cc.vel[2];
    y[3] = cc.pos[0];
    y[4] = cc.pos[1];
    y[5] = cc.pos[2];

    

    double t = cc.tinit;
    double ti=t+t1;
    //    printf("%5f\n\n",t);
    
    int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
    if (status != GSL_SUCCESS){
      printf ("error, return value=%d\n", status);
      break;
    }
    //    printf("%5f\n\n",t);

    pthread_mutex_lock (&writemutex);
    printf("%d\t%5f\t%5f\t%5f\t%5f\t%5f\t%5f\t%5f\n",
	   cc.molname,
	   t,
	   y[3],
	   y[4],
	   y[5],
	   y[0],
	   y[1],
	   y[2]
	   );                        
    pthread_mutex_unlock (&writemutex);
  }

  pthread_exit(NULL);  
  gsl_odeiv2_driver_free (d);
}
