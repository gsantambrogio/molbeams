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
#include "electrodes.h"

#define NUM_THREADS 1

#define hbar 0.063507799 //units: amu um^2/us
#define c 2.9979246e+08

pthread_mutex_t readmutex, writemutex;

//struct laser IPG;
//struct electrode StEl;

double t1;
double allconst, s, Gamma;

void *trajcalc (void *);

int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */ //Forse serve per un laser impulsato
  //  struct mol cc = *(struct mol *)params; //non so se mi serve. Forse si', per conoscere eta' delle molecole
    
  f[0] = allconst / (1.0 + s + (2.0 * (1.0 + y[0]/c) * omegal - omegaz)*(2.0 * (1.0 + y[0]/c) * omegal - omegaz)/(Gamma*Gamma) );   //a_x
  f[1] = 0.0;  //a_y
  f[2] = 0.0;  //a_z
  f[3] = y[0]; // v_x
  f[4] = y[1]; // v_y
  f[5] = y[2]; // v_z
  //    printf("ok\n");

    return GSL_SUCCESS;
}


int main (int argc, char *argv[]){
  if (argc != 2) {
    fprintf(stderr, "Usage: dipole01 <evolution t in us> <Mols> \n");
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
  s = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.s"));
  double hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.nm531.hbarkappa"));
  allconst = hbarkappa * Gamma * M_PI * s;
  
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

  gsl_odeiv2_system sys = {func, NULL, 6, &cc};
  //    gsl_odeiv2_system sys = {func, jac, 6, &cc};
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
