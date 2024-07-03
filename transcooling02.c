// gcc dipole01.c -o dipole01 dip_dynamics.c -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lpthread -lconfig

#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <pthread.h>
#include <math.h>
#include <libconfig.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "molecules.h"
#include "lasers.h"

#define NUM_THREADS 4

#define hbar 0.063507799 //units: amu um^2/us
#define c 2.9979246e+08  //units um/us = m/s

pthread_mutex_t readmutex, writemutex;
gsl_rng * r;

double t1;
struct laser las[5];
double Gamma, omegaz;

void *trajcalc (void *);

double lasacc(int lasnr, int ver, double y, double t){
  /* This calculates:
     acceleration = 1/m hbar k Gamma 0.5 s / (1 + s + (2 (delta - omegal * y * k /c)/Gamma)^2)
   */
  // Detuning: delta0 - omega * y * ver / c
  double detuning = las[lasnr].delta0+las[lasnr].delta1*t - las[lasnr].omega * y *las[lasnr].propagation[ver] / c;
  double fraction = las[lasnr].s / (1.0 + las[lasnr].s + ( 2.0 * detuning / Gamma) * ( 2.0 * detuning / Gamma)  );
  //double acc = las[lasnr].propagation[ver]*las[lasnr].constA * fraction + sqrt( las[lasnr].constB * fraction) * gsl_ran_gaussian(r, 1);
  double acc = las[lasnr].propagation[ver]*las[lasnr].constA * fraction;// + sqrt( las[lasnr].constB * fraction) * gsl_ran_gaussian(r, 1);
 
  return acc;  
}

int func (double t, const double y[], double f[], void *params){
  //(void)(t); /* avoid unused parameter warning */ //Forse serve per un laser impulsato
  //  struct mol cc = *(struct mol *)params; //non so se mi serve. Forse si', per conoscere eta' delle molecole
  f[0] = 0.0;  //a_x
  f[1] = 0.0;  //a_y
  f[2] = 0.0;  //a_z

  if ( (y[4]*y[4]+(y[5]-las[0].position[2])*(y[5]-las[0].position[2]) )<(las[0].w0*las[0].w0) ){
    f[0] += lasacc(0,0,y[0],t)+lasacc(1,0,y[0],t); //a_x
    f[1] += lasacc(0,1,y[1],t)+lasacc(1,1,y[1],t); //a_y
    f[2] += lasacc(0,2,y[2],t)+lasacc(1,2,y[2],t); //a_z
      }
  if( (y[3]*y[3]+(y[5]-las[2].position[2])*(y[5]-las[2].position[2]) )< (las[2].w0*las[2].w0) ){
    f[0] += lasacc(2,0,y[0],t)+lasacc(3,0,y[0],t); //a_x
    f[1] += lasacc(2,1,y[1],t)+lasacc(3,1,y[1],t); //a_y
    f[2] += lasacc(2,2,y[2],t)+lasacc(3,2,y[2],t); //a_z
    //    printf("%f\n",lasacc(2,0,y[1])); 
  }
  if( (y[3]*y[3]+y[4]*y[4])< (las[4].w0*las[4].w0) ){
    f[0] += lasacc(4,0,y[0],t); //a_x
    f[1] += lasacc(4,1,y[1],t); //a_y
    f[2] += lasacc(4,2,y[2],t); //a_z
    }
  f[3] = y[0]; // v_x
  f[4] = y[1]; // v_y
  f[5] = y[2]; // v_z

    return GSL_SUCCESS;
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


  las[0].propagation[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),0);
  las[0].propagation[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),1);
  las[0].propagation[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),2);
  las[0].position[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),0);
  las[0].position[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),1);
  las[0].position[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),2);
  las[0].w0 = config_setting_get_float(config_lookup(&cfg, "lasers.las1.w0"));
  las[0].lambda = config_setting_get_float(config_lookup(&cfg, "lasers.las1.lambda"));
  las[0].omega = config_setting_get_float(config_lookup(&cfg, "lasers.las1.omegal"));
  las[0].hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.las1.hbarkappa"));
  las[0].hbarkappasq = config_setting_get_float(config_lookup(&cfg, "lasers.las1.hbarkappasq"));
  las[0].s = config_setting_get_float(config_lookup(&cfg, "lasers.las1.s"));
  las[0].delta0 = config_setting_get_float(config_lookup(&cfg, "lasers.las1.delta0"));
  las[0].delta1 = config_setting_get_float(config_lookup(&cfg, "lasers.las1.delta1"));
  las[0].constA = las[0].hbarkappa * 0.5 * Gamma  / mass;
  las[0].constB = (Gamma*las[0].hbarkappa*las[0].hbarkappa)/(mass*mass*mass);

  las[1].propagation[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),0);
  las[1].propagation[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),1);
  las[1].propagation[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),2);
  las[1].position[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),0);
  las[1].position[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),1);
  las[1].position[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),2);
  las[1].w0 = config_setting_get_float(config_lookup(&cfg, "lasers.las2.w0"));
  las[1].lambda = config_setting_get_float(config_lookup(&cfg, "lasers.las2.lambda"));
  las[1].omega = config_setting_get_float(config_lookup(&cfg, "lasers.las2.omegal"));
  las[1].hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.las2.hbarkappa"));
  las[1].hbarkappasq = config_setting_get_float(config_lookup(&cfg, "lasers.las2.hbarkappasq"));
  las[1].s = config_setting_get_float(config_lookup(&cfg, "lasers.las2.s"));
  las[1].delta0 = config_setting_get_float(config_lookup(&cfg, "lasers.las2.delta0"));
  las[1].delta1 = config_setting_get_float(config_lookup(&cfg, "lasers.las2.delta1"));
  las[1].constA = las[1].hbarkappa * 0.5 * Gamma / mass;
  las[1].constB = (Gamma*las[1].hbarkappa*las[1].hbarkappa)/(mass*mass*mass);

  las[2].propagation[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),0);
  las[2].propagation[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),1);
  las[2].propagation[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),2);
  las[2].position[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),0);
  las[2].position[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),1);
  las[2].position[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),2);
  las[2].w0 = config_setting_get_float(config_lookup(&cfg, "lasers.las3.w0"));
  las[2].lambda = config_setting_get_float(config_lookup(&cfg, "lasers.las3.lambda"));
  las[2].omega = config_setting_get_float(config_lookup(&cfg, "lasers.las3.omegal"));
  las[2].hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.las3.hbarkappa"));
  las[2].hbarkappasq = config_setting_get_float(config_lookup(&cfg, "lasers.las3.hbarkappasq"));
  las[2].s = config_setting_get_float(config_lookup(&cfg, "lasers.las3.s"));
  las[2].delta0 = config_setting_get_float(config_lookup(&cfg, "lasers.las3.delta0"));
  las[2].delta1 = config_setting_get_float(config_lookup(&cfg, "lasers.las3.delta1"));
  las[2].constA = las[2].hbarkappa * 0.5 * Gamma  / mass;
  las[2].constB = (Gamma*las[2].hbarkappa*las[2].hbarkappa)/(mass*mass*mass);

  las[3].propagation[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),0);
  las[3].propagation[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),1);
  las[3].propagation[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),2);
  las[3].position[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),0);
  las[3].position[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),1);
  las[3].position[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),2);
  las[3].w0 = config_setting_get_float(config_lookup(&cfg, "lasers.las4.w0"));
  las[3].lambda = config_setting_get_float(config_lookup(&cfg, "lasers.las4.lambda"));
  las[3].omega = config_setting_get_float(config_lookup(&cfg, "lasers.las4.omegal"));
  las[3].hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.las4.hbarkappa"));
  las[3].hbarkappasq = config_setting_get_float(config_lookup(&cfg, "lasers.las4.hbarkappasq"));
  las[3].s = config_setting_get_float(config_lookup(&cfg, "lasers.las4.s"));
  las[3].delta0 = config_setting_get_float(config_lookup(&cfg, "lasers.las4.delta0"));
  las[3].delta1 = config_setting_get_float(config_lookup(&cfg, "lasers.las4.delta1"));
  las[3].constA = las[0].hbarkappa * 0.5 * Gamma  / mass;
  las[3].constB = (Gamma*las[3].hbarkappa*las[3].hbarkappa)/(mass*mass*mass);

  las[4].propagation[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),0);
  las[4].propagation[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),1);
  las[4].propagation[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),2);
  las[4].position[0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),0);
  las[4].position[1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),1);
  las[4].position[2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),2);
  las[4].w0 = config_setting_get_float(config_lookup(&cfg, "lasers.las5.w0"));
  las[4].lambda = config_setting_get_float(config_lookup(&cfg, "lasers.las5.lambda"));
  las[4].omega = config_setting_get_float(config_lookup(&cfg, "lasers.las5.omegal"));
  las[4].hbarkappa = config_setting_get_float(config_lookup(&cfg, "lasers.las5.hbarkappa"));
  las[4].hbarkappasq = config_setting_get_float(config_lookup(&cfg, "lasers.las5.hbarkappasq"));
  las[4].s = config_setting_get_float(config_lookup(&cfg, "lasers.las5.s"));
  las[4].delta0 = config_setting_get_float(config_lookup(&cfg, "lasers.las5.delta0"));
  las[4].delta1 = config_setting_get_float(config_lookup(&cfg, "lasers.las5.delta1"));
  las[4].constA = las[4].hbarkappa * 0.5 * Gamma  / mass;
  las[4].constB = (Gamma*las[4].hbarkappa*las[4].hbarkappa)/(mass*mass*mass);

  
  //  printf("%f\t%f\t%f",k[0],k[1],k[2]);
  //printf("%f\t%f\t%f",pos[0],pos[1],pos[2]);
  r = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(r, 0);

    
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

  //  gsl_odeiv2_system sys = {func, jac, 6, &cc};
  gsl_odeiv2_system sys = {func, NULL, 6, &cc};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-4, 1e-4, 0.0);
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
