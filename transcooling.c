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

#define NUM_THREADS 4

#define hbar 0.063507799 //units: amu um^2/us
#define c 2.9979246e+08  //units um/us = m/s

pthread_mutex_t readmutex, writemutex;
gsl_rng * r;

double t1;
double allconst[5], Gamma, omegaz, laser[5][13];

void *trajcalc (void *);

double lasacc(int las, int ver, double y, double t){
  /* This calculates:
     acceleration = 1/m hbar k Gamma 0.5 s / (1 + s + (2 (delta - omegal * y * k /c)/Gamma)^2)
   */
  //This was with the random broadening. Too slow
  //    double acc = laser[las][ver]*allconst[las] / (1.0 + laser[las][10] + ( ( 2.0 * (laser[las][12]*(gsl_rng_uniform(r)-0.5)+laser[las][11] - laser[las][8] * y*laser[las][ver]/c))/ Gamma) * ( ( 2.0 * (laser[las][12]*(gsl_rng_uniform(r)-0.5)+laser[las][11] - laser[las][8] * y*laser[las][ver]/c))/ Gamma)   );
  //This is with static detuning:
  //      double acc = laser[las][ver]*allconst[las] / (1.0 + laser[las][10] + ( ( 2.0 * (laser[las][11] - laser[las][8] * y*laser[las][ver]/c))/ Gamma) * ( ( 2.0 * (laser[las][11] - laser[las][8] * y*laser[las][ver]/c))/ Gamma)   );
  double acc = laser[las][ver]*allconst[las] / (1.0 + laser[las][10] + ( ( 2.0 * (laser[las][11]+laser[las][12]*t - laser[las][8] * y*laser[las][ver]/c))/ Gamma) * ( ( 2.0 * (laser[las][11]+laser[las][12]*t - laser[las][8] * y*laser[las][ver]/c))/ Gamma)   );
      

  return acc;  
}

double lasder(int las, int ver, double y){
  /* This calculates:
- k[0]*allconst *8.0 * omegal *(delta + y[0]*k[0] * omegal/c)/( c * Gamma * Gamma * (1.0 + s + 4.0*(delta + y[0]*k[0]*omegal/c)*(delta + y[0]*k[0]*omegal/c)/(Gamma*Gamma))*(1.0 + s + 4.0*(delta + y[0]*k[0]*omegal/c)*(delta + y[0]*k[0]*omegal/c)/(Gamma*Gamma)) )
  */
  double der = - laser[las][ver]*allconst[las] *8.0 * laser[las][8] *(laser[las][11] - y*laser[las][ver] * laser[las][8]/c)/( c * Gamma * Gamma * (1.0 + laser[las][10] + 4.0*(laser[las][11] - y*laser[las][ver]*laser[las][8]/c)*(laser[las][11] - y*laser[las][ver]*laser[las][8]/c)/(Gamma*Gamma))*(1.0 + laser[las][10] + 4.0*(laser[las][11] - y*laser[las][ver]*laser[las][8]/c)*(laser[las][11] - y*laser[las][ver]*laser[las][8]/c)/(Gamma*Gamma)) );
    return der;
}

int func (double t, const double y[], double f[], void *params){
  //(void)(t); /* avoid unused parameter warning */ //Forse serve per un laser impulsato
  //  struct mol cc = *(struct mol *)params; //non so se mi serve. Forse si', per conoscere eta' delle molecole
  f[0] = 0.0;  //a_x
  f[1] = 0.0;  //a_y
  f[2] = 0.0;  //a_z

  if ( (y[4]*y[4]+(y[5]-laser[0][5])*(y[5]-laser[0][5]) )<(laser[0][6]*laser[0][6]) ){
    f[0] += lasacc(0,0,y[0],t)+lasacc(1,0,y[0],t); //a_x
    f[1] += lasacc(0,1,y[1],t)+lasacc(1,1,y[1],t); //a_y
    f[2] += lasacc(0,2,y[2],t)+lasacc(1,2,y[2],t); //a_z
      }
  if( (y[3]*y[3]+(y[5]-laser[2][5])*(y[5]-laser[2][5]) )< (laser[2][6]*laser[2][6]) ){
    f[0] += lasacc(2,0,y[0],t)+lasacc(3,0,y[0],t); //a_x
    f[1] += lasacc(2,1,y[1],t)+lasacc(3,1,y[1],t); //a_y
    f[2] += lasacc(2,2,y[2],t)+lasacc(3,2,y[2],t); //a_z
    //    printf("%f\n",lasacc(2,0,y[1])); 
  }
  if( (y[3]*y[3]+y[4]*y[4])< (laser[4][6]*laser[4][6]) ){
    f[0] += lasacc(4,0,y[0],t); //a_x
    f[1] += lasacc(4,1,y[1],t); //a_y
    f[2] += lasacc(4,2,y[2],t); //a_z
    }
  f[3] = y[0]; // v_x
  f[4] = y[1]; // v_y
  f[5] = y[2]; // v_z

    return GSL_SUCCESS;
}
/*
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
  (void)(t);
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 6, 6);
  gsl_matrix *m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, lasder(0,0,y[0])+lasder(1,0,y[0]));
  gsl_matrix_set (m, 0, 1, 0.0);
  gsl_matrix_set (m, 0, 2, 0.0);
  gsl_matrix_set (m, 0, 3, 0.0);
  gsl_matrix_set (m, 0, 4, 0.0);
  gsl_matrix_set (m, 0, 5, 0.0);
  gsl_matrix_set (m, 1, 0, 0.0);
  gsl_matrix_set (m, 1, 1, lasder(0,1,y[1])+lasder(1,1,y[1]));
  gsl_matrix_set (m, 1, 2, 0.0);
  gsl_matrix_set (m, 1, 3, 0.0);
  gsl_matrix_set (m, 1, 4, 0.0);
  gsl_matrix_set (m, 1, 5, 0.0);
  gsl_matrix_set (m, 2, 0, 0.0);
  gsl_matrix_set (m, 2, 1, 0.0);
  gsl_matrix_set (m, 2, 2, lasder(0,2,y[2])+lasder(1,2,y[2]));
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
  return GSL_SUCCESS;
}
*/
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


  laser[0][0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),0);
  laser[0][1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),1);
  laser[0][2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.propagation"),2);
  laser[0][3] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),0);
  laser[0][4] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),1);
  laser[0][5] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las1.position"),2);
  laser[0][6] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.w0"));
  laser[0][7] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.lambda"));
  laser[0][8] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.omegal"));
  laser[0][9] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.hbarkappa"));
  laser[0][10] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.s"));
  laser[0][11] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.delta0"));
  laser[4][12] = config_setting_get_float(config_lookup(&cfg, "lasers.las1.delta1"));
  allconst[0] = laser[0][9] * 0.5 * Gamma *  laser[0][10] / mass;
 
  laser[1][0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),0);
  laser[1][1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),1);
  laser[1][2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.propagation"),2);
  laser[1][3] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),0);
  laser[1][4] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),1);
  laser[1][5] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las2.position"),2);
  laser[1][6] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.w0"));
  laser[1][7] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.lambda"));
  laser[1][8] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.omegal"));
  laser[1][9] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.hbarkappa"));
  laser[1][10] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.s"));
  laser[1][11] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.delta0"));
  laser[4][12] = config_setting_get_float(config_lookup(&cfg, "lasers.las2.delta1"));
  allconst[1] = laser[1][9] * 0.5 * Gamma *  laser[1][10] / mass;

  laser[2][0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),0);
  laser[2][1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),1);
  laser[2][2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.propagation"),2);
  laser[2][3] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),0);
  laser[2][4] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),1);
  laser[2][5] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las3.position"),2);
  laser[2][6] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.w0"));
  laser[2][7] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.lambda"));
  laser[2][8] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.omegal"));
  laser[2][9] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.hbarkappa"));
  laser[2][10] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.s"));
  laser[2][11] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.delta0"));
  laser[4][12] = config_setting_get_float(config_lookup(&cfg, "lasers.las3.delta1"));
  allconst[2] = laser[2][9] * 0.5 * Gamma *  laser[2][10] / mass;
 
  laser[3][0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),0);
  laser[3][1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),1);
  laser[3][2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.propagation"),2);
  laser[3][3] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),0);
  laser[3][4] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),1);
  laser[3][5] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las4.position"),2);
  laser[3][6] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.w0"));
  laser[3][7] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.lambda"));
  laser[3][8] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.omegal"));
  laser[3][9] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.hbarkappa"));
  laser[3][10] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.s"));
  laser[3][11] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.delta0"));
  laser[4][12] = config_setting_get_float(config_lookup(&cfg, "lasers.las4.delta1"));
  allconst[3] = laser[3][9] * 0.5 * Gamma *  laser[3][10] / mass;

  laser[4][0] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),0);
  laser[4][1] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),1);
  laser[4][2] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.propagation"),2);
  laser[4][3] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),0);
  laser[4][4] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),1);
  laser[4][5] = config_setting_get_float_elem(config_lookup(&cfg, "lasers.las5.position"),2);
  laser[4][6] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.w0"));
  laser[4][7] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.lambda"));
  laser[4][8] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.omegal"));
  laser[4][9] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.hbarkappa"));
  laser[4][10] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.s"));
  laser[4][11] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.delta0"));
  laser[4][12] = config_setting_get_float(config_lookup(&cfg, "lasers.las5.delta1"));
  allconst[4] = laser[4][9] * 0.5 * Gamma *  laser[4][10] / mass;
  
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
