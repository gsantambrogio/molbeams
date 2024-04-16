// gcc -o skimmer skimmer.c -lconfig
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <libconfig.h>
#include "molecules.h"


struct orifice{
  double x0;
  double y0;
  double z0;
  double radius;
};



bool inOrifice (struct mol *cc,struct orifice *or){
  //aperture has zero thickness, i.e. zw=0.0
  double t,x,y,z; //position of molecule at the aperture
  t = cc->tinit + (or->z0 - cc->pos[2])/cc->vel[2]; //time of arrival at the aperture
  x = cc->pos[0]+ (t - cc->tinit)*cc->vel[0];
  y = cc->pos[1]+ (t - cc->tinit)*cc->vel[1];
  z = or->z0;
  double  rSq = or->radius * or->radius;
  //  printf("%f %f %f %f %f\n",t,x,y,(y-or->y0)*(y-or->y0) + (x-or->x0)*(x-or->x0),or->radius);
  return ( (y-or->y0)*(y-or->y0) + (x-or->x0)*(x-or->x0) < rSq );
}

int main (int argc, char *argv[]){
  if (argc != 1) {
    fprintf(stderr, "Usage: to be written... \n");
    return 1;
  }
  struct mol cc;
  char buf[0x1000];
  struct orifice or;
  config_t cfg;

  config_init(&cfg);
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }
  or.x0=config_setting_get_float(config_lookup(&cfg, "repeller.x0"));
  or.y0=config_setting_get_float(config_lookup(&cfg, "repeller.y0"));
  or.z0=config_setting_get_float(config_lookup(&cfg, "repeller.z0"));
  or.radius=config_setting_get_float(config_lookup(&cfg, "repeller.radius"));

  for (;;){
    if (fgets(buf, sizeof(buf), stdin) == NULL){
      break;
    }
    if (buf[0] == '#') {
      fputs(buf,stdout); //puts() would produce another \n at the end of the buffer
      continue;
    }
    sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf",
	   &cc.molname, &cc.tinit, &cc.pos[0], &cc.pos[1], &cc.pos[2], &cc.vel[0], &cc.vel[1], &cc.vel[2]);
    if(inOrifice(&cc,&or)) fputs(buf,stdout);
  }
      
  return 0;
}

 
