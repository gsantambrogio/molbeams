// gcc -o apertures apertures.c -lconfig
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <libconfig.h>
#include "molecules.h"


struct aperture{
  double x0;
  double y0;
  double z0;
  double xw; // half opening, i.e. the width is 2*xw
  double yw; // half opening
};



bool inAperture (struct mol *cc,struct aperture *ap){
  //aperture has zero thickness, i.e. zw=0.0
  double t,x,y,z; //position of molecule at the aperture
  t = cc->tinit + (cc->pos[2] - ap->z0)/cc->vel[2]; //time of arrival at the aperture
  x = cc->pos[0]+ (t - cc->tinit)*cc->vel[0];
  y = cc->pos[1]+ (t - cc->tinit)*cc->vel[1];
  z = ap->z0;

  return( (fabs(x - ap->x0) < ap->xw) && (fabs(y - ap->y0) < ap->yw) );
}

int main (int argc, char *argv[]){
  if (argc != 1) {
    fprintf(stderr, "Usage: to be written... \n");
    return 1;
  }
  struct mol cc;
  char buf[0x1000];
  struct aperture ap1;
  struct aperture ap2;
  config_t cfg;

  config_init(&cfg);
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }
  ap1.x0=config_setting_get_float(config_lookup(&cfg, "apertures.ap1.x0"));
  ap1.y0=config_setting_get_float(config_lookup(&cfg, "apertures.ap1.y0"));
  ap1.z0=config_setting_get_float(config_lookup(&cfg, "apertures.ap1.z0"));
  ap1.xw=config_setting_get_float(config_lookup(&cfg, "apertures.ap1.xw"));
  ap1.yw=config_setting_get_float(config_lookup(&cfg, "apertures.ap1.yw"));
  ap2.x0=config_setting_get_float(config_lookup(&cfg, "apertures.ap2.x0"));
  ap2.y0=config_setting_get_float(config_lookup(&cfg, "apertures.ap2.y0"));
  ap2.z0=config_setting_get_float(config_lookup(&cfg, "apertures.ap2.z0"));
  ap2.xw=config_setting_get_float(config_lookup(&cfg, "apertures.ap2.xw"));
  ap2.yw=config_setting_get_float(config_lookup(&cfg, "apertures.ap2.yw"));

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
    if(inAperture(&cc,&ap1) && inAperture(&cc,&ap2)) fputs(buf,stdout);
  }
      
  return 0;
}

 
