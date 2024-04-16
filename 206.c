// gcc -o 206 206.c -lconfig
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <libconfig.h>
#include "lasers.h"
#include "molecules.h"

struct laser _206nm;

bool in206 (struct mol *cc){
  //cylindrical laser propagating along the y direction
  double x,y,z; //position of molecule when the laser shoots
  x = cc->pos[0]+(_206nm.tinit - cc->tinit)*cc->vel[0];
  //y = cc->pos[1]+(_206nm.tinit - cc->tinit)*cc->vel[1];
  z = cc->pos[2]+(_206nm.tinit - cc->tinit)*cc->vel[2];
  double  wsq = _206nm.w0*_206nm.w0;
  
  return( ((x - _206nm.x0)*(x - _206nm.x0))+((z - _206nm.z0)*(z - _206nm.z0))<wsq );
}

int main (int argc, char *argv[]){
  if (argc != 1) {
    fprintf(stderr, "Usage: to be written... \n");
    return 1;
  }
  struct mol cc;
  //  struct laser _206nm;
  char buf[0x1000];
  config_t cfg;

  config_init(&cfg);
  if(! config_read_file(&cfg, "config.cfg"))
    {
      fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
	      config_error_line(&cfg), config_error_text(&cfg));
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }
  _206nm.El0=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.El0"));
  _206nm.w0=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.w0"));
  _206nm.lambda=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.lambda0"));
  _206nm.x0=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.x0"));
  _206nm.y0=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.y0"));
  _206nm.z0=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.z0"));
  _206nm.tinit=config_setting_get_float(config_lookup(&cfg, "lasers.nm206.tinit"));

  for (;;){
    if (fgets(buf, sizeof(buf), stdin) == NULL){
      break;
    }
    if (buf[0] == '#') {
      fputs(buf,stdout); //puts() would produce another \n at the end of the buffer
      continue;
    }
    sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf %d",
	   &cc.molname, &cc.tinit, &cc.pos[0], &cc.pos[1], &cc.pos[2], &cc.vel[0], &cc.vel[1], &cc.vel[2], &cc.state);
    //    if(in206(&cc)) fputs(buf,stdout);
    if(in206(&cc)) printf("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\n", cc.molname, cc.tinit, cc.pos[0], cc.pos[1], cc.pos[2], cc.vel[0], cc.vel[1], cc.vel[2],1);
  }
      
  return 0;
}

 
