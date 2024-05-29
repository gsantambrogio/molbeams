#ifndef LASERS_H
#define LASERS_H

struct laser
{
  double El0;        //V/micron
  double power;      //mW
  double w0;         //micron
  double omega;
  double nu;
  double lambda;    //micron
  double axis;
  double position[2];
  double x0;         //micron
  double y0;         //micron
  double z0;         //micron
  double tinit;      //micrsecond
  double PulseLength;//microsecond
}; //we could introduce a propagation direction... later

//283nm>>  (along the x direction)
//const struct laser _283nm = {.El0=0.0, .power=0.0, .wzero=600, .lambda0=0.283, .x0=0.0, .y0=0.0, .z0=305000.0, .tinit=924.0, .PulseLength=0.0};

#endif
