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
  double x0;         //micron
  double y0;         //micron
  double z0;         //micron
  double tinit;      //micrsecond
  double PulseLength;//microsecond
  double propagation[3]; //norm must be 1
  double position[3]; // position on the intersection with axis z
  double hbarkappa; // in amu um/us
  double hbarkappasq; // in amu um/us
  double s;
  double delta0; // in omega = 2 pi nu; in MHz
  double delta1; //chirp
  double constA;
  double constB;
}; 

#endif
