#ifndef DIP_DYNAMICS_H
#define DIP_DYNAMICS_H

//int jac (double t, const double y[], double *dfdy,double dfdt[], void *cc);
int func (double t, const double y[], double f[], void *params);

extern struct laser IPG;

extern struct electrode StEl;

#endif /* DIP_DYNAMICS_H */
