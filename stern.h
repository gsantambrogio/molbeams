#ifndef STERN_H
#define STERN_H

int jac (double t, const double y[], double *dfdy,double dfdt[], void *cc);
int func (double t, const double y[], double f[], void *params);

#endif /* STERN_H */
