#ifndef OPTI_H
#define OPTI_H 1

int check_t1();
void update_pa();
double area_fct(Dual *pcell);
double area_fctp(Dual *pcell);
double d2E_dA2(Dual *pcell);
double rho_fct(Dual *pcell);
double rho_fctpp(Dual *pcell);
double rho_fctpa(Dual *pcell);
double rhol_fct(Bond *pbond);
double rhol_fctpl(Bond *pbond);
double d2E_dL2(Bond *pbond);
double extern_pot_fct(Dual *pcell);
vector3 extern_pot_fctp(Node *pvert);
double energy(const gsl_vector *x,void *params);
void denergy(const gsl_vector *x, void *params, gsl_vector *df);
void enerdener(const gsl_vector *x, void *params, double *f, gsl_vector *df);
int optimize(void);
void optimizet1();

#endif /* OPTI_H */
