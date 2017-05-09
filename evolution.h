
typedef struct
{
	double rho;
}Evolparam;

Evolparam *allocevolparam();

int evolfunc(double t, const double y[], double f[], void *params);
int evoljac(double t, const double y[], double *dfdy, double dfdt[], void *params);
void evol(const double delta_t);
void init_evol();
int hessian_length2d(gsl_matrix *m,int vertexindex,const double y[]);
int hessian_length3d(gsl_matrix *m,int vertexindex,const double y[]);
int hessian_area2d(gsl_matrix *m,int cellindex,const double y[]);
int hessian_area3d(gsl_matrix *m,int cellindex,const double y[]);
int test_symmetric(gsl_matrix *m);
int test_eig(gsl_matrix *m,gsl_vector *eval,gsl_matrix *evec);

void lead_contract(int *celllist,int celllist_length);
void check_cell_state(double baseline,double threshold,double dt);
