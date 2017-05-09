
void coord_mod();//modifies coords by an offset in z

void border_check();//checks to see which vertices are on the border for opti
void division_cell_slow(Dual *pcell, int ivertex);//other division algorithm to maintain planarity

void tension_modify();
void tension_modify1();
double tension_modify2(double eps);
void tension_uniform();
void tension_uniform_eraseclone();
void kagome_pre_modify();
void freeparking(int **park);
int **rnd_parking(int *oc, int *unoc, int *oclen,int *unoclen);
void kagome_nn_modify();
void kagome_nnn_modify();
void kagome_modify2();
void kagome_dblmodify();
void kagome_reg();
void init_lattice2(int *array, int num);
void import_lattice(char * filename);
void kagomet1_dd();
void kagomet1_nd();
void pulse_event(double step_size,double noise_strength);
int *division_select_once(int list_cells_to_divide[],int divlistlen,int *pnum);
void torus_division_set_parameters();
