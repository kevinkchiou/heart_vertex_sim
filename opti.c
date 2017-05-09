#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "definelat.h"
#include "lattice.h"
#include "save.h"
#include "measurements.h"
#include "pngwrite.h"
#include "compalone.h"
#include "opti.h"
#include "in_development.h"

double prefactor1 = 100.0;
double prefactor2 = 0.0;
double radius = 8.0;
/* Make a t1 process if necessary
 */
int check_t1()
{
   int i,j;
   static int nbt1=0;
   Node *pvert, *pvertn;

   for (i=0;i<nb_vertex_tot;i++){
      pvert=web[i];
      for (j=0;j<3;j++){
         pvertn=pvert->pneighb[j];
         if (dist(pvert,pvertn)<CUTOFF_INF && flagt1<4){
            printf("performing t1 transform on vertex %d and neighbor %d\n",indexvertexinweb(pvert),j);
            t1flagcheck(indexvertexinweb(pvert),j);//increments flagt1 when necessary to prevent loop
            t1transform(pvert,j);
            border_check();
            nbt1++;
            return 1;
         }
      }
   }
   if (flagt1>=4){flagt1=0;}
   return 0;
}

/* update perimeter and area
 */
void update_pa()
{
    int i;
    double carea;
    Dual *pcell;
    for (i=0;i<nb_cell_tot;i++){
	pcell=web_dual[i];
	carea=area(pcell);
	pcell->area=carea;
	pcell->sqrtarea=sqrt(carea);
	pcell->perimeter=perimeter(pcell);
	pcell->centroid=centroid(pcell);
    }
}

/* energy function of the area part
 */
    double
area_fct(Dual *pcell)
{
    double fct;
    double A,Asoll;

    A = pcell->area;
    Asoll = pcell->area_soll;

    //fct = 0.5*(A-Asoll)*(A-Asoll)/pow(A,0.3);//newer H
    fct = 0.5*(A-Asoll)*(A-Asoll)/sqrt(A);//new H

   if(TOROIDAL==0){
       if(indexcellinwebdual(pcell)==0){fct=0;}//exclude cell0 term
   }
    
	//fct=0.0;
    return fct;
}

/* derivative of energy function for the area part (wrt area)
 */
    double
area_fctp(Dual *pcell)
{
    double fct;
    double A,Asoll;

    A = pcell->area;
    Asoll = pcell->area_soll;

    //fct = (A-Asoll)/pow(A,0.3) - 0.15*(A-Asoll)*(A-Asoll)/pow(A,1.3);//newer H
    fct = (A-Asoll)/sqrt(A) - 0.25*(A-Asoll)*(A-Asoll)/sqrt(A*A*A);//new H
    
    if(TOROIDAL==0){
        if(indexcellinwebdual(pcell)==0){fct=0;}//exclude cell0 term
    }

	//fct=0.0;
    return fct;
}

double d2E_dA2(Dual *pcell)
{
	double fct;
	double A,Asoll;

	A=pcell->area;Asoll=pcell->area_soll;

	if(TOROIDAL==0){if(pcell->idx==0){fct=0;}}//exclude cell0

	fct = (1.0/sqrt(A))*(Asoll/A+3/8*(1-Asoll/A)*(1-Asoll/A));
}

/* energy function of the perimeter part
 */
    double
rho_fct(Dual *pcell)
{
    double rho,perim;
    double coeff,P0=1;
    int i,j,k;
    Node *pvert,*pnb;

	perim = pcell->perimeter;	
    //rho = pow(pcell->perimeter,RHOEXP) - RHOFCT*pcell->perimeter*cell_sizes[pcell->marker.size];
    coeff = pcell->marker.tension_index;
    //rho = coeff*(perim-P0)*(perim-P0);
	rho = coeff*perim;
	/*
    if(pcell->marker.tension_index!=1.0){
       for(i=0;i<pcell->nb_vertices;i++){
          pvert=pcell->vertexlist[i];
          if(pvert->border==1){
             for(j=0;j<3;j++){
                pnb = pvert->pneighb[j];
                if(pnb->border==1){
                   for(k=0;k<3;k++){
                      if(pnb->pncell[k]==pcell){rho+=0.5*coeff*dist(pvert,pvert->pneighb[j]);}
                   }
                }
             }
          }
       }
    }
	*/
	//if(indexcellinwebdual(pcell)==0 &&TOROIDAL==0){rho=0.0;}
	//rho=0.0;
    return rho;
}

/* derivative  of the perimeter part of the Hamiltonian
 */
    double
rho_fctpp(Dual *pcell)
{
   int i,j,k;
   double coeff;
   Node *pvert,*pnb;
   coeff = pcell->marker.tension_index;
	//coeff = 0.0;
   return coeff;
}

/* energy function of length part */

double rhol_fct(Bond *pbond)
{
	double rhol,d,kappa,L0,lambda;

	kappa = pbond->kappa;L0=pbond->L0;lambda = pbond->lambda;
	d = dist(pbond->pnvert[0],pbond->pnvert[1]);
	rhol = 0.5*kappa*(d-L0)*(d-L0)-lambda*d; 
	//printf("rhol = %f,k = %f, L0 = %f\n",rhol,kappa,L0);fflush(stdout);
	rhol=0.0;
	return rhol;
}

double rhol_fctpl(Bond *pbond)
{
	double rholl,d,kappa,L0,lambda;

	kappa = pbond->kappa;
	
	if(kappa<0.0001){return 0;}
	L0=pbond->L0;lambda = pbond->lambda;
	d = dist(pbond->pnvert[0],pbond->pnvert[1]);

	rholl = kappa*(d-L0)-lambda;
	//printf("rholl = %f,kappa = %f, L0 = %f\n",rholl,kappa,L0);fflush(stdout);
	rholl=0.;
	return rholl;
}

double d2E_dL2(Bond *pbond)
{
	double rholll,d,kappa,L0,lambda;

	kappa=pbond->kappa;

	if(kappa<0.0001){return 0;}
	L0=pbond->L0;lambda=pbond->lambda;
	d=dist(pbond->pnvert[0],pbond->pnvert[1]);

	rholll = kappa;
	return rholll;
}


/* derivative  of the perimeter part of the Hamiltonian
 */
    double
rho_fctpa(Dual *pcell)
{
    return 0;
}

double extern_pot_fct(Dual *pcell){
    double ext=0.0,x,y,z;
    int i,num = pcell->nb_vertices;
    Node *pvert;

    ext = 0.0;
    /*
       for(i=0;i<num;i++){
       pvert = pcell->vertexlist[i];
       x = pvert->x;y = pvert->y;z = pvert->z;
       ext+=0.5*prefactor1*(sqrt(x*x+y*y+z*z)-radius)*(sqrt(x*x+y*y+z*z)-radius);
       }*/
    return ext;//(distance-radius)^2 is the parameter
}

vector3 extern_pot_fctp(Node *pvert){
    vector3 dpot;
    double x,y,z;

    dpot.x=dpot.y=dpot.z=0.0;//for flat external potential case
    /*
       x=pvert->x;y=pvert->y;z=pvert->z;

       dpot.x=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
       dpot.y=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
       dpot.z=prefactor1*(sqrt(x*x+y*y+z*z)-radius);
     */
    return dpot;
}

/* ------------------------------------------------------------------------------ */

    double
energy(const gsl_vector *x,void *params)
{   /* Compute the energy */
    int i;
    double h=0;
    Dual *pcell;
	Bond *pbond;

    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(x,3*i);
	web[i]->y=gsl_vector_get(x,3*i+1);
	web[i]->z=gsl_vector_get(x,3*i+2);
    }
    update_pa();
	
    for (i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		h+=rho_fct(pcell);
		h+=((double *)params)[0]*area_fct(pcell)+extern_pot_fct(pcell);
    }
	
	for(i=0;i<nb_bond_tot;i++){
		pbond=web_bond[i];
		h+=rhol_fct(pbond);
	}

    return h;
}

    void
denergy(const gsl_vector *x, void *params, gsl_vector *df)
{
    /* Compute the gradient of the energy */
    int i,j;
    double vertx,verty,vertz;
    double dx, dy, dz, dxa, dya, dza, dxd, dyd, dzd;
    double locdist, prefact;
    double ax,ay,az,bx,by,bz;
    double xcross,ycross,zcross;
    Node *pvert, *pvert_prev, *pvert_next;
    Dual *pcell0, *pcell;
	Bond *pbond_prev,*pbond_next;
    vector3 cent,dpot;

    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(x,3*i);
	web[i]->y=gsl_vector_get(x,3*i+1);
	web[i]->z=gsl_vector_get(x,3*i+2);
    }
    update_pa();
    pcell0=web_dual[0];
    for (i=0;i<nb_vertex_tot;i++){
	pvert=web[i];
	dx=dy=dz=dxa=dya=dza=dxd=dyd=dzd=0;
	for (j=0; j<3; j++){
	    pcell    = pvert -> pncell[j];
		pbond_next = pvert->pnbond[j];
	    pvert_next = pvert ->pneighb[j];
		pbond_prev = pvert->pnbond[(j+1)%3];
	    pvert_prev = pvert -> pneighb[(j+1)%3];
	    cent = pcell->centroid;

		//check of stuff!
		//if(pbond_next->pnvert[0]!=pvert_next && pbond_next->pnvert[1]!=pvert_next){locerror("opti","bond error");}
		//if(pbond_prev->pnvert[0]!=pvert_prev && pbond_next->pnvert[1]!=pvert_prev){locerror("opti","bond error");}

	    // Displacement along x, y, z due to perimeter extension 
		locdist=dist(pvert_next,pvert);
		prefact = rho_fctpp(pcell)/locdist;
		prefact += 0.5*rhol_fctpl(pbond_next)/locdist;
		dx += distx(pvert,pvert_next)*prefact;
		dy += disty(pvert,pvert_next)*prefact;
		dz += distz(pvert,pvert_next)*prefact;

		locdist=dist(pvert_prev,pvert);
		prefact = rho_fctpp(pcell)/locdist;
		prefact += 0.5*rhol_fctpl(pbond_prev)/locdist;
		dx += distx(pvert,pvert_prev)*prefact;
		dy += disty(pvert,pvert_prev)*prefact;
		dz += distz(pvert,pvert_prev)*prefact;


		// Displacement along x and y due to area extension 

		ax=distx(pvert_next,pvert);ay=disty(pvert_next,pvert);az=distz(pvert_next,pvert); 
		bx=centdistx(cent,pvert);by=centdisty(cent,pvert);bz=centdistz(cent,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		prefact = ((double *)params)[0]*area_fctp(pcell)+rho_fctpa(pcell); //CHECK_HERE as well
		//dxa += prefact*0.5*(pvert_next -> y - pvert_prev -> y); 
		//dya += prefact*0.5*(pvert_prev -> x - pvert_next -> x);
		//dza += prefact*0.5*(pvert_prev -> z - pvert_next -> z);
		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);

		ax=centdistx(cent,pvert);ay=centdisty(cent,pvert);az=centdistz(cent,pvert); 
		bx=distx(pvert_prev,pvert);by=disty(pvert_prev,pvert);bz=distz(pvert_prev,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		

	}

	dpot=extern_pot_fctp(pvert);
	dxd+=dpot.x;
	dyd+=dpot.y;
	dzd+=dpot.z;

	gsl_vector_set(df, 3*i, dx+dxa+dxd);
	gsl_vector_set(df, 3*i+1, dy+dya+dyd);
	gsl_vector_set(df, 3*i+2, dz+dza+dzd);

   forcevector[0][i]=gsl_vector_get(df,3*i);
   forcevector[1][i]=gsl_vector_get(df,3*i+1);
   forcevector[2][i]=gsl_vector_get(df,3*i+2);
    }
}


void enerdener(const gsl_vector *x,void *params,double *f, gsl_vector *df){

	*f=energy(x,params);
	denergy(x,params,df);
}

	
/*
    void
enerdener(const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
    // Compute both the gradient and the energy 
    // Compute the gradient of the energy 
    int i,j;
    double h;
    double vertx,verty,vertz;
    double dx, dy, dz,dxa, dya, dza, dxd, dyd, dzd;
    double locdist, prefact;
    double ax,ay,az,bx,by,bz;
    double xcross,ycross,zcross;
    Node *pvert, *pvert_prev, *pvert_next;
    Dual *pcell0, *pcell;
    vector3 cent,dpot;
	Bond *pbond;

    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(x,3*i);
	web[i]->y=gsl_vector_get(x,3*i+1);
	web[i]->z=gsl_vector_get(x,3*i+2);
    }

    update_pa();
    pcell0=web_dual[0];
    for (i=0;i<nb_vertex_tot;i++){
	pvert=web[i];
	dx=dy=dz=dxa=dya=dza=dxd=dyd=dzd=0;
	for (j=0; j<3; j++){
		pbond = pvert->pnbond[j];
	    pcell    = pvert -> pncell[j];
	    pvert_next = pvert -> pneighb[j];
	    pvert_prev = pvert -> pneighb[(j+1)%3];
	    cent = pcell->centroid;

	    // Displacement along x, y, and z due to perimeter extension 
		locdist=dist(pvert_next,pvert);
		//prefact = rho_fctpp(pcell)/locdist;//CHECK_HERE
		prefact = rhol_fctpl(pbond)/locdist;
		dx += distx(pvert,pvert_next)*prefact;
		dy += disty(pvert,pvert_next)*prefact;
		dz += distz(pvert,pvert_next)*prefact;

		locdist=dist(pvert_prev,pvert);
		//prefact = rho_fctpp(pcell)/locdist;//CHECK_HERE
		prefact = rhol_fctpl(pbond)/locdist;
		dx += distx(pvert,pvert_prev)*prefact;
		dy += disty(pvert,pvert_prev)*prefact;
		dz += distz(pvert,pvert_prev)*prefact;

		// Displacement alopng x and y due to area extension

		
		//commented out for now for this project
		
		ax=distx(pvert_next,pvert);ay=disty(pvert_next,pvert);az=distz(pvert_next,pvert); 
		bx=centdistx(cent,pvert);by=centdisty(cent,pvert);bz=centdistz(cent,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		prefact = ((double *)params)[0]*area_fctp(pcell)+rho_fctpa(pcell);//CHECK_HERE
		//dxa += prefact*0.5*(pvert_next -> y - pvert_prev -> y); 
		//dya += prefact*0.5*(pvert_prev -> x - pvert_next -> x);
		//dza += prefact*0.5*(pvert_prev -> z - pvert_next -> z);

		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);

		ax=centdistx(cent,pvert);ay=centdisty(cent,pvert);az=centdistz(cent,pvert); 
		bx=distx(pvert_prev,pvert);by=disty(pvert_prev,pvert);bz=distz(pvert_prev,pvert);
		zcross = ax*by-ay*bx; xcross = ay*bz-az*by; ycross = az*bx-ax*bz;
		dxa+=prefact*0.5*(zcross*(ay-by)-ycross*(az-bz))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dya+=prefact*0.5*(xcross*(az-bz)-zcross*(ax-bx))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		dza+=prefact*0.5*(ycross*(ax-bx)-xcross*(ay-by))/sqrt(zcross*zcross+ycross*ycross+xcross*xcross);
		
	    
	}

	dpot = extern_pot_fctp(pvert);
	dxd+=dpot.x;
	dyd+=dpot.y;
	dzd+=dpot.z;

	gsl_vector_set(df, 3*i, dx+dxa+dxd);
	gsl_vector_set(df, 3*i+1, dy+dya+dyd);
	gsl_vector_set(df, 3*i+2, dz+dza+dzd);
    }
    h=0;

    for (i=0;i<nb_cell_tot;i++){
	pcell=web_dual[i];
	h+=rho_fct(pcell)+((double *)params)[0]*area_fct(pcell)+extern_pot_fct(pcell);
    }
    *f=h;
}
*/

    int
optimize(void)
{
    int i;
    int iter=0;
    int status;
    gsl_vector *startpos;
    gsl_multimin_function_fdf hamilton;
    const gsl_multimin_fdfminimizer_type *minitype;
    gsl_multimin_fdfminimizer *minimize;
    double par[1] = {(double) MUP};

    int N=3*nb_vertex_tot;
    hamilton.f      = &energy;
    hamilton.df     = &denergy;
    hamilton.fdf    = &enerdener;
    hamilton.n      = N;
    hamilton.params = &par;


	printf("nb_bond_tot = %d\n",nb_bond_tot);fflush(stdout);
    startpos = gsl_vector_alloc(N);
    for (i=0;i<nb_vertex_tot;i++){
	gsl_vector_set(startpos,3*(i),web[i]->x);
	gsl_vector_set(startpos,3*(i)+1,web[i]->y);
	gsl_vector_set(startpos,3*(i)+2,web[i]->z);
    }
    minitype=gsl_multimin_fdfminimizer_conjugate_pr;
    //   minitype=gsl_multimin_fdfminimizer_vector_bfgs;
    //   minitype=gsl_multimin_fdfminimizer_conjugate_fr;
    minimize=gsl_multimin_fdfminimizer_alloc(minitype, N);

    gsl_multimin_fdfminimizer_set(minimize, &hamilton, startpos, 0.00001, 0.00001);
    do {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate(minimize);

   //if(iter==500){printf("iter=500\n");}
   //if(iter==750){printf("iter=750\n");}
   if(iter==999){printf("iter=999, it is about to break...\n");}
	if (status) break;

	status = gsl_multimin_test_gradient(minimize->gradient, 1e-4);
    } while (status==GSL_CONTINUE && iter<1000);
    for (i=0;i<nb_vertex_tot;i++){
	web[i]->x=gsl_vector_get(minimize->x,3*(i));
	web[i]->y=gsl_vector_get(minimize->x,3*(i)+1);
	web[i]->z=gsl_vector_get(minimize->x,3*(i)+2);
    }

    gsl_multimin_fdfminimizer_free(minimize);
    gsl_vector_free(startpos);
    return 0;
}


    void
optimizet1()
{
    int icheck=1; 
    optimize();
    do {
	if (check_t1()){
	    optimize();
	    icheck=1;
	}
	else icheck=0;
    } while (icheck);
}


