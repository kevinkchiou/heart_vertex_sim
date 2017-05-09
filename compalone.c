/* Standalone program (does not need opengl) */
// export CPPFLAGS="-I/sw/include/"
// export LDFLAGS="-L/sw/lib"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "lattice.h"
#include "save.h"
#include "read.h"
#include "pngwrite.h"
#include "locerror.h"
#include "measurements.h"
#include "compalone.h"
#include "opti.h"
#include "in_development.h"
#include "dynamics.h"

	int
checkcell(int icell)
{
	int i, j; 
	Node *pvert;
	Dual *pcell0, *pcell;
	pcell0=web_dual[0];
	pcell=web_dual[icell];
	for (i=0;i<pcell->nb_vertices;i++){
		pvert=pcell->vertexlist[i];
		j=indexcell(pvert,pcell);
		if (pvert->pncell[(j+1)%3]==pcell0) return 0;
	}
	return 1;
}


	double
feedbackdiv(Dual *pcell)
{
	double feedback;

	//   feedback =1.;
	feedback=exp(-10/area(pcell));
	return feedback;
}

int divc_selective(){
	int number=0,number2=0,nbvert=0;
	int i,j,k,retval=0,rand2;
	double rand,bound,mult_factor=0.80,tot;//variables for clone cells on the boundary
	double bound2,factor2=1.2;//variables for clone cells inside boundary
	Dual *pcell;

	for(i=1;i<nb_cell_tot;i++){
		pcell = web_dual[i];
		//if(pcell->marker.border==1){number++;}
		//if(pcell->marker.tension_index!=1.0 && pcell->marker.border!=1){number2++;}
		if(pcell->marker.clone_index==1){number++;}
	}
	rand = (gsl_rng_uniform(rng));
	//if(number==0){mult_factor=0.0;/*division weight for boundary clones*/}
	//if(number2==0){factor2=0.0;/*division weight for non-boundary clones*/}
	tot = (mult_factor-1.0)*number+(factor2-1.0)*number2+nb_cell_tot; //total weight
	bound = mult_factor*number/tot;
	bound2 = factor2*number2/tot + bound;
	//following if statements test for which type of cell is dividing and then picks
	//a random one of those cells at a random vertex and calls division_cell_slow()
	if(rand<bound){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)number));
			j=0;
			for(i=0;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				//if(pcell->marker.border==1){j++;}
				if(pcell->marker.clone_index==1){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(rand>bound && rand<bound2){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)(number2)));
			j=0;
			for(i=0;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				//if(pcell->marker.border==0 && pcell->marker.tension_index!=1.0){j++;}
				if(pcell->marker.clone_index==0){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(rand>bound2){
		while(nbvert<=4){
			rand2 = ceil((gsl_rng_uniform(rng)*(double)(nb_cell_tot-number-number2-1)));
			j=0;
			for(i=1;i<nb_cell_tot;i++){
				pcell = web_dual[i];
				if(pcell->marker.border==0 && pcell->marker.tension_index==1.0){j++;}
				if(j==rand2){
					nbvert=pcell->nb_vertices;
					if(nbvert>4){division_cell_slow(pcell,floor(gsl_rng_uniform(rng)*nbvert));retval=1;break;}
				}
			}
		}
	}
	if(retval!=1){printf("number = %d, number2 = %d, rand2 = %d, rand = %lf, bound = %lf, bound2 = %lf, nb_cell_tot = %d\n",number,number2,rand2,rand,bound,bound2,nb_cell_tot);}//in case of failure for diagnostics
	return retval;
}

	int
divc_stocha()
{
	int i;
	int retval=0;
	Dual *pcell;


	i = 1 + (int) (gsl_rng_uniform(rng)*(nb_cell_tot-1));
	pcell=web_dual[i];
	division_cell_slow(pcell,(int)(gsl_rng_uniform(rng)*pcell->nb_vertices));
	retval=1;
	return retval;
}

	int
divc_feedback()
{
	int i;
	int retval=0;
	double feedback,sum,rn;
	Dual *pcell;


	sum = 0.;
	for (i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		sum+=feedbackdiv(pcell);
	}
	feedback = gsl_rng_uniform(rng);
	rn = feedback*sum;


	sum = 0.;
	i = 0;
	do {
		i ++;
		pcell=web_dual[i];
		sum+=feedbackdiv(pcell);
		//     printf("%f %f %f %i %i\n",sum,rn,feedback,i,nb_cell_tot);
	} while (i < nb_cell_tot && sum < rn);

	//   printf("%f %f %f %i %i\n",sum,rn,feedback,i,nb_cell_tot);

	pcell = web_dual[i];
	division_cell_slow(pcell,(int)(gsl_rng_uniform(rng)*pcell->nb_vertices));
	retval = 1;

	return retval;
}

unsigned int nb_killed=0;

	Dual *
compdivc(Dual *pcell)
{
	int nvert;
	int j;
	double tot=0;
	double dcell;
	Node *pvert;
	Dual *pcelln;

	for (nvert=0;nvert<pcell->nb_vertices;nvert++){
		pvert=pcell->vertexlist[nvert];
		j=indexcell(pvert,pcell);
		pcelln=pvert->pncell[(j+1)%3];
		if (pcelln!=web_dual[0]) tot+=feedbackdiv(pcelln);
	}
	dcell=tot*gsl_rng_uniform(rng);
	tot=0;
	for (nvert=0;nvert<pcell->nb_vertices;nvert++){
		pvert=pcell->vertexlist[nvert];
		j=indexcell(pvert,pcell);
		pcelln=pvert->pncell[(j+1)%3];
		if (pcelln!=web_dual[0]){
			tot+=feedbackdiv(pcelln);
			if (dcell<tot) return pcelln;
		}
	}
	return pcelln;
}

	int
killc_stocha()
{
	int i;
	int retval=0;
	Dual *pcell;

	for (i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
	}
	return retval;
}

	void
heartevol()
{
	int i, j, k;
	int nb_lat=1;
	char filename[256],filename1[256],filename2[256],filename3[256];
	FILE *pfile1, *pfile2;
	int numtimesteps=1000;

	for (i=0;i<nb_lat;i++){
		sprintf(filename,"initial");
        init_lattice();
		//init_lattice_torus();
		updatedata();
        infoprint("tested",0);
        printf("we have %d cells to play with!\n",nb_cell_tot);

		snprintf(filename1,255,"svgs/%s.svg",filename);
		snprintf(filename2,255,"dats/%s.dat",filename);
		snprintf(filename3,255,"%s",filename);

		svglattice(filename1,2);infoprint(filename3,2);

		int list[3]={1,2,3};
		init_evol(); 
		for(j=0;j<numtimesteps;j++){
			evol(0.01);
			printf("%d%%\r",(int)(100*j/numtimesteps));fflush(stdout);
			snprintf(filename1,100,"svgs/%.4i.svg",j);
			snprintf(filename2,100,"dats/%.4i.dat",j);
			snprintf(filename3,100,"%.4i",j);
			svglattice(filename1,2);infoprint(filename3,2);
			if(j%300==0){lead_contract(list,3);}
		}
		printf("\n");

		for(k=0;k<3;k++){free(forcevector[k]);}
		free(forcevector);
		forcevector=(double**)malloc(3*sizeof(double*));
		for(k=0;k<3;k++){
			forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
		}

	}

	kill_lattice();
	for(k=0;k<3;k++){free(forcevector[k]);}
	free(forcevector);
}

gsl_rng *rng;

/* Main function
 * Only call the content functions
 */
	int
main(int argc, char *argv[])
{

	rng=gsl_rng_alloc(gsl_rng_default);

	/* The following functions can be activated in compalone.h */
	heartevol();

	return 0;
}
