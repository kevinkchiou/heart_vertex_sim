#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

#include "const.h"
#include "definelat.h"
#include "lattice.h"
#include "measurements.h"
#include "save.h"
#include "compalone.h"
#include "locerror.h"
#include "in_development.h"
#include "opti.h"



//gives each vertex an offset
void coord_mod(){
    int i;
    Node *pvert;

    for(i=0;i<nb_vertex_tot;i++){
       pvert = web[i];
       pvert->z = 12.0;
    }
}

void border_check(){
   int i,j,k;
   Dual *pcell,*pcellnb;
   Node *pvert;

   for(i=0;i<nb_vertex_tot;i++){
      pvert = web[i];
      pvert->border = 0;
   }
   for(i=0;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      pcell->marker.border=0;
   }

   /*
   for(i=start;i<nb_cell_tot;i++){
      pcell = web_dual[i];
      if((pcell->marker.tension_index)>1.001){
         pcell->area_soll=0.8;
         pcell->area_soll0=pcell->area_soll;
         for(j=0;j<(pcell->nb_vertices);j++){
            pvert=pcell->vertexlist[j];
            for(k=0;k<3;k++){
               pcellnb = pvert->pncell[k];
               if((pcellnb->marker.tension_index)<1.001){
                  pvert->border=1;
                  pcell->marker.border=1;
                  pcell->area_soll=0.8;
                  pcell->area_soll0=pcell->area_soll;
               }
            }
         }
      }
   }
   */
}

void division_cell_slow(Dual *pcell, int ivertex){
   int i,max=3,icell;//max is # of steps to get to max size
   double starting_area;
   Node *pvertex;

   pvertex=pcell->vertexlist[ivertex];
   icell=indexcell(pvertex,pcell);
   if(pcell->area_soll > pcell->area_soll0){starting_area = pcell->area_soll0;}
   else{starting_area = pcell->area_soll;}
   for(i=0;i<max;i++){
      pcell->area_soll=(1.0+((double) (i+1))/((double) max))*(starting_area);
      optimizet1();
   }
   division_eq(pvertex,icell);
}

void tension_modify(){
   Dual *pcell;
   int icell;

   printf("modifying tension...\n");
   add_factor = add_factor2;
   border_check();
   for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      if(pcell->marker.tension_index!=1.0){pcell->marker.tension_index+=(add_factor2-add_factor1);}
      border_check();
   }
   optimizet1();
}

void tension_modify1(){
	Dual *pcell;
	int icell,count=0,numdel;

	printf("modifying tension...\n");
	for(icell=0;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		if(pcell->marker.clone_index==1){pcell->marker.tension_index=8.0;pcell->area_soll=13.0;count++;}
		else{pcell->marker.tension_index=10.0;pcell->area_soll=9.0;}
		border_check();
	}
	printf("we have %d delta cells, and %d notch cells!\n",count,nb_cell_tot-count);
	optimize();
}

double tension_modify2(double eps){
   Dual *pcell;
   int icell;
   
   printf("modifying tension...\n");
   border_check();
   for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      pcell->marker.tension_index=1.0+2*eps*(gsl_rng_uniform(rng)-0.5);
   }
   optimizet1();
   return eps;
}

void tension_uniform(){
    int icell;
    Dual *pcell;

    printf("modifying tensions to attempt random lattice...");fflush(stdout);
    for(icell=0;icell<nb_cell_tot;icell++){
        pcell=web_dual[icell];
        pcell->area_soll=(LX*LY-1.0)/nb_cell_tot;
        pcell->marker.tension_index=pcell->area_soll;
    }
    optimizet1();
}

void tension_uniform_eraseclone(){
    int icell;
    Dual *pcell;

    printf("modifying tensions to attempt random lattice...");fflush(stdout);
    for(icell=0;icell<nb_cell_tot;icell++){
        pcell=web_dual[icell];
        pcell->marker.clone_index=0;
        pcell->area_soll=(LX*LY-1.0)/nb_cell_tot;
        pcell->marker.tension_index=pcell->area_soll;
    }
    optimizet1();
}

void kagome_pre_modify(){
	Dual *pcell,*pncell;
	int icell,jcell,deltaflag,nbneighb;

	printf("placing delta cells in cell with low numbers of neighbors...\n");
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb<5){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb==5){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
	for(icell=1;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		if(nbneighb==6){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){deltaflag=0;}
			}
			if(deltaflag==1){pcell->marker.clone_index=1;}
		}
	}
}
int **rnd_parking(int *occupied, int *unoccupied, int *oclength, int *unoclength){
    Dual *pcell;
    int ocsize,unocsize,i,j,count=0;
    int index,cell,icell,delcount=0;
    int nbneighb;

    ocsize=*oclength;unocsize=*unoclength;
    int *tempoc1,*tempunoc1,*tempoc2,*tempunoc2;
    int *tmp;

    tempoc1=(int*) malloc(sizeof(int)*ocsize);
    if(tempoc1==NULL){locerror("rnd_parking","out of memory!");}
    tempunoc1=(int*) malloc(sizeof(int)*unocsize);
    if(tempunoc1==NULL){locerror("rnd_parking","out of memory!");}

    for(i=0;i<ocsize;i++){tempoc1[i]=occupied[i];}
    for(i=0;i<unocsize;i++){tempunoc1[i]=unoccupied[i];}

    //printf("we have ocsize = %d and unocsize = %d!\n",ocsize,unocsize);
    //sanity check!
    for(i=0;i<ocsize;i++){
        for(j=0;j<unocsize;j++){
            if(tempoc1[i]==tempunoc1[j]){locerror("rnd_parking","the vectors share a term!");}
        }
    }

    index=1+(int)(gsl_rng_uniform(rng)*(unocsize-1));
    if(TOROIDAL==1){index=(int)(gsl_rng_uniform(rng)*unocsize);}
    cell=tempunoc1[index];
    tempunoc1[index]=-1;delcount+=1;
    pcell=web_dual[cell];
    nbneighb=pcell->nb_vertices;
    for(i=0;i<nbneighb;i++){
        icell=indexcellinwebdual(pcell->celllist[i]);
        for(j=0;j<unocsize;j++){
            if(icell==tempunoc1[j]){tempunoc1[j]=-1;delcount+=1;}
        }
    }
    unocsize-=delcount;
    ocsize+=1;

    tempoc2=(int*) malloc(sizeof(int)*ocsize);
    if(tempoc2==NULL){locerror("rnd_parking","out of memory!");}
    tempunoc2=(int*) malloc(sizeof(int)*unocsize);
    if(tempunoc2==NULL){locerror("rnd_parking","out of memory!");}
    for(i=0;i<(ocsize-1);i++){tempoc2[i]=tempoc1[i];}
    tempoc2[ocsize-1]=cell;
    //printf("cell = %d, ocsize-1 = %d,tempoc2[ocsize-1]= %d",cell,ocsize-1,tempoc2[ocsize-1]);
    count=0;
    for(i=0;i<(unocsize+delcount);i++){
        if(tempunoc1[i]!=-1){tempunoc2[count]=tempunoc1[i];count++;}
    }
    if(count!=unocsize){
        printf("we have count = %d and unocsize = %d and delcount = %d\n",count,unocsize,delcount);
        locerror("rnd_parking()","something wrong with counting!");
    }

    int **result = malloc(2*sizeof(int*));
    result[0]=malloc(ocsize*sizeof(int));
    result[1]=malloc(unocsize*sizeof(int));
    for(i=0;i<ocsize;i++){result[0][i]=tempoc2[i];}
    for(i=0;i<unocsize;i++){result[1][i]=tempunoc2[i];}
    *oclength=ocsize;*unoclength=unocsize;

    free(tempoc1);free(tempoc2);
    free(tempunoc1);free(tempunoc2);

    return result;
}
void freeparking(int **park){
    int i;

    free(park[0]);
    free(park[1]);
    free(park);
}
void kagome_nn_modify(){ //places extra d-cells for cells completely surrounded by n-cells
    Dual *pcell,*pncell,*pnncell;
    int icell,jcell,kcell,deltaflag,nbneighb,count=0;
    int *occupied,*unoccupied,oclength,unoclength;
    int i,j,**result;

    //initialization in case there are cells already placed as deltas.
    oclength=0;unoclength=0;
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        nbneighb=pcell->nb_vertices;
        if(pcell->marker.clone_index==1){ //neighbors for delta clones
            oclength+=1;
            for(j=0;j<nbneighb;j++){
                pncell=pcell->celllist[j];
                if(pncell->marker.clone_index==0){pncell->marker.clone_index=2;}//set for convenience.  eliminated later
            }
        }
    }
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==0){unoclength+=1;}
    }
    occupied=(int*)malloc(oclength*sizeof(int));
    unoccupied=(int*)malloc(unoclength*sizeof(int));
    oclength=0;unoclength=0;
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==1){occupied[oclength]=i;oclength+=1;}
        if(pcell->marker.clone_index==0){unoccupied[unoclength]=i;unoclength+=1;}
    }
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(pcell->marker.clone_index==2){pcell->marker.clone_index=0;} //reset the index used for convenience
    }
    printf("placing delta cells...\n");

    while(unoclength>0){
        result=rnd_parking(occupied,unoccupied,&oclength,&unoclength);
        free(occupied);free(unoccupied);
        occupied=malloc(oclength*sizeof(int));unoccupied=malloc(unoclength*sizeof(int));
        for(i=0;i<oclength;i++){occupied[i]=result[0][i];}
        for(i=0;i<unoclength;i++){unoccupied[i]=result[1][i];}
        freeparking(result);
        /*
        printf("occupied elements: ");
        for(i=0;i<oclength;i++){printf("%d ",occupied[i]);}
        printf("\nunoccupied elements:");
        for(i=0;i<unoclength;i++){printf("%d ",unoccupied[i]);}
        printf("\n");fflush(stdout);
        */
    }
    printf("we have that %d cells are delta'd!\n",oclength);
    for(i=0;i<oclength;i++){
        pcell=web_dual[occupied[i]];
        pcell->marker.clone_index=1;
    }
}
void kagome_nnn_modify(){
   Dual *pcell,*pncell,*pnncell;
   int icell,jcell,kcell,deltaflag,nbneighb,count=0;

	printf("placing delta cells...\n");
    while(count<10000){
		deltaflag=1;icell=(int)(gsl_rng_uniform(rng)*nb_cell_tot);count++;
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		for(jcell=0;jcell<nbneighb;jcell++){
			pncell=pcell->celllist[jcell];
			if(pncell->marker.clone_index==1){deltaflag=0;}
            for(kcell=0;kcell<nbneighb;kcell++){
                pnncell=pncell->celllist[kcell];
                if(pnncell->marker.clone_index==1){deltaflag=0;}
            }
		}
		//if(deltaflag==1 && nbneighb>5){pcell->marker.clone_index=1;}
		if(deltaflag==1){pcell->marker.clone_index=1;}
	}
}
void kagome_modify2(){ //this is only for the specific toroid generated by alberto's code
    int i;
    Dual *pcell;

    printf("placing delta cells...\n");
    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        if(i==1||i==7||i==13||i==19||i==26||i==32||i==38||i==44||i==4||i==10||i==16||i==22||i==29||i==35||i==41||i==47){pcell->marker.clone_index=1;}
    }
}
void kagome_dblmodify(){
	Dual *pcell,*pncell;
	int icell,jcell,deltaflag,nbneighb,count,numdelta=0,numnotch=0;
	int numchange;

	printf("placing secondary delta cells...\n");
	for(icell=0;icell<nb_cell_tot;icell++){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		deltaflag=0;count=0;
		if(pcell->marker.clone_index==1){numdelta++;}
		if(pcell->marker.clone_index==0){numnotch++;}
		if(nbneighb>4 && pcell->marker.clone_index==0){
			deltaflag=1;
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){count++;}
			}
			if(count>1){deltaflag=0;}
		}
		if(deltaflag==1){pcell->marker.clone_index=1;numdelta++;numnotch--;}
	}
	for(icell=0;icell<nb_cell_tot;icell++){ //this is to eliminate the rare double event
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;
		count=0;
		if(pcell->marker.clone_index==1){
			for(jcell=0;jcell<nbneighb;jcell++){
				pncell=pcell->celllist[jcell];
				if(pncell->marker.clone_index==1){count++;}
			}
		}
		if(count>1){pcell->marker.clone_index=0;numdelta--;numnotch++;}
	}
	printf("we have %d delta cells and %d notch cells!\n",numdelta,numnotch+1);fflush(stdout);
	if(numdelta>numnotch/2-nb_cell_tot/80){numchange=(int)ceil(numdelta-numnotch/2+nb_cell_tot/80);}
	else{numchange=0;}
	printf("we need to change %d cells!\n",numchange);fflush(stdout);
	icell=1; //initialize
	while(numchange>0){
		pcell=web_dual[icell];
		nbneighb=pcell->nb_vertices;count=0;
		for(jcell=0;jcell<nbneighb;jcell++){
			pncell=pcell->celllist[jcell];
			if(pncell->marker.clone_index==1){count++;}
		}
		if(count>0 && pcell->marker.clone_index==1){pcell->marker.clone_index=0;numchange--;printf("changing cell %d...",icell);}
		icell++;
		if(icell>nb_cell_tot-1){break;}
	}
	printf("done!\n");fflush(stdout);
}
void kagome_reg(){
    Dual *pcell,*pncell;
    int i,j,nbneighb,delta;

    for(i=0;i<nb_cell_tot;i++){
        pcell=web_dual[i];
        delta=0;
        if(pcell->marker.clone_index==0){
            delta=1;nbneighb=pcell->nb_vertices;
            for(j=0;j<nbneighb;j++){
                pncell=pcell->celllist[j];
                if(pncell->marker.clone_index==1){delta=0;}
            }
        }
        if(delta==1){pcell->marker.clone_index=1;}
    }
}

   void
init_lattice2(int *array, int num)
{
   int icell,j,k;
   Dual *pcell;
   char filename1[100],filename2[100],fnum[4];


   printf("before coord\n");
   //coord_gen();//comment this out when using infoinput()
   printf("before shrink\n");
   sprintf(fnum,"%.4d",array[num]);
   sprintf(filename1,"info/%s.vertinfo",fnum);
   sprintf(filename2,"info/%s.cellinfo",fnum);
   infoinput(filename1,filename2);
   //cut_cell_find();
   //shrink_datafirst();//comment this out to use infoinput()
   //coord_mod();//to give the initial vertices an offset
   printf("after shrink\n");


   for (icell=0 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];
      pcell->marker.size=0;
      pcell->marker.div_rate=0;
      pcell->marker.clone_index=0;
      pcell->area_soll=1.0;
      pcell->area_soll0=pcell->area_soll;
      pcell->marker.tension_index=1.0;
   }

   forcevector=(double**)malloc(3*sizeof(double*));
   for(k=0;k<3;k++){
      forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
   }
   //svglattice("svgs/0000.svg",1);
   //latticeprint("dats/0000.dat",0);
   //cut_cell_multistep();
   /*for(icell=0;icell<nb_cell_tot;icell++){
      pcell=web_dual[icell];
      if(icell==20 || icell==23 || icell==27 || icell==24 || icell==30){pcell->marker.tension_index+=add_factor;}
      if(icell==0){pcell->area_soll=78;pcell->area_soll0=pcell->area_soll;}
      optimizet1();
   }*/
   itime = 0;
   it1process=0;
}

   void
import_lattice(char * filename)
{
   int icell,j,k;
   Dual *pcell;
   char filename1[100],filename2[100];


   sprintf(filename1,"input/%s.vertinfo",filename);
   sprintf(filename2,"input/%s.cellinfo",filename);
   infoinput(filename1,filename2,1);

   for (icell=0 ; icell<nb_cell_tot ; icell++){
      pcell=web_dual[icell];
      pcell->marker.size=0;
      pcell->marker.div_rate=0;
      pcell->area_soll=pcell->area_soll;
      pcell->area_soll0=pcell->area_soll;
      pcell->marker.tension_index=pcell->area_soll;
   }

   forcevector=(double**)malloc(3*sizeof(double*));
   for(k=0;k<3;k++){
      forcevector[k]=(double*)malloc(nb_vertex_tot*sizeof(double));
   }
   itime = 0;
   it1process=0;
}

void kagomet1_dd(){ //transforms delta cells away from each other
	int i,j,k,l,nbneighb,exec_t1;
	Dual *pcell,*pncell;
	Node *pvert,*pnvert,*povert;
	int debugcount=0,start;
	char filename[50];

    start=1;
    if(TOROIDAL==1){start=0;}
	//this first part does cells for neighboring delta'd cells
	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if(pcell->marker.clone_index==1){
			//printf("we found a cell to transform! it is cell %d\n",i);fflush(stdout);
			nbneighb=pcell->nb_vertices;
			//here we check to see if any delta'd cells are neighbors.  if so, then transform them apart. assumes at most 1 neighboring delta.
			exec_t1=0;
			for(j=0;j<nbneighb;j++){
				pncell=pcell->celllist[j];
                //printf("we are in cell %d ",indexcellinwebdual(pcell));fflush(stdout);
                //printf("and neighbor %d ",j);fflush(stdout);
                //printf("which has global index %d!\n",indexcellinwebdual(pncell));fflush(stdout);
				if(pncell->marker.clone_index==1){
					//printf("okay we found two deltas next to each other. cell %d and neighbor %d...",i,j);fflush(stdout);
					//this part just finds the index of the other vertex in the pneighb list on vertex structures
					exec_t1=1;
					pvert=pcell->vertexlist[j];pnvert=pcell->vertexlist[(j+1)%nbneighb];//transform vertices on the boundary between delta'd cells
					for(k=0;k<3;k++){
						povert=pvert->pneighb[k];
						if(povert==pnvert){l=k;}
					}
				}
			}
			if(exec_t1==1){t1transform(pvert,l);/*printf("transformed!\n");fflush(stdout);*/} //transform to separate the delta cells
		}
	}
	updatedata();
}
void kagomet1_nd(){
	int i,j,k,l,nbneighb,exec_t1;
	Dual *pcell,*pncell;
	Node *pvert,*pnvert,*povert;
	int debugcount=0,start;
	char filename[50];

    start=1;
    if(TOROIDAL==1){start=0;}
	for(i=start;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		if(pcell->marker.clone_index==1){
			//printf("we found a cell to transform! cell %d\n",i);fflush(stdout);
			//here we try to stick as many notch to the deltas as we can thru t1transforms
			nbneighb=pcell->nb_vertices;
			for(j=0;j<nbneighb;j++){
				//printf("looking at vertex %d in cell %d...",j,i);fflush(stdout);
				pvert=pcell->vertexlist[j];
				for(k=0;k<3;k++){
					if((pvert->pneighb[k]!=pcell->vertexlist[(j-1)%nbneighb]) && (pvert->pneighb[k]!=pcell->vertexlist[(j+1)%nbneighb])){pnvert=pvert->pneighb[k];l=k;}
				}
				exec_t1=1;
				for(k=0;k<3;k++){//this determines whether or not any neighbors are delta'd. if so, then don't execute t1
					pncell=pnvert->pncell[k];
					if(pncell->marker.clone_index==1){
						//printf("awww, no transform :(\n");fflush(stdout);
						exec_t1=0;
					}
				}
				if(exec_t1==1){
					//printf("performing t1!\n");fflush(stdout);
					t1transform(pvert,l);
					optimize();
					sprintf(filename,"debugsvg/%.4i.svg",debugcount);
					//svglattice(filename,1);//uncomment to see each step
					debugcount++;
				}
			}
		}
	}
	updatedata();

    Bond *pbond;
    int m;

    for(i=0;i<nb_bond_tot;i++){
        pbond=web_bond[i];m=0;l=0;k=-1;
        pvert=pbond->pnvert[0];pnvert=pbond->pnvert[1];
        for(j=0;j<3;j++){
            if(pvert->pneighb[j]==pnvert){k=j;}
        }
        if((pbond->pncell[0])->marker.clone_index==0){l+=1;}
        if((pbond->pncell[2])->marker.clone_index==0){l+=1;}
        if((pbond->pncell[1])->marker.clone_index==0){m+=1;}
        if((pbond->pncell[3])->marker.clone_index==0){m+=1;}
        if(l==2 && m==1){t1transform(pvert,k);}
    }
    updatedata();
}

void pulse_event(double step_frac,double n_str){

	int i,j;
	Dual *pcell;

	for(i=1;i<nb_cell_tot;i++){
		pcell=web_dual[i];
		pcell->area_soll=pcell->area_soll-pcell->area_soll0*step_frac+2*n_str*(gsl_rng_uniform(rng)-0.5);
	}
}

int *division_select_once(int divlist[],int divlistlen,int *pnum){ //selects one number from array and removes it

	int i,idx,len,*templist;

	len=divlistlen;
	idx=floor(gsl_rng_uniform(rng)*len);//select one cell randomly
	*pnum=divlist[idx];

	//temporarily store the remaining elements
	templist=(int*)malloc((len-1)*sizeof(int));
	for(i=0;i<idx;i++){templist[i]=divlist[i];}
	for(i=idx;i<len-1;i++){templist[i]=divlist[i+1];}

	return templist;
}

void torus_division_set_parameters(){

	int i,j;
	Dual *pcell;

	for(i=0;i<nb_cell_tot;i++){
		pcell=web_dual[i];

		pcell->marker.tension_index=3.0;
		pcell->area_soll=LX*LY/nb_cell_tot;
		pcell->area_soll0=pcell->area_soll;
	}
}
