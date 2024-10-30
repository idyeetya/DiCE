/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\
  |*                                    *|
  |*         .-------.    ______        *|
  |*        /   o   /|   /\     \       *|
  |*       /_______/o|  /o \  o  \      *|
  |*       | o     | | /   o\_____\     *|
  |*       |   o   |o/ \o   /o    /     *|
  |*       |     o |/   \ o/  o  /      *|
  |*   jgs-'-------'     \/____o/       *|
  |*                                    *|
  |*          DiCE version 0.5          *|
  |*                                    *|
  |*     written by: Aditya Gupta       *|
  |*        & Guy 'Wayyne' Dayhoff      *|
  |*                                    *|
  |*          Dr. Sameer Varma          *|
  |*   Lab of Computational Biophysics  *|  
  |*    University of South Florida     *|
  |*       Tampa, Fl - 2023-2024        *|
  |*                                    *|
  \*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

//TO-DO 
//Support for multiple files at once (hard)
//Wizard type setup (easy)

#include "ezmd.h"
#include "svm.h"

#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


// MACROS //--------------------------------------------------------------------|

#define COST 100.f 		//Max value of Lagrange multiplier (default=100)
#define GAMMA 0.4f  				//RBD Kernel width (default=0.4)

#define SYS_DATA struct system_data
#define TRJ_DATA struct trajectory_data
#define FRM_DATA struct frame_data
#define ATM_DATA struct atom_data
#define MDL_DATA struct svm_model
#define PRB_DATA struct svm_problem
#define INP_DATA struct input_data

#define XX 0
#define YY 1
#define ZZ 2

#define PDB_ENS 0
#define TRR_ENS 1
#define XTC_ENS 2
#define DCD_ENS 3
#define XYZ_ENS 4
#define TXYZ_ENS 5

// GLOBALS //------------------------------------------------------------------|

XDRFILE *xfp_out;
matrix xtc_box;
HISTOGRAM *H;

// STRUCTURES //---------------------------------------------------------------|

struct system_data {
	int 		sz;						    //# of atoms
	ATM_DATA 	*atm;			   //atomic descriptions loaded from PDB
	TRJ_DATA 	**trj;         		       //effectiely trj[2], 0=fp1, 1=fp2
	int 		status;				       //status for xdrlib calls
	int  	ens_type;
	INP_DATA *inputs;
};


struct trajectory_data{
	int 		sz;						   //# of frames
	FRM_DATA 	**frm;	   	   //the trajectory as a dynamic array of frames
};


struct frame_data {
	int 		sz;						    //# of atoms
	float   	box[3];              		   	   //box dimensions, = 3
	float  	*pos;              			    //atomic coordinates
};


struct atom_data {
	int 		ndx;							 //index
	float        eta;                              //value yielded by DiCe method
	float      avg_eta;
	char         *rid;      //resisude id, name+num for output
	char 	*rec;		   //the complete record as read in from the pdb
	MDL_DATA 	*mdl;					//libsvm model structure
	PRB_DATA 	prb;			      	      //libsvm problem structure
};

struct input_data {
	//main inputs
	char *fileNameOne; //first file name
	char *fileNameTwo; //second file name
	FILE *fileOne;     //first file
	FILE *fileTwo;     //second file

		//optional arguments
	int wizard;       //flag for the wizard
	int help;         //flag to display help
	int full;         //write out the full pdb file - single frame by default
	int avg;          //calculate an average for each eta value -- store in eta_avg.dat and *_avg.pdb
	int print;        //option to write to a new pdb file -- required input for file name
	char *printArg;    //argument containing file name for print flag
	FILE *printFile;  //file for print flag
	int threads;
};

// GLOBALS //------------------------------------------------------------------|

XDRFILE *xfp_out;
matrix xtc_box;
HISTOGRAM *H;

// MEMORY MANAGMENT //---------------------------------------------------------|
SYS_DATA *new_system(void){
	SYS_DATA *sd = malloc(sizeof*sd);
	sd->trj = malloc(sizeof*sd->trj*2);
	sd->trj[0] = malloc(sizeof**sd->trj);
	sd->trj[1] = malloc(sizeof**sd->trj);
	sd->trj[0]->sz = 0;
	sd->trj[1]->sz = 0;
	sd->trj[0]->frm = NULL;
	sd->trj[1]->frm = NULL;
	sd->atm = NULL;
	sd->sz = 0;
	sd->status = 0;

	sd->inputs = calloc(1, sizeof*sd->inputs);
	sd->inputs->fileNameOne = calloc(100, sizeof*sd->inputs->fileNameOne);
	sd->inputs->fileNameTwo = calloc(100, sizeof*sd->inputs->fileNameTwo);
	sd->inputs->printArg = calloc(100, sizeof*sd->inputs->printArg);
	sd->inputs->threads = 32;
	return sd;
}

//place a new frame into a trajectory
int new_frame(TRJ_DATA *trj){
	FRM_DATA **fd;

	if ((fd = realloc(trj->frm,sizeof*trj->frm*++trj->sz)) == NULL){
		puts("ERROR: new_frame realloc failed. Goodbye");
		exit(1);
	}
	trj->frm = fd;

	//allocate a new frame into the array
	fd[trj->sz-1] = malloc(sizeof**trj->frm);
	fd[trj->sz-1]->box[XX] = 0;
	fd[trj->sz-1]->box[YY] = 0;
	fd[trj->sz-1]->box[ZZ] = 0;
	fd[trj->sz-1]->pos = NULL;
	fd[trj->sz-1]->sz = 0;

	return 1;
}


void free_system(SYS_DATA *sd){
	if (sd->atm){
		for (int i=0;i<sd->sz;i++){
			for (int j=0;j<sd->trj[0]->sz*2;j++)
				if (sd->atm[i].prb.x[j] != NULL)
					free(sd->atm[i].prb.x[j]);

			if (sd->atm[i].prb.x != NULL)
				free(sd->atm[i].prb.x);
			if (sd->atm[i].mdl != NULL)
				free(sd->atm[i].mdl);
			if (sd->atm[i].rid != NULL)
				free(sd->atm[i].rid);
			if (sd->atm[i].rec != NULL)
				free(sd->atm[i].rec);
		}
		free(sd->atm);
	}

	for (int i=0;i<2;i++)
		if (sd->trj[i]->frm){
			for (int j=0;j<sd->trj[i]->sz;j++){
				free(sd->trj[i]->frm[j]->pos); 
				free(sd->trj[i]->frm[j]); 
			}
			free(sd->trj[i]->frm); 
		}

	free(sd->trj[0]);
	free(sd->trj[1]);
	free(sd->trj);

	free(sd->inputs->fileNameOne);
	free(sd->inputs->fileNameTwo);
	if(sd->inputs->print == 1) {
		free(sd->inputs->printArg);
		fclose(sd->inputs->printFile);
	}
	fclose(sd->inputs->fileOne);
	fclose(sd->inputs->fileTwo);
	free(sd->inputs);
}

// INPUT //-------------------------------------------------------------------|

void help() {
	puts("Welcome to DiCE! Help menu:\n");
	puts("Structure: ./DiCE [args]\n");
	puts("FULL FORM ARGUMENTS");
	puts("--fileOne [filename]: First file to be compared, must be of form pdb, xtc or etc [REQUIRED]");
	puts("--fileTwo [filename]: Second file to be comparted, must be of form pdb, xtc or etc [REQURED]");
	puts("--wizard: Enter the setup wizard for the program (not yet operational)");
	puts("--print [filename]: Print eta values in pdb structure with filename [PDB] as a template");
	puts("--full: Requires the print flag, prints out the full filled pdb file instead of just a single frame");
	puts("--avg: Creates a seperate averaged eta value per residue and prints them in accordance with the print and full flag in seperate files");
	puts("\nSHORT FORM ARGUMENTS");
	puts("-f1: fileOne");
	puts("-f2: fileTwo");
	puts("-w: wizard");
	puts("-p: print");
	puts("-f: full");
	puts("-a: avg");
	puts("\nAll arguments must be seperated by a space and filenames must be directly after their argument");
	puts("\nEXAMPLES:");
	puts("./DiCE --fileOne first.pdb --fileTwo second.pdb -f -p template.pdb -a");
	puts("./DiCE --wizard (not yet operational)");
	return;
}

void parse_inputs(int argc, char **argv, SYS_DATA *sd) {
	for(int i = 1; i < argc; i++) {

		//check for long formed arguments
		if(argv[i][0] == '-' && argv[i][1] == '-') {
			if(!strcmp(argv[i], "--help")) {
				sd->inputs->help = 1;
				help();
				exit(0);
			}

			else if(!strcmp(argv[i], "--threads")) {
				if(i+1 == argc) { //sanity check
					puts("the threads flag MUST have a following integer -- exiting");
					exit(1);
				}
				sd->inputs->threads = atoi(argv[i+1]);
				i++;
			}
			else if(!strcmp(argv[i], "--wizard")) {
				sd->inputs->wizard = 1;
			}
			else if(!strcmp(argv[i], "--full")) {
				sd->inputs->full = 1;
			}
			else if(!strcmp(argv[i], "--avg")) {
				sd->inputs->avg = 1;
			}
			else if(!strcmp(argv[i], "--print")) {
				sd->inputs->print = 1;
				if(i+1 == argc) { //sanity check
					puts("the print flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->printArg = argv[i+1];
				i++;
			}
			else if(!strcmp(argv[i], "--fileOne")) {
				if(i+1 == argc) { //sanity check
					puts("the fileOne flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->fileNameOne = argv[i+1];
				i++;
			}
			else if(!strcmp(argv[i], "--fileTwo")) {
				if(i+1 == argc) { //sanity check
					puts("the fileTwo flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->fileNameTwo = argv[i+1];
				i++;
			}
			else {
				printf("%s not recognized -- ignoring\n", argv[i]);
			}
		}


		//check for short formed arguments
		else if(argv[i][0] == '-') {
			if(!strcmp(argv[i], "-h")) {
				sd->inputs->help = 1;
				help();
				exit(0);
			}
			else if(!strcmp(argv[i], "-t")) {
				if(i+1 == argc) { //sanity check
					puts("the threads flag MUST have a following integer -- exiting");
					exit(1);
				}
				sd->inputs->threads = atoi(argv[i+1]);
				i++;
			}

			else if(!strcmp(argv[i], "-w")) {
				sd->inputs->wizard = 1;
			}
			else if(!strcmp(argv[i], "-f")) {
				sd->inputs->full = 1;
			}
			else if(!strcmp(argv[i], "-a")) {
				sd->inputs->avg = 1;
			}
			else if(!strcmp(argv[i], "-p")) {
				sd->inputs->print = 1;
				if(i+1 == argc) { //sanity check
					puts("the print flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->printArg = argv[i+1];
				i++;
			}
			else if(!strcmp(argv[i], "-f1")) {
				if(i+1 == argc) { //sanity check
					puts("the fileOne flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->fileNameOne = argv[i+1];
				i++;
			}
			else if(!strcmp(argv[i], "-f2")) {
				if(i+1 == argc) { //sanity check
					puts("the fileTwo flag MUST have a following input with a file name - exiting");
					exit(1);
				}
				sd->inputs->fileNameTwo = argv[i+1];
				i++;
			}
			else {
				printf("%s not recognized -- ignoring\n", argv[i]);
			}
		}

		else {
			printf("%s not recognized -- ignoring\n", argv[i]);
		}
	}
}

void args_checkout(SYS_DATA *sd){
	void *f1 = NULL;
	void *f2 = NULL;
/*
	if(sd->inputs->print == 0 && sd->inputs->avg == 1) {
		puts("Cannot average without reference (template) pdb. Please use the -p flag");
		exit(1);
	}
*/

	if(sd->inputs->print == 0 && sd->inputs->full == 1) {
		puts("Cannot print full without reference (template) pdb. Please use the -p flag");
		exit(1);
	}

	if(!strcmp(sd->inputs->fileNameOne, "") || !strcmp(sd->inputs->fileNameTwo, "")) {
		puts("fileOne and/or fileTwo are missing and required.");
		puts("Use --help or --wizard for help");
		exit(1);
	}

	//determine the ensemble type and try to open the ensemble files
	if (!strcmp(sd->inputs->fileNameOne+(strlen(sd->inputs->fileNameOne)-3),"trr")){
		sd->ens_type = TRR_ENS;
		if (!(f1 = (void *)xdrfile_open(sd->inputs->fileNameOne, "rb")) 
				|| !(f2 = (void *)xdrfile_open(sd->inputs->fileNameTwo,"rb"))){
			printf("Alas, one of the specified input files cannot be opened.\n");
			exit(1);
		}
	}

	else if (!strcmp(sd->inputs->fileNameOne+(strlen(sd->inputs->fileNameOne)-3),"xtc")){
		sd->ens_type = XTC_ENS;
		if (!(f1 = (void *)xdrfile_open(sd->inputs->fileNameOne, "rb")) 
				|| !(f2 = (void *)xdrfile_open(sd->inputs->fileNameTwo,"rb"))){
			printf("Alas, one of the specified input files cannot be opened.\n");
			exit(1);
		}

	}
	else if (!strcmp(sd->inputs->fileNameOne+(strlen(sd->inputs->fileNameOne)-3),"pdb")){
		sd->ens_type = PDB_ENS;
		if (!(f1 = (void *)fopen(sd->inputs->fileNameOne,"r")) 
				|| !(f2 = (void *)fopen(sd->inputs->fileNameTwo,"r"))){
			printf("Alas, one of the specified input files cannot be opened.\n");
			exit(1);
		}
	}
	else{
		puts("Alas, the trajectory type provided couldn't be validated against"
				" *supported types.\n");
		puts(" *supported trajectory types include: PDB, TRR, and XTC.");
		exit(1);
	}
	
	if(sd->inputs->print == 1) {
		//attempt to open the output file
		void *f3 = NULL;
		if (!strcmp(sd->inputs->printArg+(strlen(sd->inputs->printArg)-3),"pdb")){
			if (!(f3 = (void *)fopen(sd->inputs->printArg,"r"))) {
				printf("Alas, the specified output files cannot be opened.\n");
				exit(1);
			}
			sd->inputs->printFile = f3;
		}
		else {
			puts("Output file is not of type PDB. Goodbye.");
			exit(1);
		}
	}

	sd->inputs->fileOne = f1;
	sd->inputs->fileTwo = f2;

}

// FILE PARSING //------------------------------------------------------------|

//if the ptr is not NULL, then we have a subset of atoms we want to ignore
int parse_xtc(void *xtcFile, SYS_DATA *sd, int trj_index, void *ptr, int size){

	int nat = sd->sz, step, status = sd->status, fCnt;
	float time, prec;

	if (trj_index == 0)
		sd->sz = 0; //reset system size 

	XDRFILE *xfp = (XDRFILE *)xtcFile;
	TRJ_DATA *trj = sd->trj[trj_index]; 

	if (xfp && status == exdrOK){
		rvec *k = malloc(nat*sizeof*k);

		//when the ptr is NULL we process all atoms in the frame
		if (ptr == NULL)            
			for(fCnt=0;(status = read_xtc(xfp, nat, &step, &time, xtc_box, k, &prec)) == exdrOK;fCnt++){
				//append a new frame to the system data
				new_frame(sd->trj[trj_index]);

				if ( trj->sz == 1 && trj_index == 0){
					for (int a=0;a<nat;a++){
						//expand the atom_data array
						sd->atm = realloc(sd->atm,sizeof*sd->atm*++sd->sz);
						ATM_DATA *atom = &(sd->atm[sd->sz-1]);
						atom->ndx = a;
						atom->rid = NULL;
						atom->rec = NULL;
						atom->mdl = NULL;
					}
				}

				//loop over each atom in the frame & store their coordinates
				for (int a=0;a<nat;a++){
					FRM_DATA *f = trj->frm[trj->sz-1];
					f->pos = realloc(trj->frm[trj->sz-1]->pos,sizeof(float)*++f->sz*3);
					f->pos[(f->sz*3)-3+XX] = k[a][XX];
					f->pos[(f->sz*3)-3+YY] = k[a][YY];
					f->pos[(f->sz*3)-3+ZZ] = k[a][ZZ];
				}
			}

		//otherwise, we pull out a subset of atoms for processing
		else{
			//alloc new rvec ptr for the subset
			rvec *l = malloc(size*sizeof*l); 
			for(fCnt=0;(status = read_xtc(xfp, nat, &step, &time, xtc_box, k, &prec)) == exdrOK;){
				for (int i=0;i<size;i++){
					l[i][0] = k[(((int *)ptr)[i])][0];
					l[i][1] = k[(((int *)ptr)[i])][1];
					l[i][2] = k[(((int *)ptr)[i])][2];
				}
			}

			//housekeeping
			free(l);
		}

		//housekeeping
		xdrfile_close(xfp);
		free(k);
	}

	//be verbose
	printf("%d frames loaded containing %d atoms\n",
			sd->trj[trj_index]->sz,sd->sz);

	return fCnt;
}


//if the ptr is not NULL, then we have a subset of atoms we want to ignore
int parse_trr(void *trrFile, SYS_DATA *sd, int trj_index, void *ptr, int size){

	int nat = sd->sz, step, status = sd->status, fCnt;
	float time, lambda;

	if (trj_index == 0)
		sd->sz = 0; //reset system size 

	XDRFILE *xfp = (XDRFILE *)trrFile;
	TRJ_DATA *trj = sd->trj[trj_index]; 

	if (xfp && status == exdrOK){
		rvec *x = malloc(nat*sizeof*x);
		rvec *v = malloc(nat*sizeof*v);
		rvec *f = malloc(nat*sizeof*f);

		//when the ptr is NULL we process all atoms in the frame
		if (ptr == NULL)
			for(fCnt=0;(status = read_trr(xfp, nat, &step, &time, &lambda, xtc_box, x, v, f)) == exdrOK;fCnt++){
				//append a new frame to the system data
				new_frame(sd->trj[trj_index]);

				if ( trj->sz == 1 && trj_index == 0){
					for (int a=0;a<nat;a++){
						//expand the atom_data array
						sd->atm = realloc(sd->atm,sizeof*sd->atm*++sd->sz);
						ATM_DATA *atom = &(sd->atm[sd->sz-1]);
						atom->ndx = a;
						atom->rid = NULL;
						atom->rec = NULL;
						atom->mdl = NULL;
					}
				}

				//loop over each atom in the frame & store their coordinates
				for (int a=0;a<nat;a++){
					FRM_DATA *f = trj->frm[trj->sz-1];
					f->pos = realloc(trj->frm[trj->sz-1]->pos,sizeof(float)*++f->sz*3);
					f->pos[(f->sz*3)-3+XX] = x[a][XX];
					f->pos[(f->sz*3)-3+YY] = x[a][YY];
					f->pos[(f->sz*3)-3+ZZ] = x[a][ZZ];
				}
			}

		//otherwise, we pull out a subset of atoms for processing
		else{
			//alloc new rvec ptr for the subset
			rvec *X = malloc(size*sizeof*X);
			rvec *V = malloc(size*sizeof*V);
			rvec *F = malloc(size*sizeof*F);
			for(fCnt=0;(status = read_trr(xfp, nat, &step, &time, &lambda, xtc_box, x, v, f)) == exdrOK;){
				for (int i=0;i<size;i++){
					X[i][0] = x[(((int *)ptr)[i])][0];
					X[i][1] = x[(((int *)ptr)[i])][1];
					X[i][2] = x[(((int *)ptr)[i])][2];

					V[i][0] = v[(((int *)ptr)[i])][0];
					V[i][1] = v[(((int *)ptr)[i])][1];
					V[i][2] = v[(((int *)ptr)[i])][2];

					F[i][0] = f[(((int *)ptr)[i])][0];
					F[i][1] = f[(((int *)ptr)[i])][1];
					F[i][2] = f[(((int *)ptr)[i])][2];
				}
			}

			//housekeeping
			free(X);
			free(V);
			free(F);
		}

		//housekeeping
		xdrfile_close(xfp);
		free(x);
		free(v);
		free(f);
	}

	//be verbose
	printf("%d frames loaded containing %d atoms\n",
			sd->trj[trj_index]->sz,sd->sz);

	return fCnt;
}


//if the ptr is not NULL, then we have a subset of atoms we want to ignore
int dupe_trr(char *trrFile, void *ptr, int size){

	int nat, step, status, fCnt, output_freq = 20000;
	float time, lambda;

	status = read_trr_natoms(trrFile, &nat);
	XDRFILE *xfp = xdrfile_open(trrFile, "r");
	xfp_out = xdrfile_open("dupe.trr", "a");

	if (xfp && status == exdrOK){
		rvec *x = malloc(nat*sizeof*x);
		rvec *v = malloc(nat*sizeof*v);
		rvec *f = malloc(nat*sizeof*f);

		//when the ptr is NULL we process all atoms in the frame
		if (ptr == NULL)
			for(fCnt=0;(status = read_trr(xfp, nat, &step, &time, &lambda, xtc_box, x, v, f)) == exdrOK;){
				//process with the passed function
				if (fCnt % output_freq == 0)
					write_trr(xfp_out,nat,step,time,lambda,xtc_box,x,v,f);
				else
					write_trr(xfp_out,nat,step,time,lambda,xtc_box,x,v,NULL);
			}

		//otherwise, we pull out a subset of atoms for processing
		else{
			//alloc new rvec ptr for the subset
			rvec *X = malloc(size*sizeof*X);
			rvec *V = malloc(size*sizeof*V);
			rvec *F = malloc(size*sizeof*F);
			for(fCnt=0;(status = read_trr(xfp, nat, &step, &time, &lambda, xtc_box, x, v, f)) == exdrOK;){
				for (int i=0;i<size;i++){
					X[i][0] = x[(((int *)ptr)[i])][0];
					X[i][1] = x[(((int *)ptr)[i])][1];
					X[i][2] = x[(((int *)ptr)[i])][2];

					V[i][0] = v[(((int *)ptr)[i])][0];
					V[i][1] = v[(((int *)ptr)[i])][1];
					V[i][2] = v[(((int *)ptr)[i])][2];

					F[i][0] = f[(((int *)ptr)[i])][0];
					F[i][1] = f[(((int *)ptr)[i])][1];
					F[i][2] = f[(((int *)ptr)[i])][2];
				}
				//process with the passed function
				if (fCnt % output_freq == 0)
					write_trr(xfp_out,nat,step,time,lambda,xtc_box,X,V,F);
				else
					write_trr(xfp_out,nat,step,time,lambda,xtc_box,X,V,NULL);
			}

			//housekeeping
			free(X);
			free(V);
			free(F);
		}

		//housekeeping
		xdrfile_close(xfp_out);
		xdrfile_close(xfp);
		free(x);
		free(v);
		free(f);
	}
	else
		printf("Unable to load trajectory file (%s). Goodbye.\n\r",trrFile);

	return fCnt;
}


void parse_pdb_atom_record(char *rec, SYS_DATA *sd, int trj_index){

	char buf[1024]; float pos[3]; 
	TRJ_DATA *trj = sd->trj[trj_index]; 

	//we construct the atom_data array using the first frame of the first trj
	if ( trj->sz == 1 && trj_index == 0){

		//expand the atom_data array
		sd->atm = realloc(sd->atm,sizeof*sd->atm*++sd->sz);
		ATM_DATA *atom = &(sd->atm[sd->sz-1]);

		//allocate and initialize the residue id
		atom->rid = calloc(sizeof*atom->rid,10);

		//store a carbon copy of the input record
		atom->rec = strdup(rec);

		//pull the residue name
		char rname[10] = {'\0'};
		strncpy(rname,rec+17,3);

		//pull the residue number
		char rnumb[10] = {'\0'};
		strncpy(rnumb,rec+22,4);

		//construct the residue id
		sprintf(atom->rid,"%s_%d",rname,atoi(rnumb));

		//pull the atomic serial number (aka atom index)
		char index[10] = {'\0'};
		strncpy(index,rec+6,5);
		atom->ndx = atoi(index);
	}

	for (int j=0,i=30;i<54;i++,j++){

		buf[j+0] = rec[i]; 
		buf[j+1] = '\0';    

		switch (i){
			case (37): //X coordinates in Angstrom 
				pos[XX] = atof(buf);
				j=0; break;

			case (45): //Y coordinates in Angstrom 
				pos[YY] = atof(buf);
				j=0; break;

			case (53): //Z coordinates in Angstrom 
				pos[ZZ] = atof(buf);
				j=0; break;
			default:
				break;
		}
	}

	//append the position data to the frame
	FRM_DATA *f = trj->frm[trj->sz-1];
	f->pos = realloc(trj->frm[trj->sz-1]->pos,sizeof(float)*++f->sz*3);
	f->pos[(f->sz*3)-3+XX] = pos[XX];
	f->pos[(f->sz*3)-3+YY] = pos[YY];
	f->pos[(f->sz*3)-3+ZZ] = pos[ZZ];
}

void parse_pdb(void *trj, SYS_DATA *sd, int trj_index) {

	FILE *fp = (FILE *)trj;

	char line[1024]; int open_frame = 0;
	while(fgets(line,1024,(FILE *)fp)){

		//new frame signal
		if(strstr(line,"MODEL")==line)
			open_frame = new_frame(sd->trj[trj_index]);

		//crystal record (optional) 
		else if(strstr(line,"CYRS")==line)
			;//parse_pdb_box(line,sd);

		//atom record
		else if(strstr(line,"ATOM")==line)
			parse_pdb_atom_record(line,sd,trj_index);

		//end-of-frame signal
		else if(strstr(line,"ENDMDL")==line 
				|| strstr(line,"TER")==line)
			open_frame = 0;
	}

	//test the last line read prior to EOF
	if(strstr(line,"ENDMDL")==line 
			|| strstr(line,"TER")==line)
		open_frame = 0;

	//sanity check for potential error catching
	if (open_frame)
		printf("WARNING: pdb parsing halted within an open"
				" frame (frame #%d)\n'%s'\n",
				sd->trj[trj_index]->sz,line);
	else
		//be verbose
		printf("%d frames loaded containing %d atoms\n",
				sd->trj[trj_index]->sz,sd->sz);

	//housekeeping
	fclose((FILE *)fp);
}

// MAIN ML //-----------------------------------------------------------------|

void construct_problems(SYS_DATA *sd) {
	float scale_by = sd->ens_type == PDB_ENS ? 1.0 : 10.0;
	int fc = sd->trj[0]->sz;

	//for each atom in the atom_data array
	for(int i=0;i<sd->sz;i++) { 
		PRB_DATA *p = &(sd->atm[i].prb);

		//allocate memory for svn_node array
		p->x = malloc(sizeof*p->x*fc*2);	

		//allocate memory for training targets
		p->y = malloc(sizeof*p->y*fc*2);	

		//assign the # of training data points
		p->l = fc*2;			

		//loop over all frames in both trajectories
		for(int j=0;j<2;j++)
			for(int k=0;k<fc;k++){

				//init targets
				p->y[(fc*j)+k] = (!j) ? -1:1;

				//allocate memory an svm_node
				p->x[(fc*j)+k] = malloc(sizeof**p->x*4); 

				//loop over dims & make assignments
				for(int d=0;d<3;d++){
					p->x[(fc*j)+k][d].index = d+1;
					p->x[(fc*j)+k][d].value = sd->trj[j]->frm[k]->pos[(i*3)+d]*scale_by;
					//printf("%d %f\n", p->x[(fc*j)+k][d].index, p->x[(fc*j)+k][d].value);	
				}

				//last dim needs to be -1 (svmlib)
				p->x[(fc*j)+k][3].index = -1;
				p->x[(fc*j)+k][3].value = 0.f;
				//printf("%d %f\n", p->x[(fc*j)+k][3].index, p->x[(fc*j)+k][3].value);	
			}
	}
}


void yield_etas(SYS_DATA *sd) {
	for(int i=0;i<sd->sz;i++)
		sd->atm[i].eta = 1.0-svm_get_nr_sv(sd->atm[i].mdl) / 
			(2.0*(float)sd->trj[0]->sz);
}


void train_SVMs(SYS_DATA *sd){

	//map data in memory into svm_problem structures (libsvm)
	construct_problems(sd);
	puts("svm problem constructed successfully");

	//declare and init svm parameters, which do not change (libsvm) 
	struct svm_parameter param = { 
		C_SVC,RBF,3,GAMMA,0.f,1e+2,1e-3,COST,0,NULL,NULL,5e-1,1e-1,1,0};

	fprintf(stdout,"svm-training trajectory atoms with gamma = %f and C = %f...\n", GAMMA, COST);

#ifdef _OPENMP
	omp_set_num_threads(sd->inputs->threads);
	puts("svm training will be parallelized.");
#endif

	//train a SVM for each atom in the system
	int sz = sd->sz, i;
	ATM_DATA *atm=sd->atm;

#pragma omp parallel for schedule(dynamic) private(i) shared(sz,atm,param)
	for (i=0;i<sz;++i)
		atm[i].mdl = svm_train(&atm[i].prb,&param);
}

void avg_etas(SYS_DATA *sd) {
	if(sd->atm[0].rid == NULL) {
		puts("No residue data detected for average. Make sure that at least one input is a PDB file. Note that --print is not read for data, so will not provide the needed residues.");
		exit(1);
	}

	float avg_total = 0;
	int num_atoms = 0;
	float avg = 0;
	char *res_name = sd->atm[1].rid;
	int first = 0;

	for(int currentAtom = 0; currentAtom < sd->sz; currentAtom++) {
		if(strcmp(sd->atm[currentAtom].rid, res_name) == 0) {
			num_atoms++;
			avg_total += sd->atm[currentAtom].eta;
		}
		else {
			avg = avg_total/num_atoms;
			for(int i = first; i<currentAtom; i++) {
				sd->atm[i].avg_eta = avg;
			}
			first = currentAtom;
			num_atoms = 1;
			avg_total = sd->atm[currentAtom].eta;
			res_name = sd->atm[currentAtom].rid;
		}
	}
	avg = avg_total/num_atoms;
	for(int i = first; i<sd->sz; i++) {
		sd->atm[i].avg_eta = avg;
	}
}

// OUTPUT //------------------------------------------------------------------|

void splash(void){
	puts("\n /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*\\\n"
			" |*                                    *|\n"
			" |*         .-------.    ______        *|\n"
			" |*        /   o   /|   /\\     \\       *|\n"
			" |*       /_______/o|  /o \\  o  \\      *|\n"
			" |*       | o     | | /   o\\_____\\     *|\n"
			" |*       |   o   |o/ \\o   /o    /     *|\n"
			" |*       |     o |/   \\ o/  o  /      *|\n"
			" |*   jgs-'-------'     \\/____o/       *|\n"
			" |*                                    *|\n"
			" |*          DiCE version 0.5          *|\n"
			" |*                                    *|\n"
			" |*     written by: Aditya Gupta       *|\n"
			" |*        & Guy 'Wayyne' Dayhoff      *|\n"
			" |*                                    *|\n"
			" |*    University of South Florida     *|\n"
			" |*       Tampa, Fl - 2023-2024        *|\n"
			" |*                                    *|\n"
			" \\*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/\n");
}

void print_etas(SYS_DATA *sd) {
	FILE *fp = fopen("eta.dat","w");
	FILE *avg_fp = fopen("avg_eta.dat", "w");
	puts("\nprinting eta values to eta.dat");
	if(sd->inputs->avg) {
		puts("printing averaged eta values to avg_eta.dat");
		fprintf(avg_fp, "%5s  %7s  %s\n","atom#","residue","eta");
	}
	fprintf(fp, "%5s  %7s  %s\n","atom#","residue","eta");
	for(int i=0;i<sd->sz;i++) {
		if(sd->inputs->avg) {
			fprintf(avg_fp, "%5d  %7s  %5.3f\n",
					sd->atm[i].ndx,
					sd->atm[i].rid == NULL ? "UNKNOWN":sd->atm[i].rid,
					sd->atm[i].avg_eta);
		}
		fprintf(fp, "%5d  %7s  %5.3f\n",
				sd->atm[i].ndx,
				sd->atm[i].rid == NULL ? "UNKNOWN":sd->atm[i].rid,
				sd->atm[i].eta);
	}

	fclose(fp);
	fclose(avg_fp);
}

void write_pdb(SYS_DATA *sd) {
	char line[1024];
	char ending[15];

	FILE *outfile = fopen("writtenPDB.pdb", "w");
	puts("writing new file with etas to writtenPDB.pdb");
	FILE *avg_outfile;
	if(sd->inputs->avg) {
		puts("writing new file with averaged etas to writtenPDB_avg.pdb");
		avg_outfile = fopen("writtenPDB_avg.pdb", "w");
	}

	if(sd->inputs->full == 1) { //algorithm to copy the full file
		puts("\nf flag detected -- writing full file");
		int atom = 0;
		while(fgets(line, 1024, sd->inputs->printFile)) {
			if(strstr(line, "ATOM") != NULL) {
				//copy the end of the line after the ETA value
				for(int j = 66; j < 80; j++) {
					ending[j-66] = line[j];
				}
				if(sd->inputs->avg) {
					fprintf(avg_outfile, "%.60s %4.3f%.13s", line, sd->atm[atom%sd->sz].avg_eta, ending);
				}
				fprintf(outfile, "%.60s %4.3f%.13s", line, sd->atm[atom%sd->sz].eta, ending);
				atom++;
			}
			else {
				if(sd->inputs->avg) {
					fprintf(avg_outfile, "%s", line);
				}
				fprintf(outfile, "%s", line);
			}
		}
	}

	else {
		puts("\nno f flag -- writing one frame");
		int read = 0;
		char ending[14];
		int atom = 0;

		while(fgets(line, 1024, sd->inputs->printFile)) {
			if(strstr(line, "ATOM") == NULL) { //copy the top stuff b4 atom
				fprintf(outfile, "%s", line);
			}
			if(strstr(line, "ATOM") != NULL) {
				read = 1;
				for(int j = 66; j < 80; j++) {
					ending[j-66] = line[j];
				}
				if(sd->inputs->avg) {
					fprintf(avg_outfile, "%.60s %4.3f%.13s", line, sd->atm[atom].avg_eta, ending);
				}
				fprintf(outfile, "%.60s %4.3f%.13s", line, sd->atm[atom].eta, ending);
				atom++;

			}
			if(read == 1) {
				break;
			}
		}
		fprintf(outfile, "TER");
		fprintf(outfile, "ENDMDL");
		if(sd->inputs->avg) {
			fprintf(avg_outfile, "TER");
			fprintf(avg_outfile, "ENDMDL");
		}
	}
	fclose(outfile);
	if(sd->inputs->avg) {
		fclose(avg_outfile);
	}
}



// MAIN //--------------------------------------------------------------------|

int main(int argc, char **argv){

	//TO-DO wizard and help
	//allocate space for system data
	SYS_DATA *sd = new_system();
	parse_inputs(argc, argv, sd);
	args_checkout(sd);
	
	//populate the system data using the user-provided trajectories
	switch (sd->ens_type){

		case PDB_ENS:
			puts("PDB ENSEMBLES DETECTED");
			parse_pdb(sd->inputs->fileOne,sd,0);
			parse_pdb(sd->inputs->fileTwo,sd,1);
			break;

		case TRR_ENS:
			puts("TRR ENSEMBLES DETECTED");
			sd->status = read_trr_natoms(sd->inputs->fileNameOne, &(sd->sz));
			parse_trr(sd->inputs->fileOne,sd,0,NULL,0); 
			sd->status = read_trr_natoms(sd->inputs->fileNameTwo, &(sd->sz));
			parse_trr(sd->inputs->fileTwo,sd,1,NULL,0); 
			break;

		case XTC_ENS:
			puts("XTC ENSEMBLES DETECTED");
			sd->status = read_xtc_natoms(sd->inputs->fileNameOne, &(sd->sz));
			parse_xtc(sd->inputs->fileOne,sd,0,NULL,0); 
			sd->status = read_xtc_natoms(sd->inputs->fileNameTwo, &(sd->sz));
			parse_xtc(sd->inputs->fileTwo,sd,1,NULL,0); 
			break;

		default:
			//fall through case but args_checkout should never let us get here
			exit(1);
	}


	//do it to it
	train_SVMs(sd);
	yield_etas(sd);

	//execute on output inputs for final printing
	if(sd->inputs->avg == 1) {
		puts("\na flag detected -- etas will be copied and averaged");
		avg_etas(sd);
	}

	if(sd->inputs->print == 1) {
		puts("\np flag detected -- etas will be written in a PDB file format");
		write_pdb(sd);
	}

	print_etas(sd);
	//housekeeping
	free_system(sd);
	//good to go
	return 0;
}
