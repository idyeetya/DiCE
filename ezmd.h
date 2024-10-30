//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////                                                     
//
//                    ,----,        ____                
//     ,---,.       .'   .`|      ,'  , `.    ,---,     
//   ,'  .' |    .'   .'   ;   ,-+-,.' _ |  .'  .' `.   
// ,---.'   |  ,---, '    .',-+-. ;   , ||,---.'     ,  
// |   |   .'  |   :     ./,--.'|'   |  ;||   |  .`\  | 
// :   :  |-,  ;   | .'  /|   |  ,', |  '::   : |  '  | 
// :   |  ;/|  `---' /  ; |   | /  | |  |||   ' '  ;  : 
// |   :   .'    /  ;  /  '   | :  | :  |,'   | ;  .  | 
// |   |  |-,   ;  /  /--,;   . |  ; |--' |   | :  |  ' 
// '   :  ;/|  /  /  / .`||   : |  | ,    '   : | /  ;  
// |   |    \./__;       :|   : '  |/     |   | '` ,/   
// |   :   .'|   :     .' ;   | |`-'      ;   :  .'     
// |   | ,'  ;   |  .'    |   ;/          |   ,.'       
// `----'    `---'        '---'           '---'      
//
//        + EZMD Version 1.0
//        + Copyright 2019-2020 Guy (Wayyne) Dayhoff*
//                
//////////////////////////////////////////////////////////                                                     
///////////////// *University South Florida, Tampa Fl, USA
//////////////////////////////////////////////////////////                                                     


#ifndef _EZMD
#define _EZMD

/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"

typedef struct histogram  HISTOGRAM;
typedef struct bin        BIN;

struct bin {
  double  sum;        //running total value of the bin
  int     cnt;        //running count of bin contributions
};

struct histogram {
  double   *bins;
  float     bcnt;        //number of bins in the histogram
  float   bwidth;        //width of the bins in the histogram
  float     tcnt;        //number of accumulated tallies
};

/*booleans*/
#define bool int
#define TRUE 1
#define FALSE 0

/*global variables*/
extern matrix xtc_box;
extern HISTOGRAM *H;

/*global declarations*/
bool fread_int(int *, FILE *);
int *parse_iil(char *, int *);
void rdf_from_xtc(char **);
void process_trr_frame(int nat,rvec *x,rvec *v, rvec *f,int fIndex);
int main_2(int,char**);

#endif
