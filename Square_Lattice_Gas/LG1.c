/*
Lattice Gas according to, Hardy, Pomeau and de Pazzis (1973 and 1976)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mygraph.h>
#include <math.h>

/* The following should not be necessary, but some older compilers don't support this yet */
#define Ob0001 1
#define Ob0010 2
#define Ob0100 4
#define Ob1000 8
#define Ob0101 5
#define Ob1010 10
#define Ob1111 15
#define Ob1110 14
#define Ob1011 11

#define NV 4  // 4 velocities, square lattice, not isotropic. 
#define LX 1000
#define LY 1000

int n[LX][LY];
/* Meaning of the state: each bit of the integer implies a particle at the corresponding velocity. east: 1 north: 10 west:100 south:1000, (in binary) 
So one particle travling north and one particle travelling south would be 1010=9
 */
double R=10;
/* now some fields to be displayed */
int rhoreq=0, ureq=0,coarse=1,lx=LX,ly=LY;
double rho[LX][LY],u[LX][LY][2]; // these could be integers, but the graphics is not yet set up for that. 

/* For analysis */
double frac=0.5;
double Uxav[LY];
int Uxavreq=0;


void iterate(){
static  int nn[LX][LY];
  /* colliding the particles */
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      switch (n[x][y]){
      case Ob0101: nn[x][y]= Ob1010; break;
      case Ob1010: nn[x][y]= Ob0101; break;
      default: nn[x][y]=n[x][y]; 
      }
    }
  /* now we need to move the particles */
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      n[x][y]=
	(nn[(x+1)%LX][y]& Ob0100)
	+(nn[(x+LX-1)%LX][y]& Ob0001)
	+(nn[x][(y+1)%LY]& Ob1000)
	+(nn[x][(y+LY-1)%LY]& Ob0010);
    }  
}


void initcirc(){
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      if (pow(x-LX/2,2)+pow(y-LY/2,2)<R)
	n[x][y]=0;
      else 
	n[x][y]=rand()& Ob1111;
    }
}

void initshear(){
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      if (y<LY/2)
	n[x][y]=random()& Ob1011;
      else {
	n[x][y]=random()& Ob1110;
	if ((random()<RAND_MAX*frac)) n[x][y]=0;
      }
    }
}

void GetGraphs(){
  double *rp=&(rho[0][0]),*up=&(u[0][0][0]);
  lx=LX/coarse;
  ly=LY/coarse;
  if (rhoreq){
    rhoreq=0;
    memset(rp,0,lx*ly*sizeof(double));
    for (int x=0; x<LX;++x)
      for (int y=0; y<LY; ++y)
	rp[x/coarse*ly+y/coarse]
	  +=(n[x][y]& Ob0001)
	  +((n[x][y]& Ob0010)>>1)
	  +((n[x][y]& Ob0100)>>2)
	  +((n[x][y]& Ob1000)>>3);
  }
  if (ureq){
    ureq=0;
    memset(up,0,lx*ly*2*sizeof(double));
    for (int x=0; x<LX;++x)
      for (int y=0; y<LY; ++y){
	up[(x/coarse*ly+y/coarse)*2  ]+=
	  (n[x][y]& Ob0001)-((n[x][y]& Ob0100)>>2);    
	up[(x/coarse*ly+y/coarse)*2+1]+=
	  ((n[x][y]& Ob0010)>>1)-((n[x][y]& Ob1000)>>3);    
      }
  }
  if (Uxavreq){
    Uxavreq=0;
    memset(Uxav,0,ly*sizeof(double));
    for (int x=0; x<LX;++x)
      for (int y=0; y<LY; ++y){
	Uxav[y/coarse]+=
	  (n[x][y]& Ob0001)-((n[x][y]& Ob0100)>>2);    
      }
    for (int y=0; y<ly; ++y){
      Uxav[y]/=LX*coarse;
    }
  }
}

int main (){
  int Paused=1, Step=1, Repeat=1, done=0;

  initcirc();
  DefineGraphN_R("Ux av",&Uxav[0],&ly,&Uxavreq);
  DefineGraphNxN_R("Density",&(rho[0][0]),&lx,&ly,&rhoreq);
  DefineGraphNxN_RxR("velocity",&(u[0][0][0]),&lx,&ly,&ureq);

  StartMenu("Square Lattice Gas",1);
  DefineDouble("R^2",&R);
  DefineDouble("frac",&frac);
  DefineFunction("init circ",&initcirc);
  DefineFunction("init shear",&initshear);
  DefineInt("Coarsegrain",&coarse);
  DefineGraph(contour2d_,"Density plot");
  DefineGraph(curve2d_,"Uav graph");
  DefineInt("Repeat",&Repeat);
  DefineBool("Step",&Step);
  DefineBool("Paused",&Paused);
  DefineBool("done",&done);
  EndMenu();
  while (!done){
    Events(1);
    GetGraphs();
    DrawGraphs();
    if (!Paused || !Step){
      Step=1;
      for (int i=0;i<Repeat;i++){
	iterate();
      }
    }
    else sleep(1);
  }
}