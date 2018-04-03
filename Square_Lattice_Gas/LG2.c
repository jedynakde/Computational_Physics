/*
Lattice Gas according to Frisch Hasslacher and Pomeau
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mygraph.h>
#include <math.h>

/* The following should not be necessary, but some older compilers don't support this yet */
#define Ob000001 1
#define Ob000010 2
#define Ob000100 4
#define Ob001000 8
#define Ob010000 16
#define Ob100000 32
#define Ob001001 9
#define Ob010010 18
#define Ob100100 36
#define Ob010101 21
#define Ob101010 42
#define Ob111111 63


#define NV 6  // 6 velocities, hexagonal lattice. 
#define LX 500
#define LY 500 // This represents a parallelogram

int n[LX][LY];
/* Meaning of the state: each bit of the integer implies a particle at the corresponding velocity. east: 1 north: 10 west:100 south:1000, (in binary) 
So one particle travling north and one particle travelling south would be 1010=9
 */
double R=10,frac=0.5,frac2=0.1;
/* now some fields to be displayed */
int rhoreq=0, ureq=0,coarse=1,lx=LX,ly=LY;
double rho[LX][LY],u[LX][LY][2]; // these could be integers, but the graphics is not yet set up for that. 

/* For analysis */
double Uxav[LY];
int Uxavreq=0;

void iterate(){
  static int k=153123;
  static int nn[LX][LY];
  /* colliding the particles */
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      switch (n[x][y]){
      case Ob001001: nn[x][y]= ((++k&1)==0)? Ob100100: Ob010010; break;
      case Ob010010: nn[x][y]= ((++k&1)==0)? Ob001001: Ob100100; break;
      case Ob100100: nn[x][y]= ((++k&1)==0)? Ob010010: Ob001001; break;
      case Ob010101: nn[x][y]= Ob101010; break;
      case Ob101010: nn[x][y]= Ob010101; break;
      default: nn[x][y]=n[x][y]; 
      }
    }
  /* now we need to move the particles */
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      n[x][y]=
	(nn[(x+1)%LX][y]& Ob001000)
	|(nn[(x+LX-1)%LX][y]& Ob000001)
	|(nn[x][(y+1)%LY]& Ob010000)
	|(nn[x][(y+LY-1)%LY]& Ob000010)
	|(nn[(x+LX-1)%LX][(y+1)%LY]& Ob100000)
	|(nn[(x+1)%LX][(y+LY-1)%LY]& Ob000100);
    }  
}


void init(){
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      if (pow((x+0.5*y)-(LX/2.+LY/4.),2)+pow(sqrt(0.75)*(y-LY/2),2)<R)
	n[x][y]=0;
      else{ 
	n[x][y]=((rand()<frac*RAND_MAX)?Ob000001:0)
	  |((rand()<frac*RAND_MAX)?Ob000010:0)
	  |((rand()<frac*RAND_MAX)?Ob000100:0)
	  |((rand()<frac*RAND_MAX)?Ob001000:0)
	  |((rand()<frac*RAND_MAX)?Ob010000:0)
	  |((rand()<frac*RAND_MAX)?Ob100000:0);
      }
    }
}

void initshear(){
  for (int x=0;x<LX;++x)
    for (int y=0;y<LY;++y){
      if (y>LY/2)
	n[x][y]=((rand()<frac*RAND_MAX)?Ob000001:0)
	  |((rand()<0.25*frac*RAND_MAX)?Ob010010:0)
	  |((rand()<0.25*frac*RAND_MAX)?Ob100100:0)
	  /*|((rand()<frac*RAND_MAX)?Ob001000:0)*/
	  |((rand()<0.25*frac*RAND_MAX)?Ob010010:0)
	  |((rand()<0.25*frac*RAND_MAX)?Ob100100:0);
      else 
	n[x][y]=/*(rand()<frac*RAND_MAX)?Ob000001:0)
		   |*/((rand()<0.25*frac2*RAND_MAX)?Ob010010:0)
		   |((rand()<0.25*frac2*RAND_MAX)?Ob100100:0)
		   |((rand()<frac2*RAND_MAX)?Ob001000:0)
		   |((rand()<0.25*frac2*RAND_MAX)?Ob010010:0)
		   |((rand()<0.25*frac2*RAND_MAX)?Ob100100:0);
		   
		   
    }
}

void GetGraphs(){
  int x,y;
  double *rp=&(rho[0][0]),*up=&(u[0][0][0]);
  lx=LX/coarse;
  ly=LY/coarse;
  if (rhoreq){
    rhoreq=0;
    memset(rp,0,lx*ly*sizeof(double));
    for (x=0; x<LX;++x)
      for (y=0; y<LY; ++y)
	rp[x/coarse*ly+y/coarse]
	  +=(n[x][y]& Ob000001)
	  +((n[x][y]& Ob000010)>>1)
	  +((n[x][y]& Ob000100)>>2)
	  +((n[x][y]& Ob001000)>>3)
	  +((n[x][y]& Ob010000)>>4)
	  +((n[x][y]& Ob100000)>>5);
  }
  if (ureq){
    ureq=0;
    memset(up,0,lx*ly*2*sizeof(double));
    for (x=0; x<LX;++x)
      for (y=0; y<LY; ++y){
	up[(x/coarse*ly+y/coarse)*2  ]+=
	  (n[x][y]& Ob000001)
	  +0.5*(((n[x][y]& Ob000010)>>1)-((n[x][y]& Ob000100)>>2))
	  -((n[x][y]& Ob001000)>>3)
	  -0.5*(((n[x][y]& Ob010000)>>4)-((n[x][y]& Ob100000)>>5));    
	up[(x/coarse*ly+y/coarse)*2+1]+= sqrt(0.75)*
	  (((n[x][y]& Ob000010)>>1)+((n[x][y]& Ob000100)>>2)
	  -((n[x][y]& Ob010000)>>4)-((n[x][y]& Ob100000)>>5));    
      }
  }
  if (Uxavreq){
    memset(Uxav,0,ly*sizeof(double));
      for (y=0; y<LY; ++y){
	for (x=0; x<LX;++x)
	  Uxav[y/coarse]+=
	    (n[x][y]& Ob000001)
	    +0.5*(((n[x][y]& Ob000010)>>1)-((n[x][y]& Ob000100)>>2))
	    -((n[x][y]& Ob001000)>>3)
	    -0.5*(((n[x][y]& Ob010000)>>4)-((n[x][y]& Ob100000)>>5));    	
      }
      for (y=0;y<ly;y++){
	Uxav[y]/=LX*coarse;
      }
  }
}

int main (){
  int Paused=1, Step=1, Repeat=1, done=0;

  init();
  DefineGraphN_R("Ux av",&Uxav[0],&ly,&Uxavreq);
  DefineGraphNxN_R("Density",&(rho[0][0]),&lx,&ly,&rhoreq);
  DefineGraphNxN_RxR("velocity",&(u[0][0][0]),&lx,&ly,&ureq);

  StartMenu("Hexagonal Lattice Gas",1);
  DefineDouble("R^2",&R);
  DefineDouble("frac",&frac);
  DefineDouble("frac2",&frac2);
  DefineFunction("init circ",&init);
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