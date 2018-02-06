// PHYS370 Introduction to Computational Physics
// Homework 1
// 1-29-18
// David Jedynak
 

#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#define pi  3.14159265358979323846264338327950288419716939937
#define N 3
#define D 2

//changelog
//1. added buttons and dynamics for initial conditions for all particle velocities and positions in sub menu init
//2. added mass factor in calculating forces
//3. add traces for particles to display there trajectories.. orbits... etc
//4. add in gravitational forces
//5. Extra: adapt this system to use actual data from the solar system. i.e. initial condition for our solar system

double q[N]; // charges of the particles
double x[N][D],v[N][D]; // State of the system
double mass[N]; //masses or the particles
double initial_velocity[N][D];
double initial_position[N][D];

//initial conditions for 3 stable orbiting particles
double iv1[3][2] = {{0,1.3228},{0,0},{0,-1.3228}};
double ip1[3][2] = {{-1,0},{0,0},{1,0}};
double im1[3] = {1,1,1};
double iq1[3] = {1,-2,1};

//initial conditions 
double iv2[2][2] =  {{0,1.41},{0,0}};
double ip2[2][2] = {{0,0},{1,0}};
double im2[2] = {1,1};
double iq2[2] = {-1,1};


// parameters
double scalefac=100;
double k=1,dt=0.1,g=1;

int points=100,iterations=0;

// sets the initial conditions to 2 or 3 orbiting, stable particles,
void init2(){

if(N == 3){
	for(int n = 0;n<N;n++){
		for(int d = 0;d<D;d++){
			initial_velocity[n][d] = iv1[n][d];
			initial_position[n][d] = ip1[n][d];
		}
	q[n] = iq1[n]; // charges of the particles
	mass[n] = im1[n]; //masses or the particles
	}
}
else if(N == 2){
	for(int n = 0;n<N;n++){
		for(int d = 0;d<D;d++){
			initial_velocity[n][d] = iv2[n][d];
			initial_position[n][d] = ip2[n][d];
		}
		q[n] = iq2[n]; // charges of the particles
		mass[n] = im2[n]; //masses or the particles
		}
	}
return;
}

void F(double x[N][D], double v[N][D],double FF[N][D]){

  memset(&FF[0][0],0,N*D*sizeof(double));//sets everything in array to zero

  for (int n=0; n<N; n++)
    for (int m=n+1; m<N; m++){
      double dr[D],dR=0;
      for (int d=0; d<D; d++){

        dr[d]=x[m][d]-x[n][d];
        dR+=dr[d]*dr[d];

      }
      double Fabs=k*q[n]*q[m]/pow(dR,2) + g*mass[n]*mass[m]/pow(dR,2);
       

      for (int d=0;d<D; d++){

	//equations governing the system dynamics
        FF[n][d]-=Fabs*dr[d]/mass[n]; //divided by mass for the particle to account for the particles mass.
        FF[m][d]+=Fabs*dr[d]/mass[m]; //decided to add this using a lagrangian
      }
    }
  return;
}

void iterate(double x[N][D],double v[N][D],double dt){
  double ff[N][D];
  int n;
  F(x,v,ff);
  if (iterations==0)
    for (int n=0;n<N;n++)
      for (int d=0;d<D;d++)
        v[n][d]+=0.5*ff[n][d]*dt;
  else
    for (int n=0;n<N;n++)
      for (int d=0;d<D;d++)
        v[n][d]+=ff[n][d]*dt;
 
  for (int n=0;n<N;n++)
    for (int d=0;d<D;d++)
      x[n][d]+=v[n][d]*dt;
  iterations++;
}

void CM(){
  double cm[D],cmv[D];
  memset(&cm[0],0,D*sizeof(double));
  memset(&cmv[0],0,D*sizeof(double));
  for (int d=0; d<D; d++){
    for (int n=0; n<N; n++){
      cm[d]+=x[n][d];
      cmv[d]+=v[n][d];
    }
    cm[d]/=N;
    cmv[d]/=N;
  }
  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++){
      x[n][d]-=cm[d];
      v[n][d]-=cmv[d];
    }

}

void init(){
  memset(&x[0][0],0,N*D*sizeof(double));
  memset(&v[0][0],0,N*D*sizeof(double));
  for (int n=0; n<N; n++){
    	for (int d=0; d<D; d++){
    x[n][d]= initial_position[n][d];
    v[n][d]= initial_velocity[n][d];
    	}
}
  iterations=0;
}

void draw(int xdim, int ydim){

  for (int n=0; n<N; n++){
    int xx=(x[n][0])*scalefac;
    int yy=(x[n][1])*scalefac;
    mycircle(n+1,xdim/2+xx,ydim/2-yy,10);
  }
}

int main(){
  struct timespec ts={0,100000000};
  int cont=0;
  int sstep=0;
  int repeat=1;
  int done=0;
  int n;
 

//GUI Components
  AddFreedraw("Particles",&draw);
  StartMenu("Newton",1);
  DefineDouble("dt",&dt);
  DefineDouble("k",&k);
  DefineDouble("g",&g);
  StartMenu("init",0);
  for (int n=0; n<N; n++){
  	DefineDouble("q",&q[n]);
  	DefineDouble("mass",&mass[n]);
	for(int d = 0; d<D;d++){
	DefineDouble("init pos",&initial_position[n][d]);
	}
	for(int d = 0; d<D;d++){
	DefineDouble("init vel",&initial_velocity[n][d]);
	}
}
  DefineFunction("init",&init);
  DefineFunction("CM",&CM);
  DefineFunction("init2",&init2);
  EndMenu();
  DefineDouble("scale",&scalefac);
  DefineInt("points",&points);
  DefineGraph(freedraw_,"graph");
  DefineGraph ( curve2d_,"Graph" ); 
  DefineInt("repeat",&repeat);
  DefineBool("sstep",&sstep);
  DefineLong("NS slow",&ts.tv_nsec);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();
  while (!done){
    Events(1);
    DrawGraphs();
    if (cont||sstep){
      sstep=0;
      for (int i=0; i<repeat; i++) iterate(x,v,dt);
      nanosleep(&ts,NULL);
    }
    else sleep(1);
  }
}
