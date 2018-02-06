#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#define Nmax 1000
#define D 2
double  L=50; // size of box
double tm = 0;//total mass of system
double density = 0;
double ke = 0;//kinetic energy of the particles
double bc = 1.38064852E-23;//boltzmann constant
double pressure = 0;
int counter = 0;
//double pressure_plot_data[20000];//for collecting pressure data to plot ~ 200 seconds of data
double temperature = 0;

//typedef struct part {double xx; double yy;;} part;
//typedef struct TX{double tt; double xx;} TX;
//TX pdp[20000];

double m[Nmax]; // charges of the particles
double mass[Nmax]; //particle masses
double x[Nmax][D],v[Nmax][D]; // State of the system
int ItNo = 20000;
// parameters
double C[D],scalefac=100;
double k=1,x0[Nmax][D],v0[Nmax][D],dt=0.01,vv=1;
double mass0 = 2.001;
int N=Nmax,points=100,iterations=0;
double pdp[20000];
//int time[20000];

//function to determine the density, temperature, and pressure of the system
void density_function(){
density = 0;
tm = 0;
ke = 0;
//sum up all masses to find the total system mass
for (int n = 0; n < N;n++){
	tm = tm + m[n];
	for(int d = 0; d < D;d++){
		//calculate the total kinetic energy of the particles in the sytem
		ke = ke	+ .5*mass[n]*v[n][d]*v[n][d];	
	}
	
}
//use the box size to calculate density
density = (double) tm / (L*L);
//pressure and temperature of the system
temperature = (double)(2/(3*bc))*ke;//in Kelvin
pressure = ((2/(3.000))*N*ke/(L*L));//changed 3 to 3.000 to avoid pressur being changed to zero
//store last 200 seconds of pressure data for plotting.
pdp[counter%20000] = pressure;
//time[counter%20000] = counter;
counter = (counter + 1);

return;
}

void F(double x[N][D], double v[N][D],double FF[N][D]){

  memset(&FF[0][0],0,N*D*sizeof(double));
  for (int n=0; n<N; n++)
    for (int m=n+1; m<N; m++){
      double dr[D],dR=0;
      for (int d=0; d<D; d++){
	dr[d]=x[m][d]-(x[n][d]-L);
	double ddr;
	ddr=x[m][d]-(x[n][d]);
	if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr;
	ddr=x[m][d]-(x[n][d]+L);
	if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr;
	dR+=dr[d]*dr[d];
      }
      double dR6=dR*dR*dR;
      double dR12=dR6*dR6;
      double Fabs=12/dR12/dR-6/dR6/dR;
      for (int d=0;d<D; d++){
	FF[n][d]-=Fabs*dr[d]/mass[n];
	FF[m][d]+=Fabs*dr[d]/mass[n];
      }
    }
  return;
}

void iterate(double x[N][D],double v[N][D],double dt){
  double ff[N][D];
  F(x,v,ff);
  density_function();
  if (iterations==0)
    for (int n=0;n<N;n++)
	
      for (int d=0;d<D;d++)
	v[n][d]+=0.5*ff[n][d]/m[n]*dt;
  else
    for (int n=0;n<N;n++)
      for (int d=0;d<D;d++)
	v[n][d]+=ff[n][d]/m[n]*dt;
  
  for (int n=0;n<N;n++)
    for (int d=0;d<D;d++){
      x[n][d]+=v[n][d]*dt;
      if (x[n][d]<0) x[n][d]+=L;
      else if (x[n][d]>=L) x[n][d]-=L;
    }
  iterations++;
}

void CM(){
  double cm[D],cmv[D],M=0;
  memset(&cm[0],0,D*sizeof(double));
  memset(&cmv[0],0,D*sizeof(double));
  for (int n=0; n<N; n++) M+=m[n]; 
  for (int d=0; d<D; d++){
    for (int n=0; n<N; n++){
      cm[d]+=m[n]*x[n][d];
      cmv[d]+=m[n]*v[n][d];
    }
    cm[d]/=M;
    cmv[d]/=M;
  }
  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++){
      x[n][d]-=cm[d]+L/2;
      v[n][d]-=cmv[d]+L/2;
    }

}

void setup(){
  int M=pow(N,1./D)+1;
  for (int n=0; n<N; n++){
    m[n]=1;
    for (int d=0; d<D; d++){
      int nn=n;
      for (int dd=0; dd<d; dd++) nn/=M;
      x0[n][d]=(nn%M)*L/M;
      if (d==1){
	if (x0[n][0]<L/2)
	  v0[n][d]=vv;
	else v0[n][d]=-vv;
      }
      else v0[n][d]=0;
    }
  }
}

void init(){
  for (int n=0; n<N; n++){
    x[n][0]=x0[n][0];
    x[n][1]=x0[n][1];
    v[n][0]=v0[n][0];
    v[n][1]=v0[n][1];
//populate the mass array
	for (int n = 0; n < N;n++ ){
		mass[n] = mass0; 
	}
  }
  iterations=0;
}

void draw(int xdim, int ydim){
  int size=xdim;
  if (ydim<size) size=ydim;
  scalefac=size/L;

  mydrawline(1,0,size,size,size);
  mydrawline(1,size,0,size,size);
  for (int n=0; n<N; n++){
    int xx=x[n][0]*scalefac;
    int yy=x[n][1]*scalefac;
    myfilledcircle(n+1,xx,size-yy,0.5*scalefac);
  }
}


int main(){
  struct timespec ts={0,0};
  int cont=0;
  int sstep=0;
  int repeat=10;
  int done=0;
  char name[50],mname[N][50];

  setup();
  init();
  
  AddFreedraw("Particles",&draw);
  DefineGraphN_R("pressure",&pdp[0],&ItNo,NULL);
  StartMenu("Newton",1);
  DefineDouble("L",&L);
  DefineDouble("dt",&dt);
  DefineDouble("k",&k);
  StartMenu("init",0);
  for (int n=0; n<N; n++){
  }
  if (N<15)
    for (int n=0; n<N; n++){
      sprintf(mname[n],"Particle %i",n);
      StartMenu(mname[n],0);
      DefineDouble("m",&m[n]);
      for (int d=0; d<D; d++){
	sprintf(name,"x[%i]",d);
	DefineDouble(name,&x0[n][d]);
      }
      for (int d=0; d<D; d++){
	sprintf(name,"v[%i]",d);
	DefineDouble(name,&v0[n][d]);
      }
      EndMenu();
    }
  DefineDouble("vv",&vv);
  DefineFunction("setup",&setup);
  DefineFunction("init",&init);
//function to find and print density to terminal
  DefineFunction("density_function",&density_function);
  DefineFunction("CM",&CM);
  EndMenu();
  for (int d=0; d<D; d++){
    DefineDouble("C",&C[d]);
  }
  DefineDouble("scale",&scalefac);
  DefineDouble("density",&density);
  DefineDouble("pressure",&pressure);
  DefineDouble("tempertature",&temperature);
  DefineDouble("mass0",&mass0);
  DefineInt("points",&N);
  DefineGraph(freedraw_,"graph");
  DefineGraph ( curve2d_, "Graph" );
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

