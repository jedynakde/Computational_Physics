#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#define N 2
#define D 2

double q[N]; // charges of the particles
double x[N][D],v[N][D]; // State of the system

// parameters
double C[D],scalefac=100;
double k=1,x0=0,v0=0,dt=0.1;
int points=100,iterations=0;

void F(double x[N][D], double v[N][D],double FF[N][D]){

  memset(&FF[0][0],0,N*D*sizeof(double));
  for (int n=0; n<N; n++)
    for (int m=n+1; m<N; m++){
      double dr[D],dR=0;
      for (int d=0; d<D; d++){
        dr[d]=x[m][d]-x[n][d];
        dR+=dr[d]*dr[d];
      }
      double Fabs=k*q[n]*q[m]/pow(dR,1.5);
      for (int d=0;d<D; d++){
        FF[n][d]-=Fabs*dr[d];
        FF[m][d]+=Fabs*dr[d];
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
    x[n][0]=n;
    v[n][1]=0;
  }
  x[N-1][0]+=x0;
  v[N-1][1]+=v0;
  iterations=0;
}

void draw(int xdim, int ydim){

  for (int n=0; n<N; n++){
    int xx=(x[n][0]-C[0])*scalefac;
    int yy=(x[n][1]-C[1])*scalefac;
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
 
  AddFreedraw("Particles",&draw);
  StartMenu("Newton",1);
  DefineDouble("dt",&dt);
  DefineDouble("k",&k);
  StartMenu("init",0);
  DefineDouble("x0",&x0);
  DefineDouble("v0",&v0);
  for (int n=0; n<N; n++)
    DefineDouble("q",&q[n]);
  DefineFunction("init",&init);
  DefineFunction("CM",&CM);
  EndMenu();
  for (int d=0; d<D; d++){
    DefineDouble("C",&C[d]);
  }
  DefineDouble("scale",&scalefac);
  DefineInt("points",&points);
  DefineGraph(freedraw_,"graph");
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
