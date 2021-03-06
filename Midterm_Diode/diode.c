/*
PHYS370 intro computational physics
March 20th, 2018
Miderm Project Diode Semiconductor
David Jedynak, Jacob Connelly, Nick Zaler

Purpose:
	Model a Diode using the molecular dynamics sim to observe how changing band gap energy regions ~ (i.e. voltages) change the Voltage Current Characteristics. 

Features:
	Complete -Tested
	- Diode make from 2 fields that allow particles to travel in one direction, and 
	- Voltage source used to drive current through the diode
	- Real time Plotting of Current WRT time
	- A routine that sweeps voltage, measures current, and outputs the data to a file for VI curve plotting
	- Inputs for forward voltage gap, Reverse votlage gap, volatge sweep(start,end,dv)
	- GUI window for inputs 
	Complete - Untested
	- Resistor
	Incomplete
*/


#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#define Nmax 1000
#define D 2
#define MeasMax 100



//variables added for midterm project diode
double current_var = 0;				//global used for stoing the average current
double field_force = 0;				//for storing the force derived from diode,voltage source, or resistor
double g = 1,k = 1.1;				//gravity and electric field constants
double volt_len = 5,res_len = 5,diode_len_p = 5,diode_len_n = 5;//lengths of components
double volt_pos = 45,res_pos = 15,diode_pos_p = 22.5,diode_pos_n = 27.5;//positions of components
double resistance = 0;
double vs_v = -10,ud_v = 0.2,dd_v = -5.9;	//ideal voltages of diode band gaps
double v_f,ud_f,dd_f;				//forces for the different fields
//voltage sweep
double v_start = 5.0,v_end = -5.0,v_incr = -1.0;//voltage sweep
int v_sweep_delta = 500;			//delay for letting current settle out after a voltage change

double  L=50; 					// size of box

double mass[Nmax]; 				// masses of the particles
double x[Nmax][D],v[Nmax][D]; 			// State of the system
double q[Nmax];					// charges of the particles
double iq = 5;					//initial charge of particles

// parameters
double scalefac=100;
double x0[Nmax][D],v0[Nmax][D],dt=0.1,vv=0;	//initial conditions of system
double rho[MeasMax],Tset=0,Tmeas[MeasMax], ppnid[MeasMax],pp[MeasMax],Etot[MeasMax],Epot[MeasMax],Ekin[MeasMax],Imeas[MeasMax],IImeas[MeasMax];

int N=1000,Measlen=MeasMax,iterations=0;		//number of particles, length of values being taken for measurements

// Global variables for Isotherm
int Thermalize=10000, MeasNo=1000;

// 3d visualization
double tet=0,phi=0,tetdot=0.05,phidot=0,shift=150;

// functions defined for midterm
double getCurrent(double v[N][D]){
	double sumV = 0;
	for(int n = 0;n<N;n++){
		sumV += v[n][1];
		}
	return iq*sumV/L;
}



/*The set density function changes the value L (box length) in order to achieve a given density*/
void setdensity(){
  double fact=1/L;
  L=pow(N/rho[0],1./D);
  fact*=L;
  for (int n=0;n<N; n++)
    for (int d=0; d<D; d++)
      x[n][d]*=fact;
}

/*Calculate the temperature based off of the number of particles, kinetic energy, and the density*/
double density(){
  double r=0;
  r=N;
  for (int d=0; d<D; d++) r/=L;
  return r;
}
/*calculate the average Temperature by summing the kinetic energy, then dividing it by D dimensions and D number of particles*/
double T(double v[N][D]){
  //double
  double t=0;
  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++){
      t+=mass[n]*v[n][d]*v[n][d];
	}
  return t/N/D;
}
/*every time this is called the Temp gets set to whatever T is set to*/
void setTemp(){
  Tmeas[0]=T(v);
  double fact=sqrt(Tset/Tmeas[0]);
  for (int n=0;n<N; n++)
    for (int d=0; d<D; d++)
      v[n][d]*=fact;
}
/*Non ideal pressure*/
double Pnid(double x[N][D]){
  double p=0;
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
      double Fabs=12/dR12-6/dR6;
      p+=Fabs/D;
    }
  for (int d=0; d<D; d++) p/=L; 

  return p;
}
/*Potential Energy, calculated using the forces and the distances between each particle*/
double Ep(double x[N][D],double v[N][D]){
  double EE=0;

  for (int n=0; n<N; n++)
    for (int m=n+1; m<N; m++){
      double dr[D];
      double dR=0;
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
      EE += 1/dR12-1/dR6;
    }
  return 2*EE;
}
/*Calculates the force on each particle in each direction
-iterates though every particle calculates the distance between it and other particles*/
void F(double x[N][D], double v[N][D],double FF[N][D]){

  memset(&FF[0][0],0,N*D*sizeof(double)); //zeroize
  for (int n=0; n<N; n++) //iterate particles
    for (int m=n+1; m<N; m++){// 
      double dr[D],dR=0;
      for (int d=0; d<D; d++){ //iterate dimensions
	dr[d]=x[m][d]-(x[n][d]-L);
	double ddr;
	ddr=x[m][d]-(x[n][d]);
	if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr; // as the 
	ddr=x[m][d]-(x[n][d]+L);
	if (fabs(ddr)<fabs(dr[d])) dr[d]=ddr; //as the particle gets far away force levels out
	dR+=dr[d]*dr[d];
	
      }
      double dR6=dR*dR*dR;
      double dR12=dR6*dR6;
      double Fabs=12/dR12/dR-6/dR6/dR + k*q[n]*q[m]/pow(dR,2) + g*mass[n]*mass[m]/pow(dR,2);
      for (int d=0;d<D; d++){
	FF[n][d]-=Fabs*dr[d];
	FF[m][d]+=Fabs*dr[d];
      }
    }
  return;
}
/*primary function, updates dynamics for particles using above functions and also controls temp
-iterate through each velocity and integrate acceleration
-iterate position integrate velocity*/
void iterate(double x[N][D],double v[N][D],double dt){
  double ff[N][D];
  //itotal = 0;
  F(x,v,ff);
  if (iterations==0)
    for (int n=0;n<N;n++)
      for (int d=0;d<D;d++)
	v[n][d]+=0.5*ff[n][d]/mass[n]*dt;
  else
    for (int n=0;n<N;n++)
      for (int d=0;d<D;d++){
	//begin conditions for diode, resistor, and voltage source regions
	v[n][d]+=(ff[n][d])/mass[n]*dt;// integrate to get velocity
	
	if (d == 1){// if the dimension is in the Y dimension... vertical
		
		
		//check to see if the particle is in diode p, diode n, voltage, resistor regions, then apply the appropriate force
		if ((x[n][d] < (diode_pos_p+diode_len_p/2)) && (x[n][d] > (diode_pos_p-diode_len_p/2))&& v[n][d] < 0){//upward diode field.. if the particle has -v then force particle upwards
			field_force = ud_f - current_var;//simulate the forward diode band gap
			}
		else if((x[n][d] < (diode_pos_n+diode_len_n/2)) && (x[n][d] > (diode_pos_n-diode_len_n/2))&& (v[n][d] > 0)){//downward diode field... stops reverse bias current throgh the diode
			field_force = dd_f + current_var;//simulate the reverse diode band gap
			}
		else if((x[n][d] < (volt_pos + volt_len/2)) && (x[n][d] > (volt_pos - volt_len/2))){//voltage source field... provides electromotive force to propel particles throgh diode
			field_force = v_f;//simualate voltage source
			}
		else if((x[n][d] > (res_pos - res_len/2)) && (x[n][d] < (res_pos + res_len/2)) && abs(v[n][d]) > 0){//resistor field... just a resistor for fun
			field_force = -v[n][d]*q[n]*q[n]*resistance/(res_len*res_len);// derived from V=IR, V=E*dl, E=N/C
			}
		else{
			field_force = 0;	//particle is not in a component, then no force.	
		}
	v[n][d]+=(field_force)/mass[n]*dt;//integrate again to update velocity with the accelerations of the fields
	
	}
	else{
		field_force = 0;		
		}

	}

  current_var = getCurrent(v);		//update current value for measurements
  for (int n=0;n<N;n++)
    for (int d=0;d<D;d++){
      x[n][d]+=v[n][d]*dt;
      if (x[n][d]<0) x[n][d]+=L;
      else if (x[n][d]>=L) x[n][d]-=L;
    }
  setTemp();//added set temp to lower amout of mouse clicks for setting temperature.
  iterations++;
}

/*set initial velocities, masses*/
void setup(){
  int M=pow(N-1,1./D)+1;
  for (int n=0; n<N; n++){
    mass[n]=1;
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
void setVoltageSource(){
  //calculate forces for the given voltages
  v_f = iq*vs_v/volt_len;		//calculate the field from the set voltage
  ud_f = iq*ud_v/diode_len_p;		//calculate the field from the set upper diode voltage gap
  dd_f = iq*dd_v/diode_len_n;		//calculate the field from the lower diode votlage gap
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
    q[n] = iq;				//set initial charge for each electron
    }}

}
/*sets initial positions and velocities for particles*/
void init(){
  for (int n=0; n<N; n++){
    for (int d=0; d<D; d++){
      x[n][d]=x0[n][d];
      v[n][d]=v0[n][d];
    }
    //set initial charge for each electron
    q[n] = iq;
    
    
}
  iterations=0;
}
/*saves the current state as the initial state. for future use. pretty nice*/
void GetState(){
  for (int n=0; n<N; n++)
    for (int d=0; d<N; d++){
      x0[n][d]=x[n][d];
      v0[n][d]=v[n][d];
    }
  iterations=0;
}
/*2D graph of particles*/
void draw(int xdim, int ydim){
  int size=xdim;
  if (ydim<size) size=ydim;
  scalefac=size/L;
//code added for midterm start
//add lines for notatiting the different accelaration fields
// color, x1,y1,x2,y2

//voltage source lines
  mydrawline(2,0,scalefac*(L-volt_pos+volt_len/2),size,scalefac*(L-volt_pos+volt_len/2));
  mydrawline(2,0,scalefac*(L-volt_pos-volt_len/2),size,scalefac*(L-volt_pos-volt_len/2));

//diode field 1 lines
  mydrawline(3,0,scalefac*(L-diode_pos_p + diode_len_p/2),size,scalefac*(L-diode_pos_p + diode_len_p/2));
  mydrawline(3,0,scalefac*(L-diode_pos_p - diode_len_p/2),size,scalefac*(L-diode_pos_p - diode_len_p/2));

//diode field 2 lines
 // mydrawline(4,0,100,size,100);
  mydrawline(4,0,scalefac*(L-diode_pos_n - diode_len_n/2),size,scalefac*(L-diode_pos_n - diode_len_n/2));

//resistor field lines
  mydrawline(6,0,scalefac*(L-res_pos+res_len/2),size,scalefac*(L-res_pos+res_len/2));
  mydrawline(6,0,scalefac*(L-res_pos-res_len/2),size,scalefac*(L-res_pos-res_len/2));

//code added for midterm end

  mydrawline(1,0,size,size,size);
  mydrawline(1,size,0,size,size);
  for (int n=0; n<N; n++){
    int xx=x[n][0]*scalefac;
    int yy=x[n][1]*scalefac;
    myfilledcircle(n+1,xx,size-yy,0.5*scalefac);
  }
}
/**/
int compare(const void *x1,const void *x2){
  if (((double *) x1)[2]< ((double *)x2)[2]) return 1;
  else return -1;
}
/*3D graph of particles*/
void draw3d(int xdim, int ydim){
  int size=xdim;
  typedef struct part{
    double x[D]; // position
    int c; //color
  } part;
  part p[N];

  tet+=tetdot;								//increment the view angle
  phi+=phidot;								//increment the viewing velocity

  if (D!=3){
    printf("3d visualization requires D=3, you have D=%i!",D);
    return;
  }
  if (ydim<size) size=ydim;
  scalefac=0.7*size/L*shift;

  // center cube
  for (int n=0; n<N; n++)
    for (int d=0; d<D; d++){
      p[n].x[d]=x[n][d]-L/2;
      p[n].c=n+1; 
    }
  // rotate around y
  double c=cos(tet);
  double s=sin(tet);
  for (int n=0; n<N; n++){
    double x=c*p[n].x[0]+s*p[n].x[2];
    p[n].x[2]=-s*p[n].x[0]+c*p[n].x[2];
    p[n].x[0]=x;
  }
  // rotate around x
  c=cos(phi);
  s=sin(phi);
  for (int n=0; n<N; n++){
    double y=c*p[n].x[1]+s*p[n].x[2];
    p[n].x[2]=-s*p[n].x[1]+c*p[n].x[2];
    p[n].x[1]=y;
  }
  // shift box away from origin
  for (int n=0; n<N; n++){
    p[n].x[2]+=shift;
  }
  qsort(&p[0].x[0],N,sizeof(part),&compare);
  
  for (int n=0; n<N; n++){					//draw the circles with different sizes based on distance from observing point and colors
    int xx=p[n].x[0]/p[n].x[2]*scalefac;
    int yy=p[n].x[1]/p[n].x[2]*scalefac;
    myfilledcircle(0     ,xdim/2+xx,ydim/2-yy,0.5/p[n].x[2]*scalefac+2);
    myfilledcircle(p[n].c,xdim/2+xx,ydim/2-yy,0.5/p[n].x[2]*scalefac);
  }
}
/*store data from simulation for various dynamics*/
void Measure(){
  memmove(&rho[1],&rho[0],(MeasMax-1)*sizeof(double));
  rho[0]=density();
  memmove(&Tmeas[1],&Tmeas[0],(MeasMax-1)*sizeof(double));
  Tmeas[0]=T(v);
  //storing current values for graphing later
  memmove(&IImeas[1],&IImeas[0],(MeasMax-1)*sizeof(double));
  IImeas[0]=current_var;
  memmove(&ppnid[1],&ppnid[0],(MeasMax-1)*sizeof(double));
  ppnid[0]=Pnid(x);
  memmove(&pp[1],&pp[0],(MeasMax-1)*sizeof(double));
  pp[0]=rho[0]*Tmeas[0]+ppnid[0];
  memmove(&Ekin[1],&Ekin[0],(MeasMax-1)*sizeof(double));
  Ekin[0]=N*D*Tmeas[0];
  memmove(&Epot[1],&Epot[0],(MeasMax-1)*sizeof(double));
  Epot[0]=Ep(x,v);
  memmove(&Etot[1],&Etot[0],(MeasMax-1)*sizeof(double));
  Etot[0]=Epot[0]+Ekin[0];

}
/*isotherm file management / writing*/
void Isotherm(){
  FILE *res;
  char IsoName[100];

  sprintf(IsoName,"Iso%f_%i.dat",Tset,N);
  res=fopen(IsoName,"w");

  for (rho[0]=density(); rho[0]>0.01;rho[0]/=1.1){
    setdensity();
    // thermalize
    for (int i=0; i<Thermalize; i++){
      setTemp();
      iterate(x,v,dt);
    }
    // measure values
    double PP=0, TT=0;
    for (int i=0; i<MeasNo; i++){
      iterate(x,v,dt);
      double Tmeas=T(v);
      TT+=Tmeas;
      PP+=rho[0]*Tmeas+Pnid(x);
    }
    TT/=MeasNo;
    PP/=MeasNo;
    fprintf(res,"%e %e %e\n",rho[0],PP,TT);
    Events(1);
    DrawGraphs();
  }
  fclose(res);
}
/*isotherm routine*/
void Isotherms(){
  double LStart=L;
  for (Tset=0.05; Tset<2; Tset+=0.05){
    L=LStart;
    init();
    Isotherm();
  }
}

/*This routine sweeps the voltage from start to end with an increment, outputs a file with the voltage and current information to look at diode VI characteristics*/
void VI_Curve(){
  FILE *res;
  char IsoName[100];
  double sumI = 0;

  sprintf(IsoName,"VI_Curve_Q_%f_N_%i_ud_v_%f_dd_v_%f.dat",iq,N,ud_v,dd_v);//put info in the filename to distinguish test params
  res=fopen(IsoName,"w");
  //sweep the voltage
  if(v_incr < 0 ){
	  for (double ivolt = v_start;ivolt>v_end;ivolt = ivolt + v_incr){
		//set the voltage
		vs_v = ivolt;
		setVoltageSource();
		//loop to waste some time to get the current to settle
		for(int t = 0;t<v_sweep_delta;t++)
			{
			iterate(x,v,dt);
	       	 	Events(1);
			DrawGraphs();
			sumI += current_var;
			}
		//save the current and voltage data
		fprintf(res,"%e %f\n",sumI/v_sweep_delta,vs_v);
	   }
   }
     else{
	  for (double ivolt = v_start;ivolt<v_end;ivolt = ivolt + v_incr){
		//set the voltage
		vs_v = ivolt;
		setVoltageSource();
		//loop to waste some time to get the current to settle
		for(int t = 0;t<v_sweep_delta;t++)
			{
			iterate(x,v,dt);
	       	 	Events(1);
			DrawGraphs();
			sumI += current_var;
			}
		//save the current and voltage data
		fprintf(res,"%e %f\n",sumI/v_sweep_delta,vs_v);
	   }
   }

  fclose(res);
}



/**/
int main(){
  struct timespec ts={0,1000000};
  int cont=0;
  int sstep=0;
  int repeat=10;
  int done=0;
  char name[50],mname[N][50];
  setup();
  init();
  Measure();

  DefineGraphN_R("rho",&rho[0],&Measlen,NULL);
  DefineGraphN_R("T",&Tmeas[0],&Measlen,NULL);
  DefineGraphN_R("P",&pp[0],&Measlen,NULL);
  DefineGraphN_R("Pnid",&ppnid[0],&Measlen,NULL);
  DefineGraphN_R("E",&Etot[0],&Measlen,NULL);
  DefineGraphN_R("Epot",&Epot[0],&Measlen,NULL);
  DefineGraphN_R("Ekin",&Ekin[0],&Measlen,NULL);
  DefineGraphN_R("Average Current",&IImeas[0],&Measlen,NULL);
  
  AddFreedraw("Particles",&draw);
  AddFreedraw("Particles 3d",&draw3d);
  StartMenu("Newton",1);
  DefineMod("No part",&N,Nmax);
  DefineDouble("L",&L);
  DefineDouble("dt",&dt);
  StartMenu("measure",0);
  DefineDouble("Average Current",&IImeas[0]);
  DefineDouble("rho",&rho[0]);
  DefineFunction("setrho",setdensity);
  DefineDouble("T",&Tmeas[0]);
  DefineDouble("T set",&Tset);
  DefineFunction("set T",setTemp);
  DefineDouble("P",&pp[0]);
  DefineDouble("Pnid",&ppnid[0]);
  EndMenu();
  StartMenu("Isotherm",0);
  DefineInt("Thermalize",&Thermalize);
  DefineInt("MeasNo",&MeasNo);
  DefineFunction("Measure Isotherm",Isotherm);
  DefineFunction("Measure multiple Isotherms",Isotherms);
  EndMenu();
  StartMenu("init",0);
  for (int n=0; n<N; n++){
  }
  if (N<15)
    for (int n=0; n<N; n++){
      sprintf(mname[n],"Particle %i",n);
      StartMenu(mname[n],0);
      DefineDouble("m",&mass[n]);
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
  DefineFunction("Get State",&GetState);
  DefineFunction("init",&init);
  EndMenu();
  //menu code for diode midterm
  StartMenu("Diode",0);
  DefineDouble("Vs",&vs_v);
  DefineFunction("Update Diode Params",setVoltageSource);
  DefineDouble("forward gap V",&ud_v);
  DefineDouble("downward gap V",&dd_v);
  DefineDouble("resistance",&resistance);
  DefineDouble("charge",&iq);
  DefineDouble("avg current",&current_var);
  StartMenu("Voltage Sweep",0);
  DefineFunction("VI Curve",VI_Curve);
  DefineDouble("V start",&v_start);
  DefineDouble("V end",&v_end);
  DefineDouble("v_incr",&v_incr);
  DefineInt("v_sweep_delta",&v_sweep_delta);
  EndMenu();
  EndMenu();
  DefineGraph(curve2d_,"Measurements");
  DefineDouble("phi",&phi);
  DefineDouble("phidot",&phidot);
  DefineDouble("tet",&tet);
  DefineDouble("tetdot",&tetdot);
  DefineDouble("shift",&shift);
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
      Measure();
    }
    else       nanosleep(&ts,NULL);
  }
}


