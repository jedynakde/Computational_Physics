// Lattice gas code in 1d.
//change the initial zero point for the Ux ~ mean particle velocity

//what about dividing by total number of particles for the theoretical?
//zero points seem to be really small
//why cant we just change the probability of collisions to make the parabola fit
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <mygraph.h>
#include <string.h>

#define xdim 100
#define ydim 100
#define V 9      //number of velocities
#define MeasMax 200
int XDIM=xdim,YDIM=ydim;
int C=300;//chaneged to 300 for viscocity...
double N[xdim][ydim],NU[xdim][ydim][2];
int Nreq=0,NUreq=0;
int n[xdim][ydim][V];

// Boundary variables
#define LINKMAX 10000
#define LINKMAX_DYNAMIC 10000



double particle_flip_w = 0.032;

double viscocity = 0;

double dt = 1.0;

int x_link_var = 0;
int y_link_var = 1;
int type_link = 0;

//variables for measuring tube momentum
double last_1 = 0;
double last_2 = 0;


//particle source parameters
int src_x = 85,src_y = 85,src_len = 5,src_den = 10;
int src_x_2 = 15,src_y_2 = 15,src_len_2 = 5,src_den_2 = 10;
int var1 =10;
//variables for measuring values
//particle velocity
int particle_vx = 0,particle_vy,particle_leakage = 0;
double measure_particle_vx[MeasMax],measure_particle_vy[MeasMax];
double measure_particle_vx_filt[MeasMax],measure_particle_vy_filt[MeasMax];
int vx_m = 0,vy_m = 0;
int d1 = 0,d2 = 50;
int ux_zero_flag = 0;
double ux_zero_point[ydim];
//int ux_zero_curve_1[ydim];

//link variables
int linkcount_dynamic = 0,links_dynamic[2][LINKMAX_DYNAMIC][3];//link x y and v(for particles) added other index for right and left walls
double links_dynamic_pool[2][LINKMAX_DYNAMIC];
double links_dynamic_velocity[LINKMAX_DYNAMIC][2];//link velocities
double wall_mass = 75000.1;
int yy3 = 26,yy4 = 75,dynamic_walls_on = 1,dynamic_wall_control_on = 1;
double flow = 0;
//dynamic wall
// position
//double x_shift = 50;
double dynamic_wall_position_x = 0,dynamic_wall_position_y = 0;
//velocity
double dynamic_wall_vx = 0,dynamic_wall_vy = 0;
//momentum
int dynamic_wall_momentum_x = 0,dynamic_wall_momentum_y = 0;
//dynamic wall position for graphing
double measure_dynamic_wall_position_x[MeasMax],measure_dynamic_wall_position_y[MeasMax];
//filtered 
double measure_dynamic_wall_position_x[MeasMax],measure_dynamic_wall_position_y[MeasMax];
//
double measure_dynamic_wall_position_x_filt[MeasMax],measure_dynamic_wall_position_y_filt[MeasMax];
//dynamic wall velocity
double measure_dynamic_wall_vx[MeasMax],measure_dynamic_wall_vy[MeasMax];
//filtered 
double measure_dynamic_wall_vx_filt[MeasMax],measure_dynamic_wall_vy_filt[MeasMax];
//dynamic wall momentum
double measure_dynamic_wall_momentum_x[MeasMax],measure_dynamic_wall_momentum_y[MeasMax];
//filtered
double measure_dynamic_wall_momentum_x_filt[MeasMax],measure_dynamic_wall_momentum_y_filt[MeasMax];


double measure_particle_velocity_front[MeasMax];
double measure_particle_velocity_front_last[MeasMax];
double measure_particle_force_front[MeasMax];
double theoretical_particle_vx_front[MeasMax];
//static wall variables
int linkcount=0,links[LINKMAX][3];

//points for drawing walls
int x0=0,x1=100,yy0 = 25,yy1 = 74;

//static wall momentum
double static_wall_momentum_x = 0,static_wall_momentum_y = 0;
double measure_static_wall_momentum_x[MeasMax],measure_static_wall_momentum_y[MeasMax];
//filtered data
double measure_static_wall_momentum_x_filt[MeasMax],measure_static_wall_momentum_y_filt[MeasMax];
//double momentum_est_x[MeasMax],momentum_est_y[MeasMax];
//double momentum_est_x_filt[MeasMax],momentum_est_y_filt[MeasMax];


int MeasLen = MeasMax/2;
int range_val = 20;
int val[] = {0,0};
int filt_data[] ={0,0,0,0,0,0,0,0,0,0,0}; //dynamic walls -> x,y,px,py,vx,vy  static walls -> px py particles -> vx vy leakage


//added a running average to smooth out the graphs
void average(int range){
	//val[0] = 0;
	//val[1] = 0;
	//clear previous filtered data
	for(int i = 0; i < 11;i++){
		filt_data[i] = 0;
		}
	//update with new data
	for(int i = 0; i < range;i++){
		filt_data[0] += measure_dynamic_wall_position_x[i];
		filt_data[1] += measure_dynamic_wall_position_y[i];
		filt_data[2] += measure_dynamic_wall_momentum_x[i];
		filt_data[3] += measure_dynamic_wall_momentum_y[i];
		filt_data[4] += measure_dynamic_wall_vx[i];
		filt_data[5] += measure_dynamic_wall_vy[i];
		filt_data[6] += measure_static_wall_momentum_x[i];
		filt_data[7] += measure_static_wall_momentum_y[i];
		filt_data[8] += measure_particle_vx[i];
		filt_data[9] += measure_particle_vy[i];
		filt_data[10] += measure_particle_leakage[i];
		}
	//average values
	for(int i = 0; i < 11;i++){
		filt_data[i] = filt_data[i]/range;
		}
}

void measure_function(){
	particle_vx = 0;
	particle_vy = 0;
	//int total_particles = 0;
	int total_particles = 0;
	double shift[ydim];
	for(int i = 0; i < YDIM -1;i++){
		total_particles = 0;
		measure_particle_velocity_front[i] = 0;
		for(int x0 = 0;x0<xdim;x0++){
			measure_particle_velocity_front[i] += n[x0][i][5]+n[x0][i][2]+n[x0][i][8] - n[x0][i][3]-n[x0][i][0]-n[x0][i][6];			for(int v = 0;v<9;v++) total_particles += n[x0][i][v];	
			}
		//find average particle velocity
		if(total_particles > 0){
			measure_particle_velocity_front[i] = (double)(measure_particle_velocity_front[i]/total_particles);//+shift[i];
			theoretical_particle_vx_front[i] = (double)(particle_flip_w*6.0/9.)*(i-25)*(i-75);

			}
		else{
			measure_particle_velocity_front[i] = 0;
			theoretical_particle_vx_front[i] = 0;
			}
		//setting zero point for ux
		if(ux_zero_flag == 0)ux_zero_point[i] += measure_particle_velocity_front[i];
		if(ux_zero_flag == 1)ux_zero_point[i] = (measure_particle_velocity_front[i] + ux_zero_point[i])/2; 
		printf("zp = %f \n",ux_zero_point[i]);
			
		measure_particle_force_front[i] = (measure_particle_velocity_front[i] - measure_particle_velocity_front_last[i]);		measure_particle_velocity_front_last[i] = measure_particle_force_front[i];
		
		}


		//set flags after ux at iteration 0 and 1 have been recorded
		printf("flag = %i \n",ux_zero_flag);
		if(ux_zero_flag == 1) 
			{ux_zero_flag = 2;
			//partical_flip_w = f_partical_flip_w;
			for(int i = 0;i<ydim;i++)shift[i]=ux_zero_point[i];
			}
		if(ux_zero_flag == 0)ux_zero_flag = 1;


	particle_vx = particle_vx/100.0;
	particle_vy = particle_vy/100.0;


	
	//measure the number of particles leaking
	  
	  particle_leakage = 0;
		for(int x = 0; x < xdim; x++){
			for(int y = 0;y < ydim;y++){
				for(int v = 0; v < 9;v++){
					particle_leakage += n[x][y][v];
					}
				}
			}
}


//take measurements for plotting data
void Measure(){
 //measure the force on the particles in the tube.
  measure_function();
}

void moveParticles(){
for(int x = 0;x<xdim;x++){
	for(int y = 0;y<ydim;y++){
		//additional code for flipping some horizontal particles
		int flip_parts = particle_flip_w*((double)rand()/RAND_MAX)*n[x][y][5];
		n[x][y][3] += flip_parts;//n[x][y][3] + flip_parts;
		n[x][y][5] -= flip_parts;//n[x][y][5] - flip_parts;
		flip_parts = particle_flip_w*((double)rand()/(double)RAND_MAX)*n[x][y][8];
	
		n[x][y][6] += flip_parts;
		n[x][y][8] -= flip_parts;
	
		flip_parts = particle_flip_w*((double)rand()/(double)RAND_MAX)*n[x][y][2];
		n[x][y][0] += flip_parts;
		n[x][y][2] -= flip_parts;
		

		}
	}
}

void FindLink(){


  //horizontal walls
  for (int x=x0; x<x1; x++){
	x = ((x%XDIM)+XDIM)%XDIM;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 0;
    	linkcount++;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 1;
    	linkcount++;
    	links[linkcount][0] = x; //x-position 
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 2;
    	linkcount++;
  }
    for (int x=x0; x<x1; x++){
	x = ((x%XDIM)+XDIM)%XDIM;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy1;
    	links[linkcount][2] = 0;
    	linkcount++;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy1;
    	links[linkcount][2] = 1;
    	linkcount++;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy1;
    	links[linkcount][2] = 2;
    	linkcount++;
  }

}

void bounceback(){
  dynamic_wall_momentum_x = 0;
  dynamic_wall_momentum_y = 0;
  static_wall_momentum_x = 0;
  static_wall_momentum_y = 0;
  

  
  for (int lc=0; lc<linkcount; lc++){
    //quantity of partices
   
    int x=links[lc][0];
    int y=links[lc][1];
    int v=links[lc][2];
    int x_b = (((x)%XDIM)+XDIM)%XDIM;
    int y_b = (((y)%YDIM)+YDIM)%YDIM;
    //velocity
    int vx=v%3-1;
    int vy=1-v/3;
    int x_v_b = (((vx+x)%XDIM)+XDIM)%XDIM;
    int y_v_b = (((vy+y)%YDIM)+YDIM)%YDIM;
    int tmp = n[x_v_b][y_v_b][v];
    //summing all momemtums
    static_wall_momentum_x += -2*vx*(n[x][y][8-v]-tmp);
    static_wall_momentum_y += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[x_v_b][y_v_b][v] = n[x_b][y_b][8-v];
    n[x_b][y_b][8-v] = tmp;		
  }
}

void setrho(){
  for (int x=src_x; x<src_x+src_len; x++){
    int y=src_y;
    n[x][y][0]=src_den;
    n[x][y][1]=src_den;
    n[x][y][2]=src_den;
    n[x][y][3]=src_den;
    n[x][y][4]=src_den;
    n[x][y][5]=src_den;
    n[x][y][6]=src_den;
    n[x][y][7]=src_den;
    n[x][y][8]=src_den;
  }    
  for (int x=src_x_2; x<src_x_2+src_len_2; x++){
    int y=src_y_2;
    n[x][y][0]=src_den_2;
    n[x][y][1]=src_den_2;
    n[x][y][2]=src_den_2;
    n[x][y][3]=src_den_2;
    n[x][y][4]=src_den_2;
    n[x][y][5]=src_den_2;
    n[x][y][6]=src_den_2;
    n[x][y][7]=src_den_2;
    n[x][y][8]=src_den_2;
  }  
}

void init(){

  linkcount = 0;
  FindLink();


  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){

      dynamic_wall_position_x = 50;
     /* if(x > 25 && x <= 50 && y < 75 && y > 25){
	n[x][y][0]=d1;
	n[x][y][1]=d1;
	n[x][y][2]=d1;
	n[x][y][3]=d1;
	n[x][y][4]=d1;
	n[x][y][5]=d1;
	n[x][y][6]=d1;
	n[x][y][7]=d1;
	n[x][y][8]=d1;
	}
      else*/ if(x >= 0 && x <= 99 && y < 75 && y > 25){
	n[x][y][0]=d2;
	n[x][y][1]=d2;
	n[x][y][2]=d2;
	n[x][y][3]=d2;
	n[x][y][4]=d2;
	n[x][y][5]=d2;
	n[x][y][6]=d2;
	n[x][y][7]=d2;
	n[x][y][8]=d2;
	}
	else if(/*x <= 25 || x >= 75 ||*/ y >= 75 || y <= 25){
	n[x][y][0]=0;
	n[x][y][1]=0;
	n[x][y][2]=0;
	n[x][y][3]=0;
	n[x][y][4]=0;
	n[x][y][5]=0;
	n[x][y][6]=0;
	n[x][y][7]=0;
	n[x][y][8]=0;
	}
    }
  }
}
void initShear(){
  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){
      if (x<xdim/2){
	n[x][y][0]=1;
	n[x][y][1]=2;
	n[x][y][2]=3;
	n[x][y][3]=4;
	n[x][y][4]=5;
	n[x][y][5]=6;
	n[x][y][6]=7;
	n[x][y][7]=8;
	n[x][y][8]=9;
      }
      else {
	n[x][y][0]=8;
	n[x][y][1]=7;
	n[x][y][2]=6;
	n[x][y][3]=5;
	n[x][y][4]=4;
	n[x][y][5]=3;
	n[x][y][6]=2;
	n[x][y][7]=0;
	n[x][y][8]=0;
	}
    }
  }
}

void iterate(){
  int NI=0;
  int rate=RAND_MAX/8.;
  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){
      NI=0;
      for (int v=0; v<V; v++) NI+=n[x][y][v];
       // first implementation: random collisions
      if (NI>0){
	for (int c=0; c<C; c++){
	  //	  int r1=1.0*NI*rand()/RAND_MAX;
	  //	  int r2=1.0*NI*rand()/RAND_MAX;
	  int r1=rand()%NI;
	  int r2=rand()%NI;
	  if (r1!=r2){
	    int v1=0;
	    for (; r1>=n[x][y][v1];v1++)
	      r1-=n[x][y][v1];
	    int v2=0;
	    for (; r2>=n[x][y][v2];v2++)
	      r2-=n[x][y][v2];
	    // now we picked two particles with velocity v1 and v2.
	    int v1x=v1%3-1;
	    int v1y=v1/3-1;
	    int v2x=v2%3-1;
	    int v2y=v2/3-1;
	    // x-collision
	    if ((v1x==-1)&&(v2x==1)){
	      v1x=0;
	      v2x=0;
	    }
	    else if ((v1x==1)&&(v2x==-1)){
	      v1x=0;
	      v2x=0;
	    }
	    else if ((v1x==0)&&(v2x==0)){
	      int r=rand();
	      if (r<rate){
		if (r<rate/2){
		  v1x=-1;
		  v2x=1;
		}
		else {
		  v1x=1;
		  v2x=-1;
		}
	      }
	    }
	    else {
	      int tmp=v1x;
	      v1x=v2x;
	      v2x=tmp;
	    }
	    // y-collision
	    if ((v1y==-1)&&(v2y==1)){
	      v1y=0;
	      v2y=0;
	    }
	    else if ((v1y==1)&&(v2y==-1)){
	      v1y=0;
	      v2y=0;
	    }
	    else if ((v1y==0)&&(v2y==0)){
	      int r=rand();
	      if (r<rate){
		if (r<rate/2){
		  v1y=-1;
		  v2y=1;
		}
		else {
		  v1y=1;
		  v2y=-1;
		}
	      }
	    }
	    else {
	      int tmp=v1y;
	      v1y=v2y;
	      v2y=tmp;
	    }
	    n[x][y][v1]--;
	    n[x][y][v2]--;
	    v1=(v1y+1)*3+(v1x+1);
	    v2=(v2y+1)*3+(v2x+1);
	    n[x][y][v1]++;
	    n[x][y][v2]++;
	  }
	}
      }
    }
  }
  // move the particles in x-direction
  {
    int ntmp[ydim];
    for (int c=0; c<3; c++){
      for (int y=0; y<ydim; y++) ntmp[y]=n[xdim-1][y][c*3+2];
      for (int x=xdim-1; x>0; x--){
	measure_particle_velocity_front[x] = 0;
	for (int y=0; y<ydim; y++){
	  n[x][y][c*3+2]=n[x-1][y][c*3+2];


	}
      }
	
	//go here
      for (int y=0; y<ydim; y++) n[0][y][c*3+2]=ntmp[y];
    }
    for (int c=0; c<3; c++){
      for (int y=0; y<ydim; y++) ntmp[y]=n[0][y][c*3];
      for (int x=0; x<xdim-1; x++){
	for (int y=0; y<ydim; y++)
	  n[x][y][c*3]=n[x+1][y][c*3];
      }
      for (int y=0; y<ydim; y++) n[xdim-1][y][c*3]=ntmp[y];
    }
  }
    // move the particles in y-direction
  {
    int ntmp[xdim];
    for (int c=0; c<3; c++){
      for (int x=0; x<xdim; x++) ntmp[x]=n[x][ydim-1][c];
      for (int y=ydim-1; y>0; y--){
	for (int x=0; x<xdim; x++)
	  n[x][y][c]=n[x][y-1][c];
      }
      for (int x=0; x<xdim; x++) n[x][0][c]=ntmp[x];
    }
    for (int c=0; c<3; c++){
      for (int x=0; x<xdim; x++) ntmp[x]=n[x][0][c+6];
      for (int y=0; y<ydim-1; y++){
	for (int x=0; x<xdim; x++)
	  n[x][y][c+6]=n[x][y+1][c+6];
      }
      for (int x=0; x<xdim; x++) n[x][ydim-1][c+6]=ntmp[x];
    }
  }
  moveParticles();
  bounceback();
}

void getdata(){
  if (NUreq||Nreq){
    NUreq=0;
    Nreq=0;
    for (int x=0; x<xdim; x++)
      for (int y=0; y<ydim; y++){
	N[x][y]=0;
	NU[x][y][0]=0;
	NU[x][y][1]=0;
	for (int v=0; v<V; v++){
	  N[x][y]+=n[x][y][v];
	  NU[x][y][0]+=n[x][y][v]*(v%3-1);
	  NU[x][y][1]+=n[x][y][v]*(1-v/3);
	}
      }
  }
}

void debug(){
  int NUX=0,NUY=0,NN=0;
  for (int x=0; x<xdim; x++)
    for (int y=0; y<ydim; y++){
      for (int v=0; v<V; v++){
	NN+=n[x][y][v];
	NUX+=n[x][y][v]*(v%3-1);
	NUY+=n[x][y][v]*(1-v/3);
      }
    }
  printf("N=%i; NUX=%i; NUY=%i;\n",NN,NUX,NUY);
}

void main(){
  int done=0,sstep=0,cont=0,repeat=10;
  init();
  Measure();
  DefineGraphNxN_R("N",&N[0][0],&XDIM,&YDIM,&Nreq);
  DefineGraphNxN_RxR("NU",&NU[0][0][0],&XDIM,&YDIM,&NUreq);

  DefineGraphN_R("particle vx front",&measure_particle_velocity_front[0],&MeasLen,NULL);
  DefineGraphN_R("UX ZERO POINT",&ux_zero_point[0],&MeasLen,NULL);

  DefineGraphN_R("particle force x front",&measure_particle_force_front[0],&MeasLen,NULL);
  DefineGraphN_R("particle force x front theoretical",&theoretical_particle_vx_front[0],&MeasLen,NULL);

  StartMenu("LG",1);
  DefineFunction("init",init);
  DefineFunction("init shear",initShear);
  StartMenu("Measure",1);
  DefineInt("range_val",&range_val);
  DefineInt("ux zero point flag",&ux_zero_flag);
  DefineGraph(curve2d_,"Measurements");
  EndMenu();
  StartMenu("Wall",1);
  DefineDouble("particle flip w", &particle_flip_w);
  DefineFunction("Add Wall",FindLink);
  DefineInt("link count",&linkcount);
  EndMenu();
  DefineGraph(contour2d_,"Graph");
  DefineInt("C", &C);
  DefineInt("repeat",&repeat);
  DefineBool("sstep",&sstep);
  DefineBool("cont",&cont);
  DefineBool("done",&done);
  EndMenu();


  while (!done){
    Events(1);
    getdata();
    DrawGraphs();
    if (cont || sstep){
      sstep=0;

      for (int i=0; i<repeat; i++) {
	iterate();
	Measure();
	//setrho();
      }
    } else sleep(1);
  }
}

