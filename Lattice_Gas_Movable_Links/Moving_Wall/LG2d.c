// Lattice gas code in 1d.
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
int C=100;
double N[xdim][ydim],NU[xdim][ydim][2];
int Nreq=0,NUreq=0;
int n[xdim][ydim][V];

// Boundary variables
#define LINKMAX 10000
#define LINKMAX_DYNAMIC 10000



double dt = 1.0;



//variables for measuring tube momentum



//particle source parameters
int src_x = 85,src_y = 85,src_len = 5,src_den = 10;
int src_x_2 = 15,src_y_2 = 15,src_len_2 = 5,src_den_2 = 10;
int var1 =10;
//variables for measuring values
//particle velocity
int particle_vx = 0,particle_vy;
double measure_particle_vx[MeasMax],measure_particle_vy[MeasMax];
double measure_particle_vx_filt[MeasMax],measure_particle_vy_filt[MeasMax];
int vx_m = 0,vy_m = 0;



//link variables
int linkcount_dynamic = 0,links_dynamic[LINKMAX_DYNAMIC][3];//link x y and v(for particles)
double links_dynamic_velocity[LINKMAX_DYNAMIC][2];//link velocities
double wall_mass = 10000000000.1;
int yy3 = 26,yy4 = 75,flow = 0,dynamic_walls_on = 1;
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

//static wall variables
int linkcount=0,links[LINKMAX][3];

//points for drawing walls
int x0=25,x1=75,yy0 = 25,yy1 = 75;

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
int filt_data[] ={0,0,0,0,0,0,0,0,0,0}; //dynamic walls -> x,y,px,py,vx,vy  static walls -> px py particles -> vx vy


//added a running average to smooth out the graphs
void average(int range){
	val[0] = 0;
	val[1] = 0;
	//clear previous filtered data
	for(int i = 0; i < 10;i++){
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
		}
	//average values
	for(int i = 0; i < 10;i++){
		filt_data[i] = filt_data[i]/range;
		}
}

void measure_function(){
	particle_vx = 0;
	particle_vy = 0;
	for(int i = 0; i < YDIM -1;i++){
		particle_vx += n[50][i][2]+n[50][i][5]+n[50][i][8]-n[50][i][0]-n[50][i][3]-n[50][i][6];
		//particle_vy += n[i][50][0]+n[i][50][1]+n[i][50][2]-n[i][50][6]-n[i][50][7]-n[i][50][8];
		}


}



void Measure(){
  //store dynamic wall data
  //position
  memmove(&measure_dynamic_wall_position_x[1],&measure_dynamic_wall_position_x[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_position_x[0]=dynamic_wall_position_x;
  memmove(&measure_dynamic_wall_position_y[1],&measure_dynamic_wall_position_y[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_position_y[0]=dynamic_wall_position_y;
  //momentum
  memmove(&measure_dynamic_wall_momentum_x[1],&measure_dynamic_wall_momentum_x[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_momentum_x[0]=dynamic_wall_momentum_x;
  memmove(&measure_dynamic_wall_momentum_y[1],&measure_dynamic_wall_momentum_y[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_momentum_y[0]=dynamic_wall_momentum_y;
  //velocity
  memmove(&measure_dynamic_wall_vx[1],&measure_dynamic_wall_vx[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_vx[0]=dynamic_wall_vx;
  memmove(&measure_dynamic_wall_vy[1],&measure_dynamic_wall_vy[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_vy[0]=dynamic_wall_vy;
  
  //store static wall data
  memmove(&measure_static_wall_momentum_x[1],&measure_static_wall_momentum_x[0],(MeasMax-1)*sizeof(int));
  measure_static_wall_momentum_x[0]=static_wall_momentum_x;
  memmove(&measure_static_wall_momentum_y[1],&measure_static_wall_momentum_y[0],(MeasMax-1)*sizeof(int));
  measure_static_wall_momentum_y[0]=static_wall_momentum_y;

  //store particle data
  memmove(&measure_particle_vx[1],&measure_particle_vx[0],(MeasMax-1)*sizeof(int));
  measure_particle_vx[0]=particle_vx;
  memmove(&measure_particle_vy[1],&measure_particle_vy[0],(MeasMax-1)*sizeof(int));
  measure_particle_vy[0]=particle_vy;

  //filter data
  average(range_val);

   //position
  memmove(&measure_dynamic_wall_position_x_filt[1],&measure_dynamic_wall_position_x_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_position_x_filt[0]=filt_data[0];
  memmove(&measure_dynamic_wall_position_y_filt[1],&measure_dynamic_wall_position_y_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_position_y_filt[0]=filt_data[1];
  //momentum
  memmove(&measure_dynamic_wall_momentum_x_filt[1],&measure_dynamic_wall_momentum_x_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_momentum_x_filt[0]=filt_data[2];
  memmove(&measure_dynamic_wall_momentum_y_filt[1],&measure_dynamic_wall_momentum_y_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_momentum_y_filt[0]=filt_data[3];
  //velocity
  memmove(&measure_dynamic_wall_vx_filt[1],&measure_dynamic_wall_vx_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_vx_filt[0]=filt_data[4];
  memmove(&measure_dynamic_wall_vy_filt[1],&measure_dynamic_wall_vy_filt[0],(MeasMax-1)*sizeof(int));
  measure_dynamic_wall_vy_filt[0]=filt_data[5];
  
  //store static wall data
  memmove(&measure_static_wall_momentum_x_filt[1],&measure_static_wall_momentum_x_filt[0],(MeasMax-1)*sizeof(int));
  measure_static_wall_momentum_x_filt[0]=filt_data[6];
  memmove(&measure_static_wall_momentum_y_filt[1],&measure_static_wall_momentum_y_filt[0],(MeasMax-1)*sizeof(int));
  measure_static_wall_momentum_y_filt[0]=filt_data[7];

  //store particle data
  memmove(&measure_particle_vx_filt[1],&measure_particle_vx_filt[0],(MeasMax-1)*sizeof(int));
  measure_particle_vx_filt[0]=filt_data[8];
  memmove(&measure_particle_vy_filt[1],&measure_particle_vy_filt[0],(MeasMax-1)*sizeof(int));
  measure_particle_vy_filt[0]=filt_data[9];;
}
//re draws walls
void FindLink_Dynamic(){
	int tmp_x_shift = 0;
	

		linkcount_dynamic = 0;
		//draw vertical walls	
		for(int y = yy3; y<yy4+1;y++){

			tmp_x_shift = dynamic_wall_position_x;
			tmp_x_shift = ((tmp_x_shift%XDIM)+XDIM)%XDIM;	
			links_dynamic[linkcount_dynamic][0] = tmp_x_shift; //x-position
			links_dynamic[linkcount_dynamic][1] = y;
			links_dynamic[linkcount_dynamic][2] = 0;
			links_dynamic_velocity[linkcount_dynamic][0] = dynamic_wall_vx;//set x velocity
			links_dynamic_velocity[linkcount_dynamic][1] = 0;//set y velocity
			linkcount_dynamic++;
			links_dynamic[linkcount_dynamic][0] = tmp_x_shift; //x-position
			links_dynamic[linkcount_dynamic][1] = y;
			links_dynamic[linkcount_dynamic][2] = 3;
			links_dynamic_velocity[linkcount_dynamic][0] = dynamic_wall_vx;//set x velocity
			links_dynamic_velocity[linkcount_dynamic][1] = 0;//set y velocity
			linkcount_dynamic++;
			links_dynamic[linkcount_dynamic][0] = tmp_x_shift; //x-position
			links_dynamic[linkcount_dynamic][1] = y;
			links_dynamic[linkcount_dynamic][2] = 6;
			links_dynamic_velocity[linkcount_dynamic][0] = dynamic_wall_vx;//set x velocity
			links_dynamic_velocity[linkcount_dynamic][1] = 0;//set y velocity
			linkcount_dynamic++;	
			}
		
		dynamic_wall_position_x += dynamic_wall_vx*dt;
		//dynamic_wall_vx = momentum_est_x[0]/wall_mass;
}

void FindLink(){
	//2 links added to prevent leaking on the x0,yy2 and x2,yy0 squares 
	links[linkcount][0] = x0; //x-position
    	links[linkcount][1] = yy1;
    	links[linkcount][2] = 0;
    	linkcount++;
	links[linkcount][0] = x1; //x-position
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 0;
    	linkcount++;
  //horizontal walls
  for (int x=x0; x<x1+1; x++){
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
    for (int x=x0; x<x1+1; x++){
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
  
  //vertical walls
  for (int y=yy0; y<yy1+1; y++){
	links[linkcount][0] = x0; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 0;
	linkcount++;
	links[linkcount][0] = x0; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 3;
	linkcount++;
	links[linkcount][0] = x0; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 6;
	linkcount++;
  }

  for (int y=yy0; y<yy1+1; y++){
	links[linkcount][0] = x1; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 0;
	linkcount++;
	links[linkcount][0] = x1; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 3;
	linkcount++;
	links[linkcount][0] = x1; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 6;
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
    //velocity
    int v=links[lc][2];
    int vx=v%3-1;
    int vy=1-v/3;
    int tmp= n[x+vx][y+vy][v];
    //summing all momemtums
    static_wall_momentum_x += -2*vx*(n[x][y][8-v]-tmp);
    static_wall_momentum_y += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[x+vx][y+vy][v]= n[x][y][8-v];
    n[x][y][8-v]=tmp;		
  }
//dynamic walls
if(dynamic_walls_on == 1){
 for (int lc=0; lc<linkcount_dynamic; lc++){
    //quantity of partices
    int x=links_dynamic[lc][0];
    int y=links_dynamic[lc][1];
    //particle velocity
    int v=links_dynamic[lc][2];
    int vx=v%3-1;
    int vy=1-v/3;
    int tmp_0 = n[(((vx+x)%XDIM)+XDIM)%XDIM][((y+vy)%YDIM+YDIM)%YDIM][8-v];
    int tmp = n[(((vx+x)%XDIM)+XDIM)%XDIM][((y+vy)%YDIM+YDIM)%YDIM][v];
    //printf("switch : %i with %i \n",(((vx+x)%XDIM)+XDIM)%XDIM,x);
    //printf("tmp: %i \n",tmp);
	//find the smaller value to be the max for the random function to avoid having negative densisties
    //determining values for forward and backward particle flow
	int max_random = 1;
	if(tmp > tmp_0){
		max_random = tmp_0;
		}
	else{
		max_random = tmp;
		}
	if(max_random > 0){
		flow = (rand()%max_random)*links_dynamic_velocity[lc][0]*links_dynamic_velocity[lc][0];
		}
	else{
		flow = 0;
		}
//printf("max random = %i %i %i\n",max_random,tmp,n[x][y][8-v]);
    //summing all momemtums
    dynamic_wall_momentum_x += -2*vx*(n[x][y][8-v]-tmp+flow);
    dynamic_wall_momentum_y += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[(((vx+x)%XDIM)+XDIM)%XDIM][((y+vy)%YDIM+YDIM)%YDIM][v] = n[(((x)%XDIM)+XDIM)%XDIM][((y)%YDIM+YDIM)%YDIM][8-v] - flow;
    n[(((x)%XDIM)+XDIM)%XDIM][((y)%YDIM+YDIM)%YDIM][8-v] = tmp + flow;		
  }
    dynamic_wall_vx = dynamic_wall_momentum_x/wall_mass;
    dynamic_wall_vy = dynamic_wall_momentum_y/wall_mass;
}
  //measure routine stores values for plotting
  Measure();
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
  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){

      dynamic_wall_position_x = 50;
      if(x > 25 && x < 50 && y < 75 && y > 25){
	n[x][y][0]=50;
	n[x][y][1]=50;
	n[x][y][2]=50;
	n[x][y][3]=50;
	n[x][y][4]=50;
	n[x][y][5]=50;
	n[x][y][6]=50;
	n[x][y][7]=50;
	n[x][y][8]=50;
	}
      else if(x > 50 && x < 75 && y < 75 && y > 25){
	n[x][y][0]=60;
	n[x][y][1]=60;
	n[x][y][2]=60;
	n[x][y][3]=60;
	n[x][y][4]=60;
	n[x][y][5]=60;
	n[x][y][6]=60;
	n[x][y][7]=60;
	n[x][y][8]=60;
	}
	else if(x < 25 || x > 75 && y > 75 || y < 25){
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
      else {
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
	for (int y=0; y<ydim; y++)
	  n[x][y][c*3+2]=n[x-1][y][c*3+2];
      }
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
  if(dynamic_walls_on == 1){
	FindLink_Dynamic();	
	}
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
  DefineGraphN_R("dynamic wall position x",&measure_dynamic_wall_position_x[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall position y",&measure_dynamic_wall_position_y[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall momentum x",&measure_dynamic_wall_momentum_x[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall momentum y",&measure_dynamic_wall_momentum_y[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall vx ",&measure_dynamic_wall_vx[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall vx",&measure_dynamic_wall_vy[0],&MeasLen,NULL);
  DefineGraphN_R("static wall momentum x",&measure_static_wall_momentum_x[0],&MeasLen,NULL);
  DefineGraphN_R("static wall momentum y",&measure_static_wall_momentum_y[0],&MeasLen,NULL);
  DefineGraphN_R("particle vx",&measure_particle_vx[0],&MeasLen,NULL);
  DefineGraphN_R("particle vy",&measure_particle_vy[0],&MeasLen,NULL);
  //filtered graphs
  DefineGraphN_R("dynamic wall position x filt",&measure_dynamic_wall_position_x_filt[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall position y filt",&measure_dynamic_wall_position_y_filt[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall momentum x filt",&measure_dynamic_wall_momentum_x_filt[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall momentum y filt",&measure_dynamic_wall_momentum_y_filt[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall vx filt",&measure_dynamic_wall_vx_filt[0],&MeasLen,NULL);
  DefineGraphN_R("dynamic wall vx filt",&measure_dynamic_wall_vy_filt[0],&MeasLen,NULL);
  DefineGraphN_R("static wall momentum x filt",&measure_static_wall_momentum_x_filt[0],&MeasLen,NULL);
  DefineGraphN_R("static wall momentum y filt",&measure_static_wall_momentum_y_filt[0],&MeasLen,NULL);
  DefineGraphN_R("particle vx filt",&measure_particle_vx_filt[0],&MeasLen,NULL);
  DefineGraphN_R("particle vy filt",&measure_particle_vy_filt[0],&MeasLen,NULL);
  StartMenu("LG",1);
  DefineFunction("init",init);
  DefineFunction("init shear",initShear);
  StartMenu("Measure",0);
  DefineInt("range_val",&range_val);
  DefineGraph(curve2d_,"Measurements");
  EndMenu();
  StartMenu("Wall",0);
  DefineInt("x0", &x0);
  DefineInt("x1", &x1);
  DefineDouble("dynamic_wall_vx", &dynamic_wall_vx);
  DefineDouble("dynamic_wall_position_x", &dynamic_wall_position_x);
  DefineDouble("wall_mass", &wall_mass);
  DefineInt("yy0", &yy0);
  DefineInt("yy1", &yy1);
  DefineInt("yy3", &yy3);
  DefineInt("yy4", &yy4);
  DefineFunction("Add Wall",FindLink);
  DefineInt("link count",&linkcount);
  DefineInt("dynamic link count",&linkcount_dynamic);
  DefineInt("dynamic walls on",&dynamic_walls_on);
  DefineInt("flow",&flow);
  EndMenu();
  StartMenu("Particle Source",0);
  DefineInt("src_den",&src_den);
  DefineInt("src_x",&src_x);
  DefineInt("src_y",&src_y);
  DefineInt("src_den_2",&src_den_2);
  DefineInt("src_x_2",&src_x_2);
  DefineInt("src_y_2",&src_y_2);
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
	//setrho();
      }
    } else sleep(1);
  }
}

