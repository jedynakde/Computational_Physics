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
#define WALLMAX 10
int x0=10,x1=20,x2 = 40,y2=20,yy0 = 10,yy1 = 20, yy2 = 40,close_tube = 1;
int linkcount=0,links[LINKMAX][3];
//varaibles for dynamic walls
int linkcount_dynamic = 0,links_dynamic[LINKMAX_DYNAMIC][3];//link x y and v(for particles)
int wallcount=1;
int shape[WALLMAX][4];// shape has walls each with a 2d start and end coordinate
double shape_velocity[WALLMAX][LINKMAX_DYNAMIC][2];
double links_dynamic_velocity[LINKMAX_DYNAMIC][2];//link velocities
double dt = 1.0,vx0 = 0.0, vy0 =0.0;
int x3 = 100,yy3 = 0,yy4 = 100;
int move_period = 1,vx_dir = 1,move_count = 1,dynamic_walls_on = 1,y_shift = 0,flow = 0;
//variables for measuring tube momentum
int tot_vx=0,tot_vy=0;
int tot_vx_list[10],tot_vy_list[10];
//particle source parameters
int src_x = 12,src_y = 30,src_len = 5,src_den = 0;
int src_x_2 = 50,src_y_2 = 50,src_len_2 = 5,src_den_2 = 0;
int var1 =10;
double wall_mass = 10.1;
//variables for measuring values
int vx_m = 0,vy_m = 0;



//graphing
double x_shift = 30;
double vx_est[MeasMax],vy_est[MeasMax];
double momentum_est_x[MeasMax],momentum_est_y[MeasMax];
double momentum_est_x_filt[MeasMax],momentum_est_y_filt[MeasMax];
int MeasLen = MeasMax/2;
int range_val = 20;
int val[] = {0,0};
//added a running average to smooth out the graphs
void average(int range){
	val[0] = 0;
	val[1] = 0;
	for(int i = 0; i < range;i++){
		val[0] += momentum_est_x[i];
		val[1] += momentum_est_y[i];
		}
	val[0] = val[0]/range;
	val[1] = val[1]/range;
}

void measure_function(){
	vx_m = 0;
	vy_m = 0;
	for(int i = 0; i < YDIM -1;i++){
		vx_m += n[50][i][2]+n[50][i][5]+n[50][i][8]-n[50][i][0]-n[50][i][3]-n[50][i][6];
		vy_m += n[i][50][0]+n[i][50][1]+n[i][50][2]-n[i][50][6]-n[i][50][7]-n[i][50][8];
		}


}



void Measure(){
  memmove(&momentum_est_x[1],&momentum_est_x[0],(MeasMax-1)*sizeof(int));
  momentum_est_x[0]=tot_vx;
  memmove(&momentum_est_y[1],&momentum_est_y[0],(MeasMax-1)*sizeof(int));
  momentum_est_y[0]=tot_vy;
  //filtered data
  average(range_val);
  memmove(&momentum_est_x_filt[1],&momentum_est_x_filt[0],(MeasMax-1)*sizeof(int));
  momentum_est_x_filt[0]= val[0];
  memmove(&momentum_est_y_filt[1],&momentum_est_y_filt[0],(MeasMax-1)*sizeof(int));
  momentum_est_y_filt[0]=val[1];
  measure_function();
  memmove(&vx_est[1],&vx_est[0],(MeasMax-1)*sizeof(int));
  vx_est[0]=vx_m;
  memmove(&vy_est[1],&vy_est[0],(MeasMax-1)*sizeof(int));
  vy_est[0]=vy_m;
}
//re draws walls
void FindLink_Dynamic(){
	int tmp_x_shift = 0;
	int tmp_y_shift = 0;
	linkcount_dynamic = 0;//clear old shape walls

	for(int wall = 0; wall < wallcount; wall++){//iterate through all walls in the shape
		printf("wall: %i \n",wall);
		double ab_dx = abs(shape[wall][1]-shape[wall][0]);
		double ab_dy = abs(shape[wall][3]-shape[wall][2]);
		double start_y = 0;
		double end_y = 0;
		double start_x = 0;
		double end_x = 0;

		//determine which point to start from to allow for only positive increments
		if(shape[wall][3] > shape[wall][2]){
			start_y = shape[wall][2];
			end_y = shape[wall][3];
			}
		else{
			start_y = shape[wall][3];
			end_y = shape[wall][2];
			}
		if(shape[wall][1] > shape[wall][0]){
			start_x = shape[wall][0];
			end_x = shape[wall][1];
			}
		else{
			start_x = shape[wall][1];
			end_x = shape[wall][0];
			}
		//calculate length of the line
		double line_length = sqrt((shape[wall][1]-shape[wall][0])*(shape[wall][1]-shape[wall][0])+(shape[wall][3]-shape[wall][2])*(shape[wall][3]+shape[wall][2]));
		//iterate through all links that make up that line
		for(double l = 0; l < line_length;l = l + 1){
			for(int y_pos = start_y;y_pos < (start_y+ab_dy);y_pos++){
				tmp_x_shift = start_x;
				tmp_y_shift = start_y;
				tmp_x_shift = ((tmp_x_shift%XDIM)+XDIM)%XDIM;
				tmp_y_shift = ((tmp_y_shift%YDIM)+YDIM)%YDIM;
				//make vertical lines
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 0;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 3;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 6;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				printf("data: %i %i \n",tmp_x_shift,tmp_y_shift);
			}
			//make horizontal lines
			for(int x_pos = start_x;x_pos < (start_x+ab_dx);x_pos++){
				tmp_x_shift = start_x;
				tmp_y_shift = start_y;
				tmp_x_shift = ((tmp_x_shift%XDIM)+XDIM)%XDIM;
				tmp_y_shift = ((tmp_y_shift%YDIM)+YDIM)%YDIM;
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 0;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 1;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				shape[linkcount_dynamic][0] = tmp_x_shift; //x-position
				shape[linkcount_dynamic][1] = tmp_y_shift;
				shape[linkcount_dynamic][2] = 2;
				shape_velocity[wall][linkcount_dynamic][0] = vx0;//set x velocity
				shape_velocity[wall][linkcount_dynamic][1] = vy0;//set y velocity
				linkcount_dynamic++;
				printf("data: %i %i \n",tmp_x_shift,tmp_y_shift);
			}
			//increment the x and y positions for the links
			if(ab_dx == 0){//vertical wall
				start_y = l;
				}
			else if(ab_dy == 0){//horizontal wall
				start_x = l;
				}
			else{
				start_x += ab_dx;
				start_y += ab_dy;
				} 
			}
		

	}
	//update link positions based on velocity anf velocity based on collision momentum
	for(int wall = 0; wall < wallcount; wall++){
		for(int lc = 0; lc < linkcount_dynamic;lc++){
			shape[linkcount_dynamic][0] += shape_velocity[wall][linkcount_dynamic][0]*dt; //x-position
			shape[linkcount_dynamic][1] += shape_velocity[wall][linkcount_dynamic][1]*dt;
			shape_velocity[wall][linkcount_dynamic][0] = momentum_est_x[0]/wall_mass;
			shape_velocity[wall][linkcount_dynamic][1] = momentum_est_y[0]/wall_mass;
			}
		}

}

void FindLink(){
	//2 links added to prevent leaking on the x0,yy2 and x2,yy0 squares 
	links[linkcount][0] = x0; //x-position
    	links[linkcount][1] = yy2;
    	links[linkcount][2] = 0;
    	linkcount++;
	links[linkcount][0] = x2; //x-position
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 0;
    	linkcount++;
  //horizontal walls
  for (int x=x0; x<x1+1; x++){
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy2;
    	links[linkcount][2] = 0;
    	linkcount++;
    	links[linkcount][0] = x; //x-position
    	links[linkcount][1] = yy2;
    	links[linkcount][2] = 1;
    	linkcount++;
    	links[linkcount][0] = x; //x-position 
    	links[linkcount][1] = yy2;
    	links[linkcount][2] = 2;
    	linkcount++;
  }
    for (int x=x1; x<x2+1; x++){
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
    for (int x=x0; x<x2+1; x++){
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
  //vertical walls
if(close_tube == 1){
  for (int y=yy0; y<yy1+1; y++){
	links[linkcount][0] = x2; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 0;
	linkcount++;
	links[linkcount][0] = x2; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 3;
	linkcount++;
	links[linkcount][0] = x2; //x-position
	links[linkcount][1] = y;
	links[linkcount][2] = 6;
	linkcount++;
  }
}
  for (int y=yy1; y<yy2+1; y++){
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
  for (int y=yy0; y<yy2+1; y++){
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

}

void bounceback(){
  tot_vx =0;
  tot_vy =0;
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
    tot_vx += -2*vx*(n[x][y][8-v]-tmp);
    tot_vy += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[x+vx][y+vy][v]= n[x][y][8-v];
    n[x][y][8-v]=tmp;		
  }
//dynamic walls
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
		flow = (rand()%max_random)*links_dynamic_velocity[lc][0];
		}
	else{
		flow = 0;
		}
//printf("max random = %i %i %i\n",max_random,tmp,n[x][y][8-v]);
    //summing all momemtums
    tot_vx += -2*vx*(n[x][y][8-v]-tmp+flow);
    tot_vy += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[(((vx+x)%XDIM)+XDIM)%XDIM][((y+vy)%YDIM+YDIM)%YDIM][v] = n[(((x)%XDIM)+XDIM)%XDIM][((y)%YDIM+YDIM)%YDIM][8-v] - flow;
    n[(((x)%XDIM)+XDIM)%XDIM][((y)%YDIM+YDIM)%YDIM][8-v] = tmp + flow;		
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


  shape[0][0] = x0;
  shape[0][1] = x1;
  shape[0][2] = yy0;
  shape[0][3] = yy1;

  

  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){
      //if ((abs(xdim/2-x)<25)&&(abs(ydim/2-y)<25)){
	n[x][y][0]=10;
	n[x][y][1]=10;
	n[x][y][2]=10;
	n[x][y][3]=10;
	n[x][y][4]=10;
	n[x][y][5]=10;
	n[x][y][6]=10;
	n[x][y][7]=10;
	n[x][y][8]=10;
      

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
  DefineGraphN_R("momentum_est_x",&momentum_est_x[0],&MeasLen,NULL);
  DefineGraphN_R("momentum_est_y",&momentum_est_y[0],&MeasLen,NULL);
  DefineGraphN_R("momentum_est_x_filt",&momentum_est_x_filt[0],&MeasLen,NULL);
  DefineGraphN_R("momentum_est_y_filt",&momentum_est_y_filt[0],&MeasLen,NULL);
  DefineGraphN_R("vx_est",&vx_est[0],&MeasLen,NULL);
  DefineGraphN_R("vy_est",&vy_est[0],&MeasLen,NULL);
  StartMenu("LG",1);
  DefineFunction("init",init);
  DefineFunction("init shear",initShear);
  StartMenu("Measure",0);
  DefineInt("range_val",&range_val);
  DefineGraph(curve2d_,"Measurements");
  DefineInt("tot_vx",&tot_vx);
  DefineInt("tot_vy",&tot_vy);
  EndMenu();
  StartMenu("Wall",0);
  DefineInt("x0", &x0);
  DefineInt("x1", &x1);
  DefineInt("x2", &x2);
  DefineInt("x3", &x3);
  DefineInt("vx_dir", &vx_dir);
  DefineDouble("vx0", &vx0);
  DefineDouble("wall_mass", &wall_mass);
  DefineInt("yy0", &yy0);
  DefineInt("yy1", &yy1);
  DefineInt("yy2", &yy2);
  DefineInt("yy3", &yy3);
  DefineInt("yy4", &yy4);
  DefineInt("close_tube",&close_tube);
  DefineFunction("Add Wall",FindLink);
  DefineInt("link count",&linkcount);
  DefineInt("dynamic link count",&linkcount_dynamic);
  DefineInt("dynamic walls on",&dynamic_walls_on);
  DefineInt("move_period",&move_period);
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
