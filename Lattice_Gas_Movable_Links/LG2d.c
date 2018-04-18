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
int wall_shift = -15;

// Boundary variables
#define LINKMAX 10000
#define DYNAMIC_LINKMAX 10000
int x0=10,x1=20,x2 = 40,y2=20,yy0 = 10,yy1 = 20, yy2 = 40,close_tube = 1;
int linkcount=0,links[LINKMAX][3];
int linkcount_dynamic = 0,links_dynamic[DYNAMIC_LINKMAX][3],dynamic_walls_on = 1,auto_move = 0,wall_dx = 0;
//variables for measuring tube momentum
int tot_vx=0,tot_vy=0;
int tot_vx_list[10],tot_vy_list[10];
//particle source parameters
int src_x = 12,src_y = 30,src_len = 5,src_den = 15;
int src_x_2 = 50,src_y_2 = 50,src_len_2 = 5,src_den_2 = 0;
int var1 =10;
int link_start = 0,link_end = 0,link_flag = 0;

//graphing
double momentum_est_x[MeasMax],momentum_est_y[MeasMax];
double momentum_est_x_filt[MeasMax],momentum_est_y_filt[MeasMax];
int MeasLen = MeasMax/2;
int range_val = 20;
int val[] = {0,0};
int wall_flag = 0;
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
}
void FindLink_Dynamic(){

if(dynamic_walls_on == 1 && wall_dx != 0){
  if(link_flag == 1){
        linkcount_dynamic = link_start;
	}
else{
link_start = linkcount_dynamic;
}
        link_flag = 1;
  for (int y=yy0+2; y<yy1; y++){
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 0;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 3;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 6;
	linkcount_dynamic++;
        //additional 3 links if the wall is moving
	if(wall_dx != 0){
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 1;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 4;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 7;
	linkcount_dynamic++;
	//added more links
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 2;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 5;
	linkcount_dynamic++;
	links_dynamic[linkcount_dynamic][0] = x2+wall_shift; //x-position
	links_dynamic[linkcount_dynamic][1] = y;
	links_dynamic[linkcount_dynamic][2] = 8;
	linkcount_dynamic++;
	}
  }
 link_end = linkcount_dynamic;
}

if(wall_shift < -20){
//wall_flag = 1;
wall_dx = -wall_dx;
}
else if(wall_shift > -10){
//wall_flag = 0;
wall_dx = -wall_dx;
}

wall_shift += wall_dx;
/*
if(wall_flag == 0){
wall_shift -= 1*wall_dx;
}
else if(wall_flag == 1)
{
wall_shift += 1*wall_dx;
}
*/

printf("wall_shift = %i \n",wall_shift);

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
//for static walls
void bounceback(){
  tot_vx =0;
  tot_vy =0;
  for (int lc=0; lc<linkcount; lc++){
    //quantity of partices
    int x=links[lc][0];
    int y=links[lc][1];
    //velocity
    int v=links[lc][2];
    int vx=v%3-1;//might need to change?
    int vy=1-v/3;//might need to change?
    int tmp= n[x+vx][y+vy][v];
    //summing all momemtums
    tot_vx += -2*vx*(n[x][y][8-v]-tmp);
    tot_vy += -2*vy*(n[x][y][8-v]-tmp);
    //swapping the particles trying to enter and leave to have the effect of a wall
    n[x+vx][y+vy][v]= n[x][y][8-v];
    n[x][y][8-v]=tmp;		
  }
  //measure routine stores values for plotting
  Measure();
}

//bounceback for dynamic walls
void bounceback_dynamic(){
  tot_vx =0;
  tot_vy =0;
  for (int lc=0; lc<linkcount_dynamic; lc++){
    //quantity of partices
    int x=links_dynamic[lc][0];
    int y=links_dynamic[lc][1];
    //velocity
    int v=links_dynamic[lc][2];
    int vx=v%3-1;//might need to change?
    int vy=1-v/3;//might need to change?
    int tmp= n[x+vx][y+vy][v];		
    //summing all momemtums
//if a wall is not moving treat it like a static wall
    if(wall_dx == 0 && v == 0 || v == 3 || v == 6){

    //swapping the particles trying to enter and leave to have the effect of a wall
    n[x+vx][y+vy][v]= n[x][y][8-v];
    n[x][y][8-v]=tmp;
    tot_vx += -2*vx*(n[x][y][8-v]-tmp);
    tot_vy += -2*vy*(n[x][y][8-v]-tmp);
}
//else if the wall is moving right
else if(wall_dx < 0){
    if(v == 0 || v == 3 || v ==6){
	n[x+vx][y+vy][v] += n[x][y][v+1];	
	}
    else if(v == 1 || v == 4 || v == 7){
	n[x+vx][y+vy][v] = n[x][y][v+1];	
	}
    else{
	n[x+vx][y+vy][v] = 0;
	}
}
//else if the wall is moving left
else if(wall_dx > 0){
    if(v == 2 || v == 5 || v == 8){
	n[x+vx][y+vy][v] += n[x][y][v-1];	
	}
    else if(v == 1 || v == 4 || v == 7){
	n[x+vx][y+vy][v] = n[x][y][v-1];	
	}
    else{
	n[x+vx][y+vy][v] = 0;
	}
}
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
      if ((abs(xdim/2-x)<25)&&(abs(ydim/2-y)<25)){
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
      else
	{
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
	n[x][y][0]=20;
	n[x][y][1]=30;
	n[x][y][2]=20;
	n[x][y][3]=20;
	n[x][y][4]=10;
	n[x][y][5]=20;
	n[x][y][6]=0;
	n[x][y][7]=10;
	n[x][y][8]=0;
      }
      else {
	n[x][y][0]=0;
	n[x][y][1]=10;
	n[x][y][2]=00;
	n[x][y][3]=20;
	n[x][y][4]=10;
	n[x][y][5]=20;
	n[x][y][6]=20;
	n[x][y][7]=30;
	n[x][y][8]=20;
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
  bounceback();
if(dynamic_walls_on == 1)
{
if(wall_dx != 0)FindLink_Dynamic();
bounceback_dynamic();
}
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
  DefineInt("yy0", &yy0);
  DefineInt("yy1", &yy1);
  DefineInt("yy2", &yy2);
  DefineInt("wall_shift",&wall_shift);
  DefineInt("close_tube",&close_tube);
  DefineFunction("Add Wall",FindLink);
  DefineInt("link count",&linkcount);
  DefineInt("dynamic walls",&dynamic_walls_on);
  DefineInt("wall dx",&wall_dx);
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
	setrho();
      }
    } else sleep(1);
  }
}

