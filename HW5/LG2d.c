// Lattice gas code in 1d.
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <mygraph.h>

#define xdim 100
#define ydim 100
#define V 9      //number of velocities
int XDIM=xdim,YDIM=ydim;
int C=100;
double N[xdim][ydim],NU[xdim][ydim][2];
int Nreq=0,NUreq=0;
int n[xdim][ydim][V];

// Boundary variables
#define LINKMAX 10000
int x0=10,x1=20,x2 = 40,y2=20,yy0 = 10,yy1 = 20, yy2 = 40;
int linkcount=0,links[LINKMAX][3];
//variables for measuring tube momentum
int tot_vx=0,tot_vy=0;
int x0y1 = 0;

void FindLink(){

	links[linkcount][0] = x0; //x-position
    	links[linkcount][1] = yy2;
    	links[linkcount][2] = 0;
    	linkcount++;
	links[linkcount][0] = x2; //x-position
    	links[linkcount][1] = yy0;
    	links[linkcount][2] = 0;
    	linkcount++;
  for (int x=x0; x<x1+1; x++){
    	x0y1 = 0;
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
    	x0y1 = 0;
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
	x0y1 = 0;
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

  for (int y=yy0; y<yy1+1; y++){
	x0y1 = 1;
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
  for (int y=yy1; y<yy2+1; y++){
	x0y1 = 1;
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
	x0y1 = 1;
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
    int x=links[lc][0];
    int y=links[lc][1];
    int v=links[lc][2];
    int vx=v%3-1;
    int vy=1-v/3;
	//to find the total momentum of the system
	tot_vy += n[x][y][1] - n[x][y][7];
	tot_vx += n[x][y][3] - n[x][y][5];;
	int tmp= n[x+vx][y+vy][v];
	//if(x0y1 ==1 && x == y || y == yy2 && x == x1 || && v == 0 || v == 8){//if the link is a corner or intersect of 2 lines
	if((
	   (x == x0 && y ==yy0)||
	   (x == x1 && y ==yy1)||
	   (y == yy2 && x == x0 )||
	   ( y == yy2 && x == x1)||
	   (x == x2 && y == yy0 )||
	   (x == x2 && y == yy1)) &&
	    (v == 0) &&
	    x0y1 == 5
	    ){
			//n[x+vx][y+vy][v]= 0 ;// n[x][y][8-v];
			//n[x][y][8-v]= 0 ;//tmp;
			printf("flag %i at (%i,%i) \n",v,x,y);
		}
	else{//else if not a corner
		n[x+vx][y+vy][v]= n[x][y][8-v];
		n[x][y][8-v]=tmp;
		printf("swap %i with %i link %i at (%i,%i)\n",8 - v, v, lc,x,y);
		}
	
  }
printf("done\n");
}

void setrho(){
  for (int x=50; x<70; x++){
    int y=20;
    n[x][y][0]=10;
    n[x][y][1]=20;
    n[x][y][2]=10;
    n[x][y][3]=20;
    n[x][y][4]=10;
    n[x][y][5]=20;
    n[x][y][6]=10;
    n[x][y][7]=20;
    n[x][y][8]=10;
  }    
}

void init(){
  for (int x=0; x<xdim; x++){
    for (int y=0; y<ydim; y++){
      if ((abs(xdim/2-x)<25)&&(abs(ydim/2-y)<25)){
	n[x][y][0]=10;
	n[x][y][1]=20;
	n[x][y][2]=10;
	n[x][y][3]=20;
	n[x][y][4]=10;
	n[x][y][5]=20;
	n[x][y][6]=10;
	n[x][y][7]=20;
	n[x][y][8]=10;
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
  
  DefineGraphNxN_R("N",&N[0][0],&XDIM,&YDIM,&Nreq);
  DefineGraphNxN_RxR("NU",&NU[0][0][0],&XDIM,&YDIM,&NUreq);

  StartMenu("LG",1);
  DefineFunction("init",init);
  DefineFunction("init shear",initShear);
  StartMenu("Measure",0);
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
	setrho();
      }
    } else sleep(1);
  }
}

