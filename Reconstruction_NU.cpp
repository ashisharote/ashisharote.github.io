//		---------------------------------Reconstruction(Non_Uniform_Grid)--------------------------------------//
//		--for heart--
//		x_centre=0.5; 	y_centre=0.32;   R=0.17;
//
//		--for Torsten Stillke's flower--
//		L=2.0; H=2.0; R=0.0; with change in R you can find wiggles
//
//		_______________________________________________________________________________________________________//
// Need to find interface based on phi value not count (phi>eps || phi<1-eps)
#include<stdio.h>
#include<math.h>
int i,j,count,p,q,a,nx_I,nx_II,ny_I,ny_II;
float x[500],y[500],L,H,x_centre,y_centre,phi[500][500],xv[5],yv[5],xb[500],yb[500],R;
float nr_x[500][500],nr_y[500][500],beta[500][500],x1_int[500],y1_int[500],x2_int[500],y2_int[500];
float sl[500][500],sb[500][500],sr[500][500],st[500][500],nr_x0[500][500],nr_y0[500][500];
float xp[500],yp[500],g[500][500];
float alp_I,alp_II,dx_I[500],dx[500],dy[500],dy_II[500],sumI=0.0,sumII=0.0;
//		--------------------------------------declaring vertices--------------------------------------//
void cv()			//				Cell Vertices
				// 		4(top left)		3(top right)
 				//		1(bottom left)	 	2(bottom right)  
{
	xv[1]=xp[i];		yv[1]=yp[j];
	xv[2]=xp[i+1];		yv[2]=yp[j];
	xv[3]=xp[i+1];		yv[3]=yp[j+1];
	xv[4]=xp[i];		yv[4]=yp[j+1];
}

main()
{
//		--------------------------------domain and essentials--------------------------------------//
L=2.0;
H=2.0;
nx_I=15;
nx_II=nx_I+15;
ny_I=15;
ny_II=ny_I+15;
x_centre=1.0;
y_centre=0.64;
R=0.34;
alp_I=0.9;
alp_II=1.1;

//		------------------------------------X-Domain I-------------------------
for(i=0;i<=nx_I-1;i++)
	sumI=sumI+pow(alp_I,i);
	
dx[1]=(L/2.0)/sumI;
for(i=2;i<=nx_I;i++)
	dx[i]=alp_I*dx[i-1];
//		------------------------------------X-Domain II-------------------------
for(i=0;i<=(nx_II-1-nx_I);i++)
	sumII=sumII+pow(alp_II,i);
		
dx[nx_I+1]=(L/2.0)/sumII;
for(i=nx_I+2;i<=nx_II;i++)
	dx[i]=alp_II*dx[i-1];

//		----------------------------------PLIC X-Grid-----------------------------------------------//

xp[1]=0.0;
for(i=2;i<=nx_I+1;i++)
	xp[i]=xp[i-1]+dx[i-1];

xp[nx_I+1]=L/2.0;
for(i=nx_I+2;i<=nx_II+1;i++)
	xp[i]=xp[i-1]+dx[i-1];

sumI=0.0;
sumII=0.0;
//		------------------------------------Y-Domain I-------------------------
for(j=0;j<=ny_I-1;j++)
	sumI=sumI+pow(alp_I,j);
	
dy[1]=(H/2.0)/sumI;
for(j=2;j<=ny_I;j++)
	dy[j]=alp_I*dy[j-1];
//		------------------------------------Y-Domain II-------------------------
for(j=0;j<=(ny_II-1-ny_I);j++)
	sumII=sumII+pow(alp_II,j);
		
dy[ny_I+1]=(H/2.0)/sumII;
for(j=ny_I+2;j<=ny_II;j++)
	dy[j]=alp_II*dy[j-1];
sumI=0.0;
sumII=0.0;

//		----------------------------------PLIC Y-Grid-----------------------------------------------//

yp[1]=0.0;
for(j=2;j<=ny_I+1;j++)
	yp[j]=yp[j-1]+dy[j-1];

yp[ny_I+1]=H/2.0;
for(j=ny_I+2;j<=ny_II+1;j++)
	yp[j]=yp[j-1]+dy[j-1];

//		--------------------------------------FVM Grid--------------------------------------------//
x[1]=0.5*dx[1];
for(i=2;i<=nx_II;i++)
	x[i]=x[i-1]+(dx[i-1]+dx[i])/2.0;

y[1]=0.5*dy[1];
for(j=2;j<=ny_II;j++)
	y[j]=y[j-1]+(dy[j-1]+dy[j])/2.0;

//		-------------------------------------Patching-intializing---------------------------------//	
float xt,yt,dmx,dmy;
for(i=1;i<=nx_II;i++)
	{
	for(j=1;j<=ny_II;j++)
		{	
			count=0;
			cv();
			for(p=1;p<=4;p++)
				{
				//---------------------------Heart------------------------------------------------------------------
				pow(xv[p]-x_centre,2.0)+pow(5.0*(yv[p]-y_centre)/4.0-sqrt(fabs(xv[p]-x_centre)),2)-R<=0.0?count++:0;
				//---------------------------Circle-----------------------------------------------------------------
//				pow(xv[p]-x_centre,2.0)+pow(yv[p]-y_centre,2.0)-R*R<=0.0?count++:0;
				//---------------------------Square(rotated)--------------------------------------------------------
//				fabs(xv[p]-x_centre)+fabs(yv[p]-y_centre)-R<=0.0?count++:0;	
				//-----------------------Egg-Peanut----------------------------------------------------------------- 		
//				(pow(xv[p]-x_centre-1,2.0)+pow(yv[p]-y_centre,2.0))*(pow(xv[p]+x_centre-1,2.0)+pow(yv[p]-y_centre,2.0))-R<=0.0?count++:0;	
				//-----------------------Torsten Sillke's flower----------------------------------------------------
//				pow((pow((xv[p]-x_centre),2.0)+pow((yv[p]-y_centre),2.0)),3.0)-4.0*(pow((xv[p]-x_centre),2.0))*pow((yv[p]-y_centre),2.0)-R<=0.0?count++:0;
				}
				if(count==4)
				phi[i][j]=1.0;
			else if(count==0)
				phi[i][j]=0.0;
			else
			{	
				
				int mx=300,my=300;
				dmx=dx[i]/mx;	
				dmy=dy[j]/my;
				count=0;
					for(p=1;p<=mx;p++)
						{
						for(q=1;q<=my;q++)
							{
							xt=xv[1]+(p-0.5)*dmx;
							yt=yv[1]+(q-0.5)*dmy;
						
					//---------------------------Heart------------------------------------------------------------------
					pow(xt-x_centre,2.0)+pow(5.0*(yt-y_centre)/4.0-sqrt(fabs(xt-x_centre)),2)-R<=0.0?count++:0;
					//---------------------------------Circle-----------------------------------------------------		 
//					pow(xt-x_centre,2.0)+pow(yt-y_centre,2.0)-R*R<=0.0?count++:0;	
					//---------------------------Square(rotated)--------------------------------------------------------
//					fabs(xt-x_centre)+fabs(yt-y_centre)-R<=0.0?count++:0;		
					//-----------------------Egg-Peanut----------------------------------------------------------------- 	
//					(pow(xt-x_centre-1,2.0)+pow(yt-y_centre,2.0))*(pow(xt+x_centre-1,2.0)+pow(yt-y_centre,2.0))-R<=0.0?count++:0; 
					//-----------------------Torsten Sillke's flower----------------------------------------------------
//					pow((pow((xt-x_centre),2.0)+pow((yt-y_centre),2.0)),3.0)-4.0*(pow((xt-x_centre),2.0))*pow((yt-y_centre),2.0)-R<=0.0?count++:0;
					 		}
						}
				phi[i][j]=(count*1.0)/(mx*my); //total vol. fraction in cell [i,j]=No.of cells inside/total cells in cell [i,j]
			
			}				
		}
	}


FILE *f2;
f2=fopen("rcnst.dat","w");
//		--------------------------------------Reconstruction----------------------------------------------	
for(i=1;i<=nx_II;i++)
	{
	for(j=1;j<=ny_II;j++)
		{
			count=0;
			cv();
			for(p=1;p<=4;p++)
				//---------------------------Heart------------------------------------------------------------------
				pow(xv[p]-x_centre,2.0)+pow(5.0*(yv[p]-y_centre)/4.0-sqrt(fabs(xv[p]-x_centre)),2)-R<=0.0?count++:0;
				//---------------------------Circle-----------------------------------------------------------------
//				pow(xv[p]-x_centre,2.0)+pow(yv[p]-y_centre,2.0)-R*R<=0.0?count++:0;
				//---------------------------Square(rotated)--------------------------------------------------------
//				fabs(xv[p]-x_centre)+fabs(yv[p]-y_centre)-R<=0.0?count++:0;	
				//-----------------------Egg-Peanut----------------------------------------------------------------- 		
//				(pow(xv[p]-x_centre-1,2.0)+pow(yv[p]-y_centre,2.0))*(pow(xv[p]+x_centre-1,2.0)+pow(yv[p]-y_centre,2.0))-R<=0.0?count++:0; 
				//-----------------------Torsten Sillke's flower----------------------------------------------------
//				pow((pow((xv[p]-x_centre),2.0)+pow((yv[p]-y_centre),2.0)),3.0)-4.0*(pow((xv[p]-x_centre),2.0))*pow((yv[p]-y_centre),2.0)-R<=0.0?count++:0;
			if(count==4 || count==0)
				continue;
			else
			{
// 		-----------------------------------Normals_and_angle_with_horizontal---------------------------
				nr_x[i][j]=(1/(8.0*dx[i]))*(phi[i+1][j+1]+2.0*phi[i+1][j]+phi[i+1][j-1]-phi[i-1][j+1]-2.0*phi[i-1][j]-phi[i-1][j-1]);
				nr_y[i][j]=(1/(8.0*dy[j]))*(phi[i+1][j+1]+2.0*phi[i][j+1]+phi[i-1][j+1]-phi[i+1][j-1]-2.0*phi[i][j-1]-phi[i-1][j-1]);
//		---------------------------------------Normalized normals-----------------------------------------
			//	nr_x[i][j]=nr_x0[i][j]/sqrt(nr_x0[i][j]*nr_x0[i][j]+nr_y0[i][j]*nr_y0[i][j]);
			//	nr_y[i][j]=nr_y0[i][j]/sqrt(nr_x0[i][j]*nr_x0[i][j]+nr_y0[i][j]*nr_y0[i][j]);
//		---------------------------------------Angle of interface with the horizontal---------------------
				g[i][j]=(atan(-nr_x[i][j]/nr_y[i][j])); 	//radians
				beta[i][j]=atan((dx[i]/dy[j])*tan(g[i][j])); 	//radians
//		---------------------------------------When beta is -ve-------------------------------------------
				if((nr_x[i][j]>0.0 && nr_y[i][j]>0.0) ||(nr_x[i][j]<0.0 && nr_y[i][j]<0.0))
				beta[i][j]=M_PI/2.0+(beta[i][j]);	

//		-------------------Case Identification(Rudman Algorithm)feat. Shashwat_sir-----------------------//
	if(beta[i][j]<(M_PI/4.0))
		{
			if (phi[i][j]<=0.5*tan(beta[i][j]))
//		----------------------------------CASE I---------------------------------
				{
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					sr[i][j]=sqrt(phi[i][j]*2.0*((tan(beta[i][j]))));
					sb[i][j]=sqrt(phi[i][j]*(1.0/tan(beta[i][j]))*2.0);
					sl[i][j]=0.0;
					st[i][j]=0.0;
					x1_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
										
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
//					printf("\ncase I \t q4 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{
//					third quadrant
					sl[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					sb[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);
					sr[i][j]=0.0;
					st[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
//					printf("\ncase I \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					st[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					sl[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);	
					sr[i][j]=0.0;
					sb[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
//					printf("\ncase I \t q2 %d_%d",i,j);
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					sr[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					st[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);
					sb[i][j]=0.0;
					sl[i][j]=0.0;
					x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
//					printf("\ncase I \t q1 %d_%d",i,j);
					}					
			
				}
			else if(phi[i][j]<=(1.0-0.5*tan(beta[i][j])))	
				{
//		----------------------------------CASE II--------------------------------
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					sr[i][j]=phi[i][j]+0.5*(tan(beta[i][j]));
					sl[i][j]=phi[i][j]-0.5*(tan(beta[i][j]));
					sb[i][j]=1.0;
					st[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i];
					y1_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
										
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
//					printf("\ncase II \t q4 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{
//					third quadrant
					
					st[i][j]=phi[i][j]-0.5*(tan(beta[i][j]));
					sb[i][j]=phi[i][j]+0.5*(tan(beta[i][j]));
					sr[i][j]=0.0;
					sl[i][j]=1.0;
					x1_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y2_int[j]=y[j]+0.5*dy[j];
//					printf("\ncase II \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					sr[i][j]=phi[i][j]-0.5*(tan(beta[i][j]));
					sl[i][j]=phi[i][j]+0.5*(tan(beta[i][j]));
					sb[i][j]=0.0;
					st[i][j]=1.0;
					
					x1_int[i]=x[i]-0.5*dx[i];
					y1_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
//					printf("\ncase II \t q2 %d_%d",i,j);	
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					st[i][j]=phi[i][j]+0.5*(tan(beta[i][j]));
					sb[i][j]=phi[i][j]-0.5*(tan(beta[i][j]));
					sr[i][j]=1.0;
					sl[i][j]=0.0;
					x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y2_int[j]=y[j]-0.5*dy[j];
//					printf("\ncase II \t q1 %d_%d",i,j);	
					}
				}
			else 
				{
//		---------------------------------CASE IV---------------------------------
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					st[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sl[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					sr[i][j]=1.0;
					sb[i][j]=1.0;
					x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
//					printf("\ncase IV \t q4 %d_%d",i,j);
					}	
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{
//					third quadrant
					sr[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					st[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					sl[i][j]=1.0;
					sb[i][j]=1.0;
					x1_int[i]=x[i]+0.5*dx[i];
					y1_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y2_int[j]=y[j]+0.5*dy[j];
//					printf("\ncase IV \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					sb[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sr[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					st[i][j]=1.0;
					sl[i][j]=1.0; 
					x1_int[i]=x[i]+0.5*dx[i];
					y1_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y2_int[j]=y[j]-0.5*dy[j];
//					printf("\ncase IV \t q2 %d_%d",i,j);
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					sl[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sb[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					
					x1_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
//					printf("\ncase IV \t q1 %d_%d",i,j);
					}
				}
		}
	else 
		{
			if(phi[i][j]<=0.5*(1.0/(tan(beta[i][j]))))
//		-------------------------------CASE I------------------------------------
				{
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					sr[i][j]=sqrt(phi[i][j]*2.0*((tan(beta[i][j]))));
					sb[i][j]=sqrt(phi[i][j]*(1.0/tan(beta[i][j]))*2.0);
					sl[i][j]=0.0;
					st[i][j]=0.0;
					x1_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
										
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
//					printf("\ncase I \t q4 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{
//					third quadrant
					sl[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					sb[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);
					sr[i][j]=0.0;
					st[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
//					printf("\ncase I \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					st[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					sl[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);	
					sr[i][j]=0.0;
					sb[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
//					printf("\ncase I \t q2 %d_%d",i,j);
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					sr[i][j]=sqrt(phi[i][j]*2.0*(1.0/(tan(beta[i][j]))));
					st[i][j]=sqrt(phi[i][j]*(tan(beta[i][j]))*2.0);
					sb[i][j]=0.0;
					sl[i][j]=0.0;
					x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
//					printf("\ncase I \t q1 %d_%d",i,j);
					}					
				}
			else if (phi[i][j]<=(1.0-0.5*(1.0/(tan(beta[i][j])))))
				{
//		----------------------------CASE III-------------------------------------
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					st[i][j]=phi[i][j]-0.5*(1.0/(tan(beta[i][j])));
					sb[i][j]=phi[i][j]+0.5*(1.0/(tan(beta[i][j])));
					
					x1_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y2_int[j]=y[j]+0.5*dy[j];
//					printf("\ncase III \t q4 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{	
//					third quadrant
			 		sl[i][j]=phi[i][j]+0.5*(1.0/(tan(beta[i][j])));
					sr[i][j]=phi[i][j]-0.5*(1.0/(tan(beta[i][j])));
					st[i][j]=0.0;
					sb[i][j]=1.0;
					
					x1_int[i]=x[i]-0.5*dx[i];
					y1_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
//					printf("\ncase III \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					st[i][j]=phi[i][j]+0.5*(1.0/(tan(beta[i][j])));
					sb[i][j]=phi[i][j]-0.5*(1.0/(tan(beta[i][j])));
					sl[i][j]=1.0;
					sr[i][j]=0.0;
					x1_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y2_int[j]=y[j]+0.5*dy[j];
//					printf("\ncase III \t q2 %d_%d",i,j);
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					sl[i][j]=phi[i][j]-0.5*(1.0/(tan(beta[i][j])));
					sr[i][j]=phi[i][j]+0.5*(1.0/(tan(beta[i][j])));
					
					x1_int[i]=x[i]-0.5*dx[i];
					y1_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
					
					x2_int[i]=x[i]+0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
//					printf("\ncase III \t q1 %d_%d",i,j);
					}
				}
			else 	
				{
//		----------------------------CASE IV--------------------------------------
				if(nr_x[i][j]>0.0 && nr_y[i][j]<0.0)
					{
//					fourth quadrant
					st[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sl[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					sr[i][j]=1.0;
					sb[i][j]=1.0;
					x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];
//					printf("\ncase IV \t q4 %d_%d",i,j);
					}	
				else if(nr_x[i][j]<0.0 && nr_y[i][j]<0.0)
					{
//					third quadrant
					sr[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					st[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					sl[i][j]=1.0;
					sb[i][j]=1.0;
					x1_int[i]=x[i]+0.5*dx[i];
					y1_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y2_int[j]=y[j]+0.5*dy[j];
//					printf("\ncase IV \t q3 %d_%d",i,j);
					}
				else if(nr_x[i][j]<0.0 && nr_y[i][j]>0.0)
					{
//					second quadrant
					sb[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sr[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					st[i][j]=1.0;
					sl[i][j]=1.0; 
					x1_int[i]=x[i]+0.5*dx[i];
					y1_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y2_int[j]=y[j]-0.5*dy[j];
//					printf("\ncase IV \t q2 %d_%d",i,j);
					}
				else if(nr_x[i][j]>0.0 && nr_y[i][j]>0.0)
					{
//					first quadrant
					sl[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(1.0/(tan(beta[i][j]))));
					sb[i][j]=1.0-sqrt(2.0*(1.0-phi[i][j])*(tan(beta[i][j])));
					sr[i][j]=1.0;
					st[i][j]=1.0;
					x1_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dx[i];
					y1_int[j]=y[j]-0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i];
					y2_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
//					printf("\ncase IV \t q1 %d_%d",i,j);
					}
				}
			}			
//		----------------------------------zero normals calculations----------------------------------
	if(nr_x[i][j]==0.0 && nr_y[i][j]>0.0)
		{
			sr[i][j]=phi[i][j];
			sl[i][j]=sr[i][j];
			st[i][j]=1.0;
			sb[i][j]=0.0;
			x1_int[i]=x[i]+0.5*dx[i];
			y1_int[j]=y[j]+0.5*dy[j]-sr[i][j]*dy[j];
			
			x2_int[i]=x[i]-0.5*dx[i];
			y2_int[j]=y[j]+0.5*dy[j]-sl[i][j]*dy[j];
//			printf("\nnx=0 ny+ \t  %d_%d",i,j);	
		}else if(nr_x[i][j]==0.0 && nr_y[i][j]<0.0)
			{
				sr[i][j]=phi[i][j];
				sl[i][j]=sr[i][j];
	 			st[i][j]=0.0;
				sb[i][j]=1.0;
				x1_int[i]=x[i]+0.5*dx[i];
				y1_int[j]=y[j]-0.5*dy[j]+sr[i][j]*dy[j];
					
				x2_int[i]=x[i]-0.5*dx[i];
				y2_int[j]=y[j]-0.5*dy[j]+sl[i][j]*dy[j];	
//				printf("\nnx=0 ny- \t  %d_%d",i,j);
			}else if(nr_y[i][j]==0.0 && nr_x[i][j]<0.0)
				{
					st[i][j]=phi[i][j];
					sb[i][j]=st[i][j];
					sr[i][j]=0.0;
					sl[i][j]=1.0;
					x1_int[i]=x[i]-0.5*dx[i]+st[i][j]*dx[i];
					y1_int[j]=y[j]+0.5*dy[j];
					
					x2_int[i]=x[i]-0.5*dx[i]+sb[i][j]*dx[i];
					y2_int[j]=y[j]-0.5*dy[j];	
//					printf("\nny=0 nx- \t  %d_%d",i,j);
				}else if(nr_y[i][j]==0.0 && nr_x[i][j]>0.0)
					{
						st[i][j]=phi[i][j];
						sb[i][j]=st[i][j];
						sr[i][j]=1.0;
						sl[i][j]=0.0;
						x1_int[i]=x[i]+0.5*dx[i]-st[i][j]*dx[i];
						y1_int[j]=y[j]+0.5*dy[j];
						
						x2_int[i]=x[i]+0.5*dx[i]-sb[i][j]*dy[j];
						y2_int[j]=y[j]-0.5*dy[j];	
//						printf("\nny=0 nx+ \t  %d_%d",i,j);
					}
							
		}
		
			fprintf(f2,"%f\t%f\n%f\t%f\n\n",x1_int[i],y1_int[j],x2_int[i],y2_int[j]);
//			printf("\nsb[%d][%d]=%f\tst[%d][%d]=%f\tsr[%d][%d]=%f\tsl[%d][%d]=%f\n",i,j,sb[i][j],i,j,st[i][j],i,j,sr[i][j],i,j,sl[i][j]);
	}
}	
fclose(f2);	
//		-------------------------------------------File Writing-----------------------------------------//

FILE *f1;
f1=fopen("Plic_grid1.dat","w");
//fprintf(f1,"\nVARIABLES=""X"",""Y"", \nZONE I=%d,J=%d,ZONETYPE=ORDERED,DATAPACKING=POINT\n",nx+1,ny+1);
for(i=1;i<=nx_II+1;i++)
	{
	for(j=1;j<=ny_II+1;j++)
		{
			fprintf(f1,"\t%f\t%f\n",xp[i],yp[j]);
			
		}
		fprintf(f1,"\n");
	}
fclose(f1);

f1=fopen("Plic_grid2.dat","w");
//fprintf(f1,"\nVARIABLES=""X"",""Y"", \nZONE I=%d,J=%d,ZONETYPE=ORDERED,DATAPACKING=POINT\n",nx+1,ny+1);
for(j=1;j<=ny_II+1;j++)
	{
	for(i=1;i<=nx_II+1;i++)
		{
			fprintf(f1,"\t%f\t%f\n",xp[i],yp[j]);			
		}
		fprintf(f1,"\n");
	}
fclose(f1);

FILE *f3;
f3=fopen("Patch.dat","w");
//fprintf(f3,"\nVARIABLES=""X"",""Y"",""phi"" \nZONE I=%d,J=%d,ZONETYPE=ORDERED,DATAPACKING=POINT\n",nx,ny);
for(i=1;i<=nx_II;i++)
	{
	for(j=1;j<=ny_II;j++)
		{
			fprintf(f3,"\t%f\t%f\t%f\n",x[i],y[j],phi[i][j]);
			
		}
	//	fprintf(f1,"\n");
	}
fclose(f3);
}
