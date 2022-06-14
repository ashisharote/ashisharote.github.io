
#include<math.h>
#include<stdio.h>
#include<omp.h>
double ae,aw,as,an;
int i,j,nx,ny;
double dx,dy,L,H,dt,Re,beta,x[600],y[600],u[600][600],v[600][600],uc[600][600],vc[600][600],ut[600][600],vt[600][600],p[600][600],pc[600][600],Res,corr;
double sum,sum1,sqr,sqr1,diff,diff1,RMS,RMS1,u_in,ap,ue,uw,un,us,ve,vw,vn,vs,uwa,uea,usa,una,vsa,vea,vna,vwa,fyad,fydf,fxad,fxdf;
double s[600][600],uwp,uep,vnp,vsp,ute,utw,vts,vtn,pnew[600][600],unew[600][600],vnew[600][600];
double aw_nav,ae_nav,an_nav,as_nav,ap_nav;
double uwaq,ueaq,usaq,unaq,vwaq,vnaq,vsaq,veaq,dtmax;
double blend,umax,vmax,t,co;
double sumI,sumII,alp_I,apg;
double diffv,sqrv,RMSv,sumv;
double ue0,uw0,us0,un0,vn0,vs0,ve0,vw0,uwa0,uea0,una0,usa0,vwa0,vea0,vna0,vsa0,uwaq0,ueaq0,unaq0,usaq0,vwaq0,veaq0,vnaq0,vsaq0,fxad0,fyad0,fxdf0,fydf0,u0[600][600],v0[600][600],p0[600][600];
main()
{
omp_set_num_threads(6);
L=1.0;
H=1.0;
nx=150;
ny=150;
dt=0.0001;
dtmax=dt;
Re=7500.0;
beta=1.0;//relaxation factor

dx=L/nx;
x[1]=0.5*dx;
for(i=2;i<=nx;i++)
    {
        x[i]=x[i-1]+0.5*(dx);
    //    printf("x=%f\n",x[i]);
    }
dy=H/ny;
y[1]=0.5*dy;
for(j=2;j<=ny;j++)
    {
        y[j]=y[j-1]+0.5*(dy);
    //    printf("y=%f\n",y[j]);
    }

//________________________________________________Initialization_____________________________________________________
for(i=1;i<=nx;i++)
    {
    for(j=1;j<=ny;j++)
        {
       		u[i][j]=0.0;
       		v[i][j]=0.0;
            u0[i][j]=0.0;
       		v0[i][j]=0.0;
       		ut[i][j]=0.0;
       		vt[i][j]=0.0;
       		p[i][j]=0.0;
       		p0[i][j]=0.0;
       		pc[i][j]=0.0;
        }
    }

u_in=1.0;
int I=0;
//_______________________________________________Solution________________________________________________________
do
{
	sum=0.0;
	//Boundary conditions
	for(i=1;i<=nx;i++)
		{
			u[i][0]=-u[i][1];
            u[i][-1]=-u[i][0];
			u[i][ny+1]=2.0*u_in-u[i][ny];
            u[i][ny+2]=2.0*u_in-u[i][ny+1];
            v[i][0]=-v[i][1];
            v[i][-1]=-v[i][0];
			v[i][ny+1]=-v[i][ny];
			v[i][ny+2]=-v[i][ny+1];

            u0[i][0]=-u0[i][1];
            u0[i][-1]=-u0[i][0];
			u0[i][ny+1]=2.0*u_in-u0[i][ny];
            u0[i][ny+2]=2.0*u_in-u0[i][ny+1];
            v0[i][0]=-v0[i][1];
            v0[i][-1]=-v0[i][0];
			v0[i][ny+1]=-v0[i][ny];
			v0[i][ny+2]=-v0[i][ny+1];
		}
	for(j=1;j<=ny;j++)
		{
			u[0][j]=-u[1][j];
            u[-1][j]=-u[0][j];
			u[nx+1][j]=-u[nx][j];
			u[nx+2][j]=-u[nx+1][j];

			v[0][j]=-v[1][j];
			v[-1][j]=-v[0][j];
			v[nx+1][j]=-v[nx][j];
			v[nx+2][j]=-v[nx+1][j];

			u0[0][j]=-u0[1][j];
            u0[-1][j]=-u0[0][j];
			u0[nx+1][j]=-u0[nx][j];
			u0[nx+2][j]=-u0[nx+1][j];

			v0[0][j]=-v0[1][j];
			v0[-1][j]=-v0[0][j];
			v0[nx+1][j]=-v0[nx][j];
			v0[nx+2][j]=-v0[nx+1][j];

		}

#pragma omp parallel for schedule(dynamic) private(i,j,un,us,ve,vw,ue,uw,vn,vs,aw_nav,ae_nav,as_nav,an_nav,ap_nav,uwa,uea,una,usa,vwa,vea,vna,vsa,uwaq,ueaq,unaq,usaq,vwaq,veaq,vnaq,vsaq,fxad,fyad,fxdf,fydf,ue0,uw0,us0,un0,vn0,vs0,ve0,vw0,uwa0,uea0,una0,usa0,vwa0,vea0,vna0,vsa0,uwaq0,ueaq0,unaq0,usaq0,vwaq0,veaq0,vnaq0,vsaq0,fxad0,fyad0,fxdf0,fydf0)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{

				ue=(u[i][j]*dx+u[i+1][j]*dx)/(dx+dx);
				uw=(u[i][j]*dx+u[i-1][j]*dx)/(dx+dx);
				un=(u[i][j]*dy+u[i][j+1]*dy)/(dy+dy);
				us=(u[i][j]*dy+u[i][j-1]*dy)/(dy+dy);

				ve=(v[i][j]*dx+v[i+1][j]*dx)/(dx+dx);
				vw=(v[i][j]*dx+v[i-1][j]*dx)/(dx+dx);
				vn=(v[i][j]*dy+v[i][j+1]*dy)/(dy+dy);
				vs=(v[i][j]*dy+v[i][j-1]*dy)/(dy+dy);

				ue0=(u0[i][j]*dx+u0[i+1][j]*dx)/(dx+dx);
				uw0=(u0[i][j]*dx+u0[i-1][j]*dx)/(dx+dx);
                un0=(u0[i][j]*dy+u0[i][j+1]*dy)/(dy+dy);
				us0=(u0[i][j]*dy+u0[i][j-1]*dy)/(dy+dy);

				vn0=(v0[i][j]*dy+v0[i][j+1]*dy)/(dy+dy);
				vs0=(v0[i][j]*dy+v0[i][j-1]*dy)/(dy+dy);
                ve0=(v0[i][j]*dy+v0[i+1][j]*dy)/(dy+dy);
				vw0=(v0[i][j]*dy+v0[i-1][j]*dy)/(dy+dy);

                aw_nav=dy/(dx);
                ae_nav=aw_nav;
                as_nav=dx/(dy);
                an_nav=as_nav;

                ap_nav=-(ae_nav+aw_nav+an_nav+as_nav);
				//FOU
				if(ue>=0.0)
					{
						uea=u[i][j];
						vea=v[i][j];
					}else
						{
							uea=u[i+1][j];
							vea=v[i+1][j];
						}
				if(uw>=0.0)
					{
						uwa=u[i-1][j];
						vwa=v[i-1][j];
					}else
						{
							uwa=u[i][j];
							vwa=v[i][j];
						}
				if(vn>=0.0)
					{
						una=u[i][j];
						vna=v[i][j];
					}else
						{
							una=u[i][j+1];
							vna=v[i][j+1];
						}
				if(vs>=0.0)
					{
						usa=u[i][j-1];
						vsa=v[i][j-1];
					}else
						{
							usa=u[i][j];
							vsa=v[i][j];
						}
				if(uw>0.0)
				{
					uwaq=(3.0/8.0)*u[i][j]+(3.0/4.0)*u[i-1][j]-(1.0/8.0)*u[i-2][j];
					vwaq=(3.0/8.0)*v[i][j]+(3.0/4.0)*v[i-1][j]-(1.0/8.0)*v[i-2][j];
				}else
					{
						uwaq=(3.0/8.0)*u[i-1][j]+(3.0/4.0)*u[i][j]-(1.0/8.0)*u[i+1][j];
						vwaq=(3.0/8.0)*v[i-1][j]+(3.0/4.0)*v[i][j]-(1.0/8.0)*v[i+1][j];
					}
				if(ue>0.0)
				{
					ueaq=(3.0/8.0)*u[i+1][j]+(3.0/4.0)*u[i][j]-(1.0/8.0)*u[i-1][j];
					veaq=(3.0/8.0)*v[i+1][j]+(3.0/4.0)*v[i][j]-(1.0/8.0)*v[i-1][j];
				}else
					{
						ueaq=(3.0/8.0)*u[i][j]+(3.0/4.0)*u[i+1][j]-(1.0/8.0)*u[i+2][j];
						veaq=(3.0/8.0)*v[i][j]+(3.0/4.0)*v[i+1][j]-(1.0/8.0)*v[i+2][j];
					}
				if(vn>0.0)
				{
					unaq=(3.0/8.0)*u[i][j+1]+(3.0/4.0)*u[i][j]-(1.0/8.0)*u[i][j-1];
					vnaq=(3.0/8.0)*v[i][j+1]+(3.0/4.0)*v[i][j]-(1.0/8.0)*v[i][j-1];
				}else
					{
						unaq=(3.0/8.0)*u[i][j]+(3.0/4.0)*u[i][j+1]-(1.0/8.0)*u[i][j+2];
						vnaq=(3.0/8.0)*v[i][j]+(3.0/4.0)*v[i][j+1]-(1.0/8.0)*v[i][j+2];
					}
				if(vs>0.0)
				{
					usaq=(3.0/8.0)*u[i][j]+(3.0/4.0)*u[i][j-1]-(1.0/8.0)*u[i][j-2];
					vsaq=(3.0/8.0)*v[i][j]+(3.0/4.0)*v[i][j-1]-(1.0/8.0)*v[i][j-2];
				}else
					{
						usaq=(3.0/8.0)*u[i][j-1]+(3.0/4.0)*u[i][j]-(1.0/8.0)*u[i][j+1];
						vsaq=(3.0/8.0)*v[i][j-1]+(3.0/4.0)*v[i][j]-(1.0/8.0)*v[i][j+1];
					}

				if(ue0>=0.0)
					{
						uea0=u0[i][j];
						vea0=v0[i][j];
					}else
						{
							uea0=u0[i+1][j];
							vea0=v0[i+1][j];
						}
				if(uw0>=0.0)
					{
						uwa0=u0[i-1][j];
						vwa0=v0[i-1][j];
					}else
						{
							uwa0=u0[i][j];
							vwa0=v0[i][j];
						}
				if(vn0>=0.0)
					{
						una0=u0[i][j];
						vna0=v0[i][j];
					}else
						{
							una0=u0[i][j+1];
							vna0=v0[i][j+1];
						}
				if(vs0>=0.0)
					{
						usa0=u0[i][j-1];
						vsa0=v0[i][j-1];
					}else
						{
							usa0=u0[i][j];
							vsa0=v0[i][j];
						}
				if(uw0>0.0)
				{
					uwaq0=(3.0/8.0)*u0[i][j]+(3.0/4.0)*u0[i-1][j]-(1.0/8.0)*u0[i-2][j];
					vwaq0=(3.0/8.0)*v0[i][j]+(3.0/4.0)*v0[i-1][j]-(1.0/8.0)*v0[i-2][j];
				}else
					{
						uwaq0=(3.0/8.0)*u0[i-1][j]+(3.0/4.0)*u0[i][j]-(1.0/8.0)*u0[i+1][j];
						vwaq0=(3.0/8.0)*v0[i-1][j]+(3.0/4.0)*v0[i][j]-(1.0/8.0)*v0[i+1][j];
					}
				if(ue0>0.0)
				{
					ueaq0=(3.0/8.0)*u0[i+1][j]+(3.0/4.0)*u0[i][j]-(1.0/8.0)*u0[i-1][j];
					veaq0=(3.0/8.0)*v0[i+1][j]+(3.0/4.0)*v0[i][j]-(1.0/8.0)*v0[i-1][j];
				}else
					{
						ueaq0=(3.0/8.0)*u0[i][j]+(3.0/4.0)*u0[i+1][j]-(1.0/8.0)*u0[i+2][j];
						veaq0=(3.0/8.0)*v0[i][j]+(3.0/4.0)*v0[i+1][j]-(1.0/8.0)*v0[i+2][j];
					}
				if(vn0>0.0)
				{
					unaq0=(3.0/8.0)*u0[i][j+1]+(3.0/4.0)*u0[i][j]-(1.0/8.0)*u0[i][j-1];
					vnaq0=(3.0/8.0)*v0[i][j+1]+(3.0/4.0)*v0[i][j]-(1.0/8.0)*v0[i][j-1];
				}else
					{
						unaq0=(3.0/8.0)*u0[i][j]+(3.0/4.0)*u0[i][j+1]-(1.0/8.0)*u0[i][j+2];
						vnaq0=(3.0/8.0)*v0[i][j]+(3.0/4.0)*v0[i][j+1]-(1.0/8.0)*v0[i][j+2];
					}
				if(vs0>0.0)
				{
					usaq0=(3.0/8.0)*u0[i][j]+(3.0/4.0)*u0[i][j-1]-(1.0/8.0)*u0[i][j-2];
					vsaq0=(3.0/8.0)*v0[i][j]+(3.0/4.0)*v0[i][j-1]-(1.0/8.0)*v0[i][j-2];
				}else
					{
						usaq0=(3.0/8.0)*u0[i][j-1]+(3.0/4.0)*u0[i][j]-(1.0/8.0)*u0[i][j+1];
						vsaq0=(3.0/8.0)*v0[i][j-1]+(3.0/4.0)*v0[i][j]-(1.0/8.0)*v0[i][j+1];
					}

				blend=0.6;
                    // X-Diffusion
					fxdf=(dt/(dx*dy*Re))*(aw_nav*u[i-1][j]+ae_nav*u[i+1][j]+as_nav*u[i][j-1]+an_nav*u[i][j+1]);
					// Y-Diffusion
					fydf=(dt/(dx*dy*Re))*(aw_nav*v[i-1][j]+ae_nav*v[i+1][j]+as_nav*v[i][j-1]+an_nav*v[i][j+1]);

					fxdf0=(dt/(dx*dy*Re))*(aw_nav*u0[i-1][j]+ae_nav*u0[i+1][j]+as_nav*u0[i][j-1]+an_nav*u0[i][j+1]);
					// Y-Diffusion
					fydf0=(dt/(dx*dy*Re))*(aw_nav*v0[i-1][j]+ae_nav*v0[i+1][j]+as_nav*v0[i][j-1]+an_nav*v0[i][j+1]);


				fxad=blend*((dt/dx)*(uw*uw-ue*ue)+(dt/dy)*(us*vs-un*vn))+(1.0-blend)*((dt/dx)*(uwaq*uw-ueaq*ue)+(dt/dy)*(usaq*vs-unaq*vn));	// X-Advection
                fyad=blend*((dt/dx)*(vw*uw-ve*ue)+(dt/dy)*(vs*vs-vn*vn))+(1.0-blend)*((dt/dx)*(vwaq*uw-veaq*ue)+(dt/dy)*(vsaq*vs-vnaq*vn));	// Y-Advection


				fxad0=blend*((dt/dx)*(uw0*uw0-ue0*ue0)+(dt/dy)*(us0*vs0-un0*vn0))+(1.0-blend)*((dt/dx)*(uwaq0*uw0-ueaq0*ue0)+(dt/dy)*(usaq0*vs0-unaq0*vn0));// X-Advection previous
                fyad0=blend*((dt/dx)*(vw0*uw0-ve0*ue0)+(dt/dy)*(vs0*vs0-vn0*vn0))+(1.0-blend)*((dt/dx)*(vwaq0*uw0-veaq0*ue0)+(dt/dy)*(vsaq0*vs0-vnaq0*vn0));// Y-Advection previous

				ut[i][j]=u[i][j]*(1.0+ap_nav*dt/(Re*dx*dy))+(fxad+fxdf)*3.0/2.0-(fxad0+fxdf0)*0.5; //Adam Bashborth O(2)
                vt[i][j]=v[i][j]*(1.0+ap_nav*dt/(Re*dx*dy))+(fyad+fydf)*3.0/2.0-(fyad0+fydf0)*0.5; //Adam Bashborth O(2)
			}
		}

	//Boundary conditions
	for(i=1;i<=nx;i++)
		{
			ut[i][0]=-ut[i][1];
			ut[i][ny+1]=2.0*u_in-ut[i][ny];
			vt[i][0]=-vt[i][1];
			vt[i][ny+1]=-vt[i][ny];
			p[i][ny+1]=p[i][ny];
			p[i][0]=p[i][1];
			p0[i][ny+1]=p0[i][ny];
			p0[i][0]=p0[i][1];

		}
	for(j=1;j<=ny;j++)
		{
			ut[0][j]=-ut[1][j];
			ut[nx+1][j]=-ut[nx][j];
			vt[0][j]=-vt[1][j];
			vt[nx+1][j]=-vt[nx][j];
			p[0][j]=p[1][j];
			p[nx+1][j]=p[nx][j];
            p0[0][j]=p0[1][j];
			p0[nx+1][j]=p0[nx][j];
		}

	//Predicted Velocities
#pragma omp parallel for schedule(dynamic) private(i,j,utw,ute,vtn,vts,uwp,uep,vsp,vnp)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
                    utw=(ut[i][j]+ut[i-1][j])/2.0;
					ute=(ut[i+1][j]+ut[i][j])/2.0;
					vtn=(vt[i][j]+vt[i][j+1])/2.0;
					vts=(vt[i][j]+vt[i][j-1])/2.0;

                uwp=utw-(dt/(dx))*(p[i][j]-p[i-1][j]);
				uep=ute-(dt/(dx))*(p[i+1][j]-p[i][j]);
				vsp=vts-(dt/(dy))*(p[i][j]-p[i][j-1]);
				vnp=vtn-(dt/(dy))*(p[i][j+1]-p[i][j]);
				s[i][j]=uwp*dy-uep*dy-vnp*dx+vsp*dx;
			}
		}

int g=0;
do
{
sum1=0.0;
	//Boundary conditions
	for(i=1;i<=nx;i++)
		{
			pc[i][ny+1]=pc[i][ny];
			pc[i][0]=pc[i][1];
		}
	for(j=1;j<=ny;j++)
		{
			pc[0][j]=pc[1][j];
			pc[nx+1][j]=pc[nx][j];
		}

//ap=dt*(2.0*dy/dx+2.0*dx/dy);

#pragma omp parallel for schedule(dynamic) private(i,j,apg,sqr1,Res,corr) reduction(+:sum1)
	for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{

                apg=(dt*dy/(dx))+(dt*dy/(dx))+(dx*dt/(dy))+(dt*dx/(dy));
                Res=s[i][j]+(dt*dy/dx)*((pc[i+1][j])+(pc[i-1][j]))+(dt*dx/dy)*((pc[i][j+1])+(pc[i][j-1]))-apg*pc[i][j];
				corr=Res/(apg);
				pc[i][j]=pc[i][j]+corr;
				sqr1=Res*Res;
				sum1=sum1+sqr1;
			}
		}
RMS1=sqrt(sum1/(nx*ny));
//if(g%1000==0)
//printf("Code is here RMS=%.10f\n",RMS1);

//printf("RMS1=%f\tpc=%f\n",RMS1,pc[nx/2][ny/2]);

g++;
}while(RMS1>=0.000001);
			//			printf("Code is here\n");
//New pessure
#pragma omp parallel for schedule(dynamic) private(i,j)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
				pnew[i][j]=p[i][j]+pc[i][j];
				pc[i][j]=0.0;
			}
		}
	//Boundary conditions
	for(i=1;i<=nx;i++)
		{
			pnew[i][ny+1]=pnew[i][ny];
			uc[i][ny+1]=-uc[i][ny];
			vc[i][ny+1]=-vc[i][ny];
			pnew[i][0]=pnew[i][1];
			uc[i][0]=-uc[i][1];
			vc[i][0]=-vc[i][1];
		}
	for(j=1;j<=ny;j++)
		{
			pnew[0][j]=pnew[1][j];
			uc[0][j]=-uc[1][j];
			vc[0][j]=-vc[1][j];
			pnew[nx+1][j]=pnew[nx][j];
			uc[nx+1][j]=-uc[nx][j];
			vc[nx+1][j]=-vc[nx][j];
		}

//New velocity
sumv=0.0;

#pragma omp parallel for schedule(dynamic) private(i,j,diff,sqr,diffv,sqrv) reduction(+:sum,sumv)
for(i=1;i<=nx;i++)
		{
		for(j=1;j<=ny;j++)
			{
				uc[i][j]=0.5*(dt/dx)*((pnew[i+1][j])-(pnew[i-1][j]));
                vc[i][j]=0.5*(dt/dy)*((pnew[i][j+1])-(pnew[i][j-1]));
			//	uc[i][j]=(dt/dx)*(pnew[i+1][j]-pnew[i-1][j])/2.0;
			//	vc[i][j]=(dt/dy)*(pnew[i][j+1]-pnew[i][j-1])/2.0;

				unew[i][j]=ut[i][j]-uc[i][j];
				vnew[i][j]=vt[i][j]-vc[i][j];
				diff=unew[i][j]-u[i][j];
				sqr=diff*diff;
				sum=sum+sqr;

                diffv=vnew[i][j]-v[i][j];
				sqrv=diffv*diffv;
				sumv=sumv+sqrv;

                u0[i][j]=u[i][j];
                v0[i][j]=v[i][j];
                p0[i][j]=p[i][j];
				u[i][j]=unew[i][j];
				v[i][j]=vnew[i][j];
				p[i][j]=pnew[i][j];
			}
		}
RMS=sqrt(sum/(nx*ny));
RMSv=sqrt(sumv/(nx*ny));
if(I%200==0)
printf("\tGSSOR_iter=%d\tu_RMS=%.10f\tv_RMS=%.10f\ttime<iter=%d>=%f\n",g,RMS,RMSv,I,t);

I++;
t=t+dt;
}while(RMS>=0.000001 || RMSv>=0.000001);

FILE *f1;
f1=fopen("LDC.dat","w");
fprintf(f1,"\nVARIABLES=""X"",""Y"",U-vel"",""V-vel"",""Pressure"" \nZONE I=%d,J=%d ZONETYPE=ORDERED,DATAPACKING=POINT\n",ny,nx);
for (i=1;i<=nx;i++)
{
	for (j=1;j<=ny;j++)
	{
		fprintf(f1,"%f\t%f\t%f\t%f\t%f\n",x[i],y[j],u[i][j],v[i][j],p[i][j]);
	}
	fprintf(f1,"\n");
}
fclose(f1);
}
