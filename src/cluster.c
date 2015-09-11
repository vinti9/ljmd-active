
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gr.h"
#include "params.h"



double max(double x, double y, double z)
{
    return ( (x>=y && x>=z) ? x : ( (y>=x && y>=z)? y : z));
}



//Evaluation of Radial distribution
//Mode 0: Initialize
//Mode 1: Calculate
//Mode 2: Save
void rdf(int mode)
{

    int i,j,ig;
    double xr,yr,zr,r,r2;
    double vb,nideal;
    FILE *fp;

    if (mode==0)
    {
        nbin_gr = 200;          //No. of bins used for RDF evaluation
        //delg = max(box.xhalf,box.yhalf,box.zhalf)/nbin_gr;
        delg = rcut*2/nbin_gr;
//        delg = 0.02;
        for (i=0; i<nbin_gr; i++)
        {
            g[i] = 0;
        }
        ngr = 0;
        fp = fopen("gr.dat","w");
        fclose(fp);
    }
    else if (mode==1)
    {
    
        //Calculat the gr
        ngr = ngr+1;
        for (i=0; i<Npart-1; i++)
        {
            for(j=i+1; j<Npart; j++)
            {
                
                xr = particle[i].x-particle[j].x;
                yr = particle[i].y-particle[j].y;
                zr = particle[i].z-particle[j].z;
                if(xr>box.xhalf) xr-=box.xlen; else if (xr<-box.xhalf) xr +=box.xlen;
                if(yr>box.yhalf) yr-=box.ylen; else if (yr<-box.yhalf) yr +=box.ylen;
                if(zr>box.zhalf) zr-=box.zlen; else if (zr<-box.zhalf) zr +=box.zlen;

                r2 = xr*xr+yr*yr+zr*zr;
                r = sqrt(r2);
                
                ig = (int)(r/delg);
                if(ig<MAXRDFBIN)
                {
                    g[ig] = g[ig] + 2;
                }
            }
        }
    }
    else if (mode==2)
    {
        fp = fopen("gr.dat","w");
        fprintf(fp,"#No. of evaluations %d \n",ngr);
        for (i=0; i<nbin_gr; i++)
        {
            r = delg*(i+0.5);
            vb = (cube(i+1)-cube(i))*cube(delg);
            nideal = (4.0/3.0)*3.14159*vb*density;
            //g[i] = g[i]/(ngr*Npart*nideal);
            fprintf(fp," %3d %.3f %.3f \n",i,r,g[i]/(ngr*Npart*nideal));
        }
        fprintf(fp,"\n");
        fclose(fp);
    }
    else
    {
        printf("Error: Check RDF mode value!");
    }


}


/* Calculate the center of mass and re-reference the system by
 * adjusting the particle co-ordinates
 */
void recenter_com()
{
    int i = 0;
    double com[3];
    com[0] = 0.0;
    com[1] = 0.0;
    com[2] = 0.0;
    
    for (i=0; i<Npart; i++)
    {
        com[0] += particle[i].x;
        com[1] += particle[i].y;
        com[2] += particle[i].z;
    }
    com[0] = com[0]/Npart;
    com[1] = com[1]/Npart;
    com[2] = com[2]/Npart;
    
    //if(step%1000==0)
    //{
    //    printf("comx=%.3f comy=%.3f comz=%.3f\n",com[0],com[1],com[2]);
    //}
    
    for (i=0; i<Npart; i++)
    {
        particle[i].x -= com[0];
        particle[i].y -= com[1];
        particle[i].z -= com[2];
    }
}


/* Calculate diffusion coefficients
 * and autocorrelation functions
 * from Frenkel-Smit Algorithm 8
 * 
 * Mode 0: Initialize
 * Mode 1: Calculate
 * Mode 2: Save
 * dtime = Time diff. between the measurements
 */
void diffusion(int mode, int dtime)
{
    static int count;
    int i,j;
    
    if(mode==0) //initialize
    {
        count=0;
        vacf[0]=0.0;
        vacf[1]=0.0;
        vacf[2]=0.0;
        r2t[0]=0.0;
        r2t[1]=0.0;
        r2t[2]=0.0;
        eacf[0]=0.0;
        eacf[1]=0.0;
        eacf[2]=0.0;
        
        xinit = (double *)malloc(Npart*sizeof(double));
        yinit = (double *)malloc(Npart*sizeof(double));
        zinit = (double *)malloc(Npart*sizeof(double));
        vx0 = (double *)malloc(Npart*sizeof(double));
        vy0 = (double *)malloc(Npart*sizeof(double));
        vz0 = (double *)malloc(Npart*sizeof(double));
        ex0 = (double *)malloc(Npart*sizeof(double));
        ey0 = (double *)malloc(Npart*sizeof(double));
        ez0 = (double *)malloc(Npart*sizeof(double));
        
        for(i=0; i<Npart; i++)
        {
            xinit[i]=particle[i].x;
            yinit[i]=particle[i].y;
            zinit[i]=particle[i].z;
            vx0[i]=particle[i].vx;
            vy0[i]=particle[i].vy;
            vz0[i]=particle[i].vz;
            
            ex0[i]=particle[i].ex;
            ey0[i]=particle[i].ey;
            ez0[i]=particle[i].ez;
        }
    }
    else if(mode==1)    //Sample
    {
        count++;

        for(i=0; i<Npart; i++)
        {
            vacf[0] += particle[i].vx*vx0[i];
            vacf[1] += particle[i].vy*vy0[i];
            vacf[2] += particle[i].vz*vz0[i];
            
            r2t[0] += sqr(particle[i].x - xinit[i]);
            r2t[1] += sqr(particle[i].y - yinit[i]);
            r2t[2] += sqr(particle[i].z - zinit[i]);
            
            eacf[0] += particle[i].ex*ex0[i];
            eacf[1] += particle[i].ey*ey0[i];
            eacf[2] += particle[i].ez*ez0[i];
        }
    }
    else if(mode==2)
    {
        fprintf(diff_file,"%8d %5d %+0.3e %+0.3e %+0.3e ",step, count,r2t[0]/(Npart*count*dtime),r2t[1]/(Npart*count*dtime),r2t[2]/(Npart*count*dtime));
        fprintf(diff_file,"%+0.3e %+0.3e %+0.3e ",vacf[0]/(Npart*count*dtime),vacf[1]/(Npart*count*dtime),vacf[2]/(Npart*count*dtime));
        fprintf(diff_file,"%+0.3e %+0.3e %+0.3e\n",eacf[0]/(Npart*count*dtime),eacf[1]/(Npart*count*dtime),eacf[2]/(Npart*count*dtime));
    }
    
}




