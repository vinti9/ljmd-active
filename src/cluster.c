
#include <stdio.h>
#include <math.h>
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
        nbin_gr = 200;          //100 bins used for RDF evaluation
        delg = max(box.xhalf,box.yhalf,box.zhalf)/nbin_gr;
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
                
                if (r < box.xhalf)
                {
                    ig = (int)(r/delg);
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








