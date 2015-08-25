#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"
#include "zbin.h"


//Calculate the bin id acording to coordinates
int coords2zbin(double z)
{
    //return ((int)((z+box.zhalf)*dzbini) % nzbin);
    return ((int)(nzbin*(z+box.zhalf)*box.zleni) % nzbin);
}

//Update particle zbin
void move2zbin(int i)
{
    int zbin_new,zbin_old,count,j;
  
    zbin_new = coords2zbin(particle[i].z);  
    zbin_old = particle[i].zbin;
    if(zbin_new < 0 || zbin_new >= nzbin) 
    { 
        printf("Error: out of zbin?, i=%d, zbin_new=%d, zbin_old=%d\n",i,zbin_new,zbin_old);  
        printf("i.x=%.3f,i.y=%.3f,i.z=%.3f\n",particle[i].x,particle[i].y,particle[i].z);
        exit(1); 
    }
    if (zbin_new != zbin_old)
    {
        particle[i].zbin = zbin_new;          //Update the slab no. to particle
    }
}

//Initialize zbin parameters
void init_zbins()
{
    int i,j,k,bin;
    double zgrid;
    
    zgrid = pow(density,-0.34);
    //dzbin = zgrid/4.0;
    dzbin = zgrid/2.0;
    //dzbin = zgrid;
    nzbin = box.zlen/dzbin;
    dzbini = 1.0/dzbin;
    vzbin = box.xlen*box.ylen*dzbin;
    vzbini = 1.0/vzbin;
    
    for(j=0; j<nzbin; j++)
    {
        for(k=0; k<6; k++)
        {
            slabs_z[j].press[k] = 0.0;
            slabs_z[j].avg_press[k] = 0.0;
            slabs_z[j].press_active[k] = 0.0;
            slabs_z[j].avg_press_active[k] = 0.0;
        }
        slabs_z[j].n = 0;
        slabs_z[j].rho_z = 0.0;
        slabs_z[j].avg_rho_z = 0.0;
    }
    
    for (i=0; i<Npart; i++)
    {
        bin = coords2zbin(particle[i].z);
        particle[i].zbin = bin;
        //slabs_z[bin].part[slabs_z[bin].n] = i;      //Add the ID of particle to the slab
        slabs_z[bin].n++;                      //Increase the count of particle to the slab
    }
    printf("dzbin=%e, nzbin=%d\n",dzbin,nzbin);
}

//Measure density as per z bins
void measure_rho_z()
{
    int i;
    for(i=0; i<nzbin; i++)
    {
        slabs_z[i].rho_z = 0;       //Reinitialize to zero for next measurement
        slabs_z[i].n = 0;
    }
    
    for(i=0; i<Npart; i++)
    {
        slabs_z[particle[i].zbin].n++;
    }
    
    for(i=0; i<nzbin; i++)
    {
        slabs_z[i].rho_z = slabs_z[i].n*vzbini;
        slabs_z[i].avg_rho_z = slabs_z[i].avg_rho_z + slabs_z[i].rho_z;
    }
}

//Measure pressure as per z bins
void measure_press_z()
{
    int p,q;
    for(p=0; p<nzbin; p++)
    {
        for(q=0; q<6; q++)
        {
            //press[p][q] = Temp*rho_z[p] + press[p][q]*vzbini;
            slabs_z[p].press[q] = slabs_z[p].rho_z/beta + slabs_z[p].press[q]*vzbini;
            slabs_z[p].avg_press[q] = slabs_z[p].avg_press[q] + slabs_z[p].press[q];
        }
    }
    measure_active_press_z();
}

//Measure active pressure as per z bins
void measure_active_press_z()
{
    int p,bin,i;
    double vdotr[nzbin][3];
    
    for(p=0; p<nzbin; p++)
    {
        vdotr[p][0] = 0.0;
        vdotr[p][1] = 0.0;
        vdotr[p][2] = 0.0;
    }
    
    for(i=0; i<Npart; i++)
    {
        bin = particle[i].zbin;
        vdotr[bin][0] += particle[i].vx*particle[i].x;
        vdotr[bin][1] += particle[i].vy*particle[i].y;
        vdotr[bin][2] += particle[i].vz*particle[i].z;
    }
    
    for(p=0; p<nzbin; p++)
    {
        //Calculate inst. active pressure
        slabs_z[p].press_active[0] = friction_coeff*Fprop*vdotr[p][0]*vzbini;
        slabs_z[p].press_active[1] = friction_coeff*Fprop*vdotr[p][1]*vzbini;
        slabs_z[p].press_active[2] = friction_coeff*Fprop*vdotr[p][2]*vzbini;
        
        //Update the stored cumulative sum 
        slabs_z[p].avg_press_active[0] += slabs_z[p].press_active[0];
        slabs_z[p].avg_press_active[1] += slabs_z[p].press_active[1];
        slabs_z[p].avg_press_active[2] += slabs_z[p].press_active[2];
    }
}

//Write density profile as per z bins
void write_rho_z()
{
    int i;
    FILE *fp;
    fp = fopen("rhoz.dat","w");
    for(i=0; i<nzbin; i++)
    {
        //fprintf(fp,"%4d %+0.5e %+.5e\n",i,(dzbin*i-box.zhalf),rho_z[i]);
        fprintf(fp,"%4d %+0.5e %+.5e\n",i,(dzbin*i-box.zhalf),slabs_z[i].avg_rho_z*100/step);
    }
    fclose(fp);
}

//Write Pressure profile along z bins
//The pressure is cumulative average till the time step
void write_press_z()
{
    int bin;
    FILE *fp;
    fp = fopen("press_z.dat","w");
    for(bin=0; bin<nzbin; bin++)
    {
        fprintf(fp,"%4d %+.3e",bin,(dzbin*bin-box.zhalf));
        fprintf(fp," %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e",slabs_z[bin].avg_press[0]*100/step,slabs_z[bin].avg_press[1]*100/step,slabs_z[bin].avg_press[2]*100/step,slabs_z[bin].avg_press[3]*100/step,slabs_z[bin].avg_press[4]*100/step,slabs_z[bin].avg_press[5]*100/step);
        fprintf(fp,"\t%+.3e %+.3e %+.3e\n",slabs_z[bin].avg_press_active[0]*100/step,slabs_z[bin].avg_press_active[1]*100/step,slabs_z[bin].avg_press_active[2]*100/step);
    }
    fclose(fp);
}
