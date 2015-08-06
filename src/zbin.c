#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"
#include "zbin.h"


//Calculate the bin id acording to coordinates
int coords2zbin(double z)
{
    return ((int)((z+box.zhalf)*dzbini) % nzbin);
}

//Update particle zbin
void move2zbin(int i)
{
  int zbin_new,zbin_old;
  
  zbin_new = coords2zbin(particle[i].z);  
  zbin_old = particle[i].zbin;					
  if(zbin_new < 0 || zbin_new >= nzbin) 
  { 
      printf("Error: out of box?, i=%d, zbin_new=%d, zbin_old=%d\n",i,zbin_new,zbin_old);  
      printf("i.x=%.3f,i.y=%.3f,i.z=%.3f\n",particle[i].x,particle[i].y,particle[i].z);
      exit(1); 
  } 
  particle[i].zbin = zbin_new;
}

//Initialize zbin parameters
void init_zbins()
{
    int i,j,k,bin;
    double zgrid;
    
    zgrid = pow(density,-0.34);
    dzbin = zgrid/2.0;
    //dzbin = zgrid;
    nzbin = box.zlen/dzbin;
    dzbini = 1.0/dzbin;
    vzbin = box.xlen*box.ylen*dzbin;
    vzbini = 1.0/vzbin;
    
    for (i=0; i<Npart; i++)
    {
        bin = coords2zbin(particle[i].z);
        particle[i].zbin = bin;
    }
    
    for(j=0; j<nzbin; j++)
    {
        for(k=0; k<6; k++)
        {
            press_z[j][k] = 0.0;
            avg_press_z[j][k] = 0.0;
        }
        npart_z[j] = 0.0;
        rho_z[j] = 0.0;
    }
    printf("dzbin=%e, nzbin=%d\n",dzbin,nzbin);
}

//Measure density as per z bins
void measure_rho_z()
{
    int i;
    for(i=0; i<Npart; i++)
    {
        npart_z[particle[i].zbin]++;
    }
    
    for(i=0; i<nzbin; i++)
    {
        rho_z[i] = npart_z[i]*100/(vzbin*step);
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
            press_z[p][q] = Temp*rho_z[p] + press_z[p][q]*vzbini;
            avg_press_z[p][q] = avg_press_z[p][q] + press_z[p][q];
        }
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
        fprintf(fp,"%4d %+0.5e %+.5e\n",i,(dzbin*i-box.zhalf),rho_z[i]);
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
        fprintf(fp,"%4d %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n",bin,(dzbin*bin-box.zhalf),avg_press_z[bin][0]*100/step,avg_press_z[bin][1]*100/step,avg_press_z[bin][2]*100/step,avg_press_z[bin][3]*100/step,avg_press_z[bin][4]*100/step,avg_press_z[bin][5]*100/step);
    }
    fclose(fp);
}
