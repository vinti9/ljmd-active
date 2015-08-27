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
                slabs_z[j].avg_vdotr[k] = 0.0;
                slabs_z[j].avg_fdotv[k] = 0.0;
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
            slabs_z[p].press[q] = (slabs_z[p].rho_z/beta + slabs_z[p].press[q]*vzbini);
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
    double fdotv_z[nzbin][3];
    double v0_z[nzbin][3];
    double v02_z[nzbin][3];
    double vprop_z[nzbin][3];
    
    for(p=0; p<nzbin; p++)
    {
/*
        vdotr[p][0] = 0.0;
        vdotr[p][1] = 0.0;
        vdotr[p][2] = 0.0;
*/
///*
        fdotv_z[p][0] = 0.0;
        fdotv_z[p][1] = 0.0;
        fdotv_z[p][2] = 0.0;
        v0_z[p][0] = 0.0;
        v0_z[p][1] = 0.0;
        v0_z[p][2] = 0.0;
        v02_z[p][0] = 0.0;
        v02_z[p][1] = 0.0;
        v02_z[p][2] = 0.0;
        vprop_z[p][0] = 0.0;
        vprop_z[p][1] = 0.0;
        vprop_z[p][2] = 0.0;
//*/
    }
    
    for(i=0; i<Npart; i++)
    {
        bin = particle[i].zbin;
/*
        vdotr[bin][0] += particle[i].vx*particle[i].x;
        vdotr[bin][1] += particle[i].vy*particle[i].y;
        vdotr[bin][2] += particle[i].vz*particle[i].z;
*/
///*
        fdotv_z[bin][0] += particle[i].fx*particle[i].vx;
        fdotv_z[bin][1] += particle[i].fy*particle[i].vy;
        fdotv_z[bin][2] += particle[i].fz*particle[i].vz;
        
        //v0_z[bin][0] += particle[i].vx;
        //v0_z[bin][1] += particle[i].vy;
        //v0_z[bin][2] += particle[i].vz;
        v02_z[bin][0] += sqr(particle[i].vx);
        v02_z[bin][1] += sqr(particle[i].vy);
        v02_z[bin][2] += sqr(particle[i].vz);
//*/        
    }
    
    for(p=0; p<nzbin; p++)
    {
/*        
        //Update the stored cumulative sum of vdotr
        slabs_z[p].avg_vdotr[0] += vdotr[p][0];
        slabs_z[p].avg_vdotr[1] += vdotr[p][1];
        slabs_z[p].avg_vdotr[2] += vdotr[p][2];
        
        //Calculate inst. active pressure
        slabs_z[p].press_active[0] = ONEOVER3*friction_coeff*Fprop*slabs_z[p].avg_vdotr[0]*vzbini*100/step;
        slabs_z[p].press_active[1] = ONEOVER3*friction_coeff*Fprop*slabs_z[p].avg_vdotr[1]*vzbini*100/step;
        slabs_z[p].press_active[2] = ONEOVER3*friction_coeff*Fprop*slabs_z[p].avg_vdotr[2]*vzbini*100/step;
*/
///*
        //Update the stored cumulative sum of fdotv
        slabs_z[p].avg_fdotv[0] += fdotv_z[p][0];
        slabs_z[p].avg_fdotv[1] += fdotv_z[p][1];
        slabs_z[p].avg_fdotv[2] += fdotv_z[p][2];
        
        if(slabs_z[p].n > 0)
        {
            //vprop_z[p][0] = v0_z[p][0]*Fprop/(slabs_z[p].n*friction_coeff;
            //vprop_z[p][1] = v0_z[p][1]*Fprop/(slabs_z[p].n*friction_coeff;
            //vprop_z[p][2] = v0_z[p][2]*Fprop/(slabs_z[p].n*friction_coeff;
/*
            v0_z[p][0] = v0_z[p][0]*Fprop/(slabs_z[p].n*friction_coeff);
            v0_z[p][1] = v0_z[p][1]*Fprop/(slabs_z[p].n*friction_coeff);
            v0_z[p][2] = v0_z[p][2]*Fprop/(slabs_z[p].n*friction_coeff);
*/
///*
            v02_z[p][0] = v02_z[p][0]*sqr(Fprop/friction_coeff)/slabs_z[p].n;
            v02_z[p][1] = v02_z[p][1]*sqr(Fprop/friction_coeff)/slabs_z[p].n;
            v02_z[p][2] = v02_z[p][2]*sqr(Fprop/friction_coeff)/slabs_z[p].n;
/*
            v0_z[p][0] = sqrt(v02_z[p][0]);
            v0_z[p][1] = sqrt(v02_z[p][1]);
            v0_z[p][2] = sqrt(v02_z[p][2]);
*/        
        }

        //Calculate inst. active pressure
/*        slabs_z[p].press_active[0] = (ONEOVER6*friction_coeff*Fprop*Fprop*slabs_z[p].n + 0.5*Fprop*slabs_z[p].avg_fdotv[0]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[1] = (ONEOVER6*friction_coeff*Fprop*Fprop*slabs_z[p].n + 0.5*Fprop*slabs_z[p].avg_fdotv[1]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[2] = (ONEOVER6*friction_coeff*Fprop*Fprop*slabs_z[p].n + 0.5*Fprop*slabs_z[p].avg_fdotv[2]*100/step)*vzbini/rot_diff_coeff;
*/
/*    Using Ps = (0.5*gam*N*v0*v0+0.5*v0*avg_fdotv)/V*Dr
        slabs_z[p].press_active[0] = (0.5*friction_coeff*vprop_z[p][0]*vprop_z[p][0]*slabs_z[p].n + 0.5*vprop_z[p][0]*slabs_z[p].avg_fdotv[0]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[1] = (0.5*friction_coeff*vprop_z[p][1]*vprop_z[p][1]*slabs_z[p].n + 0.5*vprop_z[p][1]*slabs_z[p].avg_fdotv[1]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[2] = (0.5*friction_coeff*vprop_z[p][2]*vprop_z[p][2]*slabs_z[p].n + 0.5*vprop_z[p][2]*slabs_z[p].avg_fdotv[2]*100/step)*vzbini/rot_diff_coeff;
*/
/*  v8  Using Ps = 0.5*v0*(N*gam*v0 + avg_fdotv)/V*Dr
        slabs_z[p].press_active[0] = 0.5*v0_z[p][0]*(friction_coeff*v0_z[p][0]*slabs_z[p].n + slabs_z[p].avg_fdotv[0]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[1] = 0.5*v0_z[p][1]*(friction_coeff*v0_z[p][1]*slabs_z[p].n + slabs_z[p].avg_fdotv[1]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[2] = 0.5*v0_z[p][2]*(friction_coeff*v0_z[p][2]*slabs_z[p].n + slabs_z[p].avg_fdotv[2]*100/step)*vzbini/rot_diff_coeff;
*/
///* v9    Using Ps = 0.5*(N*gam*v0*v0 + v0*avg_fdotv)/V*Dr
        slabs_z[p].press_active[0] = 0.5*(friction_coeff*v02_z[p][0]*slabs_z[p].n + Fprop*slabs_z[p].avg_fdotv[0]*100/step/friction_coeff)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[1] = 0.5*(friction_coeff*v02_z[p][1]*slabs_z[p].n + Fprop*slabs_z[p].avg_fdotv[1]*100/step/friction_coeff)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[2] = 0.5*(friction_coeff*v02_z[p][2]*slabs_z[p].n + Fprop*slabs_z[p].avg_fdotv[2]*100/step/friction_coeff)*vzbini/rot_diff_coeff;
//*/
/* v10    Using Ps = 0.5*(N*gam*v0*v0 + v0*avg_fdotv)/V*Dr
        slabs_z[p].press_active[0] = 0.5*(friction_coeff*v02_z[p][0]*slabs_z[p].n + v0_z[p][0]*slabs_z[p].avg_fdotv[0]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[1] = 0.5*(friction_coeff*v02_z[p][1]*slabs_z[p].n + v0_z[p][1]*slabs_z[p].avg_fdotv[1]*100/step)*vzbini/rot_diff_coeff;
        slabs_z[p].press_active[2] = 0.5*(friction_coeff*v02_z[p][2]*slabs_z[p].n + v0_z[p][2]*slabs_z[p].avg_fdotv[2]*100/step)*vzbini/rot_diff_coeff;
*/
        //Update the stored cumulative sum of pressure
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
//        fprintf(fp," %+.3e",slabs_z[bin].avg_rho_z[0]*100/step);
        fprintf(fp,"\t%+.3e",slabs_z[bin].avg_rho_z*100/step);
        fprintf(fp,"\t%+.3e %+.3e %+.3e",slabs_z[bin].avg_press[0]*100/step,slabs_z[bin].avg_press[1]*100/step,slabs_z[bin].avg_press[2]*100/step);
        fprintf(fp,"\t%+.3e %+.3e %+.3e\n",slabs_z[bin].avg_press_active[0]*100/step,slabs_z[bin].avg_press_active[1]*100/step,slabs_z[bin].avg_press_active[2]*100/step);
//        fprintf(fp,"\t%+.3e %+.3e %+.3e\n",slabs_z[bin].avg_fdotv[0]*100/step,slabs_z[bin].avg_fdotv[1]*100/step,slabs_z[bin].avg_fdotv[2]*100/step);
    }
    fclose(fp);
}
