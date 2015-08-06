#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dSFMT.h"      //for Random Number generator
#include "params.h"



/* Function to generate a specified form of slab given the start coordinates
 * and the slab dimensions
 */
void init_lattice_slab(int lattice, double x0, double y0, double z0, double xlen, double ylen, double zlen)
{
    int i=0;
    int ix=0;
    int iy=0;
    int iz=0;
    double xgrid,ygrid,zgrid;
    
//    printf("%d,%f,%f\n",Npart,pow(Npart,0.333),box/pow(Npart,0.333));
    int p = (int)pow(Npart,0.34);
    dsfmt_t dsfmt_obj;
    dsfmt_init_gen_rand(&dsfmt_obj,12345);

    if (lattice==0)     //Random positions
    {
        for (i=0;i<Npart;i++)
        {
    
            particle[i].x = x0 + dsfmt_genrand_open_open(&dsfmt_obj)*xlen;
            particle[i].y = y0 + dsfmt_genrand_open_open(&dsfmt_obj)*ylen;
            particle[i].z = z0 + dsfmt_genrand_open_open(&dsfmt_obj)*zlen;
        }
    }
    else if(lattice ==1)        //Regular simple cubic lattice as per slab requirement
    {
        xgrid = pow(density,-0.34);
        ygrid = xgrid;
        zgrid = xgrid;
        printf("SC lattice\n");
        printf("%d, xgrid=%f, ygrid=%f, zgrid=%f\n",p,xgrid,ygrid,zgrid);
        i = 0;
        for (iz=0;(iz*zgrid)<zlen;iz++)
        {
           for(ix=0;(ix*xgrid)<xlen;ix++)
           {
               for(iy=0;(iy*ygrid<ylen);iy++)
               { 
                   //i = iz*p*p + ix*p + iy;
                   if (i<Npart)
                   {
                       particle[i].x = x0 + (ix+0.5)*xgrid;
                       particle[i].y = y0 + (iy+0.5)*ygrid;
                       particle[i].z = z0 + (iz+0.5)*zgrid;
                       i = i+1;
                   }
               }
           }
        }
    }
    else if(lattice == 2)      //Regular FCC lattice as per slab requirement
    {
        double atpos[4][3];
        int l,nx,ny,nz;
        
        xgrid = pow(4.0/density,1.0/3.0);
        ygrid = xgrid;
        zgrid = xgrid;
        nx = round(xlen/xgrid);
        ny = round(ylen/ygrid);
        nz = round(zlen/zgrid);
        printf("FCC lattice\n");
        printf("%d, xgrid=%f, ygrid=%f, zgrid=%f\n",p,xgrid,ygrid,zgrid);
        printf("nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
        
       
        atpos[0][0] = 0.;
        atpos[0][1] = 0.;
        atpos[0][2] = 0.;
    
        atpos[1][0] = 0.5;
        atpos[1][1] = 0.5;
        atpos[1][2] = 0.;
    
        atpos[2][0] = 0.;
        atpos[2][1] = 0.5;
        atpos[2][2] = 0.5;
    
        atpos[3][0] = 0.5;
        atpos[3][1] = 0.;
        atpos[3][2] = 0.5;
       
        i = 0;
        for (iz=0;(iz*zgrid)<zlen;iz++)
        {
            for(ix=0;(ix*xgrid)<xlen;ix++)
            {
                for(iy=0;(iy*ygrid<ylen);iy++)
                { 
                    for (l = 0; l<4; l++)
                    {
                         //i = iz*p*p + ix*p + iy;
                         if (i<Npart)
                         {
                             particle[i].x = x0 + (atpos[l][0] + ix + 0.25)*xgrid;
                             particle[i].y = y0 + (atpos[l][1] + iy + 0.25)*ygrid;
                             particle[i].z = z0 + (atpos[l][2] + iz + 0.25)*zgrid;
                             i = i+1;
                         }     
                    }
                }
            }
        }
    }
    else
    {
        printf("Check Lattice parameter in prm.dat\n");
        exit(1);
    }
    
    //Check if all particles allocated properly in the slab
    if (i != Npart)
    {
        printf("Cannot allocate Npart into the given slab!\n");
        printf("Check the given parameters\n");
        printf("Npart=%d, i=%d, ix=%d, iy=%d, iz=%d\n",Npart,i,ix,iy,iz);
        exit(1);
    }

    for(i=0;i<Npart;i++)
    {
        //Initialize variables with dummy values
        particle[i].xold = particle[i].x - 0.000001;
        particle[i].yold = particle[i].y - 0.000001;
        particle[i].zold = particle[i].z - 0.000001;
    }
    
    i=0;
    FILE *fplattice;
    fplattice = fopen("init_pos.dat","w");
    //fprintf(fplattice,"#Npart= %d\n",Npart);
    //fprintf(fplattice,"#ID,\tx,\ty,\tz,\tvx,\tvy,\tvz,\tfx,\tfy,\tfz\n");
    fprintf(fplattice,"%i\n",Npart);
    fprintf(fplattice,"%le %le\n",-box.xhalf,box.xhalf);
    fprintf(fplattice,"%le %le\n",-box.yhalf,box.yhalf);
    fprintf(fplattice,"%le %le\n",-box.zhalf,box.zhalf);
    for (i=0;i<Npart;i++)
    {
        fprintf(fplattice,"%5d \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e\n",i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].fx,particle[i].fy,particle[i].fz);
    }
    fclose(fplattice);

}