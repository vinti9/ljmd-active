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
//    printf("%d,%f,%f\n",Npart,pow(Npart,0.333),box/pow(Npart,0.333));
    int p = (int)pow(Npart,0.34);
    double xgrid = pow(density,-0.34);
    double ygrid = xgrid;
    double zgrid = xgrid;
    printf("%d, xgrid=%f, ygrid=%f, zgrid=%f\n",p,xgrid,ygrid,zgrid);
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt,12345);

    if (lattice==0)     //Random positions
    {
        for (i=0;i<Npart;i++)
        {
    
            particle[i].x = x0 + dsfmt_genrand_open_open(&dsfmt)*xlen;
            particle[i].y = y0 + dsfmt_genrand_open_open(&dsfmt)*ylen;
            particle[i].z = z0 + dsfmt_genrand_open_open(&dsfmt)*zlen;
        }
    }
    else if(lattice ==1)        //Regular simple cubic lattice as per slab requirement
    {
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
    //else if(lattice == 2)      //Regular FCC lattice as per slab requirement
    //{
    //       i = 0;
    //       for (iz=0;(iz*zgrid)<zlen;iz++)
    //       {
    //           for(ix=0;(ix*xgrid)<xlen;ix++)
    //           {
    //               for(iy=0;(iy*ygrid<ylen);iy++)
    //               { 
    //                   //i = iz*p*p + ix*p + iy;
    //                   if (i<Npart)
    //                   {
    //                       particle[i].x = x0 + (ix+0.5)*xgrid;
    //                       particle[i].y = y0 + (iy+0.5)*ygrid;
    //                       particle[i].z = z0 + (iz+0.5)*zgrid;
    //                       i = i+1;
    //                   }
    //               }
    //           }
    //       }
    //}
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
        particle[i].xold = particle[i].x;
        particle[i].yold = particle[i].y;
        particle[i].zold = particle[i].z;
    }
    
    i=0;
    FILE *fplattice;
    fplattice = fopen("init_pos.dat","w");
    fprintf(fplattice,"#Npart= %d\n",Npart);
    fprintf(fplattice,"#ID,\tx,\ty,\tz,\tvx,\tvy,\tvz,\tfx,\tfy,\tfz\n");
    for (i=0;i<Npart;i++)
    {
        fprintf(fplattice,"%5d \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e\n",i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].fx,particle[i].fy,particle[i].fz);
    }
    fclose(fplattice);

}