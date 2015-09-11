/*
 * This file contains the 'main subroutine' and function calls
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "dSFMT.h"      //for Random Number generator
#include "params.h"
#include "cell.h"
#include "gr.h"
#include "lattice.h"
#include "zbin.h"

//Global function Declarations
void    initialize();
void    init_lattice_pos(int);
void    init_lattice_vel();
void    force();
void    force_nlist();
void    partforce(int,int);
void    integrate_euler();
void    sample();
void    pbc(int);
void    write_pos(int);
void    write_traj();
void    init_box(double);
void    init_box_cuboid(double, double, double);
double  gen_gaussian_rand();
void    write_pressure();
void    read_conf();
void    eval_vprop();
        

//Global Structure declarations
dsfmt_t     dsfmt;
//cells_t     *cells;


//Global File pointers
FILE *prm_file;
FILE *prop_file;
FILE *pos_file;
FILE *traj_file;
FILE *press_file;
FILE *out_file;
FILE *conf_file;

//Initialization of all global parameters
void initialize()
{
    int dummy;
    int p;
    int seed;
    char str[30];
    double Lx,Ly,Lz; 
    
    prm_file = fopen("prm.dat","r");
    dummy = fscanf(prm_file,"%s %d",str,&seed);
    dummy = fscanf(prm_file,"%s %d",str,&nsteps_equil);
    dummy = fscanf(prm_file,"%s %d",str,&nsteps);
    dummy = fscanf(prm_file,"%s %lf",str,&tstep);
    dummy = fscanf(prm_file,"%s %d",str,&Npart);
    dummy = fscanf(prm_file,"%s %lf",str,&rcut);
    dummy = fscanf(prm_file,"%s %lf",str,&initTemp);
    dummy = fscanf(prm_file,"%s %lf",str,&density);
    dummy = fscanf(prm_file,"%s %lf",str,&beta);
    dummy = fscanf(prm_file,"%s %lf",str,&friction_coeff);
    dummy = fscanf(prm_file,"%s %lf",str,&epsilon);
    dummy = fscanf(prm_file,"%s %lf",str,&rot_diff_coeff);
    dummy = fscanf(prm_file,"%s %lf",str,&Fprop);
    dummy = fscanf(prm_file,"%s %d",str,&continu);
    dummy = fscanf(prm_file,"%s %d",str,&lattice);
    dummy = fscanf(prm_file,"%s %lf %lf %lf",str,&Lx,&Ly,&Lz);
    dummy = fscanf(prm_file,"%s %d",str,&print);
    fclose(prm_file);

    dsfmt_init_gen_rand(&dsfmt,seed);

    //time = 0;
    step = 0;
    tstep2 = tstep*tstep;
    tstepi = 1/tstep;
    rcut2 = rcut*rcut;
    Temp = 0.0;
    
    etot = 0.0;
    en = 0.0;
    if ( (Lx>0.0) && (Ly>0.0) && (Lz>0.0))
    {
        init_box_cuboid(Lx,Ly,Lz);
        printf("Cuboidal simulation box ! \n");
    }
    else
    {
        init_box(density);    
    }
    
//Variables for brownian dynamics
    C1 = 1.0/friction_coeff;
    sqrt_tstep = sqrt(tstep);
    C2 = pow(2.0/(beta*friction_coeff),0.5);
    C3 = pow(2.0*rot_diff_coeff,0.5);

    
    if(Npart>MAXNPART)
    {
        printf("Npart > MAXNpart, Check definition in main.c");
        exit(1);    //Exit with the error
    }
    else
    {
        //Initialize position in lattice
        particle_t particle[Npart];
    }

    if(continu == 0)
    {
//        init_lattice_pos(lattice);
        double zslab = (Npart/density)/(box.xlen*box.ylen);
        //double zslab = box.zlen;                                    //To make cubic box
        init_lattice_slab(lattice, -box.xhalf, -box.yhalf, -zslab/2.0, box.xlen, box.ylen, zslab);
        printf("Initalized positions\n");
        init_lattice_vel();
        printf("Initialized velocities\n");
    }
    else if(continu == 1)
    {
        read_conf();
        printf("Used positions from previous config\n");
        printf("Used velocities from previous config \n");
    }
//    printf("vx=%.3e,\tvy=%.3e,\tvz=%.3e\n",sumv[0],sumv[1],sumv[2]);
//    printf("vx2=%.3e,\tvy2=%.3e,\tvz2=%.3e\n",sumv2[0],sumv2[1],sumv2[2]);

#ifdef USE_NLIST
    init_cells();           //Initialize cell list
#endif

    init_zbins();           //Initialize zbin profiles

    
    ecut = 4*beta*epsilon*(1.0/(pow(rcut,12))-1.0/(pow(rcut,6)));        //calculate ecut
    //ecut = 4*epsilon*(1.0/(pow(rcut,12))-1.0/(pow(rcut,6)));        //calculate ecut

    printf("ecut=%e\n",ecut);

    //Initialize stress tensor
    for (p =0; p<6; p++)
    {
        avg_press[p] = 0.0;
        avg_fdotv[p] = 0.0;
    }
    
    //initialize running averages to zero
    for(p=0;p<6;p++)
    {
        Av[p]=0;
    }
        
    //Initialize RDF parameters
    rdf(0);
    diffusion(0,0);
  
}


//Standard Force calculation
void force()
{
    int i,j;

    for(i=0;i<Npart;i++) //Initialize forces to zero
    {
        particle[i].fx=0;
        particle[i].fy=0;
        particle[i].fz=0;
    }

    for(i=0;i<Npart;i++)
    {
        for (j=i+1;j<Npart;j++)
        {
            partforce(i,j);
        }

    }

}




//Force calculation using neighbor list
void force_nlist()
{
    int i = 0;
    int j = 0;
    int p = 0;
    int q = 0;
    int c = 0;
    int icell,inbr;
    int s;

    for(i=0;i<Npart;i++) //Initialize forces to zero
    {
        particle[i].fx=0;
        particle[i].fy=0;
        particle[i].fz=0;
    }
    
    for (s=0; s<6; s++)  //Initialize pressure tensor to zero
    {
        pressure[s]=0.0;
    }
    
    for(p=0; p<nzbin; p++)
    {
        for(q=0; q<6; q++)
        {
            slabs_z[p].press[q] = 0.0;
//            avg_press_z[p][q] = 0.0;
        }
    }    

    en = 0.0;           //Initialize potential energy to zero for the MD step

    for(i=0;i<Npart;i++)
    {
        icell = particle[i].cell;
        for (p=0; p<27; p++)
        {
            inbr = cells[icell].neighbors[p];
            for (q=0; q<cells[inbr].n; q++)
            {
                j = cells[inbr].particles[q];
                if (j>Npart)
                {
                    printf("Error!: force_nlist, step=%d, i=%d, icell=%d, p=%d, inbr=%d, q=%d, j=%d\n",step,i,icell,p,inbr,q,j);
                    exit(1);
                }
                else if (j>i)
                {
                    partforce(i,j);
                }
            }
        }

    }


}


//Function to calculate force between two particles with ID i & j
void partforce(int i, int j)
{
    double xr = 0;
    double yr = 0;
    double zr = 0;
    double r = 0;
    double r2 = 0;
    double r2i = 0;
    double r6i = 0.0;
    double ff = 0.0;
    double virial = 0.0;
    int    ibin,jbin,pz,nz,dz, low_bin, high_bin;
    double nzi,low_z, high_z;
    int    boxcross = 0;


    xr = particle[i].x-particle[j].x;
    yr = particle[i].y-particle[j].y;
    zr = particle[i].z-particle[j].z;
//  xr = xr - box*round(xr/box);
//  yr = yr - box*round(yr/box);
//  zr = zr - box*round(zr/box);
    if(xr>box.xhalf) xr-=box.xlen; else if (xr<-box.xhalf) xr +=box.xlen;
    if(yr>box.yhalf) yr-=box.ylen; else if (yr<-box.yhalf) yr +=box.ylen;
    if(zr>box.zhalf) zr-=box.zlen; else if (zr<-box.zhalf) zr +=box.zlen;
    r2 = xr*xr+yr*yr+zr*zr;
    if (r2<rcut2)
    {
        r2i = 1/r2;      
        r6i = r2i*r2i*r2i;
        ff = epsilon*48.0*r2i*r6i*(r6i-0.5);        // = (1/r)*(dU/dr)

        //Update the forces
        particle[i].fx=particle[i].fx+ff*xr;
        particle[i].fy=particle[i].fy+ff*yr;
        particle[i].fz=particle[i].fz+ff*zr;
        particle[j].fx=particle[j].fx-ff*xr;
        particle[j].fy=particle[j].fy-ff*yr;
        particle[j].fz=particle[j].fz-ff*zr;

#ifdef MEASURE_PRESSURE 
        if (mode == 1)
        {
            //Update the stress tensor components
            pressure[0] = pressure[0] + ff*xr*xr;   //Pxx
            pressure[1] = pressure[1] + ff*yr*yr;   //Pyy
            pressure[2] = pressure[2] + ff*zr*zr;   //Pzz
            pressure[3] = pressure[3] + ff*xr*yr;   //Pxy
            pressure[4] = pressure[4] + ff*yr*zr;   //Pyz
            pressure[5] = pressure[5] + ff*zr*xr;   //Pzx
              
            if(step%100 == 0)
            {
                //Pressure profile along z direction
                ibin = particle[i].zbin;
                jbin = particle[j].zbin;
                if(ibin == jbin)                //particle pair in the same slab
                {
                    slabs_z[ibin].press[0] = slabs_z[ibin].press[0] + ff*xr*xr;   //Pxx
                    slabs_z[ibin].press[1] = slabs_z[ibin].press[1] + ff*yr*yr;   //Pyy
                    slabs_z[ibin].press[2] = slabs_z[ibin].press[2] + ff*zr*zr;   //Pzz
                    //slabs_z[ibin].press[3] = slabs_z[ibin].press[3] + ff*xr*yr;   //Pxy
                    //slabs_z[ibin].press[4] = slabs_z[ibin].press[4] + ff*yr*zr;   //Pyz
                    //slabs_z[ibin].press[5] = slabs_z[ibin].press[5] + ff*zr*xr;   //Pzx
                }
                else                            //add (1/n)th contribution to all the slabs in between the particle pair
                {
                    nz = (ibin>jbin) ? (ibin-jbin+1) : (jbin-ibin+1);
                    //nz = (ibin>jbin) ? (ibin-jbin-1) : (jbin-ibin-1);
                    if(nz > nzbin/2) { boxcross = 1;}
                    //if(ibin >= nzbin-1 || ibin <1 || jbin >= nzbin-1 || jbin<1 ) goto skip;
                    //Below if-else to take care for the nearest images of particles
                    if (boxcross == 0)
                    {
                        low_bin = (ibin < jbin) ? ibin : jbin ;
                        high_bin= (ibin < jbin) ? jbin : ibin ;
                        nz      = high_bin - low_bin + 1;
                        low_z   = (ibin < jbin) ? particle[i].z : particle[j].z;
                        high_z  = (ibin < jbin) ? particle[j].z : particle[i].z;	
                    }	
                    else
                    {
                        low_bin = (ibin > jbin) ? ibin : jbin ;
                        high_bin= (ibin > jbin) ? jbin : ibin ;
                        nz      = nzbin + (high_bin - low_bin) + 1;
                        low_z   = (ibin > jbin) ? particle[i].z : particle[j].z;
                        high_z  = (ibin > jbin) ? particle[j].z : particle[i].z;		
                    }
                    //nzi = 1.0/nz;
                    for (dz=0; dz<nz ; dz++)
                    {
                        //pz = low_bin + dz;
                        //pz = (low_bin + dz)%nzbin;               //Periodic images of z-bin
                        pz = (low_bin+dz)>=nzbin ? (low_bin+dz)-nzbin : (low_bin+dz);               //Periodic images of z-bin
                        if (pz == low_bin)
                        {
                            nzi = (((pz+1)*dzbin - box.zhalf) - low_z)/zr;
                            nzi = (nzi>0) ? (nzi) : (-nzi);     //Take only positive value
                        }
                        else if (pz == high_bin)
                        {
                            nzi = (high_z - (pz*dzbin - box.zhalf))/zr;
                            nzi = (nzi>0) ? (nzi) : (-nzi);     //Take only positive value
                        }
                        else
                        {
                            nzi = dzbin/zr;
                            nzi = (nzi>0) ? (nzi) : (-nzi);     //Take only positive value
                        }

                        slabs_z[pz].press[0] = slabs_z[pz].press[0] + ff*xr*xr*nzi;   //Pxx
                        slabs_z[pz].press[1] = slabs_z[pz].press[1] + ff*yr*yr*nzi;   //Pyy
                        slabs_z[pz].press[2] = slabs_z[pz].press[2] + ff*zr*zr*nzi;   //Pzz
                        //slabs_z[pz].press[3] = slabs_z[pz].press[3] + ff*xr*yr*nzi;   //Pxy
                        //slabs_z[pz].press[4] = slabs_z[pz].press[4] + ff*yr*zr*nzi;   //Pyz
                        //slabs_z[pz].press[5] = slabs_z[pz].press[5] + ff*zr*xr*nzi;   //Pzx
                    }
                }
            }
        }
      
#endif 
        //Update the system potential energy
        en = en + 4*epsilon*beta*r6i*(r6i-1) - ecut;
        //en = en + 4*epsilon*r6i*(r6i-1) - ecut;
//Intentional mistake below to match system values with Vassilis' code
//        en = en + 4*epsilon*beta*r6i*(r6i-1) + ecut;

    }

#ifdef MEASURE_RDF
    if (mode == 1)
    {
        if (step%1000 == 0)     //Every 1000 steps, measure rdf
        {
            //r = sqrt(r2);
            //if (r<box.xhalf)
            //{
            //    g[(int)(r/delg)] += 2;
            //}
            //ngr++;
            /*rdf(1);*/
        }
    }
#endif

}


/*Function to generate simple cubic or a Random configuration
 **/
void init_lattice_pos(int lattice)
{
    int i=0;
    int ix=0;
    int iy=0;
    int iz=0;
//    printf("%d,%f,%f\n",Npart,pow(Npart,0.333),box/pow(Npart,0.333));
    int p = (int)pow(Npart,0.34);
    double xgrid = (box.xlen/pow(Npart,0.34));
    double ygrid = (box.ylen/pow(Npart,0.34));
    double zgrid = (box.zlen/pow(Npart,0.34));
    printf("%d, %f, %f, %f\n",p,xgrid,ygrid,zgrid);

    if (lattice==0)     //Random positions
    {
        for (i=0;i<Npart;i++)
        {
    
            particle[i].x = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.xhalf;
            particle[i].y = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.yhalf;
            particle[i].z = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.zhalf;
        }
    }
    else if(lattice ==1)        //Regular simple cubic lattice
    {
           for (ix=0;(ix*xgrid)<box.xlen;ix++)
           {
               for(iy=0;(iy*ygrid)<box.ylen;iy++)
               {
                   for(iz=0;(iz*zgrid<box.zlen);iz++)
                   { 
                       i = ix*p*p + iy*p + iz;
                       if (i<Npart)
                       {
                           particle[i].x = (ix+0.5)*xgrid-box.xhalf;
                           particle[i].y = (iy+0.5)*ygrid-box.yhalf;
                           particle[i].z = (iz+0.5)*zgrid-box.zhalf;
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

        for(i=0;i<Npart;i++)
        {
            //Initialize variables with dummy values
            particle[i].xold = particle[i].x;
            particle[i].yold = particle[i].y;
            particle[i].zold = particle[i].z;
        }
//  }
    
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
//    fprintf(pos_file,"\n\nCOMvx=%e,\tCOMvy=%e,\tCOMvz=%e\n",sumv[0],sumv[1],sumv[2]);


}


void init_lattice_vel()
{
    int i;
    sumv[0]=0;sumv[1]=0;sumv[2]=0;
    sumv2[0]=0;sumv2[1]=0;sumv2[2]=0;
    for (i=0;i<Npart;i++)
    {
        particle[i].vx = (dsfmt_genrand_open_open(&dsfmt)-0.5);
        particle[i].vy = (dsfmt_genrand_open_open(&dsfmt)-0.5);
        particle[i].vz = (dsfmt_genrand_open_open(&dsfmt)-0.5);
        
        particle[i].ex = 1.0;
        particle[i].ey = 0.0;
        particle[i].ez = 0.0;
    }

    for (i=0;i<Npart;i++)
    {
        sumv[0] = sumv[0] + particle[i].vx;
        sumv[1] = sumv[1] + particle[i].vy;
        sumv[2] = sumv[2] + particle[i].vz;

        sumv2[0] = sumv2[0] + particle[i].vx*particle[i].vx;
        sumv2[1] = sumv2[1] + particle[i].vy*particle[i].vy;
        sumv2[2] = sumv2[2] + particle[i].vz*particle[i].vz;
    }
//    sumv[0] = sumv[0]/Npart;
//    sumv[1] = sumv[1]/Npart;
//    sumv[2] = sumv[2]/Npart;
    Temp = (sumv2[0]+sumv2[1]+sumv2[2])/(3*Npart);

    //Scale the velocities to a temp of initTemp
    for (i=0;i<Npart;i++)
    {
        particle[i].vx = particle[i].vx*(initTemp/pow(Temp,0.5)) - sumv[0]/Npart;
        particle[i].vy = particle[i].vy*(initTemp/pow(Temp,0.5)) - sumv[1]/Npart;
        particle[i].vz = particle[i].vz*(initTemp/pow(Temp,0.5)) - sumv[2]/Npart;
    }


}


//Initialize box dimensions
void init_box(double rho)
{
    double l = pow(Npart/rho,1.0/3.0);
    box.xlen = l;
    box.ylen = l;
    box.zlen = l;
    box.xhalf = box.xlen*0.5;
    box.yhalf = box.ylen*0.5;
    box.zhalf = box.zlen*0.5;
    box.xleni = 1.0/box.xlen;
    box.yleni = 1.0/box.ylen;
    box.zleni = 1.0/box.zlen;
    volume = box.xlen*box.ylen*box.zlen;
    volumei = 1.0/volume;

    printf("Simulation Box: Xlen=%.3e, Ylen=%.3e, Zlen=%.3e\n",box.xlen,box.ylen,box.zlen);
}

//Initialize box dimensions (each length wise)
void init_box_cuboid(double lx, double ly, double lz)
{
    box.xlen = lx;
    box.ylen = ly;
    box.zlen = lz;
    box.xhalf = box.xlen*0.5;
    box.yhalf = box.ylen*0.5;
    box.zhalf = box.zlen*0.5;
    box.xleni = 1.0/box.xlen;
    box.yleni = 1.0/box.ylen;
    box.zleni = 1.0/box.zlen;
    volume = box.xlen*box.ylen*box.zlen;
    volumei = 1.0/volume;

    printf("Simulation Box: Xlen=%.3e, Ylen=%.3e, Zlen=%.3e\n",box.xlen,box.ylen,box.zlen);
}

//Integrate the equations using the Euler-Maruyama scheme
void integrate_euler()
{
    int i,s;
    double xnew=0; 
    double ynew=0; 
    double znew=0;
    double vxnew=0; 
    double vynew=0; 
    double vznew=0;
    double vv,vvi;
    double r1,r2,r3,rr,rri;
    sumv[0]=0;sumv[1]=0;sumv[2]=0;
    sumv2[0]=0;sumv2[1]=0;sumv2[2]=0;
    int p,q;
    double fdotv[6];
    fdotv[0]=0.0;fdotv[1]=0.0;fdotv[2]=0.0;
    fdotv[3]=0.0;fdotv[4]=0.0;fdotv[5]=0.0; //dummy
    double exnew,eynew,eznew,ee,eei;
    double g1,g2,g3;

    for (i=0; i<Npart ; i++)
    {

// Euler-Maruyama algorithm
        g1 = gen_gaussian_rand();
        g2 = gen_gaussian_rand();
        g3 = gen_gaussian_rand();
        //Calculate the new positions
        xnew = particle[i].x + C1*(particle[i].fx + Fprop*particle[i].ex)*tstep + C2*g1*sqrt_tstep;
        ynew = particle[i].y + C1*(particle[i].fy + Fprop*particle[i].ey)*tstep + C2*g2*sqrt_tstep;
        znew = particle[i].z + C1*(particle[i].fz + Fprop*particle[i].ez)*tstep + C2*g3*sqrt_tstep;

///* disable for passive system testing
        r1 = gen_gaussian_rand();
        r2 = gen_gaussian_rand();
        r3 = gen_gaussian_rand();

        exnew = particle[i].ex + C3*(particle[i].ey*r3-particle[i].ez*r2)*sqrt_tstep;          //direction update
        eynew = particle[i].ey + C3*(particle[i].ez*r1-particle[i].ex*r3)*sqrt_tstep;          //direction update
        eznew = particle[i].ez + C3*(particle[i].ex*r2-particle[i].ey*r1)*sqrt_tstep;          //direction update
        ee = sqrt(exnew*exnew + eynew*eynew + eznew*eznew);
        eei = (ee>0) ? 1.0/ee : 0;

        exnew = exnew*eei;                          //make unit vector
        eynew = eynew*eei;                          //make unit vector
        eznew = eznew*eei;                          //make unit vector
        
        vxnew = C1*(particle[i].fx + Fprop*exnew) + C2*g1;
        vynew = C1*(particle[i].fy + Fprop*eynew) + C2*g2;
        vznew = C1*(particle[i].fz + Fprop*eznew) + C2*g3;
//*/
        particle[i].x = xnew;
        particle[i].y = ynew;
        particle[i].z = znew;
        pbc(i);
        particle[i].vx = vxnew;
        particle[i].vy = vynew;
        particle[i].vz = vznew;
        particle[i].ex = exnew;
        particle[i].ey = eynew;
        particle[i].ez = eznew;
       
        //Update system velocity and Kinetic energy
        //sumv[0]  += vv*vxnew;
        //sumv[1]  += vv*vynew;
        //sumv[2]  += vv*vznew;
        //sumv2[0] += vv*vv*vxnew*vxnew;
        //sumv2[1] += vv*vv*vynew*vynew;
        //sumv2[2] += vv*vv*vznew*vznew;
        
        sumv[0]  += vxnew;
        sumv[1]  += vynew;
        sumv[2]  += vznew;
        sumv2[0] += vxnew*vxnew;
        sumv2[1] += vynew*vynew;
        sumv2[2] += vznew*vznew;
        
#ifdef MEASURE_PRESSURE
        fdotv[0] +=  particle[i].fx*particle[i].ex;
        fdotv[1] +=  particle[i].fy*particle[i].ey;
        fdotv[2] +=  particle[i].fz*particle[i].ez;
#endif

    }
    recenter_com();
        
    //Update Total Energy & Temp
    etot = (en + 0.5*(sumv2[0]+sumv2[1]+sumv2[2]));
    Temp = (sumv2[0]+sumv2[1]+sumv2[2])/(3*Npart);

#ifdef MEASURE_PRESSURE    
    // Pressure as per virial equation
    for (s=0; s<6; s++)
    {   
        if(s<3)
        {
            //pressure[s] = (Npart*Temp + pressure[s])*volumei;
            avg_fdotv[s] = avg_fdotv[s] + fdotv[s];                 //avg_fdotv is the cumulative sum, divide by count to take avg.
            //pressure[s] = (Npart/beta + pressure[s] + friction_coeff*Fprop*avg_fdotv[s]/step)*volumei;
            //pressure[s] = (Npart/beta + pressure[s] + 0.5*(friction_coeff*Npart*Fprop*Fprop + Fprop*avg_fdotv[s]/step)/rot_diff_coeff)*volumei;
            pressure[s] = (Npart/beta + pressure[s] + (friction_coeff*Npart*Fprop*Fprop + Fprop*avg_fdotv[s]/step)/rot_diff_coeff)*volumei;
        }
        else
        {
            //pressure[s] = (Npart/beta + pressure[s])*volumei;
            pressure[s] = pressure[s]*volumei;
        }
        avg_press[s] = avg_press[s] + pressure[s];
        //fprintf(test_file,"%5d %.3f %.3f %.3f %.3f %.3f %.3f\n",step,fdotv[0],fdotv[1],fdotv[2],avg_fdotv[0]/step,avg_fdotv[1]/step,avg_fdotv[2]/step);
    }
    
    //Pressure profile in z-direction
    if(step%100 == 0)
    {
        //for(p=0; p<nzbin; p++)
        //{
        //    for(q=0; q<6; q++)
        //    {
        //        //press_z[p][q] = (Npart*Temp + press_z[p][q])*vzbini;
        //        press_z[p][q] = (Temp*rho_z[p]*100/step + press_z[p][q])*vzbini;
        //        avg_press_z[p][q] = avg_press_z[p][q] + press_z[p][q];
        //    }
        //}
        /*measure_press_z();*/
    }
#endif
    
}

//Integration function with Verlet algorithm
void integrate()
{
    int i,s;
    double xnew=0; 
    double ynew=0; 
    double znew=0;
    double vxnew=0; 
    double vynew=0; 
    double vznew=0;
    sumv[0]=0;sumv[1]=0;sumv[2]=0;
    sumv2[0]=0;sumv2[1]=0;sumv2[2]=0;
    int p,q;

    for (i=0; i<Npart ; i++)
    {

// Simple verlet algorithm
        //Calculate the new positions and velocity using Verlet algorithm
        xnew = 2*particle[i].x - particle[i].xold + particle[i].fx*tstep2;
        ynew = 2*particle[i].y - particle[i].yold + particle[i].fy*tstep2;
        znew = 2*particle[i].z - particle[i].zold + particle[i].fz*tstep2;
 
        vxnew = (xnew - particle[i].xold)*0.5*tstepi;
        vynew = (ynew - particle[i].yold)*0.5*tstepi;
        vznew = (znew - particle[i].zold)*0.5*tstepi;

/*
        // velocity verlet algorithm
        //Calculate the new positions and velocity using Velocity Verlet algorithm
        xnew = particle[i].x + particle[i].vx*tstep + particle[i].fx*tstep2;
        ynew = particle[i].y + particle[i].vy*tstep + particle[i].fy*tstep2;
        znew = particle[i].z + particle[i].vz*tstep + particle[i].fz*tstep2;
 
        vxnew = particle[i].vx + (particle[i].fxold+particle[i].fx)*tstep/2;
        vynew = particle[i].vy + (particle[i].fyold+particle[i].fy)*tstep/2;
        vznew = particle[i].vz + (particle[i].fzold+particle[i].fz)*tstep/2;
*/
        //Update the positions and velocity
        particle[i].xold = particle[i].x;
        particle[i].yold = particle[i].y;
        particle[i].zold = particle[i].z;
        particle[i].x = xnew;
        particle[i].y = ynew;
        particle[i].z = znew;
        pbc(i);
        particle[i].vx = vxnew;
        particle[i].vy = vynew;
        particle[i].vz = vznew;
       
        //Update system velocity and Kinetic energy
        sumv[0] += vxnew;
        sumv[1] += vynew;
        sumv[2] += vznew;
        sumv2[0] += vxnew*vxnew;
        sumv2[1] += vynew*vynew;
        sumv2[2] += vznew*vznew;
    }
    //Update Total Energy & Temp
    etot = (en + 0.5*(sumv2[0]+sumv2[1]+sumv2[2]));
    Temp = (sumv2[0]+sumv2[1]+sumv2[2])/(3*Npart);

#ifdef MEASURE_PRESSURE        
    // Pressure as per virial equation
    for (s=0; s<6; s++)
    {   
        pressure[s] = (Npart*Temp + pressure[s])/volume ;
        avg_press[s] = avg_press[s] + pressure[s];
    }
        
    if(step%100 == 0)
    {
        //Pressure profile in z-direction
        //for(p=0; p<nzbin; p++)
        //{
        //    for(q=0; q<6; q++)
        //    {
        //        press_z[p][q] = (Npart*Temp + press_z[p][q])*vzbini;
        //        press_z[p][q] = (Temp*rho_z[p]*100/step + press_z[p][q])*vzbini;
        //        avg_press_z[p][q] = avg_press_z[p][q] + press_z[p][q];
        //    }
        //}
        /*measure_press_z();*/
    }
#endif


}

// Function to calculate & print average properties
void sample()
{
    int i=0;
    int freq_z = 100;            //Frequency for measuring z-profile values

    //Calculate averages
    Av[0] += (en/Npart);
    Av[1] += (etot/Npart);
    Av[2] += Temp;
    Av[3] += ONEOVER3*(avg_press[0]+avg_press[1]+avg_press[2])/step;
    //Av[4] += sqrt(sqr(sumv[0])+sqr(sumv[1])+sqr(sumv[2]));
    Av[5] += 1;             //Counter for running average

    //Print out the current average every 100 steps
    if(step%100 == 0)
    {
        //fprintf (out_file,"%5d \ten/npart= %.3e \tEtot/npart= %.3e \tTemp= %.3e Pr= %.3e \tCOMv= %+.3e\n",step,Av[0]/Av[5],Av[1]/Av[5],Av[2]/Av[5], Av[3]/Av[5],sqrt(sqr(sumv[0])+sqr(sumv[1])+sqr(sumv[2]))/Npart);
        printf ("%5d \ten/npart= %.3e \tEtot/npart= %.3e \tTemp= %.3e Pr= %.3e \tCOMv= %+.3e\n",step,Av[0]/Av[5],Av[1]/Av[5],Av[2]/Av[5], Av[3]/Av[5],sqrt(sqr(sumv[0])+sqr(sumv[1])+sqr(sumv[2]))/Npart);
        for(i=0;i<6;i++)
        {
            Av[i]=0;
        }
        write_traj(0); //Write trajectory of particle 0
//        measure_rho_z();
//#ifdef MEASURE_PRESSURE
//        measure_press_z();
//        //measure_active_press_z();
//        //Write pressure tensor to file;
//        write_pressure();
//#endif
    }

    if(mode == 1)
    {
        if(step%freq_z == 0)
        {
            measure_rho_z();
            measure_press_z();
            write_pressure();
        }
        
#ifdef MEASURE_RDF
        if(step%1000 == 0)
        {
            //ngr +=1;
            rdf(1);
            diffusion(1,1000);
        }
#endif
    
        if(step%print == 0)
        {
            fprintf (prop_file,"%7d \t%+.3e \t%+.3e \t%+.3e \n",step,(en/Npart),Temp,(pressure[0]+pressure[1]+pressure[2])/3.0);
        }
    }    
}

void write_pressure()
{
    int p = 0;
    double Pr = ONEOVER3*(avg_press[0] + avg_press[1] + avg_press[2]);
    fprintf(press_file,"%+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e ",Pr/step, avg_press[0]/step,avg_press[1]/step,avg_press[2]/step,avg_press[3]/step,avg_press[4]/step,avg_press[5]/step);
    //fprintf(press_file,"%+.3e %+.3e %+.3e %+.3e\n",ONEOVER6*friction_coeff*Npart*Fprop*Fprop*volumei/rot_diff_coeff,ONEOVER6*(Fprop*avg_fdotv[0]/step/rot_diff_coeff)*volumei,ONEOVER6*(Fprop*avg_fdotv[1]/step/rot_diff_coeff)*volumei,ONEOVER6*(Fprop*avg_fdotv[2]/step/rot_diff_coeff)*volumei);
    fprintf(press_file,"%+.3e %+.3e %+.3e %+.3e\n",friction_coeff*Npart*Fprop*Fprop*volumei/rot_diff_coeff,(Fprop*avg_fdotv[0]/step/rot_diff_coeff)*volumei,(Fprop*avg_fdotv[1]/step/rot_diff_coeff)*volumei,(Fprop*avg_fdotv[2]/step/rot_diff_coeff)*volumei);
    
    //Reinitialize avg_stress tensor to 0.0 to measure avg. of only last 's' steps
    for (p =0; p<6; p++)
    {
        avg_press[p] = 0.0;
    }
}

//function to put particle in the periodic simulation box
void pbc(int i)
{
    while(particle[i].x > box.xhalf) { particle[i].x = fmod(particle[i].x,box.xlen)-box.xlen; particle[i].xold = fmod(particle[i].xold,box.xlen)-box.xlen; }
    while(particle[i].x <-box.xhalf) { particle[i].x = fmod(particle[i].x,box.xlen)+box.xlen; particle[i].xold = fmod(particle[i].xold,box.xlen)+box.xlen; }
    while(particle[i].y > box.yhalf) { particle[i].y = fmod(particle[i].y,box.ylen)-box.ylen; particle[i].yold = fmod(particle[i].yold,box.ylen)-box.ylen; }
    while(particle[i].y <-box.yhalf) { particle[i].y = fmod(particle[i].y,box.ylen)+box.ylen; particle[i].yold = fmod(particle[i].yold,box.ylen)+box.ylen; }
    while(particle[i].z > box.zhalf) { particle[i].z = fmod(particle[i].z,box.zlen)-box.zlen; particle[i].zold = fmod(particle[i].zold,box.zlen)-box.zlen; }
    while(particle[i].z <-box.zhalf) { particle[i].z = fmod(particle[i].z,box.zlen)+box.zlen; particle[i].zold = fmod(particle[i].zold,box.zlen)+box.zlen; }
#ifdef USE_NLIST
    move2cell(i);
#endif
}


/*
//Random no.(normally distributed) generator for Brownian term
double brownian_gen()
{
    double lx,ly,lz;
    int i=0;
        lambda_tr[i] = lambda_tr[i] + (dsfmt_genrand_open_open(&dsfmt)-0.5);
    return lambda_tr;
      
}
*/

// Random no. generator wihin Normal distribution using Marsaglia Polar method 
// Refer https://en.wikipedia.org/wiki/Marsaglia_polar_method
// Variance is taken to be 1.0 for our case
double gen_gaussian_rand()
{
    static bool hasSpare = false;
    static double spare;
    static double u, v, s;
    
    if(hasSpare)
    {
        hasSpare = false;
        return spare;
    }
    
    hasSpare = true;

    do
    {
        u = 2*dsfmt_genrand_open_open(&dsfmt) - 1;
        v = 2*dsfmt_genrand_open_open(&dsfmt) - 1;
        s = u*u + v*v;
    }
    while(s >= 1 || s == 0);
    
    s = sqrt(-2.0 * log(s) / s);
    spare = v*s;
    return u*s;
}

//
//Function to evaluate mean propulsion velocity of particle
//
//
void eval_vprop()
{
    int i;
    vprop[0] = 0.0;
    vprop[1] = 0.0;
    vprop[2] = 0.0;
    double v[3];
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    
    //for(i=0; i<Npart; i++)
    //{
    //    vprop[0] += Fprop*particle[i].vx + particle[i].fx*particle[i].vx;
    //    vprop[1] += Fprop*particle[i].vy + particle[i].fy*particle[i].vy;
    //    vprop[2] += Fprop*particle[i].vz + particle[i].fz*particle[i].vz;
    //    
    //    v[0] += Fprop*particle[i].vx;
    //    v[1] += Fprop*particle[i].vy;
    //    v[2] += Fprop*particle[i].vz;
    //}
    //vprop[0] = vprop[0]/(Npart*friction_coeff);
    //vprop[1] = vprop[1]/(Npart*friction_coeff);
    //vprop[2] = vprop[2]/(Npart*friction_coeff);
        
    //fprintf(test_file," %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n",v[0]/Npart,v[1]/Npart,v[2]/Npart,vprop[0],vprop[1],vprop[2]);
}

//Read old configuration file and load position & velocities
void read_conf()
{
    int dummy = 0;
    int i = 0;
    int j = 0;
    double d1,d2,d3;
    FILE *fp;
    fp = fopen("pos.dat","r");   
    while(j<Npart && i<Npart)
    {
        dummy = fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&i,&particle[j].x,&particle[j].y,&particle[j].z,&particle[j].ex,&particle[j].ey,&particle[j].ez,&particle[j].xold,&particle[j].yold,&particle[j].zold);
        j += 1;
    }
    if(i>Npart || j!=Npart)
    {
        printf("Error!: check conf file & Npart,id=%d, j=%d\n",i,j);
    }
    fclose(fp);
    
    //for(i=0;i<Npart;i++)
    //{
    //    //Initialize variables with dummy values
    //    particle[i].xold = particle[i].x;
    //    particle[i].yold = particle[i].y;
    //    particle[i].zold = particle[i].z;
    //}
    
}



//Write the coordinates of the particles to a file
void write_pos(int mode)
{
    static int print_count = 0;
    char filename[255];
    int i=0;

    if(mode == 0)
    {
        pos_file    = fopen("pos.dat","w");
//        fprintf(pos_file,"#Npart%d\n",Npart);
//        fprintf(pos_file,"#box.x=%.3e, \tbox.y=%.3e, \tbox.z=%.3e, \tbox.xhalf=%.3e, \tbox.yhalf=%.3e, \tbox.zhalf=%.3e\n",box.xlen,box.ylen,box.zlen,box.xhalf,box.yhalf,box.zhalf);
//        fprintf(pos_file,"#ID,\tx,\ty,\tz,\tvx,\tvy,\tvz,\tfx,\tfy,\tfz\n");
        for (i=0;i<Npart;i++)
        {
            fprintf(pos_file,"%5d \t%+.8e \t%+.8e \t%+.8e \t%+.8e \t%+.8e \t%+.8e \t%+.8e \t%+.8e \t%+.8e\n",i,particle[i].x,particle[i].y,particle[i].z,particle[i].ex,particle[i].ey,particle[i].ez,particle[i].xold,particle[i].yold,particle[i].zold);
        }
//      fprintf(pos_file,"#COMvx=%e,\tCOMvy=%e,\tCOMvz=%e\n",sumv[0]/Npart,sumv[1]/Npart,sumv[2]/Npart);
//      fprintf(pos_file,"\n\nCOMvx=%e,\tCOMvy=%e,\tCOMvz=%e\n",sumv[0],sumv[1],sumv[2]);
        fclose(pos_file);
    }
    else if(mode == 1)
    {

        //Write particle configuration in standard xyz format
        //sprintf(filename,"pos_%.05i.xyz",++print_count);
        //conf_file = fopen(filename,"a");
        fprintf(conf_file,"%d\n",Npart);
        fprintf(conf_file,"VMD_XYZ: %d %d\n",++print_count,step);
        //fprintf(conf_file,"%i\n",Npart);
        //fprintf(conf_file,"%le %le\n",-box.xhalf,box.xhalf);
        //fprintf(conf_file,"%le %le\n",-box.yhalf,box.yhalf);
        //fprintf(conf_file,"%le %le\n",-box.zhalf,box.zhalf);

        for (i=0;i<Npart;i++)
        {
            fprintf(conf_file,"%5d %+.3e %+.3e %+.3e %+0.3e %+0.3e %+0.3e \t%+0.3e %+0.3e %+0.3e\n",i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].ex,particle[i].ey,particle[i].ez);
        }
        //fclose(conf_file);
    }

}

//Write trajectory of the specified particle to a file
void write_traj(int i)
{

    fprintf(traj_file,"%6d %4d \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e\n",step,i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].fx,particle[i].fy,particle[i].fz);

}


//Main MD simulation routine
void main()
{
    int p;
    clock_t start = clock(), time_diff;          //Start time of program

    traj_file   = fopen("traj.dat","w");
    press_file  = fopen("press.dat","w");
    prop_file   = fopen("prop.dat","w");
    out_file    = fopen("out.dat","w");
    test_file   = fopen("snap.dat","w");
    fclose(test_file);
    test_file   = fopen("pos.xyz","w");
    fclose(test_file);
    
    snap_file = fopen("snap.dat","a");
    conf_file = fopen("pos.xyz","a");
    diff_file = fopen("diff.dat","w");

    printf("Begin Initialization\n");
    initialize();
    printf("Initialized\n");
    printf("Density=%.3f\n",Npart/volume);
    printf("xmin=%.3f, xmax=%.3f, ymin=%.3f, ymax=%.3f, zmin=%.3f, zmax=%.3f\n",-box.xhalf,box.xhalf,-box.yhalf,box.yhalf,-box.zhalf,box.zhalf);
    printf("Temp=%.3f\n",1.0/beta);
    printf("friction_coeff=%.3f\n",friction_coeff);
    printf("beta=%.3f\n",beta);
    printf("epsilon=%.3f\n",epsilon);
    printf("Dr=%.3f\n",rot_diff_coeff);
    printf("Fp=%.3f\n",Fprop);
    printf("Begin Equilibration\n");

    //Equilibration
    mode = 0;
    while (step<nsteps_equil)
    {
        
        step += 1;
#ifdef USE_NLIST
        force_nlist();
#else
        force();
#endif
        integrate_euler();
        //integrate();
//        time = time + tstep;
//        sample();
    }

    printf("Begin MD\n");
    //Production
    mode = 1;
    step = 0;
    while (step<nsteps)
    {
        step += 1;

#ifdef USE_NLIST
        force_nlist();
#else
        force();
#endif
        
        sample();

        integrate_euler();
        //integrate();
        //time = time + tstep;

        
        if(step%1000 == 0)
        {
            write_rho_z();
            write_press_z();
            write_pos(0);       //Print last configuration for continuing the simulation run
            rdf(2);             //Print RDF values
            eval_vprop();
            diffusion(2,1000);             //Print diff_coeff values
        }
        if(step%100000 ==0)  //every 10^5 steps
        {
            write_pos(1);           //write conf in xyz format
            save_z_snap();          //Save snapshot of z-profiles of density and pressure
        }

    }
    
    printf("End program\n");
    
    fclose(prop_file);
    fclose(traj_file);
    //fclose(test_file);
    fclose(out_file);
    fclose(snap_file);
    fclose(conf_file);
    fclose(diff_file);
    
    //Compute total time taken for the simulation run
    time_diff = clock() - start;
    int msec = time_diff*1000 / CLOCKS_PER_SEC;
    printf("Time taken %d minutes %d seconds %d milliseconds\n", (msec/1000/60), (msec%60000)/1000, msec%1000);
}

