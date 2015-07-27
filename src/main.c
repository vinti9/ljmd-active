/*
 * This file contains the 'main subroutine' and function calls
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
//#include <time.h>
#include "dSFMT.h"      //for Random Number generator
#include "params.h"
#include "cell.h"

#define sqr(x) ((x)*(x))
//#define MEASURE_PRESSURE

//Global function Declarations
void    initialize();
void    init_lattice_pos();
void    init_lattice_vel();
void    force();
void    force_nlist();
void    partforce(int,int);
void    integrate_euler();
void    sample();
void    pbc(int);
void    write_pos();
void    write_traj();
void    init_box(double);
double  gen_gaussian_rand();
void    write_pressure();
void    read_conf();
        

//Global Structure declarations

//cells_t     *cells;
dsfmt_t     dsfmt;

//Global File pointers
FILE *prm_file;
FILE *prop_file;
FILE *pos_file;
FILE *traj_file;
FILE *press_file;

//Initialization of all global parameters
void initialize()
{
    int dummy;
    int p;
    int seed;
    char str[30];
    
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
    fclose(prm_file);

    dsfmt_init_gen_rand(&dsfmt,seed);

    time = 0;
    step = 0;
    tstep2 = tstep*tstep;
    tstepi = 1/tstep;
    rcut2 = rcut*rcut;
    Temp = 0.0;
    
    etot = 0.0;
    en = 0.0;
    init_box(density);

//Variables for brownian dynamics
    C1 = 1.0/friction_coeff;
    sqrt_tstep = sqrt(tstep);
    C2 = pow(2.0,0.5)/(beta*friction_coeff);
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
        init_lattice_pos();
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

    
    ecut = 4*beta*epsilon*(1.0/(pow(rcut,12))-1.0/(pow(rcut,6)));        //calculate ecut
    printf("ecut=%e\n",ecut);

    //Initialize stress tensor
    for (p =0; p<6; p++)
    {
        avg_press[p] = 0.0;
    }
  
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
    int i,j;
    int p,q,c,icell,inbr;
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

    en = 0.0;           //Initialize potential energy to zero

    for(i=0;i<Npart;i++)
    {
        icell = particle[i].cell;
        for (p=0; p<27; p++)
        {
            inbr = cells[icell].neighbors[p];
            for (q=0; q<cells[inbr].n; q++)
            {
                j = cells[inbr].particles[q];
                if (j>i)
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
        ff = epsilon*48.0*r2i*r6i*(r6i-0.5);

        //Update the forces
        particle[i].fx=particle[i].fx+ff*xr;
        particle[i].fy=particle[i].fy+ff*yr;
        particle[i].fz=particle[i].fz+ff*zr;
        particle[j].fx=particle[j].fx-ff*xr;
        particle[j].fy=particle[j].fy-ff*yr;
        particle[j].fz=particle[j].fz-ff*zr;

#ifdef MEASURE_PRESSURE 
        //Update the stress tensor components
        pressure[0] = pressure[0] + ff*xr*xr;   //Pxx
        pressure[1] = pressure[1] + ff*yr*yr;   //Pyy
        pressure[2] = pressure[2] + ff*zr*zr;   //Pzz
        pressure[3] = pressure[3] + ff*xr*yr;   //Pxy
        pressure[4] = pressure[4] + ff*yr*zr;   //Pyz
        pressure[5] = pressure[5] + ff*zr*xr;   //Pzx
#endif 

        //Update the system potential energy
        en = en + 4*epsilon*beta*r6i*(r6i-1) - ecut;

    }
}



void init_lattice_pos()
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

/* Random positions
    for (i=0;i<Npart;i++)
    {

        particle[i].x = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.xhalf;
        particle[i].y = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.yhalf;
        particle[i].z = 2*(dsfmt_genrand_open_open(&dsfmt)-0.5)*box.zhalf;
    }
*/
///*
//      Positions on a regular cubic lattice
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
//*/

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

    for (i=0; i<Npart ; i++)
    {

// Euler-Maruyama algorithm
        //Calculate the new positions
        xnew = particle[i].x + C1*(particle[i].fx + Fprop*particle[i].vx)*tstep + C2*gen_gaussian_rand()*sqrt_tstep;
        ynew = particle[i].y + C1*(particle[i].fy + Fprop*particle[i].vy)*tstep + C2*gen_gaussian_rand()*sqrt_tstep;
        znew = particle[i].z + C1*(particle[i].fz + Fprop*particle[i].vz)*tstep + C2*gen_gaussian_rand()*sqrt_tstep;

/* 
        vxnew = (xnew - particle[i].x)*tstepi;
        vynew = (ynew - particle[i].y)*tstepi;
        vznew = (znew - particle[i].z)*tstepi;
*/

        r1 = gen_gaussian_rand();
        r2 = gen_gaussian_rand();
        r3 = gen_gaussian_rand();
//        rr = sqrt(r1*r1 + r2*r2 + r3*r3);
//        rri = 1.0/sqrt(r1*r1 + r2*r2 + r3*r3);
//        r1 = r1*rri;
//        r2 = r1*rri;
//        r3 = r1*rri;
//        vxnew = C3*vv*r1*rri*sqrt_tstep;
//        vynew = C3*vv*r2*rri*sqrt_tstep;
//        vznew = C3*vv*r3*rri*sqrt_tstep;

        vxnew = C3*(particle[i].vy*r3-particle[i].vz*r2)*sqrt_tstep;          //direction update
        vynew = C3*(particle[i].vz*r1-particle[i].vx*r3)*sqrt_tstep;          //direction update
        vznew = C3*(particle[i].vx*r2-particle[i].vy*r1)*sqrt_tstep;          //direction update
        vv = sqrt(vxnew*vxnew + vynew*vynew + vznew*vznew);
        vvi = 1.0/vv;

        //Update the positions and velocity
//        particle[i].xold = particle[i].x;
//        particle[i].yold = particle[i].y;
//        particle[i].zold = particle[i].z;
        particle[i].x = xnew;
        particle[i].y = ynew;
        particle[i].z = znew;
        pbc(i);
        particle[i].vx = vxnew*vvi;                          //make unit vector & update
        particle[i].vy = vynew*vvi;                          //make unit vector & update
        particle[i].vz = vznew*vvi;                          //make unit vector & update
       
/*
        //Update system velocity and Kinetic energy
        sumv[0] += vxnew;
        sumv[1] += vynew;
        sumv[2] += vznew;
        sumv2[0] += vxnew*vxnew;
        sumv2[1] += vynew*vynew;
        sumv2[2] += vznew*vznew;
*/
        //Update Total Energy & Temp
        etot = (en + 0.5*(sumv2[0]+sumv2[1]+sumv2[2]));
        Temp = (sumv2[0]+sumv2[1]+sumv2[2])/(3*Npart);

#ifdef MEASURE_PRESSURE    
        // Pressure as per virial equation
        for (s=0; s<6; s++)
        {   
            pressure[s] = (Npart*Temp + pressure[s])*volumei ;
            avg_press[s] = avg_press[s] + pressure[s];
        }
#endif

    }
    
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
#endif
    }

}

// Function to calculate & print average properties
void sample()
{
    int i=0;
//    printf ("TotalEnergy\t E/Npart \t Temp");
//    printf ("Step=%5d,\tEtot= %.3e,\ten/npart= %.3e,\tEtot/npart= %.3e,\tTemp= %.3e,\tCOMv=%.3e\n",step,etot,(en/Npart),(etot/Npart),Temp,(sumv[0]+sumv[1]+sumv[2])/Npart);

    //Calculate averages
    Av[0] += (en/Npart);
    Av[1] += (etot/Npart);
    Av[2] += Temp;
    Av[3] += (pressure[0] + pressure[1] + pressure[2])/3.0;
    Av[5] += 1;

    //Print out the current average every 100 steps
    if(step%100 == 0)
    {
        printf ("%5d \ten/npart= %.3e \tEtot/npart= %.3e \tTemp= %.3e \tCOMv= %+.3e\n",step,Av[0]/step,Av[1]/step,Av[2]/step,(sumv[0]+sumv[1]+sumv[2])/Npart);
        Av[5] = 0;
        write_traj(0); //Write trajectory of particle 0
#ifdef MEASURE_PRESSURE
        //Write pressure tensor to file;
        write_pressure();
#endif

    }

    fprintf (prop_file,"%7d \t%+.3e \t%+.3e \t%+.3e \n",step,(en/Npart),Temp,(pressure[0]+pressure[1]+pressure[2])/3.0);
    
}

void write_pressure()
{
    double Pr = (avg_press[0] + avg_press[1] + avg_press[2])/3.0;
    fprintf(press_file,"%+.3e %+.3e %+.3e %+.3e %+.3e %+.3e %+.3e\n",Pr/step, avg_press[0]/step,avg_press[1]/step,avg_press[2]/step,avg_press[3]/step,avg_press[4]/step,avg_press[5]/step);
}

//function to put particle in the periodic simulation box
void pbc(int i)
{
    while(particle[i].x > box.xhalf) { particle[i].x = fmod(particle[i].x,box.xlen)-box.xlen; particle[i].xold = fmod(particle[i].x,box.xlen)-box.xlen; }
    while(particle[i].x <-box.xhalf) { particle[i].x = fmod(particle[i].x,box.xlen)+box.xlen; particle[i].xold = fmod(particle[i].x,box.xlen)+box.xlen; }
    while(particle[i].y > box.yhalf) { particle[i].y = fmod(particle[i].y,box.ylen)-box.ylen; particle[i].yold = fmod(particle[i].y,box.ylen)-box.ylen; }
    while(particle[i].y <-box.yhalf) { particle[i].y = fmod(particle[i].y,box.ylen)+box.ylen; particle[i].yold = fmod(particle[i].y,box.ylen)+box.ylen; }
    while(particle[i].z > box.zhalf) { particle[i].z = fmod(particle[i].z,box.zlen)-box.zlen; particle[i].zold = fmod(particle[i].z,box.zlen)-box.zlen; }
    while(particle[i].z <-box.zhalf) { particle[i].z = fmod(particle[i].z,box.zlen)+box.zlen; particle[i].zold = fmod(particle[i].z,box.zlen)+box.zlen; }
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

//Read old configuration file and load position & velocities
void read_conf()
{
    int dummy = 0;
    int i = 0;
    int j = 0;
    double d1,d2,d3;
    FILE *fp;
    fp = fopen("pos.dat","r");   
    while(dummy!=EOF && i<Npart)
    {
        dummy = fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",&i,&particle[j].x,&particle[j].y,&particle[j].z,&particle[j].vx,&particle[j].vy,&particle[j].vz,&d1,&d2,&d3);
        j += 1;
    }
    if(i>Npart || j!=Npart)
    {
        printf("Error!: check conf file & Npart,id=%d, j=%d\n",i,j);
    }
    fclose(fp);
}



//Write the coordinates of the particles to a file
void write_pos()
{
    pos_file    = fopen("pos.dat","w");
    int i=0;
//    fprintf(pos_file,"#Npart%d\n",Npart);
//    fprintf(pos_file,"#box.x=%.3e, \tbox.y=%.3e, \tbox.z=%.3e, \tbox.xhalf=%.3e, \tbox.yhalf=%.3e, \tbox.zhalf=%.3e\n",box.xlen,box.ylen,box.zlen,box.xhalf,box.yhalf,box.zhalf);
//    fprintf(pos_file,"#ID,\tx,\ty,\tz,\tvx,\tvy,\tvz,\tfx,\tfy,\tfz\n");
    for (i=0;i<Npart;i++)
    {
        fprintf(pos_file,"%5d \t%+.3e \t%+.3e \t%+.3e \t%+.3e \t%+.3e \t%+.3e \t%+.3e \t%+.3e \t%+.3e\n",i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].fx,particle[i].fy,particle[i].fz);
    }
//    fprintf(pos_file,"#COMvx=%e,\tCOMvy=%e,\tCOMvz=%e\n",sumv[0]/Npart,sumv[1]/Npart,sumv[2]/Npart);
//    fprintf(pos_file,"\n\nCOMvx=%e,\tCOMvy=%e,\tCOMvz=%e\n",sumv[0],sumv[1],sumv[2]);
    fclose(pos_file);

//Write particle configuration in standard xyz format
    FILE *conf_file;
    conf_file = fopen("pos.xyz","w");
    for (i=0;i<Npart;i++)
    {
        fprintf(conf_file,"%5d \t%+.3e \t%+.3e \t%+.3e\n",i,particle[i].x,particle[i].y,particle[i].z);
    }
    fclose(conf_file);

}

//Write trajectory of the specified particle to a file
void write_traj(int i)
{

    fprintf(traj_file,"%d %d \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e\n",step,i,particle[i].x,particle[i].y,particle[i].z,particle[i].vx,particle[i].vy,particle[i].vz,particle[i].fx,particle[i].fy,particle[i].fz);

}


//Main MD simulation routine
void main()
{
    int p;
    

    traj_file   = fopen("traj.dat","w");
    press_file  = fopen("press.dat","w");
    prop_file   = fopen("prop.dat","w");

    printf("Begin Initialization\n");
    initialize();
    printf("Initialized\n");
//    printf("%e",fmod(3.018,3.000));
    printf("Density=%.3f\n",density);
    printf("xmin=%.3f, xmax=%.3f, ymin=%.3f, ymax=%.3f, zmin=%.3f, zmax=%.3f\n",-box.xhalf,box.xhalf,-box.yhalf,box.yhalf,-box.zhalf,box.zhalf);
    printf("beta=%.3f\n",beta);
    printf("friction_coeff=%.3f\n",friction_coeff);
    printf("epsilon=%.3f\n",epsilon);
    printf("Begin Equilibration\n");

    //Equilibration
    while (step<nsteps_equil)
    {
#ifdef USE_NLIST
        force_nlist();
#else
        force();
#endif
        integrate_euler();
//        integrate();
//        time = time + tstep;
        step += 1;
//        sample();
    }

    printf("Begin MD\n");
    //Production
    step = 0;
    while (step<nsteps)
    {
#ifdef USE_NLIST
        force_nlist();
#else
        force();
#endif
        integrate_euler();
//        integrate();
        time = time + tstep;
        step += 1;
       

        sample();
    }
    
    printf("End program\n");


    write_pos();
    fclose(prop_file);
    fclose(traj_file);
}

