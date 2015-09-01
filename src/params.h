#include <stdio.h>

#define USE_NLIST       //Switch to use neighbor list 
#define MEASURE_RDF     //Switch to enable RDF measurement 
#define MEASURE_PRESSURE     //Switch to enable Pressure measurement

#define MAXNPART 5100
#define ONEOVER3 0.3333333333
#define ONEOVER6 0.1666666667

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

//Global variable Declarations
//double  time;
double  tstep;
double  tstep2;                 //sqr of timestep
double  tstepi;                 //inverse of timestep
int     step;
int     nsteps;                 //Production steps
int     nsteps_equil;           //Equilibration steps
int     Npart;
double  sumv[3];
double  sumv2[3];
double  rcut;
double  rcut2;
double  ecut;
double  en;                     //Potential energy
double  etot;                   //Total energy
double  Temp;                   //Temperature
double  initTemp;
double  beta;
double  density;
double  volume;
double  volumei;                //Inverse of volume
double  epsilon;                //interaction term for LJ potential
double  rot_diff_coeff;         //Rotational diff. coeff.
double  Fprop;                  //self-propulsion force
double  vprop[3];

double  pressure[6];
double  avg_press[6];
double  avg_fdotv[6];           //Cumulative sum of fdotv[], divide by count to get avg. value
double  Av[6];                  //Array storing average values of properties

int     continu;                //parameter indicating a continuing simulation run or starting fresh
int     lattice;                //parameter indicating type of lattice(random/cubic/bcc...etc) for initialization
int     mode;                   //Mode = 0: Equilibration/ 1: Production
int     print;                   //Frequency of printing data to prop.dat

//Variables for brownian dynamics
double  friction_coeff;                  //friction coefficient = gamma
double  C1;                     //factor for force term of Langevin equation
double  C2;                     //factor for brownian term of Langevin equation
double  C3;                     //factor for velocity propagation for rotational diffusion
double  sqrt_tstep;



//Structure for Particles
typedef struct particle_t
{
    double x;
    double y;
    double z;
    double xold;
    double yold;
    double zold;
    double vx;
    double vy;
    double vz;
    double fx;
    double fy;
    double fz;
//    double fxold;
//    double fyold;
//    double fzold;
    int    cell;
    int    zbin;
} particle_t;


//Structure for Simulation Box
typedef struct box_t
{
    double xlen;
    double ylen;
    double zlen;
    double xhalf;
    double yhalf;
    double zhalf;
    double xleni;
    double yleni;
    double zleni;
} box_t;


//Global Variables
particle_t  particle[MAXNPART];
box_t       box;

FILE *test_file;
FILE *snap_file;