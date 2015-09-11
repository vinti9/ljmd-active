#define MAXRDFBIN 200
//Variable for RDF
int nbin_gr;  //no. of bins
int     ngr;        //No. of evaluations of g(r)
float   g[MAXRDFBIN];     //Array storing g(r) graph
double  delg;       //Step size for RDF function
double r2t[3];          //Mean squared displacement for the system
double vacf[3];         //Velocity Autocorrelation
double eacf[3];         //Direction Autocorrelation 
double *xinit,*yinit,*zinit;    //Initial positions for Mean squared displacement
double *vx0,*vy0,*vz0;          //Initial velocities for Velocity autocorrelation
double *ex0,*ey0,*ez0;          //Initial direction for direction autocorrelation

//RDF function
void rdf(int);

void recenter_com();
void diffusion(int,int);
