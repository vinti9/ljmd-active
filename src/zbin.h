#define MAXZBIN 500

//Global variable declarations
int     nzbin;          //No. of bins in z-dir
double  dzbin;          //Bin width
double  vzbin;          //Volume of each bin
double  dzbini;         //dzbin inverse
double  vzbini;         //vzbin inverse


//Global function declarations
int coords2zbin(double);
void move2zbin(int);
void init_zbins();
void measure_rho_z();
void measure_press_z();
void write_rho_z();
void write_press_z();
void measure_active_press_z();

//Slab structure storing properties for bins along z direction
typedef struct
{
    double  rho_z;              //density in bin
    double  press[6];         //inst. pressure
    double  avg_rho_z;          //density avg.
    double  avg_press[6];     //avg. pressure
    double  press_active[6];  //inst. pressure due to self-propulsion
    double  avg_press_active[6];  //avg. pressure due to self-propulsion
    int     n;            //No. of particles in bin
} slabs_z_t;

//Global structure declaration
slabs_z_t   slabs_z[MAXZBIN];