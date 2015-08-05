#define MAXZBIN 5000

//Global variable declarations
int     nzbin;          //No. of bins in z-dir
double  dzbin;          //Bin width
double  vzbin;          //Volume of each bin
double  dzbini;         //dzbin inverse
double  vzbini;         //vzbin inverse
double  rho_z[MAXZBIN];         //density profile
double  press_z[MAXZBIN][6];    //inst. pressure profile
double  avg_press_z[MAXZBIN][6];    //avg. pressure profile


//Global function declarations
int coords2zbin(double);
void move2zbin(int);
void init_zbins();
void measure_rho_z();
void measure_press_z();
void write_rho_z();
void write_press_z();