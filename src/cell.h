
#define MAXNPART_CELL 60
#define MAXNCELLS 12*12*12

//Cell list structure for storing linked-list of neighbring cells
typedef struct cells_t
{
    int neighbors[27];              //stores of ID of neighboring cells
    int n;                          //Stores no. of particles in the current cell
    int particles[MAXNPART_CELL];   //Stores ID of particles in the current cell
} cell_t;


//Global structure for cells
cell_t cells[MAXNCELLS];


//Global Cell variables
int     n_cells, ncellx, ncelly, ncellz;

//Global functions
int     coords2cell(double, double, double);
void    init_cells(void);
void    move2cell(int);
