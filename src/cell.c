
#include "cell.h"
#include "params.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>



//Calculate the cell id acording to coordinates
//
int coords2cell(double x, double y, double z)
{
    int i = ((int)(ncellx*(x+box.xhalf)*box.xleni) % ncellx)*ncelly*ncellz +
            ((int)(ncelly*(y+box.yhalf)*box.yleni) % ncelly)*ncellz +
            ((int)(ncellz*(z+box.zhalf)*box.zleni) % ncellz);
    return i;
}


//Initialize the cell list

void init_cells(void)
{
  n_cells=0;							
  int x,y,z,i,j,k;
  int dx,dy,dz,a,b,c;
  ncellx = floor(0.95*box.xlen/(rcut+0.1));				//how many cells "fit" in one direction 
  ncelly = floor(0.95*box.ylen/(rcut+0.1));
  ncellz = floor(0.95*box.zlen/(rcut+0.1));
  n_cells = ncellx*ncelly*ncellz;
  printf("n_cells=%d,ncellx=%d,ncelly=%d,ncellz=%d\n",n_cells,ncellx,ncelly,ncellz);
  if (n_cells > MAXNCELLS)
  {
      printf("Error: Check MAXNCELLS"); exit(1);
  }
//  cells = (cell_t *) malloc(n_cells*sizeof(cell_t));		//Allocate the memory to the cells array

  for(x=0; x<ncellx; x++)
  {
      for(y=0; y<ncelly; y++)
      {
         for(z=0; z<ncellz; z++)
         {	  
             i = x*ncelly*ncellz + y*ncellz + z;    //Cell ID
             cells[i].n = 0;                        //initialize with 0
             k = 0;
             for(a=-1;a<2;a++) 
             {
                 for(b=-1;b<2;b++) 
                 {
                     for(c=-1;c<2;c++)
                     { 	
                         dx=a;dy=b;dz=c;
                         if(x+dx < 0) dx=ncellx-1; if(x+dx > ncellx-1) dx=-ncellx+1;
                         if(y+dy < 0) dy=ncelly-1; if(y+dy > ncelly-1) dy=-ncelly+1;
                         if(z+dz < 0) dz=ncellz-1; if(z+dz > ncellz-1) dz=-ncellz+1;
                         j = (x+dx)*ncelly*ncellz+(y+dy)*ncellz+(z+dz);
                         cells[i].neighbors[k]=j;
                         k++;
                     }
                 }
              }
         }
      }
  }

  for(i=0; i<Npart; i++){
    j = coords2cell(particle[i].x,particle[i].y,particle[i].z);
    if(cells[j].n+1 >= MAXNPART_CELL) 
    {
        printf("ERROR: Too many particles in a cell!,i=%d, cell=%d ,npart=%d\n",i,j,cells[j].n); exit(666);
    }
    cells[j].particles[cells[j].n] = i;	            //Store the particle no. i corresponding cell
    particle[i].cell = j;                           //Store the cell no. in cooresponding particle
    cells[j].n++;							        
  }
  printf("number of cells is %d\n",n_cells);
}



//Put particle in the correct cell
void move2cell(int i)
{
  int j,c_new,c_old;
  
  c_new = coords2cell(particle[i].x,particle[i].y,particle[i].z);  
  c_old = particle[i].cell;					
  if(c_new < 0 || c_new >= n_cells) 
  { 
      printf("Error: out of box?, i=%d, c_new=%d, c_old=%d\n",i,c_new,c_old);  
      printf("i.x=%.3f,i.y=%.3f,i.z=%.3f\n",particle[i].x,particle[i].y,particle[i].z);
      exit(1); 
  } 
  if(c_new != c_old)
  {						
    for(j=0; j<cells[c_old].n;j++)
    {				
      if((cells[c_old].particles[j]) == i)
      {					
        cells[c_old].particles[j] = cells[c_old].particles[cells[c_old].n-1]; 
        j = cells[c_old].n;
        break;
      }
    }
    cells[c_old].n--;					
    cells[c_new].particles[cells[c_new].n]=i;			
    cells[c_new].n++;					
    particle[i].cell = c_new;
  }							
}

