#ifndef fem_amg_preco_hpp_
#define fem_amg_preco_hpp_
#include "gslib.h"

long long * get_row_start();
long long * get_row_end();
long long * get_dof_map();
void fem_amg_setup(const sint *n_x_, const sint *n_y_, const sint *n_z_, 
                   const sint *n_elem_, const sint *n_dim_, 
                   double *x_m_, double *y_m_, double *z_m_, 
                   double *pmask_, const sint *nullspace,
                   double *param,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes
                   );
void fem_amg_solve(double*,double*);

typedef double (*Basis)(double*);
typedef void (*DBasis)(double**, double*);

void matrix_distribution();
void fem_assembly();
void quadrature_rule(double***, double**, int, int);
void mesh_connectivity(int***, int***, int, int);
void x_map(double**, double*, double**, int, Basis*);
void J_xr_map(double***, double*, double**, int, DBasis*);

double phi_2D_1(double*);
double phi_2D_2(double*);
double phi_2D_3(double*);
void dphi_2D_1(double**, double*);
void dphi_2D_2(double**, double*);
void dphi_2D_3(double**, double*);

double phi_3D_1(double*);
double phi_3D_2(double*);
double phi_3D_3(double*);
double phi_3D_4(double*);
void dphi_3D_1(double**, double*);
void dphi_3D_2(double**, double*);
void dphi_3D_3(double**, double*);
void dphi_3D_4(double**, double*);

double determinant(double**, int);
void inverse(double***, double**, int);
long long maximum(long long, long long);

int* mem_alloc_1D_int(int);
long long* mem_alloc_1D_long(int);
double* mem_alloc_1D_double(int);
int** mem_alloc_2D_int(int, int);
double** mem_alloc_2D_double(int, int);
void mem_free_1D_int(int**, int);
void mem_free_1D_long(long long**, int);
void mem_free_1D_double(double**, int);
void mem_free_2D_int(int***, int, int);
void mem_free_2D_double(double***, int, int);

struct gs_data {
  struct comm comm;
};
#endif