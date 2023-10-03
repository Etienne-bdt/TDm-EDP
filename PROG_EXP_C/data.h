#ifndef DATA_H_INCLUDED
#define DATA_H_INCLUDED

typedef struct type_donneesc
{   int nx;       // nb of cells in x-direction
    int ny;       // nb of cells in y-direction
    float Lx;     // x-length of the domain
    float Ly;     // y-length of the domain
    float U0;     // characteristic velocity
    float alpha;  // characteristic velocity
    float beta;   // characteristic velocity
    float gama;   // characteristic velocity
    float D;      // thermal/mass diffusivity
    float Ti;     // Initial temperature in the domain
    float Tg;     // Temperature at the left
    float Td;     // Temperature at the right
    float Tt;     // Temperature at the top
    float Tb;     // Temperature at the bottom
    float tf;     // final time
    int  Nout;    // number of intermediate unsaved time steps
    float CFL;    // Courant's number (advection)
    float R;      // Fourier's number (diffusion)
    int i_mesh;   // Switch for the regular/irregular mesh (0: reg, 1:irreg)
    int i_vit;    // Switch for the velocity field
    int i_solver; // Solver type (0: explicit, 1: implicit gauss, 2: SOR)
    int i_bc;     // BC type (0:uniform, 1: exponential)
    float W;      // Over-relaxation coefficient (SOR method)
    float R0;     // Relative precision (SOR method)

}type_donneesc;

#endif // DATA_H_INCLUDED
