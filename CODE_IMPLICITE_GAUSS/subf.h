struct type_donneesc read_datac();
void meshc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **vol);
void initial_conditionc(struct type_donneesc param,float **xv,float **yv, float **x,float **y,float **T0,float **U,float **V);
float calc_dtc(int nx, int ny,  float **x, float **y,float **U,float CFL,float **V, float D, float R);
void calc_flux_advc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **U, float **V, float **T0,float **Fadv);
void calc_flux_diffc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **T0,float **Fdiff);
void advance_timec(struct type_donneesc param,float dt,float **vol,float **Fadv,float **Fdiff,float **T0,float **T1);
void creation_A(struct type_donneesc param,int NA, float dt, float **x,float **y,float **xv,float **yv,float **vol,float **A);
void creation_B(struct type_donneesc param, int NA, float dt, float **x, float **y,float **xv,float **yv,float **vol, float **Fadv, float **T0, float *B);
void gaussij(int LV, float **A, float *B);
void miseajour_T(struct type_donneesc param,float **T0,float **T1,float *B);
void create_A_band(int MX,int NA,float **A, float **AB);
void SOR(int MX,int N,float **A,float *B,float R0,float W, float *X);

