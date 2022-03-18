/************************************************************************
   Program for simulating the wave-reflection analog model of subglottal 
   and supraglottal vocal tract, written by Isao Tokuda, 2 Dec. 2008.

   Physical Units are in [cm], [msec], [g] 

   Equation numbers refer to Ingo Titze's book of 
  "Myoelastic Arodynamic Theory of Phonation" (Titze, 2006)  

  Reference: I. T. Tokuda, M. Zemke, M. Kob, and H. Herzel, 
  ``Biomechanical modeling of register transitions and the role of vocal 
    tract resonators,'' Journal of the Acoustical Society of America, 
    Vol. 127, Issue 3, pp. 1528-1536 (2010)
  https://doi.org/10.1121/1.3299201
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define length       1.4000        /* length of folds in cm 
                                      (Standard Ishizaka-Flanagan: 1.4) */
#define cc              35.0       // speed of sound
#define rho          0.00113       // atmospheric density
#define P_DC         0.01000       // DC pressure in subglottis = 10 [cmH2O]

#define Num_Section_VT  44         // number of sections for vocal tract
#define Num_Section_SG  62         // number of sections for subglottis
#define dX             (17.5/44.0) // length of each section [cm]
#define dt             (dX/cc)     // integration time step [msec]

#define DIM          9

#define m_1_0        0.0090        // 1st mass [g]
#define m_2_0        0.0090        // 2nd mass [g]
#define m_3_0        0.0030        // 3rd mass [g]
#define m_b_0        0.0500        // body mass [g]
#define k_1b_0       0.0060        // Stiffness
#define k_2b_0       0.0060
#define k_3b_0       0.0020
#define k_b_0        0.0300
#define k_12_0       0.00100
#define k_23_0       0.00050
 
double m_1,m_2,m_3,m_b;
double k_1b,k_2b,k_3b,k_b,k_12,k_23;
double c_1,c_2,c_3;
double r_1b,r_2b,r_3b,r_b;
double X_[5];
double Force[5];
double XX[DIM];
double a_1;
double a_2;
double a_3;
double a_4;
double a_min;

double Ug;                         // volume velocity at glottis

double Attn_Lip;                   // attenuation factor at lip
double Area_Lip;                   // area of lip 
double Area_Input;                 // input-area from subglottis to vocal folds
double Area_VT[Num_Section_VT];    // area function of vocal tract
double Attn_VT[Num_Section_VT];    // attenuation factor of vocal tract
double Refl_VT[Num_Section_VT];    // reflection coefficient of vocal tract
double P_fw_VT[Num_Section_VT];    // foward pressure in vocal tract
double P_bk_VT[Num_Section_VT];    // backward pressure in vocal tract

double Area_SG[Num_Section_SG];    // area function of subglottis
double Attn_SG[Num_Section_SG];    // attenuation factor of subglottis
double Refl_SG[Num_Section_SG];    // reflection coefficient of subglottis
double P_fw_SG[Num_Section_SG];    // foward pressure in subglottis
double P_bk_SG[Num_Section_SG];    // backward pressure in subglottis

double P_fw_Lip;                   // foward pressure at lip
double P_bk_Lip;                   // backward pressure at lip
double p_fw_Lip_previous;          // previous state of foward pressure at lip
double P_s;                        // subglottal pressure
double P_t;                        // supraglottal pressure

int main()
{
  int    i,j,k;
  double Q;  //Q-parameter to control fundamental frequency
  void   ParameterSetting();
  void   Start();
  void   Update();
  FILE   *fp;
  
  Q = 3.0;
  ParameterSetting(Q); 
  Start(); 
  for(i=0;i<5000;i++) Update(); 

  fp=fopen("Trajectory.txt","w");
  for(j=0;j<200;j++){
    for(k=0;k<4;k++) Update();
    fprintf(fp,"%lf %lf\n",XX[8],Ug);
  }
  fclose(fp);

  return(0);
}

void Start()
{
  int k;
  /* Initialization: all pressure modulation to be zero */
  XX[0]=0.020;
  XX[2]=0.015;
  XX[4]=0.010;
  XX[6]=0.000;
  XX[1]=XX[3]=XX[5]=XX[7]=0.0;
  XX[8]=0.0;

  P_fw_Lip=P_bk_Lip=0.0;
  p_fw_Lip_previous=0.0;
  for(k=0;k<Num_Section_VT;k++) P_fw_VT[k]=P_bk_VT[k]=0.0;
  for(k=0;k<Num_Section_SG;k++) P_fw_SG[k]=P_bk_SG[k]=0.0;
}

void ParameterSetting(Q)
  double Q;
{
  int k;
  double x;
  double dratio_1=0.100;
  double dratio_2=0.400;
  double dratio_3=0.400;
  double dratio_b=0.400;
  FILE *fp;

  m_1 = m_1_0/Q;
  m_2 = m_2_0/Q;
  m_3 = m_3_0/Q;
  m_b = m_b_0/Q;
 
  k_1b = k_1b_0*Q;
  k_2b = k_2b_0*Q;
  k_3b = k_3b_0*Q;
  k_b  = k_b_0 *Q;
  k_12 = k_12_0;
  k_23 = k_23_0;

  c_1  = (3.0*k_1b);
  c_2  = (3.0*k_2b);
  c_3  = (3.0*k_3b);
 
  r_1b = (2.0*dratio_1*sqrt(m_1*k_1b));
  r_2b = (2.0*dratio_2*sqrt(m_2*k_2b));
  r_3b = (2.0*dratio_3*sqrt(m_3*k_3b));
  r_b  = (2.0*dratio_b*sqrt(m_b*k_b));
 
  X_[0] = 0.0000;
  X_[1] = 0.07500-0.0300;
  X_[2] = 0.22500-0.0300;
  X_[3] = 0.30000-0.0300;
  X_[4] = 0.32000-0.0300;

  /* Subglottal Area Function */ 
  fp = fopen("SubglottalAreaF.txt","r");
  for(k=0;k<62;k++) fscanf(fp,"%lf",&Area_SG[62-1-k]);
  fclose(fp);
  //for(k=0;k<Num_Section_SG;k++) Area_SG[k]=2.50; // [cm^2]

  for(k=0;k<Num_Section_SG;k++)
    Attn_SG[k]=1.0-0.0070*sqrt(M_PI/Area_SG[k])*dX; // Eq. (6.129)

  Area_Input=Area_SG[Num_Section_SG-1];

  for(k=0;k<(Num_Section_SG-1);k++)
    Refl_SG[k]=(Area_SG[k]-Area_SG[k+1])/(Area_SG[k]+Area_SG[k+1]); // Eq. (6.134)

  k=Num_Section_SG-1; // last section
  Refl_SG[k]=(Area_SG[k]-Area_Input)/(Area_SG[k]+Area_Input);       // Eq. (6.134)

  /* Vocal Tract Area Function */ 
  for(k=0;k<Num_Section_VT;k++) Area_VT[k]=3.00; 

  for(k=0;k<Num_Section_VT;k++)
    Attn_VT[k]=1.0-0.0070*sqrt(M_PI/Area_VT[k])*dX; // Eq. (6.129)

  Area_Lip=Area_VT[Num_Section_VT-1];
  Attn_Lip=1.0-0.0070*sqrt(M_PI/Area_Lip)*dX;       // Eq. (6.129)

  for(k=0;k<(Num_Section_VT-1);k++)
    Refl_VT[k]=(Area_VT[k]-Area_VT[k+1])/(Area_VT[k]+Area_VT[k+1]); // Eq. (6.134)

  k=Num_Section_VT-1; // last section
  Refl_VT[k]=(Area_VT[k]-Area_Lip)/(Area_VT[k]+Area_Lip);           // Eq. (6.134)
}


void Update()
{
  int    i,k;
  double kx[DIM],J();
  double R,I;                          // radiation resistance and inertance at lip
  double b1,b2,c1,c2;                  // circuit parameters for lip
  double Area_g;                       // glottal opening area
  double P_fw_VT_next[Num_Section_VT]; // next state of foward pressure in vocal tract
  double P_bk_VT_next[Num_Section_VT]; // next state of backward pressure in vocal tract
  double P_fw_SG_next[Num_Section_SG]; // next state of foward pressure in subglottis
  double P_bk_SG_next[Num_Section_SG]; // next state of backward pressure in subglottis
  double P_fw_Lip_next;                // next state of forward pressure at lip
  double P_bk_Lip_next;                // next state of backward pressure at lip
  double Aast;                         // effective vocal tract area
  double kt;                           // transglottal pressure coefficient set to be unity
  double p_fw_Lip;                     // incident wave at lip (not the same as P_fw_Lip)
  double P_t_bk,P_t_fw,P_s_bk,P_s_fw;
  void   Params();

  Params();

  kt = 1.0;
  R = 128.0/(9.0*pow(M_PI,2.0)); // Eq. (6.163)
  I = (2.0/dt)*8.0*sqrt(Area_Lip)/(3.0*pow(M_PI,1.50)*cc); // Eq. (6.165)
  c1 = -R+I-R*I; // Eq. (6.168)
  c2 = -R-I+R*I; // Eq. (6.169)
  b1 = -R+I+R*I; // Eq. (6.170)
  b2 =  R+I+R*I; // Eq. (6.171)

  P_s_fw = Attn_SG[Num_Section_SG-1]*P_fw_SG[Num_Section_SG-1]; // forward pressure from subglottis to glottis
  P_t_bk = Attn_VT[0]*P_bk_VT[0];                               // backward pressure from supraglottis to glottis

  /* Glottal area vibrates between 0 - 0.1 [cm*cm] with 120 Hz */ 
  Area_g = 0.10*0.50*(1.0+cos(2.0*M_PI*120.0*XX[8]/1000.0));
  Area_g = a_min;
  if(Area_g<0.0001) Area_g = 0.0001;

  Aast=(Area_SG[Num_Section_SG-1]*Area_VT[0])/(Area_SG[Num_Section_SG-1]+Area_VT[0]); // Eq. (6.185)

  /* source-tract coupling according to Eq. (6.184) */ 
  if((2.0*P_s_fw+P_DC) > (2.0*P_t_bk)){
    Ug=(Area_g*cc/kt)*(-(Area_g/Aast)+sqrt(pow(Area_g/Aast,2.0)+(2.0*kt/(cc*cc*rho))*(2.0*P_s_fw+P_DC-2.0*P_t_bk)));
  }
  else{
    Ug=(Area_g*cc/kt)*(-(Area_g/Aast)-sqrt(pow(Area_g/Aast,2.0)+(2.0*kt/(cc*cc*rho))*(2.0*P_s_fw+P_DC-2.0*P_t_bk)));
  }

  P_s_bk = P_s_fw-(rho*cc*Ug/Area_SG[Num_Section_SG-1]); // backward pressure from glottis to subglottis
  P_t_fw = P_t_bk+(rho*cc*Ug/Area_VT[0]);                // forward pressure from glottis to supraglottis
  P_s    = P_DC+P_s_fw+P_s_bk; // subglottal pressure
  P_t    = P_t_fw+P_t_bk;      // supraglottal pressure

  for(i=0;i<8;i++){
    Params();
    for(k=0;k<DIM;k++) kx[k]=(dt/8.0)*J(k);
    for(k=0;k<DIM;k++) XX[k]+=kx[k];
  }

  P_fw_SG_next[0] = 0.0;    // zero input from lung

  /* Computation of next state according to Eqs. (6.139) and (6.140) */ 
  for(k=0;k<(Num_Section_SG-1);k++){
    P_bk_SG_next[k]  =(1.0-Refl_SG[k])*Attn_SG[k+1]*P_bk_SG[k+1]+Refl_SG[k]*Attn_SG[k]*P_fw_SG[k];
    P_fw_SG_next[k+1]=(1.0+Refl_SG[k])*Attn_SG[k]  *P_fw_SG[k]  -Refl_SG[k]*Attn_SG[k+1]*P_bk_SG[k+1];
  }

  P_bk_SG_next[Num_Section_SG-1]=P_s_bk; // backward pressure from glottis

  P_fw_VT_next[0]=P_t_fw; // forward pressure from glottis

  /* Computation of next state according to Eqs. (6.139) and (6.140) */ 
  for(k=0;k<(Num_Section_VT-1);k++){
    P_bk_VT_next[k]  =(1.0-Refl_VT[k])*Attn_VT[k+1]*P_bk_VT[k+1]+Refl_VT[k]*Attn_VT[k]*P_fw_VT[k];
    P_fw_VT_next[k+1]=(1.0+Refl_VT[k])*Attn_VT[k]  *P_fw_VT[k]  -Refl_VT[k]*Attn_VT[k+1]*P_bk_VT[k+1];
  }

  k=Num_Section_VT-1; // last section
  P_bk_VT_next[k]=(1.0-Refl_VT[k])*Attn_Lip*P_bk_Lip+Refl_VT[k]*Attn_VT[k]*P_fw_VT[k];

  p_fw_Lip = Attn_VT[Num_Section_VT-1]*P_fw_VT[Num_Section_VT-1];              // Eq. (6.178)
  P_bk_Lip_next = (c1*p_fw_Lip_previous+c2*p_fw_Lip+b1*P_bk_Lip)/b2;           // Eq. (6.176)
  P_fw_Lip_next = ((c1-b1)*p_fw_Lip_previous+(c2+b2)*p_fw_Lip+b1*P_fw_Lip)/b2; // Eq. (6.177)

  /* Update of state */
  for(k=0;k<Num_Section_SG;k++){
    P_bk_SG[k]=P_bk_SG_next[k]; P_fw_SG[k]=P_fw_SG_next[k];
  }
  for(k=0;k<Num_Section_VT;k++){
    P_fw_VT[k]=P_fw_VT_next[k]; P_bk_VT[k]=P_bk_VT_next[k];
  }

  P_bk_Lip = P_bk_Lip_next; 
  P_fw_Lip = P_fw_Lip_next;
  p_fw_Lip_previous = p_fw_Lip;
}


double J(j)
{
  double y;
  double Theta();
  if     (j==0) y = XX[1];
  else if(j==1)
    y = (1.0/m_1)*(Force[1]*length-r_1b*(XX[1]-XX[7])-
                   k_1b*(XX[0]-XX[6])-
                   Theta(-a_1)*c_1*(a_1/(2.0*length))-
                   k_12*(XX[0]-XX[2]));
  else if(j==2) y = XX[3];
  else if(j==3)
    y = (1.0/m_2)*(Force[2]*length-r_2b*(XX[3]-XX[7])-
                   k_2b*(XX[2]-XX[6])-
                   Theta(-a_2)*c_2*(a_2/(2.0*length))-
                   k_12*(XX[2]-XX[0])-k_23*(XX[2]-XX[4]));

  else if(j==4) y = XX[5];
  else if(j==5)
    y = (1.0/m_3)*(Force[3]*length-r_3b*(XX[5]-XX[7])-
                   k_3b*(XX[4]-XX[6])-
                   Theta(-a_3)*c_3*(a_3/(2.0*length))-
                   k_23*(XX[4]-XX[2]));

  else if(j==6) y = XX[7];
  else if(j==7)
    y = (1.0/m_b)*(-r_b*XX[7]-
                   k_b *XX[6]-
                   k_1b*(XX[6]-XX[0])-
                   k_2b*(XX[6]-XX[2])-
                   k_3b*(XX[6]-XX[4])-
                   r_1b*(XX[7]-XX[1])-
                   r_2b*(XX[7]-XX[3])-
                   r_3b*(XX[7]-XX[5]));

  else if(j==8) y = 1.0;

  return (y);
}

void Params()
{
  int    i,min;
  double Theta();
  double P0,P1,H1,H2,Integral1,Integral2;
  double d,dc;
  double integralmax,hmax;
  double x01,x12,x23,x34;
  double x,y;
  double h[5];
  double Fr[5];
  double Fl[5];
  double a_0;
  double a_01,a_02,a_03;
  double Boundary_Layer;

  a_0 = length*2.0*1.800;
  a_4 = length*2.0*1.800;

  a_01 = length*2.0*0.0180;
  a_02 = length*2.0*0.0180;
  a_03 = length*2.0*0.0180;
 
  a_1 = a_01+2.0*length*XX[0];
  a_2 = a_02+2.0*length*XX[2];
  a_3 = a_03+2.0*length*XX[4];

  if     (a_1<=a_2 && a_1<=a_3){ a_min=a_1; min=1; }
  else if(a_2<=a_1 && a_2<=a_3){ a_min=a_2; min=2; }
  else                         { a_min=a_3; min=3; }

  h[0]=a_0/length;
  h[1]=a_1/length;
  h[2]=a_2/length;
  h[3]=a_3/length;
  h[4]=a_4/length;

  Boundary_Layer = 2.0*0.00500*length;

  if(0 < a_min){
    if(a_min > Boundary_Layer){
      for(i=1;i<5;i++){

        H1 = (h[i]-h[i-1])/(X_[i]-X_[i-1]);
        H2 = h[i-1]-H1*X_[i-1];

        Integral1 = -1.0*((1.0/h[i])-(1.0/h[i-1]))/H1;
        Integral1 = (X_[i]-X_[i-1])/(h[i]*h[i-1]);
        Integral2 = (log(h[i]/h[i-1])+H2*((1.0/h[i])-(1.0/h[i-1])))/pow(H1,2.0);
        P0 = P_s;
        P1 = P_s*pow(a_min/length,2.0);

        Fr[i-1]=(X_[i]/(X_[i]-X_[i-1]))*((P0*(X_[i]-X_[i-1]))-(P1*Integral1))
               -(1.0/(X_[i]-X_[i-1]))*(
                ((P0/2.0)*(pow(X_[i],2.0)-pow(X_[i-1],2.0)))-(P1*Integral2));

        x=Fl[i]=(1.0/(X_[i]-X_[i-1]))*(
              ((P0/2.0)*(pow(X_[i],2.0)-pow(X_[i-1],2.0)))-(P1*Integral2))
             -(X_[i-1]/(X_[i]-X_[i-1]))*((P0*(X_[i]-X_[i-1]))-(P1*Integral1));

        y=Fl[i]=P0*(X_[i]-X_[i-1])-P1*Integral1-Fr[i-1];

        hmax = h[i-1];
        if(hmax < h[i]) hmax = h[i];
        integralmax = (X_[i]-X_[i-1])*(P0-(P1/pow(hmax,2.0)));

        if((Fr[i-1]+Fl[i]) > 1.0*integralmax){ 
          Fr[i-1] = Fl[i] = 0.5*integralmax;
        }

        if(min < i) 
          Fl[i]=Fr[i-1]=0.50*P_t*(X_[i]-X_[i-1]);
      }
    }
    else{
      if(a_1 > Boundary_Layer && a_2 > Boundary_Layer && a_3 <= Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.50*P_s*(X_[2]-X_[1]);
        Fl[3]=Fr[2]=0.50*P_s*(X_[3]-X_[2]);
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
        Fl[4]=Fr[3]=0.00;
      }
      else if(a_1 > Boundary_Layer && a_2 <= Boundary_Layer && a_3 > Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.50*P_s*(X_[2]-X_[1]);
        Fl[3]=Fr[2]=0.50*P_t*(X_[3]-X_[2]);
        Fl[3]=Fr[2]=0.00;
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
      }
      else if(a_1 <= Boundary_Layer && a_2 > Boundary_Layer && a_3 > Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.50*P_t*(X_[2]-X_[1]);
        Fl[2]=Fr[1]=0.00;
        Fl[3]=Fr[2]=0.50*P_t*(X_[3]-X_[2]);
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
      }
      else if(a_1 > Boundary_Layer && a_2 <= Boundary_Layer && a_3 <= Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.50*P_s*(X_[2]-X_[1]);
        Fl[3]=Fr[2]=0.00;
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
        Fl[4]=Fr[3]=0.00;
      }
      else if(a_1 <= Boundary_Layer && a_2 <= Boundary_Layer && a_3 > Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.00;
        Fl[3]=Fr[2]=0.50*P_t*(X_[3]-X_[2]);
        Fl[3]=Fr[2]=0.0;
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
      }
      else if(a_1 <= Boundary_Layer && a_2 > Boundary_Layer && a_3 <= Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.00;
        Fl[3]=Fr[2]=0.00;
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
        Fl[4]=Fr[3]=0.00;

      }
      else if(a_1 <= Boundary_Layer && a_2 <= Boundary_Layer && a_3 <= Boundary_Layer){ 
        Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
        Fl[2]=Fr[1]=0.00;
        Fl[3]=Fr[2]=0.00;
        Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
        Fl[4]=Fr[3]=0.00;
      }
    }
  }
  else{
    if(0 < a_1 && 0 < a_2 && a_3 <= 0){ 
      Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);
      Fl[2]=Fr[1]=0.50*P_s*(X_[2]-X_[1]);

      x23 = X_[2]+((X_[3]-X_[2])*h[2]/(h[2]-h[3]));
      if(x23 < X_[2]) x23 = X_[2];
      if(x23 > X_[3]) x23 = X_[3];
      Fl[3]=0.50*P_s*pow(x23-X_[2],2.0)/(X_[3]-X_[2]);
      Fr[2]=P_s*(x23-X_[2])-Fl[3];

      x34 = X_[3]+((X_[4]-X_[3])*(-h[3])/(h[4]-h[3]));
      if(x34 < X_[3]) x34 = X_[3];
      if(x34 > X_[4]) x34 = X_[4];
      Fr[3]=0.50*P_t*pow(X_[4]-x34,2.0)/(X_[4]-X_[3]);
      Fl[4]=P_t*(X_[4]-x34)-Fr[3];
    }
    else if(0 < a_1 && a_2 <= 0 && a_3 <= 0){ 
      Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);

      x12 = X_[1]+((X_[2]-X_[1])*h[1]/(h[1]-h[2]));
      if(x12 < X_[1]) x12 = X_[1];
      if(x12 > X_[2]) x12 = X_[2];
      Fl[2]=0.50*P_s*pow(x12-X_[1],2.0)/(X_[2]-X_[1]);
      Fr[1]=P_s*(x12-X_[1])-Fl[2];

      Fl[3]=Fr[2]=0.0;

      x34 = X_[3]+((X_[4]-X_[3])*(-h[3])/(h[4]-h[3]));
      if(x34 < X_[3]) x34 = X_[3];
      if(x34 > X_[4]) x34 = X_[4];
      Fr[3]=0.50*P_t*pow(X_[4]-x34,2.0)/(X_[4]-X_[3]);
      Fl[4]=P_t*(X_[4]-x34)-Fr[3];
    }
    else if(a_1 <= 0 && a_3 <= 0){ 
      x01 = X_[0]+((X_[1]-X_[0])*h[0]/(h[0]-h[1]));
      if(x01 < X_[0]) x01 = X_[0];
      if(x01 > X_[1]) x01 = X_[1];
      Fl[1]=0.50*P_s*pow(x01-X_[0],2.0)/(X_[1]-X_[0]);
      Fr[0]=P_s*(x01-X_[0])-Fl[1];

      Fl[2]=Fr[1]=0.0;
      Fl[3]=Fr[2]=0.0;

      x34 = X_[3]+((X_[4]-X_[3])*(-h[3])/(h[4]-h[3]));
      if(x34 < X_[3]) x34 = X_[3];
      if(x34 > X_[4]) x34 = X_[4];
      Fr[3]=0.50*P_t*pow(X_[4]-x34,2.0)/(X_[4]-X_[3]);
      Fl[4]=P_t*(X_[4]-x34)-Fr[3];
    }
    else if(a_1 <= 0 && 0 < a_2 && 0 < a_3){ 
      x01 = X_[0]+((X_[1]-X_[0])*h[0]/(h[0]-h[1]));
      if(x01 < X_[0]) x01 = X_[0];
      if(x01 > X_[1]) x01 = X_[1];
      Fl[1]=0.50*P_s*pow(x01-X_[0],2.0)/(X_[1]-X_[0]);
      Fr[0]=P_s*(x01-X_[0])-Fl[1];

      x12 = X_[1]+((X_[2]-X_[1])*(-h[1])/(h[2]-h[1]));
      if(x12 < X_[1]) x12 = X_[1];
      if(x12 > X_[2]) x12 = X_[2];
      Fr[1]=0.50*P_t*pow(X_[2]-x12,2.0)/(X_[2]-X_[1]);
      Fl[2]=P_t*(X_[2]-x12)-Fr[1];

      Fl[3]=Fr[2]=0.50*P_t*(X_[3]-X_[2]);
      Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
    }
    else if(a_1 <= 0 && a_2 <= 0 && 0 < a_3){ 
      x01 = X_[0]+((X_[1]-X_[0])*h[0]/(h[0]-h[1]));
      if(x01 < X_[0]) x01 = X_[0];
      if(x01 > X_[1]) x01 = X_[1];
      Fl[1]=0.50*P_s*pow(x01-X_[0],2.0)/(X_[1]-X_[0]);
      Fr[0]=P_s*(x01-X_[0])-Fl[1];

      Fl[2]=Fr[1]=0.0;

      x23 = X_[2]+((X_[3]-X_[2])*(-h[2])/(h[3]-h[2]));
      if(x23 < X_[2]) x23 = X_[2];
      if(x23 > X_[3]) x23 = X_[3];
      Fr[2]=0.50*P_t*pow(X_[3]-x23,2.0)/(X_[3]-X_[2]);
      Fl[3]=P_t*(X_[3]-x23)-Fr[2];

      Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
    }
    else if(0 < a_1 && a_2 <= 0 && 0 < a_3){ 
      Fl[1]=Fr[0]=0.50*P_s*(X_[1]-X_[0]);

      x12 = X_[1]+((X_[2]-X_[1])*h[1]/(h[1]-h[2]));
      if(x12 < X_[1]) x12 = X_[1];
      if(x12 > X_[2]) x12 = X_[2];
      Fl[2]=0.50*P_s*pow(x12-X_[1],2.0)/(X_[2]-X_[1]);
      Fr[1]=P_s*(x12-X_[1])-Fl[2];

      x23 = X_[2]+((X_[3]-X_[2])*(-h[2])/(h[3]-h[2]));
      if(x23 < X_[2]) x23 = X_[2];
      if(x23 > X_[3]) x23 = X_[3];
      Fr[2]=0.50*P_t*pow(X_[3]-x23,2.0)/(X_[3]-X_[2]);
      Fl[3]=P_t*(X_[3]-x23)-Fr[2];

      Fl[4]=Fr[3]=0.50*P_t*(X_[4]-X_[3]);
    }
  }

  for(i=1;i<4;i++) Force[i]=Fr[i]+Fl[i];
}

double Theta(x)
  double x;
{
  if (x < 0.0) return(0.0);
  else         return(1.0);
}
