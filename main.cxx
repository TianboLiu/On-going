#include <omp.h>/*{{{*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "TSystem.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"
#include "Math/GSLIntegrator.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_gamma.h"
#include <Math/SpecFuncMathCore.h>
#include "TLatex.h"
#include "hoppet_v1.h"
//#include "/work/halla/solid/yez/TMDs/hoppet-1.1.5/include/hoppet_v1.h"
#include <ctime>
#include "TStyle.h" 
#define Log log
#define Power pow
using namespace std;/*}}}*/

/*Extern-C and predefine{{{*/
extern "C" {/*{{{*/
    void setct10_(int * iset);
}

extern "C" {
    double ct10pdf_(int * iparton, double * x, double * q);
}/*}}}*/

double pdf ( int  iparton , double  x , double q2  ) /*{{{*/
{ int Iset = 100; 
    setct10_(&Iset);
    double q = sqrt(q2);
    return ct10pdf_(&iparton,&x,&q);
} /*}}}*/

const int num_x=1000;
double list_x[num_x];
double list_pdf[5][num_x];

typedef struct { int FINI; } mycommonF77;

extern mycommonF77 fraginid_;

extern "C" {/*{{{*/
    void fdss_(int*IH,int*IC,int*IO,double*X,double*Q2,double*U,double*UB,
            double*D,double*DB,double*S,double*SB,double*C,double*B,double*GL);
}/*}}}*/

double ffu (int hadron,int charge ,double x,double q2) /*{{{*/
{
    int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return U/x;
} /*}}}*/

double ffub (int hadron,int charge ,double x,double q2) /*{{{*/
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return UB/x;
} 

double ffd (int hadron,int charge ,double x,double q2) 
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return D/x;
} 

double ffdb (int hadron,int charge ,double x,double q2) 
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return DB/x;
} 

double ffs (int hadron,int charge ,double x,double q2) 
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return S/x;
} 

double ffsb (int hadron,int charge ,double x,double q2) 
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return SB/x;
} 

double ffg (int hadron,int charge ,double x,double q2) 
{int IO=1;
    double U, UB, D, DB, S, SB, C, B, GL;
    fdss_(&hadron,&charge,&IO,&x,&q2,&U,&UB,&D,&DB,&S,&SB,&C,&B,&GL);
    return GL/x;
} /*}}}*/

extern "C" {/*{{{*/
    void dssvfit_(double*X,double*Q2,double*DUV,double*DDV,double*DUBar,double*DDBar,double*DSTR,double*DGLU);
}

extern "C" {
    void dssvini_(int * iset);
}/*}}}*/

double gpu (double x,double q2) /*{{{*/
{int Iset = 0; 
    dssvini_(&Iset);

    double U, UB, D, DB, S, GL;
    dssvfit_(&x,&q2,&U,&D,&UB,&DB,&S,&GL);
    return (U+UB)/x;
} 

double gpd (double x,double q2) 
{int Iset = 0; 
    dssvini_(&Iset);

    double U, UB, D, DB, S, GL;
    dssvfit_(&x,&q2,&U,&D,&UB,&DB,&S,&GL);
    return (D+DB)/x;
} /*}}}*/

double list_ff[5][num_x];

void list_ff_f()/*{{{*/
{
    for( int j=0 ;j< num_x ;++j)
    {
        fraginid_.FINI=0;
        list_ff[1+1][j]=ffu (1,1 ,list_x[j],2.4);
        fraginid_.FINI=0;
        list_ff[1-1][j]=ffu (1,-1 ,list_x[j],2.4);
        fraginid_.FINI=0;
        list_ff[1+3][j]=ffs (1, 1 ,list_x[j],2.4);

    }
}/*}}}*/

void list_pdf_f()/*{{{*/
{


    for (int i=0; i<num_x ;++i )
    {
        list_x[i]=1.0/pow(1.02,num_x-1)*pow(1.02,i);
        //    cout<<list_x[i]<<endl;
    }

    list_x[num_x-1]=1; 

    double value_u;
    double value_d;


    for( int j=0 ;j< num_x ;++j)
    {

        if(list_x[j]<0.00001   )
        { value_u=gpu(0.00001,2.4);  
            value_d=gpd(0.00001,2.4);
        }else{
            value_u=gpu(list_x[j],2.4 );
            value_d=gpd(list_x[j],2.4 );

        }

        list_pdf[2+1][j]=0.5*(pdf(1,list_x[j],2.4)+value_u);
        list_pdf[2-1][j]=pdf(-1,list_x[j],2.4);
        list_pdf[2+2][j]=0.5*(pdf(2,list_x[j],2.4)+value_d);
        list_pdf[2-2][j]=pdf(-2,list_x[j],2.4);    

    }
}/*}}}*/

double pdf_n( int  iparton , double  x   ) /*{{{*/
{
    ROOT::Math::Interpolator intera_x (num_x,ROOT::Math::Interpolation::kAKIMA);
    //  ROOT::Math::Interpolator intera_x (num_x,ROOT::Math::Interpolation::kAKIMA);
    intera_x.SetData(num_x,list_x,list_pdf[2+iparton]);

    return intera_x.Eval(x) ;


}/*}}}*/

double ff_nn( int  iparton , double  x   ) /*{{{*/
{
    ROOT::Math::Interpolator intera_y (num_x,ROOT::Math::Interpolation::kAKIMA);
    //  ROOT::Math::Interpolator intera_x (num_x,ROOT::Math::Interpolation::kAKIMA);
    intera_y.SetData(num_x,list_x,list_ff[1+iparton]);

    return intera_y.Eval(x) ;
}/*}}}*/

const double pi= TMath::Pi(); 
double S;
const double compass=160.0*2*0.94;
const double hermas=27.6*2*0.94;
const double alf=1.0/137;
const double C1=2.0*exp(-TMath::EulerGamma());
//const double nf=4.0;
const double CA=3.0;
const double CF=4.0/3;
//const double Bb0=(11.0-(2.0/3)*nf);
//const double Mp=0.93;
const double A1=CF;
const double B1=-(3.0/2)*CF;
//const double Q02=2.4;
const double bm=1.5;
const double LQCD2_5=0.225*0.225;
double b_l=0.00001;
double b_u=10; 

double as ( double X )/*{{{*/
{ 
    double  Bb0f_5=  (11.0-(2.0/3.0)*5.0);
    double  Bb1f_5=  102.0-38.0*5.0/3.0;
    double LQCD2f_5=LQCD2_5;
    // double X2=X;
    //   return 4*pi*(1/(Bb0f(X2)*Log(X/LQCD2f(X2))) - (Bb1f(X2)*Log(Log(X/LQCD2f(X2))))/(Power(Bb0f(X2),3)*Power(Log(X/LQCD2f(X2)),2)));
    return 4*pi*(1/(Bb0f_5*Log(X/LQCD2f_5)) - (Bb1f_5*Log(Log(X/LQCD2f_5)))/(Power(Bb0f_5,3)*Power(Log(X/LQCD2f_5),2)));
}/*}}}*/

const double re3=1.202;
const double Nc=3.0;

double Sudakovkernel( double mu2 , double Q2  ,double LQCD2,double nf   )/*{{{*/
{
    double Bb0=(11.0-(2.0/3)*nf);
    double Bb1=102.0-38.0*nf/3;
    double A1=CF;
    double A2= CF*((67.0/36.0-pi*pi/12.0)*CA-5.0/18.0*nf);
    double B1=-(3.0/2.0)*CF;
    // double B2=CF*CF*(pi*pi/4.0-3.0/16.0-3*re3)+CF*Nc*(11.0/36.0*pi*pi-193.0/48.0+3.0/2.0*re3)+CF*nf*(-1.0/18.0*pi*pi+17.0/24.0);
    double B2=0;
    return
        (2*(54*A1*Power(Bb0,5)*Power(Log(mu2/LQCD2),4)*(-1 + Log(Log(mu2/LQCD2))) - 
            8*Power(Bb1,2)*(B2 + A2*Log(Q2/mu2))*(2 + 6*Log(Log(mu2/LQCD2)) + 9*Power(Log(Log(mu2/LQCD2)),2)) + 
            27*Power(Bb0,3)*Power(Log(mu2/LQCD2),3)*(-8*A2*Bb0 + 2*A1*Bb1 + 2*(-4*A2*Bb0 + B1*Power(Bb0,2) + A1*Bb1
                    + A1*Power(Bb0,2)*Log(Q2/mu2))*Log(Log(mu2/LQCD2)) + 
                A1*Bb1*Power(Log(Log(mu2/LQCD2)),2)) + 2*Bb1*Log(mu2/LQCD2)*
            (54*B2*Power(Bb0,2) + 19*A2*Bb1 + 6*(18*B2*Power(Bb0,2) + 5*A2*Bb1)*Log(Log(mu2/LQCD2)) + 18*A2*Bb1*Power(Log(Log(mu2/LQCD2)),2) + 
             54*A2*Power(Bb0,2)*Log(Q2/mu2)*(1 + 2*Log(Log(mu2/LQCD2)))) - 
            54*Power(Bb0,2)*Power(Log(mu2/LQCD2),2)*(4*B2*Power(Bb0,2) + 6*A2*Bb1 - B1*Bb0*Bb1 + (4*A2 - B1*Bb0)*Bb1*Log(Log(mu2/LQCD2)) + 
                Bb0*Log(Q2/mu2)*(4*A2*Bb0 - A1*Bb1
                    - A1*Bb1*Log(Log(mu2/LQCD2))))))/(27.*Power(Bb0,6)*Power(Log(mu2/LQCD2),3));
}/*}}}*/

double Sudakov1( double Q2 , double b )/*{{{*/
{
    double bs2= b*b/(1+b*b/(bm*bm) )/C1/C1;
    double upL1= Sudakovkernel( Q2 ,  Q2  , LQCD2_5 ,5.0  );
    double downL1= Sudakovkernel( 1/bs2 , Q2 ,   LQCD2_5 ,5.0    );
    return exp(-(upL1-downL1));
}/*}}}*/

double Sudakov2( double Q2 ,double b )/*{{{*/
{
    double g2=0.842;

    double Q12=2.4;

    double bs2= b*b/(1+b*b/(bm*bm) );

    return exp(-(g2*log(Q2/Q12)/2)*log(b*b/bs2)/2 );
}/*}}}*/

const  double BNLYd= 0.042;
const double euq2=4.0/9;
const double edq2=1.0/9;
const double CB=0.3894E9;
double phtit,ptit;
double BNLYit;
double Q2it;
double zit;
double xit;

#pragma omp threadprivate(zit)
#pragma omp threadprivate(xit)
#pragma omp threadprivate(ptit)
#pragma omp threadprivate(phtit)
#pragma omp threadprivate(Q2it)

//Parameters in transversity and Collins
double Nuit, Ndit, Nuit_t, Nubit_t ,Ndit_t  ,Ndbit_t ,
       auit, adit,   auit_t, aubit_t,  adit_t, adbit_t ,
       buit, bdit,  buit_t, bdit_t, bubit_t, bdbit_t   ;

void  heralhc_init(const double & x,/*{{{*/
        const double & Q,
        double * pdf_a) { 

    // Nuit=       1.00000e+01  ;
    // Ndit=      -1.53931e+00  ;
    // auit=       7.96969e+00  ;
    // adit=       1.43211e+00  ;
    // buit=       1.18983e+00  ;
    // bdit=       7.20196e-09  ;

    double  uv  =x* Nuit * pow(x,auit) * pow(1-x,buit)   * ff_nn(1 ,x) ;
    double  ub  =x* Ndit * pow(x,adit) * pow(1-x,bdit)   * ff_nn(-1,x) ;
    double  sv  =x* Ndit * pow(x,adit) * pow(1-x,bdit)   * ff_nn(3,x) ;

    // double  uv  =x* Nuit * pow(x,auit) * pow(1-x,buit)   * ff_nn (1 ,x) ;
    // double  ub  =x* Ndit * pow(x,adit) * pow(1-x,bdit)   * ff_nn (-1 ,x) ;
    // double  dv  =x* Ndit * pow(x,adit) * pow(1-x,bdit)   * ffd (1,1 ,x,2.4) ;
    // double  db  =x* Nuit * pow(x,auit) * pow(1-x,buit)   * ffd (1,-1 ,x,2.4) ;

    double  uv_t  =pow(auit_t+buit_t,auit_t+buit_t)/pow(auit_t,auit_t)/pow(buit_t,buit_t)*x* Nuit_t * pow(x,auit_t) * pow(1-x,buit_t)   * pdf_n (1 ,x) ;
    // double  ub_t  =x* Nubit_t * pow(x,aubit_t) * pow(1-x,bubit_t)   * pdf_n (-1 ,x) ;
    double  dv_t  =pow(adit_t+bdit_t,adit_t+bdit_t)/pow(adit_t,adit_t)/pow(bdit_t,bdit_t)*x* Ndit_t * pow(x,adit_t) * pow(1-x,bdit_t)   * pdf_n (2 ,x) ;
    // double  db_t  =x* Ndbit_t * pow(x,adbit_t) * pow(1-x,bdbit_t)   * pdf_n (-2 ,x) ;

    // double  uv_t  =0 ;
    // double  ub_t  =0 ;
    // double  dv_t  =0 ;
    // double  db_t  =0 ;

    pdf_a[ 0+6] = 0;
    pdf_a[2+6] = uv_t;
    pdf_a[-2+6] = dv_t;
    pdf_a[ 3+6] = sv;
    pdf_a[-3+6] = 0;
    pdf_a[ 1+6] = uv;
    pdf_a[-1+6] = ub;
    pdf_a[ 4+6] = 0;
    pdf_a[-4+6] = 0;
    pdf_a[ 5+6] = 0;
    pdf_a[ 6+6] = 0; 
    pdf_a[-5+6] = 0;
    pdf_a[-6+6] = 0;
}/*}}}*/

double pof_s_pi(  double x,  double Q2  )/*{{{*/
{
    double Q= sqrt(Q2);
    double pdf_a[13];
    hoppetEval(x, Q, pdf_a);
    return pdf_a[6+3]/x;
}/*}}}*/

double pof_u_pi(  double x,  double Q2  )/*{{{*/
{
    double Q= sqrt(Q2);
    double pdf_a[13];
    hoppetEval(x, Q, pdf_a);
    return pdf_a[6+1]/x;
}/*}}}*/

double pof_ub_pi(  double x,  double Q2  )/*{{{*/
{
    double Q= sqrt(Q2);
    double pdf_a[13];
    hoppetEval(x, Q, pdf_a);
    return pdf_a[6-1]/x;
}/*}}}*/

double pof_u_mpi( double x,  double Q2  )/*{{{*/
{ 
    return pof_ub_pi(  x,   Q2  );
}/*}}}*/

double pof_ub_mpi(  double x,  double Q2  )/*{{{*/
{
    return pof_u_pi(   x,   Q2  );
}/*}}}*/

double pof_d_pi(   double x,  double Q2  )/*{{{*/
{
    return  pof_ub_pi(   x,   Q2  );
}/*}}}*/

double pof_db_pi(   double x,  double Q2  )/*{{{*/
{
    return  pof_u_pi(    x,   Q2  );
}/*}}}*/

double pof_d_mpi(   double x,  double Q2  )/*{{{*/
{
    return  pof_ub_mpi(   x,   Q2  );
}/*}}}*/

double pof_db_mpi(   double x,  double Q2  )/*{{{*/
{
    return  pof_u_mpi(   x,   Q2  );
}/*}}}*/

double Sit;
const  double S_compass=160*0.94*2;
const  double S_hermes=27.6*0.94*2;
const  double S_Jlab=5.9*0.94*2;
/*}}}*/

////////////////////////////////////////////////////////////////////////
////////////////////// BEGINNING /////////////////*{{{*/
//const  double S_Jlab12=11.*0.94*2;

double Singma0( double Q2 , double Si ,double x )/*{{{*/
{
    return (2*pi*alf*alf)*(1+(1-Q2/(x*Si))*(1-Q2/(x*Si)))/(Q2*Q2);
}/*}}}*/

double Singma0_collins( double Q2 , double Si ,double x )/*{{{*/
{

    double y=Q2/(Si*x);

    double val= (2*pi*alf*alf)*(1+(1-Q2/(x*Si))*(1-Q2/(x*Si)))/(Q2*Q2);

    double cof=(2*(1-y)/(1+(1-y)*(1-y)));

    return  val*cof;
}/*}}}*/

const  double BNLYf= 0.212/2.0;
const  double  g3=0;

double Sudakov3unpoDIS( double b , double x, double z )/*{{{*/
{
    return exp(-BNLYf*b*b-g3*pow(0.01/x,0.2)-BNLYd*b*b/(z*z));
}/*}}}*/

double Sudakov3poDIS(double adterm , double b , double x, double z )/*{{{*/
{
    return exp(-BNLYf*b*b-g3*pow(0.01/x,0.2)-BNLYd*b*b/(z*z)+adterm*b*b/(z*z));
}/*}}}*/

double unpopdf_dis_N(double x1, double  x2 , double xi1a, double xi2a ,double b ,int fh  ,int ch   )/*{{{*/
{ 
    //fraginid_.FINI=0;
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(2,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi*pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffu(fh,ch,x2,mu2)/(1-x2) 
             +as(mu2)/pi*ffu(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 

            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-2,x1,mu2)/(1-x1)+as(mu2)/pi* pdf(-2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffub(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffub(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(1,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffd(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffd(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-1,x1,mu2)/(1-x1)
             +as(mu2)/pi* pdf(-1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffdb(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffdb(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    return pdfuf+pdfdf;
}/*}}}*/

double unpopdf_dis_P(double x1, double  x2 , double xi1a, double xi2a ,double b ,int fh  ,int ch   )/*{{{*/
{ 
    //fraginid_.FINI=0;
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(1,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi*pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffu(fh,ch,x2,mu2)/(1-x2) 
             +as(mu2)/pi*ffu(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 

            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-1,x1,mu2)/(1-x1)+as(mu2)/pi* pdf(-1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffub(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffub(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(2,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffd(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffd(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-2,x1,mu2)/(1-x1)
             +as(mu2)/pi* pdf(-2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffdb(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffdb(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    return pdfuf+pdfdf;
}/*}}}*/

double unpopdf_dis_D(double x1, double  x2 , double xi1a, double xi2a ,double b ,int fh  ,int ch   )/*{{{*/
{ 
    //fraginid_.FINI=0;
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf1=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(1,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffu(fh,ch,x2,mu2)/(1-x2) 
             +as(mu2)/pi*ffu(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 

            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-1,x1,mu2)/(1-x1)+as(mu2)/pi* pdf(-1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffub(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffub(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    double pdfuf2=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(2,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffu(fh,ch,x2,mu2)/(1-x2) 
             +as(mu2)/pi*ffu(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 

            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-2,x1,mu2)/(1-x1)+as(mu2)/pi* pdf(-2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffub(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffub(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    double pdfdf1=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(2,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffd(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffd(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-2,x1,mu2)/(1-x1)
             +as(mu2)/pi* pdf(-2,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffdb(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffdb(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );
    double pdfdf2=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(1,x1,mu2)/(1-x1) +as(mu2)/pi*pdf(1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffd(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffd(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) + 
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pdf(-1,x1,mu2)/(1-x1)
             +as(mu2)/pi* pdf(-1,xi1,mu2)/xi1 *2.0/3.0*(1-z1)
             +as(mu2)/pi* pdf(0,xi1,mu2)/xi1 *1.0/2.0*(1-z1)*z1  ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*ffdb(fh,ch,x2,mu2)/(1-x2)
             +as(mu2)/pi*ffdb(fh,ch,xi2,mu2)/xi2*(2.0/3.0*(1-z2)+4.0/3.0*(1+z2*z2)/(1-z2)*log(z2))
             +as(mu2)/pi*ffg(fh,ch,xi2,mu2)/xi2 *(2.0/3.0*z2 +4.0/3.0*(1+(1-z2)*(1-z2))/z2*log(z2)   )  ) );

    return pdfuf1+pdfdf1+pdfuf2+pdfdf2;
}/*}}}*/

double ker1_fun_unpo_N(double b,double Q2, double pt, double z ,double x , double xi1a, double zi2a ,int cn   )/*{{{*/
{ 

    double funpo=unpopdf_dis_N( x, z ,  xi1a, zi2a , b ,1  ,  cn   );
    double crsunpo=Singma0(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*pi))*(TMath::BesselJ0(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3unpoDIS(b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_unpo_P(double b,double Q2, double pt, double z ,double x , double xi1a, double zi2a ,int cn   )/*{{{*/
{ 
    double funpo=unpopdf_dis_P( x, z ,  xi1a, zi2a , b ,1  ,  cn   );
    double crsunpo=Singma0(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*pi))*(TMath::BesselJ0(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3unpoDIS(b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_unpo_D(double b,double Q2, double pt, double z ,double x , double xi1a, double zi2a ,int cn   )/*{{{*/
{ 
    double funpo=unpopdf_dis_D( x, z ,  xi1a, zi2a , b ,1  ,  cn   );
    double crsunpo=Singma0(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*pi))*(TMath::BesselJ0(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3unpoDIS(b,x,z);
    return crsunpo ;
} /*}}}*/

double ker2_fun_unpo_pip_xz_N(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_N(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],1  );
}/*}}}*/

double ker2_fun_unpo_pim_xz_N(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_N(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],-1  );
}/*}}}*/

double ker2_fun_unpo_pip_xz_P(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_P(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],1  );
}/*}}}*/

double ker2_fun_unpo_pim_xz_P(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_P(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],-1  );
}/*}}}*/

double ker2_fun_unpo_pip_xz_D(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_D(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],1  );
}/*}}}*/

double ker2_fun_unpo_pim_xz_D(const double * xi)/*{{{*/
{
    return ker1_fun_unpo_D(xi[0],Q2it,phtit/zit,zit,xit, xi[1],xi[2],-1  );
}/*}}}*/

double UNPOLARISED_N_PIP (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit = zi;
    xit = xi;
    Q2it = Q2i;
    phtit = phti;    
    ROOT::Math::Functor wfunpo(&ker2_fun_unpo_pip_xz_N,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);    
    return igunpo.Integral(a,b);
}/*}}}*/

double UNPOLARISED_N_PIM (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_unpo_pim_xz_N,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double UNPOLARISED_P_PIP (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_unpo_pip_xz_P,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double UNPOLARISED_P_PIM (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };    
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_unpo_pim_xz_P,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double pop_u(  double x,  double Q2  )/*{{{*/
{
    double Q= sqrt(Q2);
    double pdf_a[13];
    hoppetEval(x, Q, pdf_a);
    return pdf_a[6+2]/x;
    //  return ffu (1,1 ,x,Q2) ;
}/*}}}*/

double pop_d(   double x,  double Q2  )/*{{{*/
{
    double Q= sqrt(Q2);
    double pdf_a[13];
    hoppetEval(x, Q, pdf_a);
    return pdf_a[6-2]/x;
    //    return ffd (1,1 ,x,Q2);
}/*}}}*/

double popdf_dis_N_pi(double x1, double  x2  ,double b , double xi1a, double xi2a     )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    //double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    return pdfuf+pdfdf;
    //  return pdfuf+pdfdf+pdfubf+pdfdbf;
}/*}}}*/

double popdf_dis_N_pim(double x1, double  x2  ,double b  , double xi1a, double xi2a   )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    // double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    return pdfuf+pdfdf;
    // return pdfuf+pdfdf+pdfubf+pdfdbf;
}/*}}}*/

double popdf_dis_P_pi(double x1, double  x2  ,double b , double xi1a, double xi2a     )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    //double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_pi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                          ) ;

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_pi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    return pdfuf+pdfdf;
    //  return pdfuf+pdfdf+pdfubf+pdfdbf;
}/*}}}*/

double popdf_dis_P_pim(double x1, double  x2  ,double b  , double xi1a, double xi2a   )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    //double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_mpi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    double pdfdf=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_mpi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    return pdfuf+pdfdf;
    // return pdfuf+pdfdf+pdfubf+pdfdbf;
}/*}}}*/

double popdf_dis_D_pi(double x1, double  x2  ,double b  , double xi1a, double xi2a    )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    //double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;

    double pdfuf1=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf1=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_pi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;
    double pdfuf2=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf2=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_pi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;
    double pdfdf1=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf1=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_pi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    double pdfdf2=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_pi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf2=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_pi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_pi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    return pdfuf1+pdfdf1+pdfuf2+pdfdf2;
    // return pdfuf1+pdfdf1+pdfubf1+pdfdbf1+pdfuf2+pdfdf2+pdfubf2+pdfdbf2;
}/*}}}*/

double popdf_dis_D_pim(double x1, double  x2  ,double b  , double xi1a, double xi2a    )/*{{{*/
{
    double bm2=1.5;
    double bs2= b*b/(1+b*b/(bm2*bm2) )/C1/C1;
    double mu2=1.0/bs2;
    //double xi1=(1-x1)*xi1a+x1;
    double xi2=(1-x2)*xi2a+x2;
    //double z1=x1/xi1;
    double z2=x2/xi2;
    double pdfuf1=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf1=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_mpi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    double pdfuf2=(1-x1)*(1-x2)*euq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_u_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_u_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )                          
            ) ;

    // double pdfubf2=(1-x1)*(1-x2)*euq2*(
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				     ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_ub_mpi(x2,mu2)/(1-x2) 
    // 				                     +as(mu2)/pi*pof_ub_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;
    double pdfdf1=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_d(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf1=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_db(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_mpi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;

    double pdfdf2=(1-x1)*(1-x2)*edq2*(
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_u(x1,mu2)/(1-x1)   ) *
            ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_d_mpi(x2,mu2)/(1-x2) 
             +as(mu2)/pi*pof_d_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
            )
            ) ;

    // double pdfdbf2=(1-x1)*(1-x2)*edq2*(
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pop_ub(x1,mu2)/(1-x1)   ) *
    // 				      ((1+as(mu2)/pi/3.0*(-8.0 ) )*pof_db_mpi(x2,mu2)/(1-x2) 
    //                                                     +as(mu2)/pi*pof_db_mpi(xi2,mu2)/(xi2)*(4.0/3.0*(2.0*z2)/(1-z2)*log(z2))
    // 		                                      )
    //                                                                                           ) ;


    return pdfuf1+pdfdf1+pdfuf2+pdfdf2;
    // return pdfuf1+pdfdf1+pdfubf1+pdfdbf1+pdfuf2+pdfdf2+pdfubf2+pdfdbf2;
}/*}}}*/

double ker1_fun_po_N_pip(double b,double Q2, double pt, double z ,double x  ,  double BNLY , double xi1a, double xi2a  )/*{{{*/
{ 
    double funpo=popdf_dis_N_pi( x, z  , b ,  xi1a,  xi2a  );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_po_N_pim(double b,double Q2, double pt, double z ,double x   ,  double BNLY , double xi1a, double xi2a  )/*{{{*/
{ 
    double funpo=popdf_dis_N_pim( x, z  , b,  xi1a,  xi2a   );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_po_D_pip(double b,double Q2, double pt, double z ,double x   ,  double BNLY , double xi1a, double xi2a )/*{{{*/
{ 
    double funpo=popdf_dis_D_pi( x, z  , b ,  xi1a,  xi2a  );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_po_D_pim(double b,double Q2, double pt, double z ,double x   ,  double BNLY , double xi1a, double xi2a )/*{{{*/
{ 
    double funpo=popdf_dis_D_pim( x, z  , b ,  xi1a,  xi2a   );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_po_P_pip(double b,double Q2, double pt, double z ,double x  ,  double BNLY , double xi1a, double xi2a  )/*{{{*/
{ 
    double funpo=popdf_dis_P_pi( x, z  , b ,  xi1a,  xi2a  );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker1_fun_po_P_pim(double b,double Q2, double pt, double z ,double x   ,  double BNLY , double xi1a, double xi2a  )/*{{{*/
{ 
    double funpo=popdf_dis_P_pim( x, z  , b,  xi1a,  xi2a   );
    double crsunpo=Singma0_collins(Q2,Sit,x)*(funpo)*(2*pi*pt/z)*(b/(2.0*z))*(b/(2.0*pi))*(TMath::BesselJ1(pt*b))*
        Sudakov1(Q2,b)*Sudakov2(Q2,b)*Sudakov3poDIS(BNLY,b,x,z);
    return crsunpo ;
} /*}}}*/

double ker2_fun_po_pip_xz_N(const double * xi)/*{{{*/
{
    return ker1_fun_po_N_pip(xi[0],Q2it,phtit/zit,zit,xit , BNLYit ,xi[1] ,xi[2]   );
}/*}}}*/

double ker2_fun_po_pim_xz_N(const double * xi)/*{{{*/
{
    return ker1_fun_po_N_pim(xi[0],Q2it,phtit/zit,zit,xit  , BNLYit ,xi[1] ,xi[2]  );
}/*}}}*/

double ker2_fun_po_pip_xz_P(const double * xi)/*{{{*/
{
    return ker1_fun_po_P_pip(xi[0],Q2it,phtit/zit,zit,xit , BNLYit ,xi[1] ,xi[2]   );
}/*}}}*/

double ker2_fun_po_pim_xz_P(const double * xi)/*{{{*/
{
    return ker1_fun_po_P_pim(xi[0],Q2it,phtit/zit,zit,xit  , BNLYit ,xi[1] ,xi[2]  );
}/*}}}*/

double ker2_fun_po_pip_xz_D(const double * xi)/*{{{*/
{
    return ker1_fun_po_D_pip(xi[0],Q2it,phtit/zit,zit,xit  , BNLYit ,xi[1] ,xi[2]  );
}/*}}}*/

double ker2_fun_po_pim_xz_D(const double * xi)/*{{{*/
{
    return ker1_fun_po_D_pim(xi[0],Q2it,phtit/zit,zit,xit  , BNLYit ,xi[1] ,xi[2]  );
}/*}}}*/

double POLARISED_N_PIP (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_po_pip_xz_N,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double POLARISED_N_PIM (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_po_pim_xz_N,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double POLARISED_P_PIP (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_po_pip_xz_P,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

double POLARISED_P_PIM (double Eb, double zi, double xi, double Q2i, double phti) /*{{{*/
{ 
    Sit=Eb*0.94*2;
    double a[3] =  {b_l   ,0  ,0 };
    double b[3] =  {b_u   ,1  ,1 };
    zit=zi;
    xit=xi;
    Q2it=Q2i;
    phtit=phti;
    ROOT::Math::Functor wfunpo(&ker2_fun_po_pim_xz_P,3); 
    ROOT::Math::IntegratorMultiDim igunpo(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    igunpo.SetFunction(wfunpo);
    return igunpo.Integral(a,b);
}/*}}}*/

// tensor charge at Q2 = 10 GeV^2 //////////////////

double inttensorup( double x, void * t){/*{{{*/
    double Q2 = 10;
    return pop_u(x,Q2);  
}/*}}}*/

double inttensordown( double x, void * t){/*{{{*/
    double Q2 = 10;
    return pop_d(x,Q2);  
}/*}}}*/

double tensorcharge_up( double x_min, double x_max){/*{{{*/
    // Create the Integrator
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);

    // Set parameters of the integration
    //ROOT::Math::Functor tens(&inttensorup,1); 
    ig.SetFunction(&inttensorup);
    ig.SetRelTolerance(0.001);

    return ig.Integral(x_min, x_max);
}/*}}}*/

double tensorcharge_down( double x_min, double x_max){/*{{{*/
    // Create the Integrator
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
    // Set parameters of the integration
    ig.SetFunction(&inttensordown);
    ig.SetRelTolerance(0.001);

    return ig.Integral(x_min, x_max);
}/*}}}*/

/*}}}*/
////////////////////// END ///////////////////////////////////////

///// Main ///////////////////////////////////////////////////////////
/*{{{*/
double   ymax  = 12.;        //      ! max value of ln 1/x
double   dy    = 0.1;        //     ! the internal grid spacing (smaller->higher accuarcy)
// ! 0.1 should provide at least 10^{-3} accuracy 
double   Qmin  = 0.5;  //   ! smallest Q value in tabulation
double   Qmax  = 1e+4;       // ! largest Q value in tabulation
double   dlnlnQ = dy/4.;     // ! tabulation spacing in dlnlnQ (dy/4 recommended)
int   nloop  = 1;            // ! the number of loops to initialise (max=3!)
int   order  = -6;           // ! numerical interpolation order (-6 is a good choice)
int   scheme = 5;            // ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar, 
// ! 4=frag, 5=TransMsbar

bool VERBOSE =kTRUE; // print everything or not, level-0, print only necessary
bool VERBOSE1 =kFALSE; // print everything or not, level-1, print more
bool VERBOSE2 =kFALSE;// print everything or not, level-2, print much more

bool CALCULATE_ChiSQ = true;// calculate ChiSQ 
bool GENERATE_TRANS = false;//

//========================================================= main
int main(int argc, char **argv)
{
    // SETTINGS FOR THE PROGRAM BEGIN  
    TString EXPERIMENT; 
    TString TARGET; 
    TString PARTICLE; 
    TString ERROR;//STAT, TOTAL
    
    // FIRST argument : EXPERIMENT
    if(strcmp(argv[1],"-SOLID")== 0)
      EXPERIMENT = "solid";
    else if(strcmp(argv[1],"-CLAS12")== 0)
      EXPERIMENT = "clas";
    else {
        cerr << "Experiment " << argv[1] << " is not implemented";
        return 1;
    }

    // SECOND argument : TARGET
    if(strcmp(argv[2],"-NEUTRON")== 0)
      TARGET = "_neutron";
    else if(strcmp(argv[2],"-PROTON")== 0)
      TARGET = "_proton";
    else {
        cerr << "Target " << argv[2] << " is not implemented";
        return 1;
    }

    // THIRD argument : PRODUCED PARTICLE
    if(strcmp(argv[3],"-PI+")== 0)
      PARTICLE = "_pip";
    else if(strcmp(argv[3],"-PI-")== 0)
      PARTICLE = "_pim";
    else {
        cerr << "Particle " << argv[3] << " is not implemented";
        return 1;
    }

    // FOURTH argument : ERROR TYPE
    if(strcmp(argv[4], "-STAT") == 0)
      ERROR = "stat";
    else if(strcmp(argv[4], "-TOTAL") == 0)
      ERROR = "total";
    else if(strcmp(argv[4], "-SYST") == 0)
      ERROR = "syst";
    else {
      cerr << "Error " << argv[4] << " is not implemented";
      return 1;
    } 

    // Sixth: Read in how many segments 
    int segs= int(strtod(argv[5], NULL));

    // Seventh: Read in which segment to use 
    int seg_used= int(strtod(argv[6], NULL));

    //Sixth: READ THE GENERATED 13-PARAMETER SETS BEGIN 
    int N = 6;
    BNLYit=strtod(argv[N+1], NULL);
    Nuit=strtod(argv[N+2], NULL);
    Ndit=strtod(argv[N+3], NULL);
    auit=strtod(argv[N+4], NULL);
    adit=strtod(argv[N+5], NULL);
    buit=strtod(argv[N+6], NULL);
    bdit=strtod(argv[N+7], NULL);

    Nuit_t=strtod(argv[N+8], NULL);
    Ndit_t=strtod(argv[N+9], NULL);
    auit_t=strtod(argv[N+10], NULL);
    adit_t=strtod(argv[N+11], NULL);
    buit_t=strtod(argv[N+12], NULL);
    bdit_t=strtod(argv[N+13], NULL);
    /*}}}*/

    if(VERBOSE){
       cout<<endl<<"======================================"<<endl;
       cout<<"=== EXP = "<<EXPERIMENT.Data()<<",  TARGET="<<TARGET.Data()<<",   PID="<<PARTICLE.Data()<<endl;
       cout<<Form("=== 13 Parameters from Inputs : ===")<<endl;
       cout<<Form("===  B=%f, Nu=%f, Nd=%f, Au=%f, Ad=%f, Bu=%f, Bd=%f",
                BNLYit, Nuit, Ndit, auit, adit, buit, bdit)<<endl;
       cout<<Form("===  Nu_t=%f, Nd_t=%f, Au_t=%f, Ad_t=%f, Bu_t=%f, Bd_t=%f",
                Nuit_t, Ndit_t, auit_t, adit_t, buit_t, bdit_t)<<endl;
       cout<<"======================================"<<endl<<endl;
    }
    
    // START HOPPET/*{{{*/
    hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order, scheme);
    list_pdf_f();
    list_ff_f();
    double asQ0 = as(2.4), Q0=sqrt(2.4);
    hoppetEvolve(asQ0, Q0, nloop, 1.0, heralhc_init, Q0);  /*}}}*/


    double ChiSQ = 0.0;
    if (CALCULATE_ChiSQ){
      TString datafile = "./data/"+EXPERIMENT+TARGET+PARTICLE+".root";
      TFile * fd = new TFile(datafile, "r");
      if (fd->IsOpen() && VERBOSE){
	cout << ">>>> Reading data from file: " << datafile.Data() << endl;
	cout << ">>>> Divided into " << segs << " segments and the " << seg_used << " is used." << endl;
      }
      else {
	cout << "!!! Fail to open data file: " << datafile.Data() << endl;
	return 1;
      }
      TTree * Td = (TTree *) fd->Get("data");
      int Nt = (int) Td->GetEntries();
      double Ebeam, z, x, Q2, Pt, unpol, Asym, Estat, Erel, Eabs;
      Td->SetBranchAddress("Ebeam", &Ebeam);
      Td->SetBranchAddress("z", &z);
      Td->SetBranchAddress("x", &x);
      Td->SetBranchAddress("Q2", &Q2);
      Td->SetBranchAddress("Pt", &Pt);
      Td->SetBranchAddress("unpol", &unpol);
      Td->SetBranchAddress("Asym", &Asym);
      Td->SetBranchAddress("Estat", &Estat);
      Td->SetBranchAddress("Erel", &Erel);
      Td->SetBranchAddress("Eabs", &Eabs);
      int Nseg = Nt / segs;
      int Nstart = Nseg * seg_used;
      int Nend = Nseg * (seg_used + 1);
      if (seg_used == segs - 1) Nend = Nt;
      if (VERBOSE){
	cout << "<<<< Number of data points in this segment is " << Nend - Nstart << endl;
	cout << "<<<< " << ERROR.Data() << " error is used." << endl;
	cout << "---- CALCULATING ChiSQ ----" << endl;
      }
      double newpol, Err2;
      for (int i = Nstart; i < Nend; i++){
	Td->GetEntry(i);
	if (TARGET.EqualTo("_neutron") && PARTICLE.EqualTo("_pip"))
	  newpol = POLARISED_N_PIP(Ebeam, z, x, Q2, Pt);
	else if (TARGET.EqualTo("_neutron") && PARTICLE.EqualTo("_pim"))
	  newpol = POLARISED_N_PIM(Ebeam, z, x, Q2, Pt);
	else if (TARGET.EqualTo("_proton") && PARTICLE.EqualTo("_pip"))
	  newpol = POLARISED_P_PIP(Ebeam, z, x, Q2, Pt);
	else if (TARGET.EqualTo("_proton") && PARTICLE.EqualTo("_pim"))
	  newpol = POLARISED_P_PIM(Ebeam, z, x, Q2, Pt);
	else {
	  cout << "!!! Fail to match target and particle implemented!" << endl;
	  return 1;
	}
	if (ERROR.EqualTo("stat"))
	  Err2 = Estat * Estat;
	else if (ERROR.EqualTo("total"))
	  Err2 = Estat * Estat + Eabs * Eabs + Erel * Erel * Asym * Asym;
	else {
	  cout << "!!! Fail to match error type!" << endl;
	  return 1;
	}
	ChiSQ = ChiSQ + pow( newpol/unpol - Asym, 2) / Err2;
      }
      if (VERBOSE){
	cout << "---- End Calculating ChiSQ ----" << endl;
	cout << "--- The ChiSQ for this replica is " << ChiSQ << endl;
      }
    }


    // CALCULATE TENSOR CHARGE/*{{{*/
    if (GENERATE_TRANS){
      double x_min_tens = 1.e-5;
      double x_max_tens = 0.9999;
      double x_min_jlab = 0.05;
      double x_max_jlab = 0.65;
      double tensor_up        = tensorcharge_up(x_min_tens, x_max_tens);
      double tensor_down      = tensorcharge_down(x_min_tens, x_max_tens);
      double tensor_up_jlab   = tensorcharge_up(x_min_jlab, x_max_jlab);
      double tensor_down_jlab = tensorcharge_down(x_min_jlab, x_max_jlab);/*}}}*/

      TString filename;
      filename.Form("data/ChiSQ_%s%s%s_%d_%d.dat",EXPERIMENT.Data(),TARGET.Data(),PARTICLE.Data(),segs, seg_used);
      ofstream outfile2(filename);	
      outfile2<<Form(" %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e",
		     BNLYit, Nuit, Ndit, auit, adit, buit, bdit);
      outfile2<<Form("  %12.4e %12.4e  %12.4e  %12.4e  %12.4e  %12.4e  %12.4e",
		     Nuit_t, Ndit_t, auit_t, adit_t, buit_t, bdit_t, ChiSQ)<<endl;
      outfile2<<Form("%12.4e  %12.4e  %12.4e  %12.4e", tensor_up, tensor_up_jlab, tensor_down, tensor_down_jlab)<<endl;

    // write files with transversity at 10 GeV^2 
      double Q2 = 10;//TL: I prefer 2.4 because 10 is out of JLab range and it is hard to determine the imcomplete tensor charge measured in JLab x range due to the evolution effect
      double xmin = 1.e-4;
      double xmax = 0.999;
      int npoints = 100;
      double step = (xmax-xmin)/(npoints-1);
    
      for( int i = 0; i < npoints; i++){
	double x = xmin + i * step;
	outfile2 << x << " " << x * pop_u(x,Q2) << " " << x * pop_d(x,Q2) << endl;
      }    
      outfile2.close();
    }
    if(VERBOSE) cout<<"done"<<endl;
    return ChiSQ;
}/*}}}*/
