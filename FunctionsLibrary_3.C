#include <iostream>
#include "TMath.h"
#include "TMath.h"
#include "TF1.h"
#include "TROOT.h"

Double_t scaling_factor = 1.05154; //factor introduced to pass from the sigma of Jpsi to the sigma of Psi(2S) -> WARNING! TO SUBSTITUTE
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
Double_t Func_VWG(Double_t *x, Double_t *par){
  Double_t sigma = par[2] + par[3]*((x[0] - par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0] - par[1])*(x[0] - par[1])/(2.*sigma*sigma));
}
//------------------------------------------------------------------------------
Double_t Func_POL4EXP(Double_t *x, Double_t *par){
  return par[0]*TMath::Exp(x[0]/par[1])*(1 + par[2]*x[0] + par[3]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0]*x[0]);
}
//==============================================================================
// J/psi
//==============================================================================
Double_t Func_Jpsi_CB2_VWG(Double_t *x, Double_t *par){
  Double_t t = (x[0] - par[5])/par[6];
  if (par[7] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);
  if (t >= -absAlpha && t < absAlpha2){return par[4]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha) {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;
    return par[4]*(a/TMath::Power(b - t, par[8]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;
   return  par[4]*(c/TMath::Power(d + t, par[10]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Jpsi_CB2_POL4EXP(Double_t *x, Double_t *par){
  Double_t t = (x[0] - par[7])/par[8];
  if (par[9] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[9]);
  Double_t absAlpha2 = fabs((Double_t)par[11]);
  if (t >= -absAlpha && t < absAlpha2){return par[6]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha) {
    Double_t a =  TMath::Power(par[10]/absAlpha,par[10])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[10]/absAlpha - absAlpha;
    return par[6]*(a/TMath::Power(b - t, par[10]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[12]/absAlpha2,par[12])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[12]/absAlpha2 - absAlpha2;
   return  par[6]*(c/TMath::Power(d + t, par[12]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Jpsi_CB2_fix(Double_t *x, Double_t *par){
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  if (t >= -absAlpha && t < absAlpha2){return par[0]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha){
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;
   return  par[0]*(c/TMath::Power(d + t, par[6]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Jpsi_NA60_VWG(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[13];
  Double_t t2 = par[14];
  Double_t t = (x[0]-par[5])/par[6];

  if (t >= t1 && t < t2){
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[7]*(t1-t),par[8]-par[9]*TMath::Sqrt(t1-t)));
  } else if (t >= t2){
    t0 = (1+TMath::Power(par[10]*(t-t2),par[11]-par[12]*TMath::Sqrt(t-t2)));
  }
  return par[4]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//------------------------------------------------------------------------------
Double_t Func_Jpsi_NA60_POL4EXP(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[15];
  Double_t t2 = par[16];
  Double_t t = (x[0]-par[7])/par[8];

  if (t >= t1 && t < t2){
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[9]*(t1-t),par[10]-par[11]*TMath::Sqrt(t1-t)));
  } else if (t >= t2){
    t0 = (1+TMath::Power(par[12]*(t-t2),par[13]-par[14]*TMath::Sqrt(t-t2)));
  }
  return par[6]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//------------------------------------------------------------------------------
Double_t Func_Jpsi_NA60_fix(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[9];
  Double_t t2 = par[10];
  Double_t t = (x[0]-par[1])/par[2];

  if (t >= t1 && t < t2) {
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[3]*(t1-t),par[4]-par[5]*TMath::Sqrt(t1-t)));
  } else if (t >= t2) {
    t0 = (1+TMath::Power(par[6]*(t-t2),par[7]-par[8]*TMath::Sqrt(t-t2)));
  }
  return par[0]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//==============================================================================
// Psi(2S)
//==============================================================================
Double_t Func_Psi2S_CB2_VWG(Double_t *x, Double_t *par){
  Double_t t = (x[0] - (par[5]+(3.686097-3.0969)))/(par[6]*1.05154);
  if (par[7] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);
  if (t >= -absAlpha && t < absAlpha2){return par[4]*par[11]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha) {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;
    return par[4]*par[11]*(a/TMath::Power(b - t, par[8]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;
   return  par[4]*par[11]*(c/TMath::Power(d + t, par[10]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Psi2S_CB2_POL4EXP(Double_t *x, Double_t *par){
  Double_t t = (x[0] - (par[7]+(3.686097-3.0969)))/(par[8]*1.05154);
  if (par[9] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[9]);
  Double_t absAlpha2 = fabs((Double_t)par[11]);
  if (t >= -absAlpha && t < absAlpha2){return par[6]*par[13]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha) {
    Double_t a =  TMath::Power(par[10]/absAlpha,par[10])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[10]/absAlpha - absAlpha;
    return par[6]*par[13]*(a/TMath::Power(b - t, par[10]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[12]/absAlpha2,par[12])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[12]/absAlpha2 - absAlpha2;
   return  par[6]*par[13]*(c/TMath::Power(d + t, par[12]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Psi2S_NA60_VWG(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[13];
  Double_t t2 = par[14];
  Double_t t = (x[0] - (par[5]+(3.686097-3.0969)))/(par[6]*1.05154);

  if (t >= t1 && t < t2){
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[7]*(t1-t),par[8]-par[9]*TMath::Sqrt(t1-t)));
  } else if (t >= t2){
    t0 = (1+TMath::Power(par[10]*(t-t2),par[11]-par[12]*TMath::Sqrt(t-t2)));
  }
  return par[4]*par[15]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//------------------------------------------------------------------------------
Double_t Func_Psi2S_NA60_POL4EXP(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[15];
  Double_t t2 = par[16];
  Double_t t = (x[0] - (par[7]+(3.686097-3.0969)))/(par[8]*1.05154);

  if (t >= t1 && t < t2){
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[9]*(t1-t),par[10]-par[11]*TMath::Sqrt(t1-t)));
  } else if (t >= t2){
    t0 = (1+TMath::Power(par[12]*(t-t2),par[13]-par[14]*TMath::Sqrt(t-t2)));
  }
  return par[6]*par[17]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//------------------------------------------------------------------------------
Double_t Func_Psi2S_CB2_fix(Double_t *x, Double_t *par){
  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0){t = -t;}

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);
  if (t >= -absAlpha && t < absAlpha2){return par[0]*par[7]*(exp(-0.5*t*t));}

  //left tail
  if (t < -absAlpha){
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*par[7]*(a/TMath::Power(b - t, par[4]));
  }
  //right tail
  if (t >= absAlpha2){
   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;
   return  par[0]*par[7]*(c/TMath::Power(d + t, par[6]));
  }
  return 0. ;
}
//------------------------------------------------------------------------------
Double_t Func_Psi2S_NA60_fix(Double_t *x, Double_t *par){
  Double_t t0 = 0.;
  Double_t t1 = par[9];
  Double_t t2 = par[10];
  Double_t t = (x[0]-par[1])/par[2];

  if (t >= t1 && t < t2) {
    t0 = 1;
  } else if (t < t1) {
    t0 = (1+TMath::Power(par[3]*(t1-t),par[4]-par[5]*TMath::Sqrt(t1-t)));
  } else if (t >= t2) {
    t0 = (1+TMath::Power(par[6]*(t-t2),par[7]-par[8]*TMath::Sqrt(t-t2)));
  }
  return par[0]*par[11]*TMath::Exp(-0.5*((t/t0)*(t/t0)));
}
//==============================================================================
// Total functions
//==============================================================================
Double_t Func_tot_CB2_VWG(Double_t *x, Double_t *par){
  return Func_VWG(x,par) + Func_Jpsi_CB2_VWG(x,par) + Func_Psi2S_CB2_VWG(x,par);
}
//------------------------------------------------------------------------------
Double_t Func_tot_CB2_POL4EXP(Double_t *x, Double_t *par){
  return Func_POL4EXP(x,par) + Func_Psi2S_CB2_POL4EXP(x,par) + Func_Jpsi_CB2_POL4EXP(x,par);
}
//------------------------------------------------------------------------------
Double_t Func_tot_NA60_VWG(Double_t *x, Double_t *par){
  return Func_VWG(x,par) + Func_Jpsi_NA60_VWG(x,par) + Func_Psi2S_NA60_VWG(x,par);
}
//------------------------------------------------------------------------------
Double_t Func_tot_NA60_POL4EXP(Double_t *x, Double_t *par){
  return Func_POL4EXP(x,par) + Func_Jpsi_NA60_POL4EXP(x,par) + Func_Psi2S_NA60_POL4EXP(x,par);
}
