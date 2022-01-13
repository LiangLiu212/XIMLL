#ifndef ANALYSIS_H
#define ANALYSIS_H
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <numeric>
#include <TMath.h>
#include <math.h>
#include "TF1.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TRandom.h>
using namespace std;
namespace jj
{
  TLorentzVector ep, em, collision; 
  TLorentzVector sgm, sgmbar, n_vx, nbar_vx; 
  std::vector<TLorentzVector> fourvec;
  std::vector<TLorentzVector> vec;
std::vector<double> HelicityAngles(std::vector<TLorentzVector> fourvec){
  std::vector<double> HelicityAngles;
  HelicityAngles.resize(0);
  const double m_Jpsi    = 3.0969;
  em.SetPxPyPzE(0.011*m_Jpsi/2. ,0. , sqrt((pow(m_Jpsi,2)/4.)-(0.000511*0.000511)),( m_Jpsi/( sqrt(1-0.011*0.011) ) )/2. ); 
  ep.SetPxPyPzE(0.011*m_Jpsi/2. ,0. ,-sqrt((pow(m_Jpsi,2)/4.)-(0.000511*0.000511)),( m_Jpsi/( sqrt(1-0.011*0.011) ) )/2. );
  collision = ep + em;
  TVector3 bCMS = (collision).Vect(); bCMS *= 1/((collision).E());
  double sgmth, sgmphi, sgmbth, sgmbphi;
  TLorentzVector sgm, sgmbar, n, nbar;
  sgm = fourvec.at(0); sgmbar = fourvec.at(1); n = fourvec.at(2); nbar = fourvec.at(3); 
  //TVector3 bCMS = (sgm+sgmbar).Vect(); bCMS *= 1/((sgm+sgmbar).E());
  // Boost Xis in CM frame
  sgm.Boost(-bCMS);    n.Boost(-bCMS);
  sgmbar.Boost(-bCMS);  nbar.Boost(-bCMS);
  HelicityAngles.push_back(sgm.Theta());
  sgmth = sgm.Theta(); sgmphi = sgm.Phi();
   // Velocity of Xim and Xip
  TVector3 bSgm(sgm.Vect());    bSgm*=1/(sgm.E());
  n.Boost(-bSgm);   n.RotateZ(-sgmphi);   n.RotateY(-sgmth); 
  HelicityAngles.push_back(n.Theta()); HelicityAngles.push_back(n.Phi());
  // Doing the same for Xibar
  sgmbth = sgmbar.Theta(); sgmbphi = sgmbar.Phi();
  TVector3 bSgmbar(sgmbar.Vect()); bSgmbar*=1/(sgmbar.E());
  nbar.Boost(-bSgmbar);
  nbar.RotateZ(-sgmbphi);   nbar.RotateY(-sgmbth);
  HelicityAngles.push_back(nbar.Theta()); HelicityAngles.push_back(nbar.Phi());
  return HelicityAngles;
}
}  
#endif