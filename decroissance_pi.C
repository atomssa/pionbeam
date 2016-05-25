#include <iostream>
#include <cstdio>
#include <math>
//#include <vector>
#include <cstdlib>

#include "TGraph.h"
#include "TRandom3.h"

using namespace std;

const float m_nu = 0;
const float m_mu = 0.1057;
const float m_pi = 0.1396;

//const float _Ppi = 1.0;

const float Ecm_mu = ((m_pi*m_pi) + (m_mu*m_mu) ) / (2*m_pi);
const float Ecm_nu = ((m_pi*m_pi) - (m_mu*m_mu) ) / (2*m_pi);

const float Pcm_mu = sqrt(Ecm_mu*Ecm_mu - m_mu*m_mu);
const float Pcm_nu = sqrt(Ecm_nu*Ecm_nu - m_nu*m_nu);

float E_p(float p,float m){
	return sqrt(p*p + m*m);
}

float beta_e(float E,float m){
	return sqrt(1 - ((m*m)/(E*E)));
}
float beta_p(float p,float m){
	return 1./sqrt(1 + ((m*m)/(p*p)));
}
float gamma_b(float beta){
	return 1./sqrt(1-beta*beta);
}
float gamma_p(float p,float m){
	return gamma_b(beta_p(p, m));
}
float gamma_e(float E,float m){
	return gamma_b(beta_e(E, m));
}

float gamma_beta_p(float p, float m){
  return gamma_p(p,m)*beta_p(p,m);
}

float PLab_mu(float p_pi, float pcm, float Ecm){
	return gamma_p(p_pi, m_pi) * ( pcm + (beta_p(p_pi, m_pi) * Ecm));
}

float ELab_mu(float p_pi,float pcm, float Ecm){
	return gamma_p(p_pi, m_pi) * ( Ecm + (beta_p(p_pi, m_pi) * pcm));
}

void ThetaLab_mu(float p_pi, float theta_cm_mu, double &_th, double &_plab, bool do_print=false){
  float Pzcm_mu = Pcm_mu * TMath::Cos(theta_cm_mu);
  float Ptcm_mu = Pcm_mu * TMath::Sin(theta_cm_mu);
  float PzLab_mu = PLab_mu(p_pi, Pzcm_mu, Ecm_mu);
  float PtLab_mu = Ptcm_mu;
  _plab = TMath::Hypot(PzLab_mu, PtLab_mu);
  _th = TMath::ACos(PzLab_mu/_plab);
  if (do_print){
    cout << "Pcm_mu=" << Pcm_mu <<  "Ecm_mu=" << Ecm_mu << "PzLab_mu=" << PzLab_mu << "PtLab_mu=" << PtLab_mu << "Plab_mu=" << _plab << "th=" << _th <<   endl;
  }
}


/*int decroissance_pi(float _p_pi, float _theta_pi, float _phi_pi,
		float &_p_mu, float &_theta_mu, float &_phi_mu) {

	double _theta_cm_mu = rand.Uniform(TMath::Pi());
	double _phi_cm_mu = rand.Uniform(2*TMath::Pi());

	double _theta_lab_mu = 0.0, _p_lab_mu = 0.0;
	ThetaLab_mu(_p_pi, _theta_cm_mu, _theta_lab_mu, _p_lab_mu);

	//_theta_mu = _theta_pi + _theta_lab_mu;
	//_phi_mu = _phi_pi;

	//_theta_mu = _theta_pi;
	//_phi_mu = _phi_pi + _theta_lab_mu;

	_p_mu = _p_lab_mu;

	_theta_mu = _theta_pi + _theta_lab_mu*TMath::Cos(_phi_pi + _phi_cm_mu);
	_phi_mu = _phi_pi + _theta_lab_mu*TMath::Sin(_phi_pi + _phi_cm_mu);
}*/

int decroissance_pi(float _p_pi = 1.0){

	//	random.seed(60934386)
	TRandom3 rand;
	rand.SetSeed();
	//float _p_pi;
	//cout << "Entrez l'impulsion des pions (en GeV) : p = ";
	//cin >> _p_pi;

	double x[1000] = {0.0};
	double y[1000] = {0.0};

	for (int i=0;i<1000;i++){
	  //double _theta_cm_mu = rand.Uniform(TMath::Pi());
	  double _theta_cm_mu = TMath::ACos(gRandom->Uniform(2.0)-1.0);
	  double _phi_cm_mu = rand.Uniform(2*TMath::Pi());
	  if (i < 10){
	    cout << "i: " << i <<  " th_cm_mu = " << _theta_cm_mu << endl;
	  }
	  double _theta_lab_mu = 0.0, _p_lab_mu = 0.0;
	  ThetaLab_mu(_p_pi, _theta_cm_mu, _theta_lab_mu, _p_lab_mu);
	  cout << "th_lab_mu= " << _theta_lab_mu << " _p_lab_mu= " << _p_lab_mu << endl;
	  x[i] = 1000*_theta_lab_mu;
	  y[i] = _p_lab_mu;
	}

	TGraph *graph = new TGraph(1000, x, y);
	graph->SetTitle("Theta vs p_{LAB}");

	graph->Draw("A*");

	return 0;

}
