#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TStyle.h"
#include "TParameter.h"
#include "TGaxis.h"

void set_style_ana(TH1* h, int col, int rebin, bool sumw2) {
  return;
  if (rebin>1)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetLabelOffset(1.2);
  h->GetXaxis()->SetTitleOffset(0.9);

  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(1.2);
  h->GetYaxis()->SetTitleOffset(0.9);

  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);

  if (col>0) {
    h->SetLineWidth(2);
    h->SetLineColor(col);
  }
}

void log() {
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
}

void plot_mu(TString fin="output/ms00_ds00_xy0.5mm_dth10_dph10_p1.30.root"){

  gStyle->SetOptStat("e");
  //gStyle->SetStatStyle(0);
  TGaxis::SetMaxDigits(3);

  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetNdivisions(507,"X");
  gStyle->SetNdivisions(507,"Y");
  gStyle->SetNdivisions(507,"Z");

  //  history[0].fName= beam
  //  history[1].fName= det1
  //  history[2].fName= det2
  //  history[3].fName= plane
  //  history[4].fName= diamond
  //  history[5].fName= hades
  //  history[6].fName= target
  //  ------------
  //  vms_history[0].fName= ORIG_BEAM_PROFILE
  //  ------------
  //  vpidec_history[0].fName= pidec_PRE_PIDEC
  //  vpidec_history[1].fName= pidec_POST_PIDEC
  //  vpidec_history[2].fName= pidec_POST_PIDEC_BEAM_PROF
  //  vpidec_history[3].fName= pidec_POST_PIDEC_RETRACE
  //  ------------
  //  vsolution[0].fName= TDR_SOLUTION_BEAM_INITIAL
  //  vsolution[1].fName= TDR_SOLUTION_DIAM
  //  vsolution[2].fName= TDR_SOLUTION_HAD
  //  vsolution[3].fName= MIN_SOLUTION_BEAM_INITIAL
  //  vsolution[4].fName= MIN_SOLUTION_DIAM
  //  vsolution[5].fName= MIN_SOLUTION_HAD
  //  vsolution[6].fName= DIGI_POS_det1
  //  vsolution[7].fName= DIGI_POS_det2
  //  vsolution[8].fName= DIGI_POS_diamond

  TFile *f = TFile::Open(fin);
  TTree *t = (TTree*) f->Get("t");

  // infer reference momentum
  double p_ref = -9999;
  for (int i=0; i<t->GetUserInfo()->GetEntries(); ++i)  {
    TParameter<double> *ref_mom = (TParameter<double>*) t->GetUserInfo()->At(0);
    const char * name = ref_mom->GetName();
    cout << "TParameter " << i << ": name= " << name << " val= " << ref_mom->GetVal() << endl;
    if (strcmp(name, "ref_beam_mom")==0) {
      p_ref = ref_mom->GetVal();
      break;
    }
  }
  if (p_ref<0) {
    cout << "Reference momentum not found in ntuple... stop! " << endl;
    return;
  }

  t->SetAlias("pbeam","p[0]");
  t->SetAlias("thbeam","th[0]");
  t->SetAlias("phbeam","ph[0]");

  t->SetAlias("ppi","ppidec[0]");
  t->SetAlias("dppi",Form("100*(ppi-%3.5f)/%3.5f",p_ref,p_ref));

  t->SetAlias("pmu", "ppidec[1]");
  t->SetAlias("dpmu",Form("100*(pmu-%3.5f)/%3.5f",p_ref,p_ref));

  t->SetAlias("pmub","ppidec[2]");
  t->SetAlias("pmur","ppidec[3]");

  t->SetAlias("thmu","thpidec[1]");
  t->SetAlias("thmub","thpidec[2]");
  t->SetAlias("thmur","thpidec[3]");
  t->SetAlias("thmur_er","thmur-thmu");

  t->SetAlias("phmu","phpidec[1]");
  t->SetAlias("phmub","phpidec[2]");
  t->SetAlias("phmur","phpidec[3]");
  t->SetAlias("phmur_er","phmur-phmu");

  t->SetAlias("theta_mu_pi", "TMath::Hypot(thpidec[1]-thpidec[0],phpidec[1]-phpidec[0])");

  t->SetAlias("del0",Form("(abs((ppi-%3.5f)/%3.5f)-1.0)<0.001",p_ref,p_ref));
  for (int idel=-6; idel<7; ++idel) {
    TString alias = "del0";
    if (idel>0) alias = Form("del_p%d",idel);
    if (idel<0) alias = Form("del_m%d",idel);
    TString formula = "abs(dppi)<0.1";
    if (idel>=0) formula = Form("abs(dppi-%3.1f)<0.1",(float)idel);
    if (idel<0) formula = Form("abs(dppi+%3.1f)<0.1",(float)-idel);
    cout << "Setting alias " << alias << " for formula: " << formula << endl;
    t->SetAlias(alias,formula);
  }

  float w = 100000./t->GetEntries();
  cout <<  "nent= " << t->GetEntries() << " weight = " << w << endl;

  //t->Draw("ppi-pmu", "del_p2");
  //t->Draw("dppi", "del_p2");
  //t->Draw("pmu:theta_mu_pi", "del_p5");

  //t->Draw("pmu:theta_mu_pi", "ipidec<33&&acc");

  TH1F* h1=new TH1F("h1","num. #mu^{-} decay;element id; number of #pi^{-} decay b/n element i-1 and i",35,-0.5,34.5);
  TH1F* h2=new TH1F("h2","num. #mu^{-} in acceptance;element id; numbr of accepted #mu^{-};",35,-0.5,34.5);
  set_style_ana(h1,4,0,0);
  set_style_ana(h2,4,0,0);
  t->Project("h1","ipidec","dpidec<0");
  t->Project("h2","ipidec","dpidec<0&&acc","same");
  //t->Project("h1","ipidec","");
  //t->Project("h2","ipidec","acc","same");
  TH1F* hr = (TH1F*)h2->Clone("hr");
  hr->SetTitle("Percentage of #mu^{-} accepted;element id;acceptance");
  hr->Divide(h1);
  hr->Scale(100);
  set_style_ana(hr,4,0,0);

  TCanvas *tc = new TCanvas("tc","tc");
  tc->Divide(3,2);
  tc->cd(2);
  h1->Draw();
  //h2->Draw("same");
  //h2->SetLineColor(2);
  //log();
  tc->cd(3);
  h2->Draw();
  log();
  tc->cd(4);
  hr->Draw();
  //log();
  tc->cd(1);
  t->Draw("pmu:theta_mu_pi>>h2d","","box");

  TH2F *h2d = (TH2F*) gDirectory->Get("h2d");
  h2d->SetTitle("Decay kinematics;decay polar angle (#pi^{-} rest frame)[mrad]; momentum [GeV/c]");
  cout << "h2d_ptr= " << h2d << " axX= " << h2d->GetXaxis() << endl;
  h2d->GetYaxis()->SetRangeUser(0,p_ref*1.2);

  tc->cd(5);
  t->Draw("pmu>>hpmu","dpidec<0&&acc");
  TH1F* hpmu = (TH1F*) gDirectory->Get("hpmu");
  hpmu->SetTitle("Mom. distr. of accepted #mu^{-}; Momentum [GeV/c]");

  tc->cd(6);
  t->Draw("ppi>>hppi","dpidec>0&&acc");
  TH1F* hppi = (TH1F*) gDirectory->Get("hppi");
  hppi->SetTitle("Mom. distr. of accepted #pi^{-}; Momentum [GeV/c]");

  float r1 = h1->GetEntries()/t->GetEntries();
  float r2 = h2->GetEntries()/h1->GetEntries();
  float r3 = h2->GetEntries()/(t->GetEntries()-h1->GetEntries());
  float r4 = h2->GetEntries()/(hppi->GetEntries());

  cout << "R1 = " << r1 << endl;
  cout << "R2 = " << r2 << endl;
  cout << "R3(=contam) = " << r3 << endl;
  cout << "R4(=contam) = " << r4 << endl;

  //TCanvas *tc = new TCanvas("tc","tc");
  //tc->Divide(3,2);
  //tc->cd(1);
  //h1->Draw();
  ////h2->Draw("same");
  ////h2->SetLineColor(2);
  ////log();
  //tc->cd(2);
  //h2->Draw();
  //log();
  //tc->cd(3);
  //hr->Draw();
  ////log();
  //tc->cd(4);
  //t->Draw("pmu:theta_mu_pi>>h2d","","box");
  //gDirectory->ls();
  //TH2F *h2d = (TH2F*) gDirectory->Get("h2d");
  //cout << "h2d_ptr= " << h2d << " axX= " << h2d->GetXaxis() << endl;
  ////h2d->GetYaxis()->SetRangeUser(0,1.2);
  //
  //tc->cd(5);
  //t->Draw("pmu","dpidec<0&&acc");

  //t->Draw("ipidec","ipidec<33");
  //t->Draw("ipidec","ipidec<33&&acc","same");

  //t->Draw("pmur-pmu","del_p5&&dpidec<0");

  //new TCanvas();
  //t->Draw("phmur_er","dpidec<0&&abs(phmur_er)<10");
  //new TCanvas();
  //t->Draw("thmur_er","dpidec<0&&abs(thmur_er)<10");

  //new TCanvas();
  //t->Draw("phbeam:thbeam>>h_beam_ph_vs_th(200,-15,15,200,-15,15)","");
  //h_beam_ph_vs_th->Draw();
  //new TCanvas();
  //t->Draw("phmub:thmub>>h_beam_back_prop_ph_vs_th(200,-500,500,200,-500,500)","dpidec<0");
  //h_beam_back_prop_ph_vs_th->Draw();
  //new TCanvas();
  //t->Draw("phmur:thmur>>h_beam_retrace_ph_vs_th(200,-40,40,200,-40,40)","dpidec<0");
  //h_beam_retrace_ph_vs_th->Draw();

}
