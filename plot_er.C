#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

void set_style_ana(TH1* h, int col, int rebin, bool sumw2) {
  if (rebin>1)h->Rebin(rebin);
  if (sumw2)h->Sumw2();
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetLabelOffset(0.);
  h->GetXaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelOffset(0.02);
  h->GetYaxis()->SetTitleOffset(0.8);

  h->SetMarkerStyle(20);
  h->SetMarkerColor(col);
  h->SetMarkerSize(0.5);

  if (col>0) {
    h->SetLineWidth(1);
    h->SetLineColor(col);
  }
}

void log() {
  gPad->SetLogy();
  gPad->SetGridy();
  gPad->SetGridx();
}

void plot_er(TString fin="output/ms11_ds11_xy0.5mm_dth10_dph50_dp-6.0_p1.70.root"){

  //gStyle->SetOptStat(1111111);
  gStyle->SetStatStyle(0);

  //------------
  //history[0].fName= beam
  //history[1].fName= det1_POST_MS
  //history[2].fName= det2_POST_MS
  //history[3].fName= plane
  //history[4].fName= diamond
  //history[5].fName= hades
  //history[6].fName= target
  //------------
  //vms_history[0].fName= ORIG_BEAM_PROFILE
  //vms_history[1].fName= det1_PRE_MS
  //vms_history[2].fName= det1_POST_MS_BEAM_PROF
  //vms_history[3].fName= det1_POST_MS_RETRACE
  //vms_history[4].fName= det2_PRE_MS
  //vms_history[5].fName= det2_POST_MS_BEAM_PROF
  //vms_history[6].fName= det2_POST_MS_RETRACE
  //------------
  //vpd_history[0].fName= pidec_PRE_PD
  //vpd_history[1].fName= pidec_POST_PD
  //vpd_history[2].fName= pidec_POST_PD_BEAM_PROF
  //vpd_history[3].fName= pidec_POST_PD_RETRACE
  //------------
  //vsolution[0].fName= TDR_SOLUTION_BEAM_INITIAL
  //vsolution[1].fName= TDR_SOLUTION_DIAM
  //vsolution[2].fName= TDR_SOLUTION_HAD
  //vsolution[3].fName= MIN_SOLUTION_BEAM_INITIAL
  //vsolution[4].fName= MIN_SOLUTION_DIAM
  //vsolution[5].fName= MIN_SOLUTION_HAD
  //vsolution[6].fName= DIGI_POS_det1
  //vsolution[7].fName= DIGI_POS_det2
  //vsolution[8].fName= DIGI_POS_diamond

  TFile *f = TFile::Open(fin);
  TTree *t = (TTree*) f->Get("t");

  float p_ref = 0.65;

  t->SetAlias("pbeam","p[0]");
  t->SetAlias("thbeam","th[0]");
  t->SetAlias("phbeam","ph[0]");

  t->SetAlias("ppi","p[0]");
  t->SetAlias("dppi","dp[0]");

  //t->SetAlias("ppi","ppd[0]");
  //t->SetAlias("dppi",Form("100*(ppi-%3.5f)/%3.5f",p_ref,p_ref));

  t->SetAlias("x0", "x[0]");
  t->SetAlias("th0", "th[0]");
  t->SetAlias("y0", "y[0]");
  t->SetAlias("ph0", "ph[0]");

  t->SetAlias("x1", "xms[1]");
  t->SetAlias("y1", "yms[1]");

  t->SetAlias("x0er_d1ms", "x[0]-xms[2]");
  t->SetAlias("th0er_d1ms", "th[0]-thms[2]");
  t->SetAlias("y0er_d1ms", "y[0]-yms[2]");
  t->SetAlias("ph0er_d1ms", "ph[0]-phms[2]");

  //t->SetAlias("pmu", "ppd[1]");
  //t->SetAlias("dpmu",Form("100*(pmu-%3.5f)/%3.5f",p_ref,p_ref));
  //
  //t->SetAlias("pmub","ppd[2]");
  //t->SetAlias("pmur","ppd[3]");
  //
  //t->SetAlias("thmu","thpd[1]");
  //t->SetAlias("thmub","thpd[2]");
  //t->SetAlias("thmur","thpd[3]");
  //t->SetAlias("thmur_er","thmur-thmu");
  //
  //t->SetAlias("phmu","phpd[1]");
  //t->SetAlias("phmub","phpd[2]");
  //t->SetAlias("phmur","phpd[3]");
  //t->SetAlias("phmur_er","phmur-phmu");
  //
  //t->SetAlias("theta_mu_pi", "TMath::Hypot(thpd[1]-thpd[0],phpd[1]-phpd[0])");

  t->SetAlias("del0",Form("(abs((ppi-%3.5f)/%3.5f)-1.0)<0.001",p_ref,p_ref));
  for (int idel=-6; idel<7; ++idel) {
    TString alias = "del_p0";
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

  TCanvas *tc_d1ms_x1 = new TCanvas("tc_d1ms_x1","tc_d1ms_x1",1000,1000);
  tc_d1ms_x1->Divide(2,2);
  tc_d1ms_x1->cd(1);
  t->Draw("x1:x0er_d1ms","acc");
  tc_d1ms_x1->cd(2);
  t->Draw("x1:th0er_d1ms","acc");
  tc_d1ms_x1->cd(3);
  t->Draw("x1:y0er_d1ms","acc");
  tc_d1ms_x1->cd(4);
  t->Draw("x1:ph0er_d1ms","acc");

  TCanvas *tc_d1ms = new TCanvas("tc_d1ms","tc_d1ms",1000,1000);
  tc_d1ms->Divide(2,2);
  tc_d1ms->cd(1);
  t->Draw("y1:x0er_d1ms","acc");
  tc_d1ms->cd(2);
  t->Draw("y1:th0er_d1ms","acc");
  tc_d1ms->cd(3);
  t->Draw("y1:y0er_d1ms","acc");
  tc_d1ms->cd(4);
  t->Draw("y1:ph0er_d1ms","acc");

  //t->Draw("ppi-pmu", "del_p2");
  //t->Draw("dppi", "del_p2");
  //t->Draw("pmu:theta_mu_pi", "del_p5");

  //t->Draw("pmu:theta_mu_pi", "ipidec<33&&acc");


  //TH1F* h1=new TH1F("h1","h1",35,-0.5,34.5);
  //TH1F* h2=new TH1F("h2","h2",35,-0.5,34.5);
  //set_style_ana(h1,4,0,0);
  //set_style_ana(h2,4,0,0);
  //t->Project("h1","ipidec","dpidec<0");
  //t->Project("h2","ipidec","dpidec<0&&acc","same");
  ////t->Project("h1","ipidec","");
  ////t->Project("h2","ipidec","acc","same");
  //TH1F* hr = h2->Clone("hr");
  //hr->Divide(h1);
  //hr->Scale(100);
  //set_style_ana(hr,4,0,0);

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
  //log();
  //tc->cd(4);
  //t->Draw("pmu:theta_mu_pi>>h2d","","box");
  //gDirectory->ls();
  //TH2F *h2d = (TH2F*) gDirectory->Get("h2d");
  //cout << "h2d_ptr= " << h2d << " axX= " << h2d->GetXaxis() << endl;
  //h2d->GetYaxis()->SetRangeUser(0,1.2);

  //tc->cd(5);
  //t->Draw("pmu","dpidec<0&&acc");
  //
  //float r1 = h1->GetEntries()/t->GetEntries();
  //float r2 = h2->GetEntries()/h1->GetEntries();
  //
  //cout << "R1 = " << r1 << endl;
  //cout << "R2 = " << r2 << endl;
  //
  //tc->cd(6);
  //t->Draw("ppi","dpidec>0&&acc");

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
