#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "string.h"

#include "utils.C"

void plot_retrace(TString fin="output/pidec0_ms001_ds11_xy0.5mm_dth10_dph50_dp6.0_p0.70.root"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *f = TFile::Open(fin);
  TTree *t = (TTree*) f->Get("t");

  // infer reference momentum
  double p_ref = get_ref_mom(t);
  cout << "p_ref= " << p_ref << endl;
  cout << "PiTr MS Width= " << get_ms_width_mrad(0.03, 9.36, p_ref) << endl;
  cout << "Diam MS Width= " << get_ms_width_mrad(0.03, 18.8, p_ref) << endl;

  int ims = fin.Index("ms") + 2;
  bool bmsd1 = fin(ims,1)(0) == '1';
  bool bmsd2 = fin(ims+1,1)(0) == '1';
  bool bmsdd = fin(ims+2,1)(0) == '1';
  int imsd1 = bmsd1?bmsd1*3:100;
  int imsd2 = bmsd2?bmsd1*3+bmsd2*3:100;
  int imsdd = bmsdd?bmsd1*3+bmsd2*3+bmsdd*3:100;

  cout << "Substr: " << fin(ims,3) << endl;
  cout << "ms: " << bmsd1 << " "<< bmsd2 << " "<< bmsdd << " " << endl;
  cout << "ims: " << imsd1 << " "<< imsd2 << " "<< imsdd << " " << endl;

  const char* st[4] = {"x", "y", "th", "ph"};
  const char* sd[3] = {"d1", "d2", "dd"};
  const char* sdf[3] = {"Pion Tracker 1", "Pion Tracker 2", "Diamond Detector"};
  t->SetAlias("del","dp[0]");
  for (int ii=0; ii < 4; ++ii) {
    t->SetAlias(Form("%sd1",st[ii]),Form("%s[1]",st[ii]));
    t->SetAlias(Form("%sd2",st[ii]),Form("%s[2]",st[ii]));
    t->SetAlias(Form("%sdd",st[ii]),Form("%s[3]",st[ii]));
    if (bmsd1) t->SetAlias(Form("%sd1_ret",st[ii]),Form("%sms[%d]",st[ii],imsd1));
    if (bmsd2) t->SetAlias(Form("%sd2_ret",st[ii]),Form("%sms[%d]",st[ii],imsd2));
    if (bmsdd) t->SetAlias(Form("%sdd_ret",st[ii]),Form("%sms[%d]",st[ii],imsdd));
    if (bmsd1) t->SetAlias(Form("%sd1_pre",st[ii]),Form("%sms[%d]",st[ii],imsd1-2));
    if (bmsd2) t->SetAlias(Form("%sd2_pre",st[ii]),Form("%sms[%d]",st[ii],imsd2-2));
    if (bmsdd) t->SetAlias(Form("%sdd_pre",st[ii]),Form("%sms[%d]",st[ii],imsdd-2));
  }

  TCanvas *tc_ret[3];
  TCanvas *tc_ret2d[3];
  TCanvas *tc_pre[3];
  TH1F *hmsw[3][4];

  for (int ii=0; ii < 3; ++ii) {

    if (ii==0&&!bmsd1) continue;
    if (ii==1&&!bmsd2) continue;
    if (ii==2&&!bmsdd) continue;

    tc_ret[ii] = new TCanvas(Form("tc_ret_%s",sd[ii]),Form("Retrace Error: %s",sdf[ii]));
    tc_ret[ii]->Divide(2,2);
    for (int jj=0; jj < 4; ++jj) {
      tc_ret[ii]->cd(1+jj);
      t->Draw(Form("%s%s-%s%s_ret",st[jj],sd[ii],st[jj],sd[ii]));
      gPad->SetLogy();
    }

    tc_ret2d[ii] = new TCanvas(Form("tc_ret2d_%s",sd[ii]),Form("Retrace Error: %s",sdf[ii]));
    tc_ret2d[ii]->Divide(2,2);
    for (int jj=0; jj < 4; ++jj) {
      tc_ret2d[ii]->cd(1+jj);
      t->Draw(Form("del:%s%s-%s%s_ret",st[jj],sd[ii],st[jj],sd[ii]));
    }

    tc_pre[ii] = new TCanvas(Form("tc_pre_%s",sd[ii]),Form("MS width: %s",sdf[ii]), 750, 10, 700, 500);
    tc_pre[ii]->Divide(2,2);
    for (int jj=0; jj < 4; ++jj) {
      tc_pre[ii]->cd(1+jj);
      t->Draw(Form("%s%s-%s%s_pre>>hmsw_%s_%s",st[jj],sd[ii],st[jj],sd[ii],st[jj],sd[ii]));
      if (jj>1) {
	hmsw[ii][jj] = (TH1F*)gDirectory->Get(Form("hmsw_%s_%s",st[jj],sd[ii]));
	hmsw[ii][jj]->Fit("gaus");
      }
    }

  }

}
