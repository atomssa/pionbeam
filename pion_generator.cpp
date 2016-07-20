//#include "hphysicsconstants.h"
#include "HBeam.h"
//#include "hhistmap.h"

#include "TRandom.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"

#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <assert.h>

#include <stdlib.h>

using namespace std;

int main( int argc, const char **argv )
{

  if (argc!=7&&argc!=10) {
    cout << "Need minimum either 5 or 7 arguments. " << argc-1 << " given" << endl;
    cout << "./gen Nevt[10000] dphi[50mrad] dth[10mrad] xyr[0.5mm] p[GeV] dp[6\%] [Digitzie?{00,01,10/11} [MS{00,01,10/11}] [PIDEC{0,1}]" << endl;
    return -1;
  }

  cout << "argc= " << argc << " nEvents= " << atoi(argv[1]) << endl;

  const int nEvents=atoi(argv[1]);
  const double dph = atof(argv[2]);
  const double dth = atof(argv[3]);
  const double xyr = atof(argv[4]);
  const double pbeam = atof(argv[5]);
  const double dp_beam = atof(argv[6]);

  TParameter<double> tpar_beam_mom("ref_beam_mom", pbeam);

  const bool digi_pitrk = argc!=10 ? false : (atoi(argv[7])/10 == 1);
  const bool digi_diam = argc!=10 ? false : (atoi(argv[7])%10 == 1);
  const double diam_seg = 2.45;

  const bool do_ms_diam = true;
  const bool do_ms_pitrk1 = argc!=10?false:(atoi(argv[8])/10 == 1);
  const bool do_ms_pitrk2 = argc!=10?false:(atoi(argv[8])%10 == 1);
  cout << "do_ms_pitrk1 = " << (do_ms_pitrk1?1:0) << "  do_ms_pitrk2 = " << (do_ms_pitrk2?1:0) << endl;

  const bool do_pidec = argc!=10?false:(atoi(argv[9])==1);

  TString outfile=Form("output/pidec%d_ms%d%d%d_ds%d%d%s_xy%3.1fmm_dth%d_dph%d_dp%3.1f_p%4.2f.evt",
		       do_pidec?1:0, do_ms_pitrk1?1:0, do_ms_pitrk2?1:0, do_ms_diam?1:0, digi_pitrk?1:0, digi_diam?1:0,
		       (argc==11?Form("_seg%3.1f",diam_seg):""), xyr, (int)dth, (int)dph, dp_beam, pbeam);

  cout << "Output file name: " << outfile << endl;

  map<Int_t,TString> elementNames;

  elementNames[0]  = "in_Q1";
  elementNames[1]  = "out_Q1";
  elementNames[2]  = "in_Q2";
  elementNames[3]  = "out_Q2";
  elementNames[4]  = "out_Q2_+_0.6_m";
  elementNames[5]  = "in_FOPI";
  elementNames[6]  = "in_FOPI";
  elementNames[7]  = "out_FOPI";
  elementNames[8]  = "out_FOPI";
  elementNames[9]  = "in_dip_1";
  elementNames[10] = "out_dip_1";
  elementNames[11] = "in_Q3";
  elementNames[12] = "out_Q3";
  elementNames[13] = "in_Q4";
  elementNames[14] = "out_Q4";
  elementNames[15] = "inter_focal_plane";
  elementNames[16] = "detector_I"; // Position of Det1
  elementNames[17] = "in_Q5";
  elementNames[18] = "out_Q5";
  elementNames[19] = "in_Q6";
  elementNames[20] = "out_Q6";
  elementNames[21] = "in_dip_2";
  elementNames[22] = "out_dip_2";
  elementNames[23] = "in_Q7";
  elementNames[24] = "out_Q7";
  elementNames[25] = "inter_point_det_2"; // Position of Det1
  elementNames[26] = "in_Q8";
  elementNames[27] = "out_Q8";
  elementNames[28] = "in_Q9";
  elementNames[29] = "out_Q9";
  elementNames[30] = "inter_point";
  elementNames[31] = "HADES_target";
  elementNames[32] = "target";

  Int_t cuttup[] = {      // 0 = no, 1=radial, 2=box
    1,1,1,1,2, //5
    2,2,2,2,2, //10
    2,1,1,1,1, //15
    1,2,1,1,1, //20
    1,2,2,1,1, //25
    2,1,1,1,1, //30
    1,1,1
  };

  Double_t xcut[] =  {   // [mm]
    60.,60.,60.,60.,70.,   //5
    70.,70.,70.,70.,90.,   //10
    90.,60.,60.,60.,60.,   //15
    60.,50.,60.,60.,60.,   //20
    60.,90.,90.,60.,60.,   //25
    50.,60.,60.,60.,60.,   //30
    60.,60.,60.
  };

  Double_t ycut[] =        // [mm]
    {
      60.,60.,60.,60.,35., //5
      35.,35.,35.,35.,35., //10
      35.,60.,60.,60.,60., //15
      60.,50.,60.,60.,60., //20
      60.,35.,35.,60.,60., //25
      50.,60.,60.,60.,60., //30
      60.,60.,60.
    };

  HBeam pionbeam;
  //pionbeam.setBeam           (HPhysicsConstants::pid("pi-"),pbeam,60,60,0.0,0.0); // id, totl mom [GeV], beamtube x and y size, xoff,yoff
  // http://web-docs.gsi.de/~halo/docs/hydra/classDocumentation/dev/src/HPhysicsConstants.cxx.html#f0lQqE

  pionbeam.setBeam           (9, pbeam, 60, 60, 0.0, 0.0); // id, totl mom [GeV], beamtube x and y size, xoff,yoff
  pionbeam.setBeamProfile    (xyr,0.0);                // sigma [mm], flatradius [mm]
  pionbeam.setBeamResolution (dth/1000.,dph/1000.,dp_beam/100);          // dpx [rad],dpy [rad] ,dpz [relative]  [+-]

  if(!pionbeam.initBeamLine  ("par_files/pibeam_250414.data",32)) return 1;

  // Locations along z compatible with the parameter file
  pionbeam.addDetector("det1",    do_ms_pitrk1?0.03:0.0, 9.36, digi_pitrk, 0.78,     -17209, 2, 50., 50.);  // 17 du fichier
  pionbeam.addDetector("det2",    do_ms_pitrk2?0.03:0.0, 9.36, digi_pitrk, 0.78,     -5442,  2, 50., 50.); //
  //pionbeam.addDetector("plane",   0.0,                  1.0,  false,      0.0,      -1300.0,  1, 60., 60.);
  pionbeam.addDetector("diamond", do_ms_diam?0.03:0.0,    18.8, digi_diam,  3.0,      -170.0,   2, 7.1, 7.1);
  pionbeam.addDetector("pe_targ",   0.0,                  1.0,  false,      0.0,      0.0,      1, 6.0, 6.0);
  pionbeam.addDetector("pidec",  -1.0,                  1.0,  false,      0,        100,      1, 1.e9, 1.e9 /*no acceptacne cut*/);

  pionbeam.set_detector_names("det1", "det2", "diamond", "pe_targ");
  pionbeam.set_pion_decay_plane_name("pidec", do_pidec);

  vector<HBeamElement>& elements  = pionbeam.getElements();
  vector<HBeamElement>& detectors = pionbeam.getDetectors();

  if(elementNames.size() != elements.size()) {
    cout<<"Number of elements differs from name map! nelements = "<<elements.size()<<", map size = "<<elementNames.size()<<endl;
    return 1;
  }

  for(UInt_t i = 0 ; i < elements.size(); i++){
    elements[i].setElement(elementNames[i],cuttup[i],xcut[i],ycut[i]);
  }

  pionbeam.print_elements();
  pionbeam.print_detectors();

  //pionbeam.printBeamLine(kTRUE);   // kTRUE : print transform matrices in addition to name and distance, kFALSE : don't
  //pionbeam.printDetectors();
  //pionbeam.printBeamProfile();

  Bool_t bWriteDetectors = kFALSE; // true : write det1,det2 as beam particles

  TString histfile = outfile;
  histfile.ReplaceAll(".evt","");
  histfile.ReplaceAll(".root","");
  histfile += ".root";
  TFile *fout = TFile::Open(histfile,"RECREATE");
  fout->mkdir("xy");
  fout->mkdir("xth");
  fout->mkdir("yph");
  fout->mkdir("mom");
  fout->mkdir("del");
  fout->mkdir("mom_x");
  fout->mkdir("mom_y");
  fout->mkdir("mom_z");
  fout->mkdir("acceptance");
  fout->mkdir("dir");

  vector<HBeamParticle>& vhistory = pionbeam.newParticle();
  const unsigned int ndet = vhistory.size();
  vector<HBeamParticle>& vms_history = pionbeam.get_ms_history();
  const unsigned int nmspt = vms_history.size();
  vector<HBeamParticle>& vpidec_history = pionbeam.get_pidec_history();
  const unsigned int npidecpt = vpidec_history.size();
  vector<HBeamParticle>& vsolution = pionbeam.get_solution();
  const unsigned int nrec = vsolution.size();

  cout << "------------" << endl;
  for (unsigned int kk=0; kk<vhistory.size(); ++kk) cout << "history["<< kk << "].fName= " << vhistory[kk].fName << endl;
  cout << "------------" << endl;
  for (unsigned int kk=0; kk<vms_history.size(); ++kk) cout << "vms_history["<< kk << "].fName= " << vms_history[kk].fName << endl;
  cout << "------------" << endl;
  for (unsigned int kk=0; kk<vpidec_history.size(); ++kk) cout << "vpidec_history["<< kk << "].fName= " << vpidec_history[kk].fName << endl;
  cout << "------------" << endl;
  for (unsigned int kk=0; kk<vsolution.size(); ++kk) cout << "vsolution["<< kk << "].fName= " << vsolution[kk].fName << endl;

  TTree *tt = new TTree("t","t");

  int _acc=0;      tt->Branch("acc",&_acc,"acc/I");
  int _accDiam=0;      tt->Branch("accDiam",&_accDiam,"accDiam/I");
  int _accTarg=0;      tt->Branch("accTarg",&_accTarg,"accTarg/I");

  int _nel = (int) elements.size();     tt->Branch("nel",&_nel,"nel/I");
  int _acce[_nel];   tt->Branch("accElem[nel]",&_acce,"accElem[nel]/I");

  int _ndet=ndet;    tt->Branch("ndet",&_ndet,"ndet/I");
  float _p[ndet];    tt->Branch("p[ndet]",&_p,"p[ndet]/F");
  float _dp[ndet];   tt->Branch("dp[ndet]",&_dp,"dp[ndet]/F");
  float _x[ndet];    tt->Branch("x[ndet]",&_x,"x[ndet]/F");
  float _y[ndet];    tt->Branch("y[ndet]",&_y,"y[ndet]/F");
  float _th[ndet];   tt->Branch("th[ndet]",&_th,"th[ndet]/F");
  float _ph[ndet];   tt->Branch("ph[ndet]",&_ph,"ph[ndet]/F");

  int _nrec=nrec;     tt->Branch("nrec",&_nrec,"nrec/I");
  float _pr[nrec];    tt->Branch("pr[nrec]",&_pr,"pr[nrec]/F");
  float _xr[nrec];    tt->Branch("xr[nrec]",&_xr,"xr[nrec]/F");
  float _yr[nrec];    tt->Branch("yr[nrec]",&_yr,"yr[nrec]/F");
  float _thr[nrec];   tt->Branch("thr[nrec]",&_thr,"thr[nrec]/F");
  float _phr[nrec];   tt->Branch("phr[nrec]",&_phr,"phr[nrec]/F");
  int _rstat[nrec];   tt->Branch("rstat[nrec]",&_rstat,"rstat[nrec]/F");

  // msi == multiple scattering initial state
  int _nmspt = nmspt;  tt->Branch("nmspt",&_nmspt,"nmspt/I");
  float _xms[nmspt];  tt->Branch("xms[nmspt]",&_xms,"xms[nmspt]/F");
  float _thms[nmspt]; tt->Branch("thms[nmspt]",&_thms,"thms[nmspt]/F");
  float _yms[nmspt];  tt->Branch("yms[nmspt]",&_yms,"yms[nmspt]/F");
  float _phms[nmspt]; tt->Branch("phms[nmspt]",&_phms,"phms[nmspt]/F");
  float _pms[nmspt]; tt->Branch("pms[nmspt]",&_pms,"pms[nmspt]/F");

  // msi == multiple scattering initial state
  int _npidecpt = npidecpt;  tt->Branch("npidecpt",&_npidecpt,"npidecpt/I");
  int _i_pi_dec=0;     tt->Branch("ipidec",&_i_pi_dec,"ipidec/I");
  float _d_pi_dec;     tt->Branch("dpidec",&_d_pi_dec,"dpidec/F");
  float _xpidec[npidecpt];  tt->Branch("xpidec[npidecpt]",&_xpidec,"xpidec[npidecpt]/F");
  float _thpidec[npidecpt]; tt->Branch("thpidec[npidecpt]",&_thpidec,"thpidec[npidecpt]/F");
  float _ypidec[npidecpt];  tt->Branch("ypidec[npidecpt]",&_ypidec,"ypidec[npidecpt]/F");
  float _phpidec[npidecpt]; tt->Branch("phpidec[npidecpt]",&_phpidec,"phpidec[npidecpt]/F");
  float _ppidec[npidecpt]; tt->Branch("ppidec[npidecpt]",&_ppidec,"ppidec[npidecpt]/F");

  TH2F* h_xy[ndet], *h_xyAcc[ndet];
  TH2F* h_xth[ndet], *h_xthAcc[ndet];
  TH2F* h_yph[ndet], *h_yphAcc[ndet];
  TH1F* h_mom[ndet], *h_momAcc[ndet];
  TH1F* h_mom_x[ndet], *h_momAcc_x[ndet];
  TH1F* h_mom_y[ndet], *h_momAcc_y[ndet];
  TH1F* h_mom_z[ndet], *h_momAcc_z[ndet];
  TH1F* h_delta0, *h_delta0Acc, *h_delta0AccDiam, *h_delta0AccDiamTarg;

  int indexOut = 0;
  int indexDet[2] = {0};
  double r_pos = 100.0;
  double r_ang = 20;
  for (UInt_t idet=0; idet< ndet; ++idet) {
    const char *det_name = vhistory[idet].fName.Data();
    if (strcmp(det_name, "plane")==0) indexOut = idet;
    if (strcmp(det_name, "det1")==0) indexDet[0] = idet;
    if (strcmp(det_name, "det2")==0) indexDet[1] = idet;
    cout << Form("hxy_%s",det_name) << " " << Form("hxyAcc_%s",det_name) << endl;

    h_xy[idet] = new TH2F( Form("hxy_%s",det_name), Form("x vs. y [%s]; x[mm]; y[mm]; counts",det_name), 2000, -r_pos, r_pos, 2000, -r_pos, r_pos );
    h_xyAcc[idet] = new TH2F( Form("hxyAcc_%s",det_name), Form("x vs. y, accepted [%s]; x[mm]; y[mm]",det_name), 2000, -r_pos, r_pos, 2000, -r_pos, r_pos );

    h_xth[idet] = new TH2F( Form("hxth_%s",det_name), Form("#theta vs. x [%s]; x[mm]; #theta[mrad]; counts",det_name), 2000, -r_pos, r_pos, 2000, -r_ang, r_ang );
    h_xthAcc[idet] = new TH2F( Form("hxthAcc_%s",det_name), Form("#theta vs. x, accepted [%s]; x[mm]; #theta[mrad]; counts",det_name), 2000, -r_pos, r_pos, 2000, -r_ang, r_ang );

    h_yph[idet] = new TH2F( Form("hyph_%s",det_name), Form("#phi vs. x [%s]; y[mm]; #phi[mrad]; counts",det_name), 2000, -r_pos, r_pos, 2000, -r_ang, r_ang );
    h_yphAcc[idet] = new TH2F( Form("hyphAcc_%s",det_name), Form("#phi vs. y, accepted [%s]; y[mm]; #phi[mrad]; counts",det_name), 2000, -r_pos, r_pos, 2000, -r_ang, r_ang );

    h_mom[idet] = new TH1F( Form("hmom_%s",det_name), Form("momentum [%s]; p [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_momAcc[idet] = new TH1F( Form("hmomAcc_%s",det_name), Form("momentum, accepted [%s]; p [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_mom_x[idet] = new TH1F( Form("hmom_x_%s",det_name), Form("p_{x} [%s]; p_{x} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_momAcc_x[idet] = new TH1F( Form("hmomAcc_x_%s",det_name), Form("p_{x}, accepted [%s]; p_{x} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_mom_y[idet] = new TH1F( Form("hmom_y_%s",det_name), Form("p_{y} [%s]; p_{y} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_momAcc_y[idet] = new TH1F( Form("hmomAcc_y_%s",det_name), Form("p_{y}, accepted [%s]; p_{y} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_mom_z[idet] = new TH1F( Form("hmom_z_%s",det_name), Form("p_{z} [%s]; p_{z} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
    h_momAcc_z[idet] = new TH1F( Form("hmomAcc_z_%s",det_name), Form("p_{z}, accepted [%s]; p_{z} [GeV/c]; counts",det_name), 500, -0.10, 0.10 );
  }

  h_delta0 = new TH1F("hdelta0", "generated momentum offset ; delta [%]; counts", 500, -10., 10. );
  h_delta0Acc = new TH1F("hdelta0Acc", "generated momentum offset accepted in Q1-Q9; delta [%]; counts", 500, -10., 10. );
  h_delta0AccDiam = new TH1F("hdelta0AccDiam", "generated momentum offset accepted in Q1-Q9 + Diamond; delta [%]; counts", 500, -10., 10. );
  h_delta0AccDiamTarg = new TH1F("hdelta0AccDiamTarg", "generated momentum offset accepted in Q1-Q9 + Diamond + Target; delta [%]; counts", 500, -10., 10. );

  // These indexes are a bit different than the "history" vector indexes
  int indexDiam = 0;
  int indexTarg = 0;
  for (int idet=0; idet < (int)ndet; ++idet) {
    const char *det_name = detectors[idet].fName.Data();
    if (strcmp(det_name, "diamond")==0) indexDiam = idet;
    if (strcmp(det_name, "pe_targ")==0) indexTarg = idet;
  }

  const int nelt = elements.size();
  TH1F* hAccElmnt = new TH1F("hAcceptanceElement", "hAcceptanceElement", nelt, 0, nelt );
  TH1F* hAccCumul = new TH1F("hAcceptanceAccumulated", "hAcceptanceAccumulated", nelt, 0, nelt );
  TH2F* hxElement = new TH2F("hxElement","X Enveloppe vs. Element position;; x pos [mm]", nelt, 0, nelt, 200, -r_pos, r_pos);
  TH2F* hyElement = new TH2F("hyElement","Y Enveloppe vs. Element position;; y pos [mm]", nelt, 0, nelt, 200, -r_pos, r_pos);
  TH2F* hxElementTotalAcc = new TH2F("hxElementTotalAcc","hxElementTotalAcc;; x pos [mm]", nelt, 0, nelt, 200, -r_pos, r_pos);
  TH2F* hyElementTotalAcc = new TH2F("hyElementTotalAcc","hyElementTotalAcc;; y pos [mm]", nelt, 0, nelt, 200, -r_pos, r_pos);
  TH1F* hxDir = new TH1F("hxDir","hxDir; x dir [mrad]", 100, -100, 100);
  TH1F* hyDir = new TH1F("hyDir","hyDir; y dir [mrad]", 100, -100, 100);

  for (int ielt=0; ielt<nelt; ++ielt) {
    hAccElmnt->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
    hAccCumul->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
    hxElement->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
    hyElement->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
    hxElementTotalAcc->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
    hyElementTotalAcc->GetXaxis()->SetBinLabel(1+ielt, elementNames[ielt].Data());
  }

  hAccElmnt->GetXaxis()->LabelsOption("v");
  hAccCumul->GetXaxis()->LabelsOption("v");
  hxElement->GetXaxis()->LabelsOption("v");
  hyElement->GetXaxis()->LabelsOption("v");
  hxElementTotalAcc->GetXaxis()->LabelsOption("v");
  hyElementTotalAcc->GetXaxis()->LabelsOption("v");

  Int_t    flag               =  4 ; // if( getVERTEX && writeINDEX == 0)
  Double_t blast              =  0 ; // not available
  Int_t    nParticle          =  1 ;
  if(bWriteDetectors) nParticle = ndet-2; // no beam, no target
  Float_t  event_impact_param = -1.; // not available
  Double_t weight             =  1.; // no gen weights
  Int_t    parentID           = -1 ; // no parent
  Int_t    sourceID           = -1 ; // no source

  FILE* asciiFile = 0;
  // do not write ASCII file
  //asciiFile = fopen(outfile.Data(),"w");

  Int_t ctEvents   = 0;
  Int_t ctTotalTry = 0;
  Int_t maxTry     = 300; // max number of allowed try to prevent endless loop
  Bool_t Accepted  = kTRUE;
  //-----------------------------------------------------------
  for(Int_t ctEvt = 1; ctEvt <= nEvents; ctEvt++)
    {

      if (ctEvt % 1000 == 0) cout << "Event " << ctEvt << endl;

      if(asciiFile) fprintf(asciiFile," %i %i %f %f %i\n",ctEvt,nParticle,blast,event_impact_param,flag);
      ctEvents ++;

      if (ctTotalTry>=nEvents) break;

      Int_t nPart             = 0;
      Int_t ctTryParticle = 0;
      while (nPart < nParticle && ctTryParticle <= maxTry)
	{

	  if (ctTotalTry>=nEvents) break;

	  vector<TLorentzVector> vPion;

	  //-----------------------------------------------------------
	  // create particles
	  Bool_t reDo = kTRUE;
	  Int_t ctTry = 0;
	  while(reDo && ctTry < maxTry) {  // try as long a particle is accepted by the beam line

	    if (ctTotalTry>=nEvents) break;

	    _acc = 0;

	    ctTotalTry ++;
	    ctTry++;

	    vector<HBeamParticle>& vhistory = pionbeam.newParticle();
	    vector<HBeamParticle>& vms_history = pionbeam.get_ms_history();
	    vector<HBeamParticle>& vpidec_history = pionbeam.get_pidec_history();
	    vector<HBeamParticle>& vsolution = pionbeam.get_solution();

	    // make sure that the vectors have the same entries from event to event
	    // They should also have the entries in the same order, but that's more cotly to verify for every event
	    assert(ndet==vhistory.size());
	    assert(nmspt==vms_history.size());
	    assert(npidecpt==vpidec_history.size());
	    assert(nrec==vsolution.size());

	    vPion.clear();
	    for(UInt_t i = 0; i < vhistory.size(); i++){
	      TLorentzVector vpion;
	      vpion.SetXYZM(vhistory[i].fP.X(),vhistory[i].fP.Y(),vhistory[i].fP.Z(), 139.56995*0.001);
	      vPion.push_back(vpion);
	    }

	    //-----------------------------------------------------------
	    // qa plots

	    hxDir->Fill(TMath::ATan2(vPion[0].Px(),vPion[0].Pz())*1000);
	    hyDir->Fill(TMath::ATan2(vPion[0].Py(),vPion[0].Pz())*1000);
	    for(UInt_t i=0; i< vhistory.size();i++){
	      const double the_mrad = TMath::ATan2(vPion[i].Px(),vPion[i].Pz())*1000;
	      const double phi_mrad = TMath::ATan2(vPion[i].Py(),vPion[i].Pz())*1000;
	      _p[i] = vPion[i].P();
	      _dp[i] = (vPion[i].P() - pbeam)*100./pbeam;
	      _x[i] = vhistory[i].fPos.X();
	      _y[i] = vhistory[i].fPos.Y();
	      _th[i] = the_mrad;
	      _ph[i] = phi_mrad;
	      h_xy[i]->Fill(vhistory[i].fPos.X(),vhistory[i].fPos.Y());
	      h_xth[i]->Fill(vhistory[i].fPos.X(), the_mrad);
	      h_yph[i]->Fill(vhistory[i].fPos.Y(), phi_mrad);
	      h_mom[i]->Fill(vPion   [i].P());
	      h_mom_x[i]->Fill(vPion   [i].Px());
	      h_mom_y[i]->Fill(vPion   [i].Py());
	      h_mom_z[i]->Fill(vPion   [i].Pz());
	    }
	    h_delta0->Fill(_dp[0]);

	    for (unsigned int j=0; j<vms_history.size() && j<(unsigned int)nmspt; ++j) {
	      TLorentzVector ms_pion;
	      ms_pion.SetXYZM(vms_history[j].fP.X(),vms_history[j].fP.Y(),vms_history[j].fP.Z(), 139.56995*0.001);
	      const double ms_the_mrad = TMath::ATan2(ms_pion.Px(),ms_pion.Pz())*1000;
	      const double ms_phi_mrad = TMath::ATan2(ms_pion.Py(),ms_pion.Pz())*1000;
	      //_dpms[j] = (ms_pion.P() - pbeam)*100./pbeam;
	      _pms[j] = ms_pion.P();
	      _xms[j] = vms_history[j].fPos.X();
	      _thms[j] = ms_the_mrad;
	      _yms[j] = vms_history[j].fPos.Y();
	      _phms[j] = ms_phi_mrad;
	    }

	    for (unsigned int j=0; j<vpidec_history.size() && j<(unsigned int)_npidecpt; ++j) {
	      TLorentzVector pidec_muon;
	      pidec_muon.SetXYZM(vpidec_history[j].fP.X(),vpidec_history[j].fP.Y(),vpidec_history[j].fP.Z(), 139.56995*0.001);
	      const double pidec_the_mrad = TMath::ATan2(pidec_muon.Px(),pidec_muon.Pz())*1000;
	      const double pidec_phi_mrad = TMath::ATan2(pidec_muon.Py(),pidec_muon.Pz())*1000;
	      //_dppidec[j] = (pidec_pion.P() - pbeam)*100./pbeam;
	      _d_pi_dec = vpidec_history[j].fPos.Z();
	      _ppidec[j] = pidec_muon.P();
	      _xpidec[j] = vpidec_history[j].fPos.X();
	      _thpidec[j] = pidec_the_mrad;
	      _ypidec[j] = vpidec_history[j].fPos.Y();
	      _phpidec[j] = pidec_phi_mrad;
	    }

	    // retreive the element index where pion decay took place
	    for (_i_pi_dec=0; _i_pi_dec<nelt; ++_i_pi_dec) {
	      if (elements[_i_pi_dec].fPid!=9) break;
	    }

	    for (unsigned int j=0; j<vsolution.size() && j<(unsigned int)_nrec; ++j) {
	      TLorentzVector rec_pion;
	      rec_pion.SetXYZM(vsolution[j].fP.X(),vsolution[j].fP.Y(),vsolution[j].fP.Z(), 139.56995*0.001);
	      const double rec_the_mrad = TMath::ATan2(rec_pion.Px(),rec_pion.Pz())*1000;
	      const double rec_phi_mrad = TMath::ATan2(rec_pion.Py(),rec_pion.Pz())*1000;
	      //_pr[j] = (rec_pion.P() - pbeam)*100./pbeam;
	      _pr[j] = rec_pion.P();
	      _xr[j] = vsolution[j].fPos.X();
	      _thr[j] = rec_the_mrad;
	      _yr[j] = vsolution[j].fPos.Y();
	      _phr[j] = rec_phi_mrad;
	      _rstat[j] = vsolution[j].fStatus;
	    }

	    Accepted = kTRUE;
	    for(UInt_t i = 0 ; i < elements.size(); i++){
	      hxElement->Fill(i,elements[i].fout[0]);
	      hyElement->Fill(i,elements[i].fout[2]);
	      _acce[i] = elements[i].fAccepted;
	      if(!elements[i].fAccepted) Accepted = kFALSE;
	      if( i < (elements.size()-3) && Accepted) hAccCumul->Fill(i);
	    }
	    _acc = Accepted?1:0;

	    _accDiam = detectors[indexDiam].fAccepted;
	    _accTarg = detectors[indexTarg].fAccepted;

	    if(Accepted){
	      reDo = kFALSE;
	    }

	    tt->Fill();


	    //-----------------------------------------------------------
	    // qa plots
	    if(Accepted){
	      for(UInt_t i=0; i< vhistory.size();i++){

		const double the_mrad = TMath::ATan2(vPion[i].Px(),vPion[i].Pz())*1000;
		const double phi_mrad = TMath::ATan2(vPion[i].Py(),vPion[i].Pz())*1000;
		h_xyAcc[i]->Fill(vhistory[i].fPos.X(),vhistory[i].fPos.Y());
		h_xthAcc[i]->Fill(vhistory[i].fPos.X(), the_mrad);
		h_yphAcc[i]->Fill(vhistory[i].fPos.Y(), phi_mrad);
		h_momAcc[i]->Fill(vPion   [i].P());
		h_momAcc_x[i]->Fill(vPion   [i].Px());
		h_momAcc_y[i]->Fill(vPion   [i].Py());
		h_momAcc_z[i]->Fill(vPion   [i].Pz());
	      }

	      h_delta0Acc->Fill(_dp[0]);
	      if (_accDiam) h_delta0AccDiam->Fill(_dp[0]);
	      if (_accDiam && _accTarg) h_delta0AccDiamTarg->Fill(_dp[0]);

	      for(UInt_t i = 0 ; i < elements.size(); i++){
		hxElementTotalAcc->Fill(i,elements[i].fout[0]);
		hyElementTotalAcc->Fill(i,elements[i].fout[2]);
	      }
	    }
	    //-----------------------------------------------------------

	  } // end while redo

	  if(ctTry == maxTry && !_acc){ //
	    ctTryParticle++;
	    cout<<"evt "<<ctEvt<<" no accepted beam particle after "<<ctTry<<" attempts! redo one particle "<<ctTryParticle<<endl;

	    if(ctTryParticle == maxTry)  {
	      cout<<"ERROR : Too many attempts. Outputfile will be corrupted (empty or to small event)"<<endl;
	      if(asciiFile) fclose(asciiFile);
	      return 1;
	    }

	    continue;

	  }

	  if(vhistory.size() < 3) {
	    cout<<"no history available!"<<endl;
	    continue;
	  }

	  nPart++; // succesfull

	  //-----------------------------------------------------------
	  //-----------------------------------------------------------
	  // new particle
	  if(asciiFile )
	    {
	      //target first
	      UInt_t itarget = indexOut ;
	      fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %e\n",                 // 12 vars
		      vPion[itarget].E(),vPion[itarget].Px(),vPion[itarget].Py(),vPion[itarget].Pz(),
		      0.       ,vhistory[itarget].fPos.X(),vhistory[itarget].fPos.Y(),vhistory[itarget].fPos.Z(),
		      vhistory[itarget].fPid , sourceID , parentID, weight);
	      if(bWriteDetectors){
		UInt_t ndet = sizeof(indexDet)/sizeof(Int_t);
		for(UInt_t j = 0; j < ndet; j ++){   // write out wanted beam detector hits
		  Int_t i = j;
		  fprintf(asciiFile," %e %e %e %e %e %e %e %e %i %i %i %e\n",                 // 12 vars
			  vPion[i].E(),vPion[i].Px(),vPion[i].Py(),vPion[i].Pz(),
			  0.       ,vhistory[i].fPos.X(),vhistory[i].fPos.Y(),vhistory[i].fPos.Z(),
			  vhistory[i].fPid , sourceID , parentID, weight);
		}
	      }
	    }
	  //-----------------------------------------------------------

	} // end particle loop

    } // end event loop
  //-----------------------------------------------------------

  // Note: Since the last three elements change depending on the settings, we do not
  // fill their acceptance in the acceptance histogram to avoid confusion. They have to
  // be caluclated by using the variables addDiam and accTarg in the ntuple.
  for(UInt_t i = 0 ; i < elements.size()-3; i++){
    hAccElmnt->SetBinContent(i+1,elements[i].fCtAll > 0 ? 100 - (elements[i].fCtFail/(Double_t)elements[i].fCtAll)*100 : 100);
  }
  hAccCumul->Scale(1./(elements[0].fCtAll/100.));

  if(asciiFile) fclose(asciiFile);

  //pionbeam.printBeamLine(kFALSE);   // print acceptance statistic for detecors and elements
  //pionbeam.printDetectors();

  for (UInt_t idet=0; idet< ndet; ++idet) {
    fout->cd("xy");
    h_xy[idet]->Write();
    h_xyAcc[idet]->Write();
    fout->cd("xth");
    h_xth[idet]->Write();
    h_xthAcc[idet]->Write();
    fout->cd("yph");
    h_yph[idet]->Write();
    h_yphAcc[idet]->Write();
    fout->cd("mom");
    h_mom[idet]->Write();
    h_momAcc[idet]->Write();
    fout->cd("mom_x");
    h_mom_x[idet]->Write();
    h_momAcc_x[idet]->Write();
    fout->cd("mom_y");
    h_mom_y[idet]->Write();
    h_momAcc_y[idet]->Write();
    fout->cd("mom_z");
    h_mom_z[idet]->Write();
    h_momAcc_z[idet]->Write();
  }

  fout->cd("del");
  h_delta0->Write();
  h_delta0Acc->Write();
  h_delta0AccDiam->Write();
  h_delta0AccDiamTarg->Write();

  fout->cd("acceptance");
  hAccElmnt->Write();
  hAccCumul->Write();
  hxElement->Write();
  hyElement->Write();
  hxElementTotalAcc->Write();
  hyElementTotalAcc->Write();

  fout->cd("dir");
  hxDir->Write();
  hyDir->Write();

  fout->cd();
  tt->GetUserInfo()->Add(&tpar_beam_mom);
  tt->Write();

  fout->Close();

  cout<<"Total Acceptance : "<< (ctTotalTry > 0 ?  (ctEvents/(Double_t)ctTotalTry)*100 : 100)<<" %"<<endl;

  return 0;
}
