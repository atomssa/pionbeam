#include "TTree.h"
#include "TParameter.h"
#include "TLorentzVector.h"

double get_ref_mom(TTree *t) {
  double p_ref = -9999.;
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
  }
  return p_ref;
}

double get_ms_width_mrad(double t, double x0, double mom) {
  TLorentzVector p = TLorentzVector(TVector3(0,0,mom),TMath::Hypot(mom,0.134));
  return 13.6*TMath::Sqrt(t/x0)/p.Vect().Mag()/p.Beta();
}
