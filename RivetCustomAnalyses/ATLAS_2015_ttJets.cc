// -*- C++ -*-
// ATLAS 13 TeV tt+jets  differential  analysis in dilepton channel
// !!!!! Run with asetup 19.2.5 !!!!!

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Cuts.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

#include <bitset>

namespace Rivet {

  class ATLAS_2015_ttJets : public Analysis {
  public:
    
    ATLAS_2015_ttJets() : Analysis("ATLAS_2015_ttJets"){ }
    
    void init() {

      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT >= 1.0*MeV;
      Cut eta_lep = Cuts::abseta < 2.5;

      // Collect final state particles
      FinalState FS(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(FS);
      photons.acceptIdPair(PID::PHOTON);
      addProjection(photons, "photons");

      // Projection to find the electrons
      IdentifiedFinalState el_id(FS);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(false);
      addProjection(electrons, "electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1, Cuts::abseta< 2.5 && Cuts::pT>= 25.0*GeV, true, true);
      addProjection(dressedelectrons, "dressedelectrons");
      DressedLeptons vetodressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT >= 15.0*GeV, true, true);
      addProjection(vetodressedelectrons, "vetodressedelectrons");
      DressedLeptons fulldressedelectrons(photons, electrons, 0.1, eta_full, true, true);
      addProjection(fulldressedelectrons, "fulldressedelectrons");
      
      // Projection to find the muons
      IdentifiedFinalState mu_id(FS);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(false);
      addProjection(muons, "muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25.0*GeV, true, true);
      addProjection(dressedmuons, "dressedmuons");
      DressedLeptons vetodressedmuons(photons, muons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 15.0*GeV, true, true);
      addProjection(vetodressedmuons, "vetodressedmuons");
      DressedLeptons fulldressedmuons(photons, muons, 0.1, eta_full, true, true);
      addProjection(fulldressedmuons, "fulldressedmuons");

      // Projection to find neutrinos to exclude from jets
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(false);

      // Jet clustering
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(fulldressedelectrons);
      vfs.addVetoOnThisFinalState(fulldressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);			
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      addProjection(jets, "jets");

      eventsAll = 0.;
      eventsM1 = 0.;
      eventsM2 = 0.;
      eventsM3 = 0.;
      eventsM4 = 0.;

      eventsBefore = 0;
      eventsEMU = 0;
      eventsEMU2B = 0;
      eventsEMU2B1J = 0;

      // Book Histograms
      initializeHistos();

    }

    void initializeHistos() {
      // ========== Event Count ===============
      _h_Nevents           = bookHisto1D ("Nevents", 3, 0, 3);  

      _h_bjet_pt           = bookHisto1D (1,1,1);
      _h_2bjet_pt          = bookHisto1D (2,1,1);
      _h_ljet_pt           = bookHisto1D (3,1,1);

      // ========== Jet multiplicity ========== 
      _h_njet25 	  = bookHisto1D ("Mult25", 5, -0.5, 4.5);
      _h_njet40 	  = bookHisto1D ("Mult40", 4, -0.5, 3.5);
      _h_njet60 	  = bookHisto1D ("Mult60", 4, -0.5, 3.5);
      _h_njet80 	  = bookHisto1D ("Mult80", 4, -0.5, 3.5);

      // ========== Gap Fraction===============
      const double q0[] = {25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,310.};
      const std::vector<double> binEdgesQ0(q0,q0+sizeof(q0)/sizeof(q0[0]));

      const double qsum[] = {25.,35.,45.,55.,65.,75.,85.,95.,110.,130.,150.,170.,190.,210.,230.,250.,270.,300.,340.,380.,420.,500.,510.};
      const std::vector<double> binEdgesQSum(qsum,qsum+sizeof(qsum)/sizeof(qsum[0]));

      _h_pTQ01		= bookHisto1D ("pTQ01", binEdgesQ0);
      _h_pTQ02		= bookHisto1D ("pTQ02", binEdgesQ0);
      _h_pTQ03		= bookHisto1D ("pTQ03", binEdgesQ0);
      _h_pTQ04		= bookHisto1D ("pTQ04", binEdgesQ0);
      _h_pTQsum1	= bookHisto1D ("pTQsum1", binEdgesQSum);
      _h_pTQsum2	= bookHisto1D ("pTQsum2", binEdgesQSum);
      _h_pTQsum3	= bookHisto1D ("pTQsum3", binEdgesQSum);
      _h_pTQsum4	= bookHisto1D ("pTQsum4", binEdgesQSum);

      _h_pTMQ01		= bookHisto1D ("pTMQ01", binEdgesQ0);
      _h_pTMQ02		= bookHisto1D ("pTMQ02", binEdgesQ0);
      _h_pTMQ03		= bookHisto1D ("pTMQ03", binEdgesQ0);
      _h_pTMQ04		= bookHisto1D ("pTMQ04", binEdgesQ0);
      _h_pTMQsum1	= bookHisto1D ("pTMQsum1", binEdgesQSum);
      _h_pTMQsum2	= bookHisto1D ("pTMQsum2", binEdgesQSum);
      _h_pTMQsum3	= bookHisto1D ("pTMQsum3", binEdgesQSum);
      _h_pTMQsum4	= bookHisto1D ("pTMQsum4", binEdgesQSum);

      _h_gapFracQ01	= bookScatter2D ("GapFracQ01", true);
      _h_gapFracQ02	= bookScatter2D ("GapFracQ02", true);
      _h_gapFracQ03	= bookScatter2D ("GapFracQ03", true);
      _h_gapFracQ04	= bookScatter2D ("GapFracQ04", true);
      _h_gapFracQsum1	= bookScatter2D ("GapFracQsum1", true);
      _h_gapFracQsum2	= bookScatter2D ("GapFracQsum2", true);
      _h_gapFracQsum3	= bookScatter2D ("GapFracQsum3", true);
      _h_gapFracQsum4	= bookScatter2D ("GapFracQsum4", true);

      _h_gapFracMQ01	= bookScatter2D ("GapFracMQ01", true);
      _h_gapFracMQ02	= bookScatter2D ("GapFracMQ02", true);
      _h_gapFracMQ03	= bookScatter2D ("GapFracMQ03", true);
      _h_gapFracMQ04	= bookScatter2D ("GapFracMQ04", true);
      _h_gapFracMQsum1	= bookScatter2D ("GapFracMQsum1", true);
      _h_gapFracMQsum2	= bookScatter2D ("GapFracMQsum2", true);
      _h_gapFracMQsum3	= bookScatter2D ("GapFracMQsum3", true);
      _h_gapFracMQsum4	= bookScatter2D ("GapFracMQsum4", true);

    }

    void analyze(const Event& event) 
    {

      const double weight = event.weight();

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = applyProjection<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons = applyProjection<DressedLeptons>(event, "dressedmuons").dressedLeptons();
 
      Jets all_jets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
    
      // Event Selection: Remove electrons and muons if they overlap with jets
      vector<DressedLepton> _electrons, _muons;
	
      foreach (const DressedLepton& elec, electrons) {                                                                                                         
        bool hasOverlap = false;
        for (unsigned int j = 0; j < all_jets.size(); j++) {                                                                                                   
          const Jet& jet = all_jets[j]; 
          if (deltaR(jet, elec) < 0.4) hasOverlap = true;
        }
        if(!hasOverlap) _electrons.push_back(elec);
      }
	
      foreach (const DressedLepton& muon, muons) {
        bool hasOverlap = false;
        for (unsigned int j = 0; j < all_jets.size(); j++) {
          const Jet& jet = all_jets[j];
          if (deltaR(jet, muon) < 0.4) hasOverlap = true;
        }
        if(!hasOverlap) _muons.push_back(muon);
      }

      eventsBefore  += weight;
      _h_Nevents->fill(0.5, weight);

      // Event Selection: Require exactly one electron and exactly one muon
      if (_electrons.size() != 1 || _muons.size() != 1) return;  
      
      // Event Selection: Require opposite charge
      if (_electrons[0].charge() == _muons[0].charge()) return;     
  
      eventsEMU += weight;

      // Event Selection: Require at least two b-tagged jets. Collect vector of two highest pT b-jets
      Jets bjets, extrajets;
      for (unsigned int i = 0; i < all_jets.size(); i++) {
        Particles bTags = all_jets[i].particles();
        int count_hadrons = 0;
        foreach( Particle b, bTags ){
          if(b.containsBottom() && b.pT() > 5*GeV && deltaR(all_jets[i],b)<0.3)
            count_hadrons++;
        }
        bool b_tagged = false;
        if(count_hadrons > 0)
          b_tagged = true;

        if (bjets.size() < 2 && b_tagged) bjets += all_jets[i]; 
        else extrajets += all_jets[i];
      }
      if (bjets.size() < 2)  return;
      eventsAll += weight;
      eventsEMU2B += weight;
      _h_Nevents->fill(1.5, weight);

      // ========== pT Spectra ================
      double ljetpt = 0.;
      double bjetpt = 0.;
      double b2jetpt = 0.;

      if(extrajets.size()>0){
        _h_Nevents->fill(2.5, weight);
        eventsEMU2B1J += weight;
        if (extrajets[0].pt()>250.*GeV) ljetpt = 275*GeV;
        else ljetpt = extrajets[0].pt();
        _h_ljet_pt->fill(ljetpt, weight);  
      }
      if(bjets[0].pt()>250.*GeV) bjetpt = 275*GeV;
      else bjetpt = bjets[0].pt();

      if(bjets[1].pt()>150.*GeV) b2jetpt = 175*GeV;
      else b2jetpt = bjets[1].pt();

      _h_bjet_pt->fill(bjetpt, weight); 
      _h_2bjet_pt->fill(b2jetpt, weight); 

      // ========== Jet multiplicity ========== 
      unsigned int njet25 = 0;
      unsigned int njet40 = 0;
      unsigned int njet60 = 0;
      unsigned int njet80 = 0;

      for (unsigned int j = 0; j < extrajets.size(); j++) {
        if (extrajets[j].pt() > 25*GeV) njet25++;
	if (extrajets[j].pt() > 40*GeV) njet40++;
	if (extrajets[j].pt() > 60*GeV) njet60++;
	if (extrajets[j].pt() > 80*GeV) njet80++;	
      }

      // overflow hack for now...
      if(njet25>4) njet25 = 4;
      if(njet40>3) njet40 = 3;
      if(njet60>3) njet60 = 3;
      if(njet80>3) njet80 = 3;

      _h_njet25->fill(njet25, weight); 
      _h_njet40->fill(njet40, weight);
      _h_njet60->fill(njet60, weight);
      _h_njet80->fill(njet80, weight);

      // ========== Gap Fraction===============
      // extrajets is array of extra jets from 0 to n, highest pT first

      double pt1 = 0.;
      double pt2 = 0.;
      double pt3 = 0.;
      double pt4 = 0.;
      double ptsum1 = 0.;
      double ptsum2 = 0.;
      double ptsum3 = 0.;
      double ptsum4 = 0.;

      for (unsigned int i = 0; i < extrajets.size(); i++) {
        if(extrajets[i].absrap()<0.8 && extrajets[i].pt()>pt1) pt1 = extrajets[i].pt(); 
        else if(extrajets[i].absrap()>0.8 && extrajets[i].absrap()<1.5 && extrajets[i].pt()>pt2) pt2 = extrajets[i].pt(); 
        else if(extrajets[i].absrap()>1.5 && extrajets[i].absrap()<2.1 && extrajets[i].pt()>pt3) pt3 = extrajets[i].pt();
        if(extrajets[i].absrap()<2.1 && extrajets[i].pt()>pt4) pt4 = extrajets[i].pt(); 

        if(extrajets[i].absrap()<0.8) ptsum1 = ptsum1 + extrajets[i].pt(); 
        else if(extrajets[i].absrap()>0.8 && extrajets[i].absrap()<1.5) ptsum2 = ptsum2 + extrajets[i].pt(); 
        else if(extrajets[i].absrap()>1.5 && extrajets[i].absrap()<2.1) ptsum3 = ptsum3 + extrajets[i].pt();
        if(extrajets[i].absrap()<2.1) ptsum4 = ptsum4 + extrajets[i].pt(); 
      }
      if(pt1>305.*GeV) pt1 = 305.*GeV;
      if(pt2>305.*GeV) pt2 = 305.*GeV;
      if(pt3>305.*GeV) pt3 = 305.*GeV;
      if(pt4>305.*GeV) pt4 = 305.*GeV;

      if(ptsum1>505.*GeV) ptsum1 = 505.*GeV;
      if(ptsum2>505.*GeV) ptsum2 = 505.*GeV;
      if(ptsum3>505.*GeV) ptsum3 = 505.*GeV;
      if(ptsum4>505.*GeV) ptsum4 = 505.*GeV;

      _h_pTQ01->fill(pt1,weight);
      _h_pTQ02->fill(pt2,weight);
      _h_pTQ03->fill(pt3,weight);
      _h_pTQ04->fill(pt4,weight);
      _h_pTQsum1->fill(ptsum1,weight);
      _h_pTQsum2->fill(ptsum2,weight);
      _h_pTQsum3->fill(ptsum3,weight);
      _h_pTQsum4->fill(ptsum4,weight);

      // Invariant mass calculation / regions
      double Memubb = 0.;
      double Etot = _electrons[0].E()+_muons[0].E()+bjets[0].E()+bjets[1].E();
      double pe2 = pow(_electrons[0].pt()*cosh(_electrons[0].eta()),2);
      double pmu2 = pow(_muons[0].pt()*cosh(_muons[0].eta()),2);
      double pb12 = pow(bjets[0].pt()*cosh(bjets[0].eta()),2);
      double pb22 = pow(bjets[1].pt()*cosh(bjets[1].eta()),2);

      double pex = _electrons[0].pt()*cos(_electrons[0].phi());
      double pey = _electrons[0].pt()*sin(_electrons[0].phi());
      double pez = _electrons[0].pt()*sinh(_electrons[0].eta());
      double pmux = _muons[0].pt()*cos(_muons[0].phi());
      double pmuy = _muons[0].pt()*sin(_muons[0].phi());
      double pmuz = _muons[0].pt()*sinh(_muons[0].eta());
      double pb1x = bjets[0].pt()*cos(bjets[0].phi());
      double pb1y = bjets[0].pt()*sin(bjets[0].phi());
      double pb1z = bjets[0].pt()*sinh(bjets[0].eta());
      double pb2x = bjets[1].pt()*cos(bjets[1].phi());
      double pb2y = bjets[1].pt()*sin(bjets[1].phi());
      double pb2z = bjets[1].pt()*sinh(bjets[1].eta());

      double pepmu = pex*pmux + pey*pmuy + pez*pmuz;
      double pepb1 = pex*pb1x + pey*pb1y + pez*pb1z;
      double pepb2 = pex*pb2x + pey*pb2y + pez*pb2z;
      double pmupb1 = pmux*pb1x + pmuy*pb1y + pmuz*pb1z;
      double pmupb2 = pmux*pb2x + pmuy*pb2y + pmuz*pb2z;
      double pb1pb2 = pb1x*pb2x + pb1y*pb2y + pb1z*pb2z;

      Memubb = pow(pow(Etot,2)-(pe2+pmu2+pb12+pb22+2*(pepmu+pepb1+pepb2+pmupb1+pmupb2+pb1pb2)),0.5);
      if(Memubb<300.*GeV){
        eventsM1 += weight;
        _h_pTMQ01->fill(pt4,weight);
        _h_pTMQsum1->fill(ptsum4,weight);
      }
      if(Memubb>300.*GeV && Memubb<425.*GeV){
        eventsM2 += weight;
        _h_pTMQ02->fill(pt4,weight);
        _h_pTMQsum2->fill(ptsum4,weight);
      }
      if(Memubb>425.*GeV && Memubb<600.*GeV){
        eventsM3 += weight;
        _h_pTMQ03->fill(pt4,weight);
        _h_pTMQsum3->fill(ptsum4,weight);
      }
      if(Memubb>600.*GeV){
        eventsM4 += weight;
        _h_pTMQ04->fill(pt4,weight);
        _h_pTMQsum4->fill(ptsum4,weight);
      }
    }

    void finalize() {

      // Normalize to cross-section

      std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "        crossSection() = " << crossSection() << std::endl;
      std::cout	<< "        sumOfWeights() = " << sumOfWeights() << std::endl;
      std::cout	<< "+++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

      normalize(_h_ljet_pt);  
      normalize(_h_bjet_pt);  
      normalize(_h_2bjet_pt);  
      normalize(_h_njet25);  
      normalize(_h_njet40);  
      normalize(_h_njet60);  
      normalize(_h_njet80);  

      // Build gap fraction plots
      for (unsigned int j = 0; j < 22; j++) {
        if(j<18){
          double jetEvents1 = _h_pTQ01->integralRange(j,17)+_h_pTQ01->bin(17).sumW();
          double gap1 = (eventsAll - jetEvents1)/eventsAll;
          double uncert1 = sqrt(eventsAll*gap1*(1-gap1))/eventsAll;
          _h_gapFracQ01->point(j).setY(gap1,uncert1);
          double jetEvents2 = _h_pTQ02->integralRange(j,17)+_h_pTQ02->bin(17).sumW();
          double gap2 = (eventsAll - jetEvents2)/eventsAll;
          double uncert2 = sqrt(eventsAll*gap2*(1-gap2))/eventsAll;
          _h_gapFracQ02->point(j).setY(gap2,uncert2);
          double jetEvents3 = _h_pTQ03->integralRange(j,17)+_h_pTQ03->bin(17).sumW();
          double gap3 = (eventsAll - jetEvents3)/eventsAll;
          double uncert3 = sqrt(eventsAll*gap3*(1-gap3))/eventsAll;
          _h_gapFracQ03->point(j).setY(gap3,uncert3);
          double jetEvents4 = _h_pTQ04->integralRange(j,17)+_h_pTQ04->bin(17).sumW();
          double gap4 = (eventsAll - jetEvents4)/eventsAll;
          double uncert4 = sqrt(eventsAll*gap4*(1-gap4))/eventsAll;
          _h_gapFracQ04->point(j).setY(gap4,uncert4);

          double jetEventsM1 = _h_pTMQ01->integralRange(j,17)+_h_pTMQ01->bin(17).sumW();
          double gapM1 = (eventsM1 - jetEventsM1)/eventsM1;
          double uncertM1 = sqrt(eventsAll*gapM1*(1-gapM1))/eventsAll;
          _h_gapFracMQ01->point(j).setY(gapM1,uncertM1);
          double jetEventsM2 = _h_pTMQ02->integralRange(j,17)+_h_pTMQ02->bin(17).sumW();
          double gapM2 = (eventsM2 - jetEventsM2)/eventsM2;
          double uncertM2 = sqrt(eventsAll*gapM2*(1-gapM2))/eventsAll;
          _h_gapFracMQ02->point(j).setY(gapM2,uncertM2);
          double jetEventsM3 = _h_pTMQ03->integralRange(j,17)+_h_pTMQ03->bin(17).sumW();
          double gapM3 = (eventsM3 - jetEventsM3)/eventsM3;
          double uncertM3 = sqrt(eventsAll*gapM3*(1-gapM3))/eventsAll;
          _h_gapFracMQ03->point(j).setY(gapM3,uncertM3);
          double jetEventsM4 = _h_pTMQ04->integralRange(j,17)+_h_pTMQ04->bin(17).sumW();
          double gapM4 = (eventsM4 - jetEventsM4)/eventsM4;
          double uncertM4 = sqrt(eventsAll*gapM4*(1-gapM4))/eventsAll;
          _h_gapFracMQ04->point(j).setY(gapM4,uncertM4);
        }
        double jetEventsSum1 = _h_pTQsum1->integralRange(j,21)+_h_pTQsum1->bin(21).sumW();
        double gapSum1 = (eventsAll - jetEventsSum1)/eventsAll;
        double uncertSum1 = sqrt(eventsAll*gapSum1*(1-gapSum1))/eventsAll;
        _h_gapFracQsum1->point(j).setY(gapSum1,uncertSum1);
        double jetEventsSum2 = _h_pTQsum2->integralRange(j,21)+_h_pTQsum2->bin(21).sumW();
        double gapSum2 = (eventsAll - jetEventsSum2)/eventsAll;
        double uncertSum2 = sqrt(eventsAll*gapSum2*(1-gapSum2))/eventsAll;
        _h_gapFracQsum2->point(j).setY(gapSum2,uncertSum2);
        double jetEventsSum3 = _h_pTQsum3->integralRange(j,21)+_h_pTQsum3->bin(21).sumW();
        double gapSum3 = (eventsAll - jetEventsSum3)/eventsAll;
        double uncertSum3 = sqrt(eventsAll*gapSum3*(1-gapSum3))/eventsAll;
        _h_gapFracQsum3->point(j).setY(gapSum3,uncertSum3);
        double jetEventsSum4 = _h_pTQsum4->integralRange(j,21)+_h_pTQsum4->bin(21).sumW();
        double gapSum4 = (eventsAll - jetEventsSum4)/eventsAll;
        double uncertSum4 = sqrt(eventsAll*gapSum4*(1-gapSum4))/eventsAll;
        _h_gapFracQsum4->point(j).setY(gapSum4,uncertSum4);

        double jetEventsSumM1 = _h_pTMQsum1->integralRange(j,21)+_h_pTMQsum1->bin(21).sumW();
        double gapSumM1 = (eventsM1 - jetEventsSumM1)/eventsM1;
        double uncertSumM1 = sqrt(eventsM1*gapSumM1*(1-gapSumM1))/eventsM1;
        _h_gapFracMQsum1->point(j).setY(gapSumM1,uncertSumM1);
        double jetEventsSumM2 = _h_pTMQsum2->integralRange(j,21)+_h_pTMQsum2->bin(21).sumW();
        double gapSumM2 = (eventsM2 - jetEventsSumM2)/eventsM2;
        double uncertSumM2 = sqrt(eventsM2*gapSumM2*(1-gapSumM2))/eventsM2;
        _h_gapFracMQsum2->point(j).setY(gapSumM2,uncertSumM2);
        double jetEventsSumM3 = _h_pTMQsum3->integralRange(j,21)+_h_pTMQsum3->bin(21).sumW();
        double gapSumM3 = (eventsM3 - jetEventsSumM3)/eventsM3;
        double uncertSumM3 = sqrt(eventsM3*gapSumM3*(1-gapSumM3))/eventsM3;
        _h_gapFracMQsum3->point(j).setY(gapSumM3,uncertSumM3);
        double jetEventsSumM4 = _h_pTMQsum4->integralRange(j,21)+_h_pTMQsum4->bin(21).sumW();
        double gapSumM4 = (eventsM4 - jetEventsSumM4)/eventsM4;
        double uncertSumM4 = sqrt(eventsM4*gapSumM4*(1-gapSumM4))/eventsM4;
        _h_gapFracMQsum4->point(j).setY(gapSumM4,uncertSumM4);
      }
    }
  private:
   
    /// @name Histogram helper functions
    Histo1DPtr _h_Nevents;  

    Histo1DPtr _h_ljet_pt;
    Histo1DPtr _h_bjet_pt;
    Histo1DPtr _h_2bjet_pt;
    
    Histo1DPtr _h_njet25;
    Histo1DPtr _h_njet40;
    Histo1DPtr _h_njet60;
    Histo1DPtr _h_njet80;

    Histo1DPtr _h_pTQ01;
    Histo1DPtr _h_pTQ02;
    Histo1DPtr _h_pTQ03;
    Histo1DPtr _h_pTQ04;
    Histo1DPtr _h_pTQsum1;
    Histo1DPtr _h_pTQsum2;
    Histo1DPtr _h_pTQsum3;
    Histo1DPtr _h_pTQsum4;

    Histo1DPtr _h_pTMQ01;
    Histo1DPtr _h_pTMQ02;
    Histo1DPtr _h_pTMQ03;
    Histo1DPtr _h_pTMQ04;
    Histo1DPtr _h_pTMQsum1;
    Histo1DPtr _h_pTMQsum2;
    Histo1DPtr _h_pTMQsum3;
    Histo1DPtr _h_pTMQsum4;

    Scatter2DPtr _h_gapFracQ01;
    Scatter2DPtr _h_gapFracQ02;
    Scatter2DPtr _h_gapFracQ03;
    Scatter2DPtr _h_gapFracQ04;
    Scatter2DPtr _h_gapFracQsum1;
    Scatter2DPtr _h_gapFracQsum2;
    Scatter2DPtr _h_gapFracQsum3;
    Scatter2DPtr _h_gapFracQsum4;

    Scatter2DPtr _h_gapFracMQ01;
    Scatter2DPtr _h_gapFracMQ02;
    Scatter2DPtr _h_gapFracMQ03;
    Scatter2DPtr _h_gapFracMQ04;
    Scatter2DPtr _h_gapFracMQsum1;
    Scatter2DPtr _h_gapFracMQsum2;
    Scatter2DPtr _h_gapFracMQsum3;
    Scatter2DPtr _h_gapFracMQsum4;

    double eventsAll;
    double eventsM1;
    double eventsM2;
    double eventsM3;
    double eventsM4;

    // Cut flow check
    double eventsBefore;
    double eventsEMU;
    double eventsEMU2B;
    double eventsEMU2B1J;

  };

  // Declare the class as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_ttJets);

}

