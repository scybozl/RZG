// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/Vector4.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  struct ATLAS_sherpa_parton_Plots { };



  /// Top pair production with central jet veto
  class ATLAS_sherpa_parton : public Analysis {

  #include "NLOHisto1D.cc"


  public:

    /// Constructor
    ATLAS_sherpa_parton()
      : Analysis("ATLAS_sherpa_parton"),

    _jet_ntag(0)
    //_hMap(),
    //_histLimit(6)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      addProjection(fs, "ALL_FS");

      addProjection(MissingMomentum(fs), "MissingET");

/*
      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      elec_fs.acceptIdPair(PID::ELECTRON);
      addProjection(elec_fs, "ELEC_FS");


      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      addProjection(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 5.0);
      neutrino_fs.acceptNeutrinos();
      addProjection(neutrino_fs, "NEUTRINO_FS");
*/
      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
    //  VetoedFinalState jet_input(fs);
   //   jet_input.vetoNeutrinos();
  //    jet_input.addVetoPairId(PID::MUON);
 //     jet_input.addVetoPairId(PID::ELECTRON);
//      addProjection(jet_input, "JET_INPUT");

      // Get the jets
//      FastJets jets(jet_input, FastJets::ANTIKT, 0.5);
      FastJets jets(fs, FastJets::ANTIKT, 0.5);
      addProjection(jets, "JETS");
/*
      for (unsigned int ihist = 0; ihist < _histLimit ; ihist++) {
        const unsigned int threshLimit = _thresholdLimit(ihist);
        for (unsigned int ithres = 0; ithres < threshLimit; ithres++) {
          _histogram(ihist, ithres); // Create all histograms
        }
      }
*/

	// Create histograms
	_histPt1 = bookNLOHisto1D("Pt1", 25, 0, 400);
	_histEta1 = bookNLOHisto1D("Eta1", 20, -2.5, 2.5);
	_histPhiemu = bookNLOHisto1D("Phiemu", 20, 0, 2*PI);
	_histdeltaRb = bookNLOHisto1D("deltaRb", 20, 0, 5);
	_histdeltaRl = bookNLOHisto1D("deltaRl", 20, 0, 5);
	_histPtMiss = bookNLOHisto1D("PtMiss", 25, 0, 400);
	_histHT = bookNLOHisto1D("HT", 20, 0, 1200);
	_histMlb = bookNLOHisto1D("Mlb", 25, 0, 200);
	_histWWphi = bookNLOHisto1D("WWphi", 25, 0, 2*PI);
        _histWPt = bookNLOHisto1D("WPt", 50, 0, 1000);
	_histWEta = bookNLOHisto1D("WEta", 25, -6, 6);
	_histWPhi = bookNLOHisto1D("WPhi", 25, 0, 2*PI);
	_histMlb2 = bookNLOHisto1D("Mlb2", 25, 0, 200);
	_histMlbBare = bookNLOHisto1D("MlbBare", 25, 0, 200);
	_histHT2 = bookNLOHisto1D("HT2", 20, 0, 1200);
	_histHT3 = bookNLOHisto1D("HT3", 20, 0, 1200);

	_histScales = bookHisto2D("mu1mu2", 60, 200, 800, 60, 200, 800);
	_histmu1PtT = bookScatter2D("mu1PtT");
	_histmu2PtT = bookScatter2D("mu2PtT");
	_histmu1MT = bookScatter2D("mu1MT");
	_histmu2MT = bookScatter2D("mu2MT");
	_histttbarpt = bookHisto1D("ttbarpt", 
	_histttbarmass = bookHisto1D("ttbarmass",


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double weight = event.weight();
//      const double normWW = 1.0/3.9465;
//      weight *= normWW;
/*
      /// Get the various sets of final state particles
      Particles elecFS = applyProjection<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      Particles muonFS = applyProjection<IdentifiedFinalState>(event, "MUON_FS").particlesByPt(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);

      const Particles& neutrinoFS = applyProjection<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();
*/
      // Get all jets with pT > 30 GeV
      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt();

      vector<const Jet*> hardJets;
      foreach (const Jet& j, jets) {
	if (j.pT()>30.0*GeV) hardJets.push_back(&j);
      }
      if (hardJets.size()<2) cout << "H";

      // Keep any jets that pass the initial rapidity cut
      vector<const Jet*> central_jets;
      foreach(const Jet* j, hardJets) {
        if (j->absrap() < 2.5) central_jets.push_back(j);
//        else cout << " 0 ";
      }
      if (central_jets.size()<2) cout << "R";
//      cout << " " << central_jets.size() << " ";


/*
      // Get b hadrons with pT > 5 GeV
      /// @todo This is a hack -- replace with UnstableFinalState
      vector<GenParticle const *> B_hadrons;
      vector<GenParticle const *> allParticles = particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        const GenParticle* p = allParticles[i];
        if (!PID::isHadron(p->pdg_id()) || !PID::hasBottom(p->pdg_id())) continue;
        //if (p->momentum().perp() < 5*GeV) continue;
        B_hadrons.push_back(p);
      }
*/
      // For each of the good jets, check whether any are b-jets (via dR matching)
      vector<const Jet*> b_jets;
      foreach (const Jet* j, central_jets) {
        //bool isbJet = false;
        //foreach (const GenParticle* b, B_hadrons) {
          //if (deltaR(j->momentum(), FourMomentum(b->momentum())) < 0.35) isbJet = true;
        //}
	//if(isbJet) b_jets.push_back(j);
//        b_jets.push_back(j);
//        if (j->containsBottom()) b_jets.push_back(j);
	  int nb = 0;
	  foreach (const Particle& p, j->particles()) {
		if (p.pdgId() == 5) ++nb;
		if (p.pdgId() == -5) --nb;
	  }
          if (nb != 0) b_jets.push_back(j);
//        else if (b_jets.size() <= 1) cout << " 1 ";
      }
      if (b_jets.size()<2) cout << "B";
//	foreach (const Jet* j, b_jets) cout << j->pT() << " ";

/*
      // Get the MET by taking the vector sum of all neutrinos
      /// @todo Use MissingMomentum instead?
      double MET = 0;
      FourMomentum p_MET;
      foreach (const Particle& p, neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();
*/
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingET");
      double MET = met.missingMomentum().pT();
//      cout << MET << " ";
      bool passed_emu = false;
      bool passed_Mlb = false;

      ParticleVector diLeptons;
      ParticleVector bBares;
      vector<GenParticle const *> allParticles = particles(event.genEvent());
      for (size_t i = 0; i < allParticles.size(); i++) {
        const GenParticle* p = allParticles[i];
        if (abs(p->pdg_id())==11 || abs(p->pdg_id())==13) diLeptons.push_back(Particle(*p));
	if (abs(p->pdg_id())==5) bBares.push_back(Particle(*p));
        //if (p->momentum().perp() < 5*GeV) continue;
        //if (Particle(*p).hasAncestor(5) && Particle(*p).pT()>maxpTHad){
        // B_hadrons = *p;
        //maxpTHad = Particle(*p).pT();
        //foundB=true;
      }

      if (diLeptons.size() < 2) cout << "Less than 2 leptons found!\n";

      ParticleVector Top;
      ParticleVector antiTop;
      bool foundTop=false;
      bool foundAntiTop=false;
      for (size_t i = 0; i < allParticles.size() && !(foundTop && foundAntiTop); i++) {
	const GenParticle* p = allParticles[i];
        if (p->pdg_id()==6 && foundTop==false) {Top.push_back(Particle(*p)); foundTop=true;}
	if (p->pdg_id()==-6 && foundAntiTop==false) {antiTop.push_back(Particle(*p)); foundAntiTop=true;}
      }
      if (Top.size() + antiTop.size() != 2) cout << "Less than 2 tops found!\n";

      // Get the electrons and muons coming from the W bosons
/*      ParticleVector diLeptons;
      ParticleVector diBosons;
      bool diLeptonsFound = false;
      bool elecFound = false;
      bool muonFound = false;
      bool diBosonsFound = false;
      bool WPFound = false; 
      bool WMFound = false;
      // Generate the complete event
      const HepMC::GenEvent *ge = event.genEvent();

      // Loop over all vertices, and find the W bosons
      for (HepMC::GenEvent::vertex_const_iterator iv = ge->vertices_begin(); iv != ge->vertices_end() && (diLeptonsFound ==  false || diBosonsFound == false); ++iv) {
         for (HepMC::GenVertex::particles_in_const_iterator ip1 = (*iv)->particles_in_const_begin(); ip1 != (*iv)->particles_in_const_end(); ++ip1) {

		// If one of the particles coming in the vertex is a W boson, match the leptons
		//if (abs((*ip1)->pdg_id()) == 24) {
			for (HepMC::GenVertex::particles_out_const_iterator ip2 = (*iv)->particles_out_const_begin(); ip2 != (*iv)->particles_out_const_end(); ++ip2) {
				if (((*ip2)->pdg_id() == -11 && !elecFound) || ((*ip2)->pdg_id() == 13 && !muonFound)) {
				  diLeptons.push_back(Particle(*ip2));
				  if ((*ip2)->pdg_id() == -11) elecFound = true;
				  if ((*ip2)->pdg_id() == 13) muonFound = true;
				//GenVertex* prodVtx = Particle(*ip2).production_vertex();
				//foreach (const GenParticle* ancestor, particles(*iv, HepMC::ancestors)) {
     					// cout << " " << ancestor->pdg_id() << " ";
    				//}
			//	cout << "\n";
				  if (diLeptons.size() == 2) diLeptonsFound = true;
				}
				if (((*ip2)->pdg_id() == 24 && !WPFound) || ((*ip2)->pdg_id() == -24 && !WMFound)) {
				  diBosons.push_back(Particle(*ip2));
				  if ((*ip2)->pdg_id() == 24) WPFound = true;
				  if ((*ip2)->pdg_id() == -24) WMFound = true;
				  if (diBosons.size() == 2) diBosonsFound = true;
				}
			}
		//}
	}
      }
	cout << diBosonsFound;
//      cout << diLeptons[0].pT() << " " << diLeptons[1].pT() << " ";
*/
/*
      //(event.genEvent())->print();
      if (elecFS.empty() || muonFS.empty()) {
         foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
	   //cout << p->pdg_id() << " ";
           if (p->pdg_id() != -11 && p->pdg_id() != 13) continue;
           //const GenVertex* pv = p->production_vertex();
           bool passed = true;
           //if (pv) {
           //  foreach (const GenParticle* pp, particles_in(pv)) {
           //    if ( p->pdg_id() == pp->pdg_id() ) {
           //      passed = false;
           //      break;
           //    }
           //  }
           //}
      	if (passed && p->pdg_id() == -11 && Particle(*p).hasAncestor(24)) {
	 if (!elecFS.empty()) {
		if (Particle(*p).pT() > elecFS[0].pT()) elecFS[0] = Particle(*p);
	 }
	 else elecFS.push_back(Particle(*p));
	}
      	if (passed && p->pdg_id() == 13 && Particle(*p).hasAncestor(-24)) {
	 if (!muonFS.empty()) {
		if (Particle(*p).pT() > muonFS[0].pT()) muonFS[0] = Particle(*p);
	 }
	 else muonFS.push_back(Particle(*p));
	}
      }
      }
*/
      //if (elecFS.size() == 1 && muonFS.size() == 1) {
        // With the desired charge signs
        //if (charge(elecFS[0]) >= 0 && charge(muonFS[0]) <= 0) {
          // Calculate HT: scalar sum of the pTs of the leptons and all good jets
          double HT = 0;
	  double HT2 = 0;
	  double HT3 = 0;
	  if(diLeptons.size() == 2) {
          HT += diLeptons[0].pT();
          HT += diLeptons[1].pT();
	  HT2 += diLeptons[0].pT();
          HT2 += diLeptons[1].pT();
	  HT3 += diLeptons[0].pT();
          HT3 += diLeptons[1].pT();
	  }
          foreach (const Jet* j, central_jets) HT += j->pT();
	  foreach (const Jet* j, hardJets) HT2+= j->pT();
	  foreach (const Jet* j, b_jets) HT3 += j->pT();
//          HT += MET;

//	  for (size_t i = 0; i < allParticles.size(); i++) {
 //           const GenParticle* p = allParticles[i];
//	    HT += fabs(Particle(*p).pT());
//	  }
//	  for (size_t i = 0; i < allParticles.size(); i++) {
//            const GenParticle* p = allParticles[i];
//            HT2 += Particle(*p).pT();
//          }
//	  if(diLeptons.size() ==2) {
//          HT3 += diLeptons[0].pT();
//          HT3 += diLeptons[1].pT();
//	  }
//          foreach (const Jet* j, central_jets) HT3 += fabs(j->pT());
//	  HT3 += MET;

	  double mu1=0;
	  double mu2=0;

	  double ttbarmass = 0;
	  double ttbarpt = 0;
 
	  if (Top.size() + antiTop.size() == 2){
	    mu1=(Top[0].momentum()+antiTop[0].momentum()).mass();
            double mT2 = Top[0].mass2()+Top[0].pT2() + antiTop[0].mass2()+antiTop[0].pT2();
	    if (mT2 >=0) mu2=sqrt(mT2);
	    else cout << "mT2 < 0!\n";

	    ttbarmass = (Top[0].momentum()+antiTop[0].momentum()).mass();
	    ttbarpt = (Top[0].momentum()+antiTop[0].momentum()).pT();
	  }

          // Keep events with HT > 130 GeV
//          if (HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            if (diLeptons.size() == 2) {
		 if (MET > 20.0*GeV) {
			if(diLeptons[0].pT() > 20.*GeV && diLeptons[1].pT() > 20.*GeV && fabs(diLeptons[0].eta()) < 2.5 && fabs(diLeptons[1].eta()) < 2.5) {
			  if (b_jets.size() > 1) {
			  passed_emu = true;

			// Additional requirements for the invariant mass event selection
			if(diLeptons.size() == 2) {

			if(deltaR(b_jets[0]->momentum(), diLeptons[0].momentum()) >= 0.4 && deltaR(b_jets[0]->momentum(), diLeptons[1].momentum()) >= 0.4
			&& deltaR(b_jets[1]->momentum(), diLeptons[0].momentum()) >= 0.4 && deltaR(b_jets[1]->momentum(), diLeptons[1].momentum()) >= 0.4) {
				if(HT3>=130*GeV) {
					passed_Mlb = true;
				}
			}
			}
			}
			}
			else cout << "RAP";
		}
		else cout << "MET";
            }
//	    else cout << "L";
//	}
//	else cout << "HT";
	//}
	//}

      if (Top.size() + antiTop.size() == 2 ){
	 _histScales->fill(mu1, mu2, 1);
	 _histmu1Pt->fill(mu1

      if (passed_emu == true) {

/*
	double pTcut[4] = {30.,40.,60.,80.};
	// Count the jet multiplicity for 30, 40, 60 and 80GeV
      	unsigned int ithres, jet_n[4] = {0,0,0,0};
      	for (ithres = 0; ithres < 4; ithres++) {
		foreach(const Jet* j, central_jets) {
        		if (j->pT() > pTcut[ithres]) jet_n[ithres]++;
		}
      	unsigned int ncutoff[4] = {8,7,6,5};
      	if (jet_n[ithres] > ncutoff[ithres]) jet_n[ithres] = ncutoff[ithres];
      }


      // Fill histograms
      for (unsigned int ihist = 0; ihist < 6; ihist++) {
        unsigned int threshLimit = _thresholdLimit(ihist);
        for (ithres = 0; ithres < threshLimit; ithres++) {
          if (jet_n[ithres] < 2) continue; // 2 or more jets for ljets
          // Fill
          if (ihist == 0) _histogram(ihist, ithres)->fill(jet_n[ithres], event); // njets
          else if (ihist == 1) _histogram(ihist, ithres)->fill(central_jets[0]->pT(), event); // leading jet pT
          else if (ihist == 2) _histogram(ihist, ithres)->fill(central_jets[1]->pT(), event); // 2nd jet pT
          else if (ihist == 3 && jet_n[ithres] >= 3) _histogram(ihist, ithres)->fill(central_jets[2]->pT(), event); // 3rd jet pT
          else if (ihist == 4 && jet_n[ithres] >= 4) _histogram(ihist, ithres)->fill(central_jets[3]->pT(), event); // 4th jet pT
          else if (ihist == 5 && jet_n[ithres] >= 5) _histogram(ihist, ithres)->fill(central_jets[4]->pT(), event); // 5th jet pT
        }
      }
*/

	_histPt1->fill(b_jets[0]->pT(), event);
	_histEta1->fill(b_jets[0]->eta(), event);
	_histdeltaRb->fill(deltaR(b_jets[0]->momentum(), b_jets[1]->momentum()), event);
	if(diLeptons.size() == 2) {
		//_histPhiemu->fill(deltaPhi(elecFS[0], muonFS[0]), event);
		double emuPhi = diLeptons[0].momentum().phi() - diLeptons[1].momentum().phi();
		if (emuPhi > 2*PI) emuPhi -= 2*PI;
		if (emuPhi < 0) emuPhi += 2*PI;
		_histPhiemu->fill(emuPhi, event);
		//_histPhiemu->fill(mapAngle0To2Pi(diLeptons[0].momentum().phi() - diLeptons[1].momentum().phi()), event);
		_histdeltaRl->fill(deltaR(diLeptons[0], diLeptons[1]), event);
	}
	_histPtMiss->fill(MET, event);
	_histHT->fill(HT, event);
	_histHT2->fill(HT2, event);
	_histHT3->fill(HT3, event);

	// Without the additional requirement from DDKP

        if((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass()
                > (diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass()) {
                        _histMlb2->fill(((diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass())/2, event);
                }
                else _histMlb2->fill(((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass())/2, event);

	// Old definition

	if(passed_Mlb == true) {
		if((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass()
		> (diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass()) {
			_histMlb->fill(((diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass())/2, event);
		}
		else _histMlb->fill(((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass())/2, event);
	}

	// Definition without jets, simply with the bare bs code

	if((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass()
                > (diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass()) {
                        _histMlbBare->fill(((diLeptons[0].momentum() + b_jets[1]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[0]->momentum()).mass())/2, event);
                }
                else _histMlbBare->fill(((diLeptons[0].momentum() + b_jets[0]->momentum()).mass() + (diLeptons[1].momentum() + b_jets[1]->momentum()).mass())/2, event);

      }
      
/*      if(diBosons.size() == 2) {
		double WWPhi = diBosons[0].momentum().phi() - diBosons[1].momentum().phi();
		if(WWPhi > 2*PI) WWPhi -= 2*PI;
		if(WWPhi < 0) WWPhi += 2*PI;
		_histWWphi->fill(WWPhi, event);
		_histWPt->fill(diBosons[0].pT(), event);
		_histWPt->fill(diBosons[1].pT(), event);
		_histWPhi->fill(diBosons[0].phi(), event);
		_histWPhi->fill(diBosons[1].phi(), event);
		_histWEta->fill(diBosons[0].eta(), event);
		_histWEta->fill(diBosons[1].eta(), event);
      }
*/    }

    /// Normalise histograms etc., after the run
    void finalize() {
//	normalize(_histPt1);
//	normalize(_histEta1);
//	normalize(_histPhiemu);
//	normalize(_histdeltaRb);
//	normalize(_histdeltaRl);
//	normalize(_histPtMiss);
//	normalize(_histHT);
//	normalize(_histMlb);

	_histPt1->finalize();
	_histEta1->finalize();
	_histPhiemu->finalize();
	_histdeltaRb->finalize();
	_histdeltaRl->finalize();
	_histPtMiss->finalize();
	_histHT->finalize();
	_histHT2->finalize();
	_histHT3->finalize();
	_histMlb->finalize();
	_histWWphi->finalize();
	_histWPt->finalize();
	_histWPhi->finalize();
	_histWEta->finalize();
	_histMlb2->finalize();
	_histMlbBare->finalize();

	const double norm = crossSection()/sumOfWeights();
	scale(_histPt1, norm);
	scale(_histEta1, norm);
	scale(_histPhiemu, norm);
	scale(_histdeltaRb, norm);
	scale(_histdeltaRl, norm);
	scale(_histPtMiss, norm);
	scale(_histHT, norm);
	scale(_histHT2, norm);
	scale(_histHT3, norm);
	scale(_histMlb, norm);
	scale(_histWWphi, norm);
	scale(_histWPt, norm);
	scale(_histWPhi, norm);
	scale(_histWEta, norm);
	scale(_histMlb2, norm);
	scale(_histMlbBare, norm);

      //const double norm = crossSection()/sumOfWeights();
      //typedef map<unsigned int, NLOHisto1DPtr>::value_type IDtoNLOHisto1DPtr; ///< @todo Remove when C++11 allowed
      //foreach (IDtoNLOHisto1DPtr ihpair, _hMap) scale(ihpair.second, norm); ///< @todo Use normalize(ihpair.second, crossSection())
    }



  private:

    /// @name Histogram helper functions
    //@{

	NLOHisto1DPtr _histPt1;
	NLOHisto1DPtr _histEta1;
	NLOHisto1DPtr _histPhiemu;
	NLOHisto1DPtr _histdeltaRb;
	NLOHisto1DPtr _histdeltaRl;
	NLOHisto1DPtr _histPtMiss;
	NLOHisto1DPtr _histHT;
	NLOHisto1DPtr _histHT2;
	NLOHisto1DPtr _histHT3;
	NLOHisto1DPtr _histMlb;
	NLOHisto1DPtr _histWWphi;
	NLOHisto1DPtr _histWPt;
	NLOHisto1DPtr _histWPhi;
	NLOHisto1DPtr _histWEta;
	NLOHisto1DPtr _histMlb2;
	NLOHisto1DPtr _histMlbBare;
	Histo2DPtr _histScales;

/*
    unsigned int _thresholdLimit(unsigned int histId) {
      if (histId == 0) return 4;
      return 1;
    }

    NLOHisto1DPtr _histogram(unsigned int histId, unsigned int thresholdId) {
      assert(histId < _histLimit);
      assert(thresholdId < _thresholdLimit(histId));

      const unsigned int hInd = (histId == 0) ? thresholdId : (_thresholdLimit(0) + (histId-1) + thresholdId);
      if (_hMap.find(hInd) != _hMap.end()) return _hMap[hInd];

      if (histId == 0) _hMap.insert(make_pair(hInd,bookNLOSHisto1D(1,thresholdId+1,1)));
      else _hMap.insert(make_pair(hInd,bookNLOSHisto1D(2,histId,1)));
      return _hMap[hInd];
    }
*/

private:

    unsigned int _jet_ntag;

    //map<unsigned int, NLOHisto1DPtr> _hMap;
    //unsigned int _histLimit;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_sherpa_parton);

}
