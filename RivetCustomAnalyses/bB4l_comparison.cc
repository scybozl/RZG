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

namespace Rivet {


  struct bB4l_comparison_Plots { };



  /// Top pair production with central jet veto
  class bB4l_comparison : public Analysis {
  public:

    /// Constructor
    bB4l_comparison()
      : Analysis("bB4l_comparison"),

    _jet_ntag(0)
    //_hMap(),
    //_histLimit(6)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(Cuts::abseta < 100);
      addProjection(fs, "ALL_FS");

      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(Cuts::abseta < 100);
      elec_fs.acceptIdPair(PID::ELECTRON);
      addProjection(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(Cuts::abseta < 100);
      muon_fs.acceptIdPair(PID::MUON);
      addProjection(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 100);
      neutrino_fs.acceptNeutrinos();
      addProjection(neutrino_fs, "NEUTRINO_FS");

      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(PID::MUON);
      addProjection(jet_input, "JET_INPUT");

      // Get the jets
      FastJets jets(jet_input, FastJets::ANTIKT, 0.5);
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
	_histPt1 = bookHisto1D("Pt1", 25, 0, 400);
	_histEta1 = bookHisto1D("Eta1", 20, -2.5, 2.5);
	_histPhiemu = bookHisto1D("Phiemu", 20, 0, 2*PI);
	_histdeltaRb = bookHisto1D("deltaRb", 20, 0, 5);
	_histdeltaRl = bookHisto1D("deltaRl", 20, 0, 5);
	_histPtMiss = bookHisto1D("PtMiss", 25, 0, 400);
	_histHT = bookHisto1D("HT", 20, 0, 1200);
	_histMlb = bookHisto1D("Mlb", 25, 0, 200);

	_histmWjB = bookHisto1D("mWjB", 20, 150, 200);
	_histmljB = bookHisto1D("mljB", 15, 0, 350);
	_histmjB  = bookHisto1D("mjB", 20, 5, 50);
	_histdeltaR = bookHisto1D("deltaR", 100, 0, 2);
	_histxB   = bookHisto1D("xB", 20, 0, 1);
	_histpTbDec = bookHisto1D("pTbDec", 15, 0, 30);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double weight = event.weight();
//      const double normWW = 1.0/3.9465;
//      weight *= normWW;

      /// Get the various sets of final state particles
      Particles elecFS = applyProjection<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt();
      Particles muonFS = applyProjection<IdentifiedFinalState>(event, "MUON_FS").particlesByPt();
      const Particles& neutrinoFS = applyProjection<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();

      // Get all jets with pT > 30 GeV
      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt();

      // Keep any jets that pass the initial rapidity cut
      vector<const Jet*> central_jets;
      foreach(const Jet& j, jets) {
        if (j.absrap() < 2.5) central_jets.push_back(&j);
      }


      // Get b hadrons with pT > 5 GeV
      /// @todo This is a hack -- replace with UnstableFinalState
      GenParticle B_hadrons;
      GenParticle B_bbar_hadrons;
      vector<GenParticle const *> allParticles = particles(event.genEvent());
      double maxpTHad = 0;
      double maxpTHadBar = 0;
      for (size_t i = 0; i < allParticles.size(); i++) {
        const GenParticle* p = allParticles[i];
        if (!PID::isHadron(p->pdg_id()) || !PID::hasBottom(p->pdg_id())) continue;
        //if (p->momentum().perp() < 5*GeV) continue;
	if (p->pdg_id()>0 && Particle(*p).pT()>maxpTHad){
	  B_hadrons = *p;
	  maxpTHad = Particle(*p).pT();
	}
	if (p->pdg_id()<0 && Particle(*p).pT()>maxpTHadBar){
	  B_bbar_hadrons = *p;
	  maxpTHadBar = Particle(*p).pT();
	}

//	cout << (p->pdg_id()>0 && Particle(*p).pT()>maxpTHad) << "\n";
//	cout << i << ": " << Particle(*p).pT() << "\n";
      }

      vector<const Jet*> b_jets;
      foreach (const Jet* j, central_jets) {
//        foreach (const GenParticle* b, B_hadrons) {
        if (j->containsParticle(B_hadrons)) b_jets.push_back(j);
//        }
      }

      vector<const Jet*> b_bar_jets;
      foreach (const Jet* j, central_jets) {
	if (j->containsParticle(B_bbar_hadrons)) b_bar_jets.push_back(j);
      }
 

      GenParticle WPlus;
      GenParticle WMinus;
      bool found = false;
      foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
           //cout << p->pdg_id() << " ";
          if (p->pdg_id() == 24){ WPlus = *p; found = true;}
	  if (p->pdg_id() == -24){ WMinus = *p; found = true;}
           //const GenVertex* pv = p->production_vertex();
      }

      // Get the MET by taking the vector sum of all neutrinos
      /// @todo Use MissingMomentum instead?
      double MET = 0;
      FourMomentum p_MET;
      foreach (const Particle& p, neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();

/*
      // Get the electrons and muons if there aren't any in the final state
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
      // Finally, the same again with the emu channel
      //if (elecFS.size() == 1 && muonFS.size() == 1) {
        // With the desired charge signs
        //if (charge(elecFS[0]) >= 0 && charge(muonFS[0]) <= 0) {
          // Calculate HT: scalar sum of the pTs of the leptons and all good jets
/*          double HT = 0;
	  if(elecFS.size() >= 1 && muonFS.size() >= 1) {
          HT += elecFS[0].pT();
          HT += muonFS[0].pT();
	  }
          foreach (const Jet* j, central_jets)
            HT += fabs(j->pT());
          // Keep events with HT > 130 GeV
          if (HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            if (b_jets.size() > 1 && elecFS.size() >= 1 && muonFS.size() >= 1) {
		 if (MET >= 20.0*GeV) {
              		passed_emu = true;

			// Additional requirements for the invariant mass event selection
			if(elecFS.size() >= 1 && muonFS.size() >= 1) {
			if(elecFS[0].pT() >= 20.*GeV && muonFS[0].pT() >= 20.*GeV && fabs(elecFS[0].eta()) < 2.5 && fabs(muonFS[0].eta()) < 2.5) {

			if(deltaR(b_jets[0]->momentum(), elecFS[0].momentum()) >= 0.4 && deltaR(b_jets[0]->momentum(), muonFS[0].momentum()) >= 0.4
			&& deltaR(b_jets[1]->momentum(), elecFS[0].momentum()) >= 0.4 && deltaR(b_jets[1]->momentum(), muonFS[0].momentum()) >= 0.4) {
				if(MET >= 60*GeV && (elecFS[0].momentum() + muonFS[0].momentum()).mass() >= 15*GeV) {
					if(abs((elecFS[0].momentum() + muonFS[0].momentum()).mass() - 91*GeV) >= 10*GeV) {
						passed_Mlb = true;
					}
				}
			}
			}
			}
		}
            }
	}
	//}
	//}

      if (passed_emu == true) {
*/
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
          if (ihist == 0) _histogram(ihist, ithres)->fill(jet_n[ithres], weight); // njets
          else if (ihist == 1) _histogram(ihist, ithres)->fill(central_jets[0]->pT(), weight); // leading jet pT
          else if (ihist == 2) _histogram(ihist, ithres)->fill(central_jets[1]->pT(), weight); // 2nd jet pT
          else if (ihist == 3 && jet_n[ithres] >= 3) _histogram(ihist, ithres)->fill(central_jets[2]->pT(), weight); // 3rd jet pT
          else if (ihist == 4 && jet_n[ithres] >= 4) _histogram(ihist, ithres)->fill(central_jets[3]->pT(), weight); // 4th jet pT
          else if (ihist == 5 && jet_n[ithres] >= 5) _histogram(ihist, ithres)->fill(central_jets[4]->pT(), weight); // 5th jet pT
        }
      }
*/
    if(b_jets.size()>=1 && b_bar_jets.size()>=1) {
	_histPt1->fill(0.5*(b_jets[0]->pT()+b_bar_jets[0]->pT()), weight);
	_histEta1->fill(0.5*(b_jets[0]->eta()+b_bar_jets[0]->eta()), weight);
	_histdeltaRb->fill(deltaR(b_jets[0]->momentum(), b_bar_jets[0]->momentum()), weight);
		//_histPhiemu->fill(deltaPhi(elecFS[0], muonFS[0]), weight);
		//_histPhiemu->fill(mapAngle0To2Pi(elecFS[0].momentum().phi() - muonFS[0].momentum().phi()), weight);
		//_histdeltaRl->fill(deltaR(elecFS[0], muonFS[0]), weight);
	_histmWjB->fill(0.5*(((Particle(WPlus).momentum()+b_jets[0]->momentum()).mass())+((Particle(WMinus).momentum()+b_bar_jets[0]->momentum()).mass())), weight);	
     }
	_histPtMiss->fill(MET, weight);
	//_histHT->fill(HT, weight);
	//if(passed_Mlb == true) {
	//	if((elecFS[0].momentum() + b_jets[0]->momentum()).mass() + (muonFS[0].momentum() + b_jets[1]->momentum()).mass()
	//	> (elecFS[0].momentum() + b_jets[1]->momentum()).mass() + (muonFS[0].momentum() + b_jets[0]->momentum()).mass()) {
	//		_histMlb->fill(((elecFS[0].momentum() + b_jets[1]->momentum()).mass() + (muonFS[0].momentum() + b_jets[0]->momentum()).mass())/2, weight);
	//	}
	//	else _histMlb->fill(((elecFS[0].momentum() + b_jets[0]->momentum()).mass() + (muonFS[0].momentum() + b_jets[1]->momentum()).mass())/2, weight);
	//}

      
    }


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

	const double norm = crossSection()/sumOfWeights();
	scale(_histPt1, norm);
	scale(_histEta1, norm);
	scale(_histPhiemu, norm);
	scale(_histdeltaRb, norm);
	scale(_histdeltaRl, norm);
	scale(_histPtMiss, norm);
	scale(_histHT, norm);
	scale(_histMlb, norm);
	scale(_histmWjB, norm);
      //const double norm = crossSection()/sumOfWeights();
      //typedef map<unsigned int, Histo1DPtr>::value_type IDtoHisto1DPtr; ///< @todo Remove when C++11 allowed
      //foreach (IDtoHisto1DPtr ihpair, _hMap) scale(ihpair.second, norm); ///< @todo Use normalize(ihpair.second, crossSection())
    }



  private:

    /// @name Histogram helper functions
    //@{

	Histo1DPtr _histPt1;
	Histo1DPtr _histEta1;
	Histo1DPtr _histPhiemu;
	Histo1DPtr _histdeltaRb;
	Histo1DPtr _histdeltaRl;
	Histo1DPtr _histPtMiss;
	Histo1DPtr _histHT;
	Histo1DPtr _histMlb;

	Histo1DPtr _histmWjB;
        Histo1DPtr _histmljB;
        Histo1DPtr _histmjB;
        Histo1DPtr _histdeltaR;
        Histo1DPtr _histxB;
        Histo1DPtr _histpTbDec;
/*
    unsigned int _thresholdLimit(unsigned int histId) {
      if (histId == 0) return 4;
      return 1;
    }

    Histo1DPtr _histogram(unsigned int histId, unsigned int thresholdId) {
      assert(histId < _histLimit);
      assert(thresholdId < _thresholdLimit(histId));

      const unsigned int hInd = (histId == 0) ? thresholdId : (_thresholdLimit(0) + (histId-1) + thresholdId);
      if (_hMap.find(hInd) != _hMap.end()) return _hMap[hInd];

      if (histId == 0) _hMap.insert(make_pair(hInd,bookHisto1D(1,thresholdId+1,1)));
      else _hMap.insert(make_pair(hInd,bookHisto1D(2,histId,1)));
      return _hMap[hInd];
    }
*/

private:

    unsigned int _jet_ntag;

    //map<unsigned int, Histo1DPtr> _hMap;
    //unsigned int _histLimit;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(bB4l_comparison);

}
