/**
 * RESOLVEDDIFFXS
 * Calculation of the boosted ttbar production differential cross-section @ 13 TeV
 * Author: Steffen Henkelmann
 * steffen.henkelmann@cern.ch
 * based on <ATLAS_2015_I1404878> 8 TeV
 **/


// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/Math.hh"
#include "Rivet/Tools/Utils.hh"

#include "HepMC/PdfInfo.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/GenVertex.h"

namespace Rivet {
  class Resolved13TeVljets : public Analysis {
    public:

      /// Constructor
      Resolved13TeVljets()
        : Analysis("Resolved13TeVljets")
      {
        setNeedsCrossSection(true);
      }


      void init() {
        // Eta ranges
	//        Cut eta_full = (Cuts::abseta < 4.2) & (Cuts::pT >= 1.0*MeV);
	Cut eta_full = (Cuts::abseta < 5.0);
        Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);

        // All final state particles
        FinalState fs(eta_full);

	// Get photons to dress leptons
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);

        // Projection to find the electrons
        IdentifiedFinalState el_id(fs);
        el_id.acceptIdPair(PID::ELECTRON);

        PromptFinalState electrons(el_id);
        electrons.acceptTauDecays(true);
        addProjection(electrons, "electrons");

        DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true, true);
        addProjection(dressedelectrons, "dressedelectrons");

	DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true, true);
        addProjection(ewdressedelectrons, "ewdressedelectrons");

        // Projection to find the muons
        IdentifiedFinalState mu_id(fs);
        mu_id.acceptIdPair(PID::MUON);

        PromptFinalState muons(mu_id);
        muons.acceptTauDecays(true);
        addProjection(muons, "muons");

        DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true, true);
        addProjection(dressedmuons, "dressedmuons");

        DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true, true);
        addProjection(ewdressedmuons, "ewdressedmuons");

        // Projection to find neutrinos and produce MET
        IdentifiedFinalState nu_id(fs);
        nu_id.acceptNeutrinos();
	addProjection(nu_id, "neutrinos");
        PromptFinalState neutrinos(nu_id);
        neutrinos.acceptTauDecays(true);
	//	addProjection(neutrinos, "neutrinos");

        // Jet clustering.
	//	VetoedFinalState vfs(FinalState(Cuts::pT > 0.5GeV))
        VetoedFinalState vfs(fs);
        vfs.addVetoOnThisFinalState(ewdressedelectrons);
        vfs.addVetoOnThisFinalState(ewdressedmuons);
        vfs.addVetoOnThisFinalState(neutrinos);
        FastJets jets(vfs, FastJets::ANTIKT, 0.4);
        jets.useInvisibles(true);
        addProjection(jets, "jets");

        //pseudotop
	_h["ptpseudotophadron"]        = bookHisto1D(1, 1, 1);
	_h["ptttbar"]                  = bookHisto1D(1, 1, 2);
        _h["absrappseudotophadron"]    = bookHisto1D(1, 1, 3);
        _h["absrapttbar"]              = bookHisto1D(1, 1, 4);
        _h["massttbar"]                = bookHisto1D(1, 1, 5);
	_h["absPout"]                  = bookHisto1D(1, 1, 6);
        _h["dPhittbar"]                = bookHisto1D(1, 1, 8);
        //just to test the lepton reconstruction
        _h_ptlepton = bookHisto1D("pt_lepton",10,0.0,400.0);
        _h_selected_events = bookHisto1D("h_selected_events", 5, 0.0, 5.0);
	_h_cutflow = bookHisto1D("h_cutflow", 20, 0.0, 20.0);
        _n_events = 0;
        _n_passed_events = 0;
        _n_events_single_lepton = 0;
        _n_events_jets = 0;
        _n_events_bjets = 0;

      }

      void analyze(const Event& event) {
	const double weight = event.weight();
        // Get the selected objects, using the projections.
        _dressedelectrons     = applyProjection<DressedLeptons>(  event, "dressedelectrons").dressedLeptons();
        _dressedmuons         = applyProjection<DressedLeptons>(  event, "dressedmuons").dressedLeptons();
        _neutrinos            = applyProjection<IdentifiedFinalState>(event, "neutrinos").particlesByPt();
        const Jets& all_jets  = applyProjection<FastJets>(        event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);


        // Calculate the missing ET, using the prompt neutrinos only

        FourMomentum met;
        foreach (const Particle& p, _neutrinos)  met += p.momentum();


        Jets bjets, lightjets;

        vector<DressedLepton> _electrons;
        vector<DressedLepton> _muons;





        // overlap removal: if a lepton is near a jet, it will be excluded.

        foreach (const DressedLepton& el, _dressedelectrons) {
            //bool overlap = false;

            //foreach (const Jet& jet, all_jets)

            //    overlap |= deltaR(jet, el) < 0.4;

            //if (!overlap)

            _electrons += el;
        }
        foreach (const DressedLepton& mu, _dressedmuons) {
            //bool overlap = false;
            //foreach (const Jet& jet, all_jets)

            //    overlap |= deltaR(jet, mu) < 0.4;

            //if (!overlap)

            _muons += mu;
        }


	//        foreach (const Jet& jet, all_jets){
	//          HepMC::GenEvent::particle_const_iterator itrParton     = event.genEvent()->particles_begin();
	//          HepMC::GenEvent::particle_const_iterator itrPartonLast = event.genEvent()->particles_end();
	//          bool b_tagged = false;

	//          double maxEnergy = 0.;
	//          int flavourPdgId = 0;
	//          for( ; itrParton != itrPartonLast ; ++itrParton ) {
	//              HepMC::GenParticle * p = *itrParton;
	//              if(p->status() == 1) continue;
	//              if( p->pdg_id() == 22 ) continue; // skip photons
	//              if(abs(p->pdg_id()) > 5 && p->pdg_id() != 21 && p->pdg_id() != 9) continue;
	//              const FourVector& particle = p->momentum();
	//              const double dR = Rivet::deltaR( particle, jet.momentum() );
	//              if( dR >= 0.4 ) continue;

	//              const double energy = p->momentum().e();
	//              if( energy > maxEnergy ) {
	//                  flavourPdgId = p->pdg_id();
	//                  maxEnergy    = energy;
	//              }

	//          } /* loop over partons */

	//          // check if the jet is b-tagged
	//          b_tagged |= abs(flavourPdgId) == 5;
	//         if (b_tagged && bjets.size() < 2) {
	//            bjets += jet;
	//          }
	//          else {
	//            lightjets += jet;
	//          }
	//
	//        }  /* loop over jets */

	Jets jets;
	foreach (const Jet& jet, all_jets) {
	  bool keep = true;
	  if((_electrons.size() == 1 && _muons.size() == 0)){
	    foreach( DressedLepton lep, _dressedelectrons )  keep |= deltaR(jet, lep) < 0.4;
	  }else if(_electrons.size()==0 && _muons.size() ==1){
	    foreach( DressedLepton lep, _dressedmuons )  keep |= deltaR(jet, lep) < 0.4;
	  }
	  if (keep)  jets += jet;
	}

	foreach( Jet jet, jets) {
	  bool b_tagged = jet.bTags(Cuts::pT > 5*GeV).size();
	  if ( b_tagged && bjets.size() < 2)  bjets     +=jet;
	  else             lightjets += jet;
	}

	bool debug = false;
	if(debug){
	  if((_electrons.size() == 1 && _muons.size() == 0)){/*  ||
								 (_electrons.size()==0 && _muons.size() ==1)){*/
      	if(_electrons.size()>0)
      	  std::cout << " DEBUG :: electron -- (" << _electrons.size() << ")"
      		    << " pt :: > " << _electrons.at(0).momentum().pt()
      		    << " eta :: > "<< _electrons.at(0).momentum().eta()
      		    << " phi :: > "<< _electrons.at(0).momentum().phi()
      		    << std::endl;

      	if(_muons.size()>0)
      	  std::cout << " DEBUG :: muon -- (" << _muons.size() << ")"
      		    << " pt :: > " << _muons.at(0).momentum().pt()
      		    << " eta :: > "<< _muons.at(0).momentum().eta()
      		    << " phi :: > "<< _muons.at(0).momentum().phi()
      		    << std::endl;


      	if(all_jets.size()>0)
      	  foreach(const Jet& j, all_jets)
      	    std::cout << " DEBUG :: jets (R=0.4) -- (" << all_jets.size() <<")"
      		      << " pt :: > " << j.momentum().pt()
      		      << " eta :: > "<< j.momentum().eta()
      		      << " phi :: > "<< j.momentum().phi()
      		      << std::endl;



	if(bjets.size()>0)
      	  foreach(const Jet& i, bjets)
      	    std::cout << " DEBUG :: bjets  -- (" << bjets.size() <<")"
      		      << " pt :: > " << i.momentum().pt()
      		      << " eta :: > "<< i.momentum().eta()
      		      << " phi :: > "<< i.momentum().phi()
      		      << std::endl;
	std::cout << " DEBUG :: neutrino sum -- "
		  << " MET :: > "<< met.pt()
	  //		  << " MET :: > "<<  neutrino_sum.Et() //* neutrino_sum.vector3().setZ(0.0).unit()
		  << std::endl;

	}
	}//debug


	unsigned int icut = 0;
	//	double transmass = 0 ;
	CountEvent(_h_cutflow, icut, weight);
	if( _electrons.size() >= 1 ){
	  CountEvent(_h_cutflow, icut, weight);
	  if(_electrons.size() == 1){
	    CountEvent(_h_cutflow, icut, weight);
	    if( _muons.size() == 0){
	      CountEvent(_h_cutflow, icut, weight);
	      if( all_jets.size() >=1 ){
		CountEvent(_h_cutflow, icut, weight);
		if( all_jets.size() >=2){
		  CountEvent(_h_cutflow, icut, weight);
		  if( all_jets.size() >=3){
		    CountEvent(_h_cutflow, icut, weight);
		    if( all_jets.size() >=4){
		      CountEvent(_h_cutflow, icut, weight);
		      if(bjets.size()>=1){
			CountEvent(_h_cutflow, icut, weight);
			if(bjets.size()>=2){
			  CountEvent(_h_cutflow, icut, weight);
			  if(met.pt() > 20.){
			    CountEvent(_h_cutflow, icut, weight);
			    double transmass = TransMass( _electrons.at(0).momentum().pt(),
							  _electrons.at(0).momentum().phi(),
						met.pt(),
						met.phi());
			if(transmass + met.pt() > 60.) CountEvent(_h_cutflow, icut, weight);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	}
      }



        _n_events++;
        // Evaluate basic event selection



        bool single_electron=(_electrons.size() == 1) && (_muons.empty());
        bool single_muon=(_muons.size() == 1) && (_electrons.empty());



        DressedLepton *lepton = NULL;
        if (single_electron)       lepton = &_electrons[0];
        else if (single_muon)      lepton = &_muons[0];

        if(!single_electron && !single_muon) vetoEvent;
        else _n_events_single_lepton++;
        // just to have a distribution
       _h_ptlepton->fill( lepton->pt(),1.0);


        if (all_jets.size() < 4) vetoEvent;
        else _n_events_jets++;
        if (bjets.size() <2) vetoEvent;
        else _n_events_bjets++;

        FourMomentum pbjet1; //Momentum of bjet1
        FourMomentum pbjet2; //Momentum of bjet

        if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
          pbjet1 = bjets[0].momentum();
          pbjet2 = bjets[1].momentum();
        } else {
          pbjet1 = bjets[1].momentum();
          pbjet2 = bjets[0].momentum();
        }

        //if (bjets.size() < 2 || lightjets.size() < 2)  vetoEvent;

        double bestWmass = 1000.0*TeV;
        double mWPDG = 80.399*GeV;
        int Wj1index = -1, Wj2index = -1;
        for (unsigned int i = 0; i < (lightjets.size() - 1); ++i) {
          for (unsigned int j = i + 1; j < lightjets.size(); ++j) {
            double wmass = (lightjets[i].momentum() + lightjets[j].momentum()).mass();
            if (fabs(wmass - mWPDG) < fabs(bestWmass - mWPDG)) {
              bestWmass = wmass;
              Wj1index = i;
              Wj2index = j;
            }
          }
        }


        FourMomentum pjet1 = lightjets[Wj1index].momentum();
        FourMomentum pjet2 = lightjets[Wj2index].momentum();

	// compute hadronic W boson
	FourMomentum pWhadron = pjet1 + pjet2;
        double pz = computeneutrinoz(lepton->momentum(), met);
        FourMomentum ppseudoneutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz)), met.px(), met.py(), pz);



        //compute leptonic, hadronic, combined pseudo-top
        FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
        FourMomentum ppseudotophadron = pbjet2 + pWhadron;
        FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

	Vector3 z_versor(0,0,1);
	Vector3 vpseudotophadron = ppseudotophadron.vector3();
	Vector3 vpseudotoplepton = ppseudotoplepton.vector3();
        // Variables
        double ystar = (ppseudotophadron.pt() > ppseudotoplepton.pt()) ? 0.5 * (ppseudotophadron.rap()-ppseudotoplepton.rap()) : 0.5*(ppseudotoplepton.rap()-ppseudotophadron.rap());
        double chi_ttbar = exp(2 * fabs(ystar));
	double deltaPhi_ttbar = deltaPhi(ppseudotoplepton,ppseudotophadron);
	double HT_ttbar = ppseudotophadron.pt() + ppseudotoplepton.pt();
	double Yboost = 0.5 * (fabs(ppseudotophadron.rapidity() + ppseudotoplepton.rapidity()));
	double R_Wt = pWhadron.pt() / ppseudotophadron.pt();
	double absPout = fabs(vpseudotophadron.dot((vpseudotoplepton.cross(z_versor))/(vpseudotoplepton.cross(z_versor).mod())));

        // Fill histograms
	//        const double weight = event.weight();

        //pseudotop hadrons and leptons fill histogram
        _h["ptpseudotophadron"]->fill(ppseudotophadron.pt(), weight); //pT of pseudo top hadron
        _h["ptttbar"]->fill(          pttbar.pt(),           weight); //fill pT of ttbar in combined channel
	_h["absrappseudotophadron"]->fill(ppseudotophadron.absrap(), weight);
	_h["absrapttbar"]->fill(pttbar.absrap(), weight);
	_h["massttbar"]->fill(pttbar.mass(), weight);
	_h["absPout"]->fill(absPout, weight);


      }

      void finalize() {
        // Normalize to cross-section
        const double sf = (crossSection() / sumOfWeights());
        for (map<string, Histo1DPtr>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
          scale(hit->second, sf);
        }
        _h_selected_events->fill (0+0.5,_n_events);
        _h_selected_events->fill (1+0.5,_n_events_single_lepton);
        _h_selected_events->fill (2+0.5,_n_events_jets);
        _h_selected_events->fill (3+0.5,_n_events_bjets);
        _h_selected_events->fill (4+0.5,_n_passed_events);
      }

    private:
  double TransMass(double ptLep, double phiLep, double met, double phiMet) {
      return std::sqrt(2.0*ptLep*met*( 1 - std::cos( phiLep-phiMet ) ) );
    }


    void CountEvent(Histo1DPtr h, unsigned int& icut, double weight) {
      h->fill(icut+0.5, weight);
      ++icut;
    }


      double computeneutrinoz(const FourMomentum& lepton, FourMomentum& met) const {
        //computing z component of neutrino momentum given lepton and met
        double pzneutrino;
        double m_W = 80.399; // in GeV, given in the paper
        double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.px() + lepton.py() * met.py());
        double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
        double b = -2*k*lepton.pz();
        double c = sqr( lepton.E() ) * sqr( met.pT() ) - sqr( k );
        double discriminant = sqr(b) - 4 * a * c;
        double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
        if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
        else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
          double absquad[2];
          for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
          if (absquad[0] < absquad[1])  pzneutrino = quad[0];
          else                          pzneutrino = quad[1];
        }
        //if ( !std::isfinite(pzneutrino) )   pzneutrino = 1;

        return pzneutrino;
      }



      /// @name Objects that are used by the event selection decisions
      vector<DressedLepton> _dressedelectrons, _dressedmuons;
      Particles _neutrinos;
      map<string, Histo1DPtr> _h;
      Histo1DPtr _h_selected_events;
      Histo1DPtr _h_ptlepton;
    Histo1DPtr _h_cutflow;
      //fastjet::AreaDefinition* _area_def;
      double _n_events;
      double _n_passed_events;
      double _n_events_single_lepton;
      double _n_events_jets;
      double _n_events_bjets;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(Resolved13TeVljets);

}
