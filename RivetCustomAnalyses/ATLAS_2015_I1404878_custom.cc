// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief ATLAS 8 TeV pseudo-top resolved analysis
  ///
  /// @author F. La Ruffa <francesco.la.ruffa@cern.ch>
  ///
  class ATLAS_2015_I1404878_custom : public Analysis {
    public:

      /// Constructor
      ATLAS_2015_I1404878_custom()
        : Analysis("ATLAS_2015_I1404878_custom")
      {
        setNeedsCrossSection(true);
      }


      void init() {
        // Eta ranges
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
        IdentifiedFinalState nu_id;
        nu_id.acceptNeutrinos();
        PromptFinalState neutrinos(nu_id);
        neutrinos.acceptTauDecays(true);
        addProjection(neutrinos, "neutrinos");

        // Jet clustering.
        VetoedFinalState vfs;
        vfs.addVetoOnThisFinalState(ewdressedelectrons);
        vfs.addVetoOnThisFinalState(ewdressedmuons);
        vfs.addVetoOnThisFinalState(neutrinos);
        FastJets jets(vfs, FastJets::ANTIKT, 0.4);
        jets.useInvisibles();
        addProjection(jets, "jets");

        //pseudotop 
        _h["ptpseudotophadron"]        = bookHisto1D(7, 1, 1);
        _h["ptttbar"]                  = bookHisto1D(3, 1, 1);
	_h["absrappseudotophadron"]    = bookHisto1D(9, 1, 1);
	_h["absrapttbar"]              = bookHisto1D(5, 1, 1);
	_h["massttbar"]                = bookHisto1D(1, 1, 1);
	_h["absPout"]                  = bookHisto1D(11, 1, 1);
	_h["chittbar"]                 = bookHisto1D(19, 1, 1);
	_h["dPhittbar"]                = bookHisto1D(13, 1, 1);
	_h["HTttbar"]                  = bookHisto1D(15, 1, 1);
	_h["Yboost"]                   = bookHisto1D(17, 1, 1);
	_h["RWt"]                      = bookHisto1D(21, 1, 1);
        //just to test the lepton reconstruction
        _h_ptlepton = bookHisto1D("pt_lepton",10,0.0,400.0);
        _h_selected_events = bookHisto1D("h_selected_events", 5, 0.0, 5.0);
        _n_events = 0;
        _n_passed_events = 0;
        _n_events_single_lepton = 0;
        _n_events_jets = 0;
        _n_events_bjets = 0;

      }

      void analyze(const Event& event) {

        // Get the selected objects, using the projections.
        _dressedelectrons     = applyProjection<DressedLeptons>(  event, "dressedelectrons").dressedLeptons();
        _dressedmuons         = applyProjection<DressedLeptons>(  event, "dressedmuons").dressedLeptons();
        _neutrinos            = applyProjection<PromptFinalState>(event, "neutrinos").particlesByPt();
        const Jets& all_jets  = applyProjection<FastJets>(        event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);


        // Calculate the missing ET, using the prompt neutrinos only

        FourMomentum met;
        foreach (const Particle& p, _neutrinos)  met += p.momentum();


        Jets bjets, lightjets;
        //Jets goodjets;
        vector<DressedLepton> _electrons;
        vector<DressedLepton> _muons; 
        //vector<unsigned int> pos_bad_jets;

/*	foreach (const DressedLepton& el, _dressedelectrons) {
           double dR_min=0;
	   for (unsigned int i = 0;i<all_jets.size(); ++i) {
	      double dR = deltaR(el,all_jets[i]);
	      if(i==0)   dR_min = dR;
	      else {
                 if(dR < dR_min)   dR_min = dR;
              }
           }
           if (dR_min < 0.2) {
              for (unsigned int j = 0; j < all_jets.size(); ++j) {
                 if(deltaR(el,all_jets[j]) == dR_min)   pos_bad_jets.push_back(j);   
	      }
	   }
	}

        for (unsigned int k = 0; k < all_jets.size(); ++k) {
	   bool remove = false;
           for (unsigned int l = 0; l<pos_bad_jets.size(); ++l){
	      if(k == pos_bad_jets[l])   remove = true;
           }
	   if(!remove)   goodjets += all_jets[k];
        }*/
   
        // overlap removal: if a lepton is near a jet, it will be excluded.
        foreach (const DressedLepton& el, _dressedelectrons) {
            //bool overlap = false;
            //foreach (const Jet& jet, goodjets)
            //    overlap |= deltaR(jet, el) < 0.4;
            //if (!overlap)   
            _electrons += el;
        }
        foreach (const DressedLepton& mu, _dressedmuons) {
            //bool overlap = false;
            //foreach (const Jet& jet, goodjets)
            //    overlap |= deltaR(jet, mu) < 0.4;
            //if (!overlap)   
            _muons += mu;
        }
         
        //// Count the number of b-tags
        foreach (const Jet& jet, all_jets){
          bool b_tagged = false;                                      
	  if (jet.containsBottom()) b_tagged = true;
          if ( b_tagged && bjets.size() < 2  )  bjets += jet;
          else lightjets += jet;
        }
       
        // Evaluate basic event selection
        bool pass_eljets = (_electrons.size() == 1) &&
                           (_muons.empty()) &&
                           (all_jets.size() >= 4) &&
                           (bjets.size() >= 2);
        bool pass_mujets = (_muons.size() == 1) &&
                           (_electrons.empty()) &&
                           (all_jets.size() >= 4) &&
                           (bjets.size() >= 2);
        _n_events++;        
        
        bool single_lepton = ((_electrons.size() == 1) && (_muons.empty())) || ((_muons.size() == 1) && (_electrons.empty()));
        if (single_lepton)  _n_events_single_lepton++;
        if (all_jets.size() >=4) _n_events_jets++;
        if (bjets.size() >=2) _n_events_bjets++;  
        // basic event selection requirements
        if (!pass_eljets && !pass_mujets) vetoEvent;

        _n_passed_events++;
    
        DressedLepton *lepton = NULL;
        if (pass_eljets)       lepton = &_electrons[0];
        else if (pass_mujets)  lepton = &_muons[0];

        // just to have a distribution
       _h_ptlepton->fill( lepton->pt(),1.0); 

        //
        //If I comment here till the end, it works on Athena and produce ptlepton histogram
        //
        FourMomentum pbjet1; //Momentum of bjet1
        FourMomentum pbjet2; //Momentum of bjet2
        if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
          pbjet1 = bjets[0].momentum();
          pbjet2 = bjets[1].momentum();
        } else {
          pbjet1 = bjets[1].momentum();
          pbjet2 = bjets[0].momentum();
        } 

        //if (bjets.size() < 2 || lightjets.size() < 2)  vetoEvent;

        double bestWmass = 999.e9;
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
        const double weight = event.weight();

        //pseudotop hadrons and leptons fill histogram
        _h["ptpseudotophadron"]->fill(ppseudotophadron.pt(), weight); //pT of pseudo top hadron
        _h["ptttbar"]->fill(          pttbar.pt(),           weight); //fill pT of ttbar in combined channel
	_h["absrappseudotophadron"]->fill(ppseudotophadron.absrap(), weight);
	_h["absrapttbar"]->fill(pttbar.absrap(), weight);
	_h["massttbar"]->fill(pttbar.mass(), weight);
	_h["absPout"]->fill(absPout, weight);
	_h["chittbar"]->fill(chi_ttbar, weight);
	_h["dPhittbar"]->fill(deltaPhi_ttbar, weight);
	_h["HTttbar"]->fill(HT_ttbar, weight);
	_h["Yboost"]->fill(Yboost, weight);
	_h["RWt"]->fill(R_Wt, weight);
	
        
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
        if ( !std::isfinite(pzneutrino) )   pzneutrino = 1;  
           
        return pzneutrino;
      }

      double _mT(const FourMomentum &l, FourMomentum &nu) const {
        return sqrt( 2 * l.pT() * nu.pT() * (1 - cos(deltaPhi(l, nu))) );
      }

      /// @name Objects that are used by the event selection decisions
      vector<DressedLepton> _dressedelectrons, _dressedmuons;
      Particles _neutrinos;
      map<string, Histo1DPtr> _h;
      Histo1DPtr _h_selected_events;
      Histo1DPtr _h_ptlepton;
      double _n_events;
      double _n_passed_events;
      double _n_events_single_lepton;
      double _n_events_jets;
      double _n_events_bjets;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1404878_custom);

}
