// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Math/Math.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  class ATLAS_2015_I1404878_Parton : public Analysis {
  public:

    /// Constructor
    ATLAS_2015_I1404878_Parton()
      : Analysis("ATLAS_2015_I1404878_Parton")
    {    }


  public:

    void init() {

 

      _h_ttbar_m = bookHisto1D(23,1,1);
      _h_ttbar_m_norm = bookHisto1D(24,1,1);

      _h_ttbar_pt = bookHisto1D(25,1,1);
      _h_ttbar_pt_norm = bookHisto1D(26,1,1);

      _h_ttbar_rap = bookHisto1D(27,1,1);
      _h_ttbar_rap_norm = bookHisto1D(28,1,1);

      _h_top_pt = bookHisto1D(29,1,1);
      _h_top_pt_norm = bookHisto1D(30,1,1);     

      _h_top_rap = bookHisto1D(31,1,1);
      _h_top_rap_norm = bookHisto1D(32,1,1);     

      _h_ttbar_ptout = bookHisto1D(33,1,1);
      _h_ttbar_ptout_norm = bookHisto1D(34,1,1);

      _h_ttbar_deltaphi = bookHisto1D(35,1,1);
      _h_ttbar_deltaphi_norm = bookHisto1D(36,1,1);

      _h_ttbar_ht = bookHisto1D(37,1,1);
      _h_ttbar_ht_norm = bookHisto1D(38,1,1);

      _h_ttbar_rapboost = bookHisto1D(39,1,1);
      _h_ttbar_rapboost_norm = bookHisto1D(40,1,1);

      _h_ttbar_chi = bookHisto1D(41,1,1);
      _h_ttbar_chi_norm = bookHisto1D(42,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


   
      const double weight = event.weight();

      // Vectors to hold any status=3 (anti-)top quarks (Pythia 6)
      vector<const GenParticle*> v_status3_top;
      vector<const GenParticle*> v_status3_antitop;

      // Vectors to hold any status=155 (anti-)top quarks (Herwig 6)
      vector<const GenParticle*> v_status155_top;
      vector<const GenParticle*> v_status155_antitop;

      // Vectors to hold any status=11 (anti-)top quarks for Herwig++  
      vector<const GenParticle*> v_status11_top;
      vector<const GenParticle*> v_status11_antitop;
      
      // Vectors to hold any status=22 (anti-)top quarks
      vector<const GenParticle*> v_statusOther_top;
      vector<const GenParticle*> v_statusOther_antitop;



      // Loop over all the HepMC GenParticles and keep the status 3 or status 155 tops.
      foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {

        if(p->status()==3) {
          if(p->pdg_id()==6){ v_status3_top.push_back(p);}
          if(p->pdg_id()==-6){ v_status3_antitop.push_back(p);}
        }

        if(p->status() == 155) {
          if(p->pdg_id()==6) v_status155_top.push_back(p);
          if(p->pdg_id()==-6) v_status155_antitop.push_back(p);
        }

	// for Herwig++: take only the tops that decay into Wb!!!
	if(p->status() == 11){
	  
	  if(p->pdg_id()==6) {
	  
	    GenVertex *vertex = p -> end_vertex();
	    if(vertex -> particles_out_size() == 2)
	      v_status11_top.push_back(p);
	    
	  }

	  if(p->pdg_id()==-6) {

            GenVertex *vertex = p -> end_vertex();
            if(vertex -> particles_out_size() == 2)
              v_status11_antitop.push_back(p);

          }
	}
	
	// modif. Andrea: take only the last top for Py8 and Herwig++ !!!
	if(p->pdg_id()==6) {
	  v_statusOther_top.push_back(p);
	}
	if(p->pdg_id()==-6){	  
	  v_statusOther_antitop.push_back(p);
	}
	
      }
      
      //      std::cout << "\t" << std::endl;

      // If there are more than 1 status 3 tops or anti-tops, only keep the last one put into the vector
      if(v_status3_top.size()>1) {
        MSG_DEBUG("Found more than one status 3 top! Keeping only the last one.");
        v_status3_top = vector<const GenParticle*>(v_status3_top.end() - 1, v_status3_top.end());
      }

      if(v_status3_antitop.size()>1) {
        MSG_DEBUG("Found more than one status 3 antitop! Keeping only the last one");
        v_status3_antitop = vector<const GenParticle*>(v_status3_antitop.end() - 1, v_status3_antitop.end());
      }

      // If there are more than 1 status 11 tops or anti-tops, only keep the last one put into the vector                                                             
      if(v_status3_top.size()>1) {
        MSG_DEBUG("Found more than one status 11 top! Keeping only the last one.");
        v_status11_top = vector<const GenParticle*>(v_status11_top.end() - 1, v_status11_top.end());
      }

      if(v_status11_antitop.size()>1) {
        MSG_DEBUG("Found more than one status 11 antitop! Keeping only the last one");
        v_status11_antitop = vector<const GenParticle*>(v_status11_antitop.end() - 1, v_status11_antitop.end());
      }


      // =======================================
      // Rach: check for Pythia 8 as well
      // If there are more than 1 status 3 tops or anti-tops, only keep the last one put into the vector
      if(v_statusOther_top.size()>1) {
        MSG_DEBUG("Found more than one status 22 top! Keeping only the last one.");
        v_statusOther_top = vector<const GenParticle*>(v_statusOther_top.end() - 1, v_statusOther_top.end());
      }

      if(v_statusOther_antitop.size()>1) {
        MSG_DEBUG("Found more than one status 22 antitop! Keeping only the last one");
        v_statusOther_antitop = vector<const GenParticle*>(v_statusOther_antitop.end() - 1, v_statusOther_antitop.end());
      }
      // =======================================

      Particle* top=0;
      Particle* antitop=0;

      // If there are status 3 tops and no status 155 tops this is probably a Pythia event, so used the status 3s.
      if(v_status3_top.size() == 1 && v_status3_antitop.size()==1 && v_status155_top.size() == 0 && v_status155_antitop.size()==0) {
        top     = new Particle(v_status3_top[0]);
        antitop = new Particle(v_status3_antitop[0]);
      }

      // If there are status 155 tops this must be a Herwig event, so use the status 155s.
      if( v_status155_top.size() == 1 && v_status155_antitop.size()==1 && v_status3_top.size() == 0 && v_status3_antitop.size()==0) {
        top     = new Particle(v_status155_top[0]);
        antitop = new Particle(v_status155_antitop[0]);
      }

       // If there are tops with other status this must be a Pythia8 event, so use them.
      if( v_statusOther_top.size() == 1 && v_statusOther_antitop.size()==1 && v_status155_top.size() == 0 && v_status155_antitop.size()==0 && v_status3_top.size()==0 && v_status3_antitop.size()==0) {
        top     = new Particle(v_statusOther_top[0]);
        antitop = new Particle(v_statusOther_antitop[0]);
      }

      // If there are status 155 tops this must be a Herwig event, so use the status 155s.                                                                           
      if( v_status11_top.size() == 1 && v_status11_antitop.size()==1 && v_status3_top.size() == 0 && v_status3_antitop.size()==0) {
        top     = new Particle(v_status11_top[0]);
        antitop = new Particle(v_status11_antitop[0]);
      }


 
	
      if (top && antitop) 
      {

	bool isTopLep    = hasLeptonicDecay(top);
	bool isAntiTopLep = hasLeptonicDecay(antitop);
	FourMomentum toplep;
	FourMomentum tophad;
	if (isTopLep) {
	  toplep = top->momentum();
	  tophad = antitop->momentum();
	} else {
	  toplep = antitop->momentum();
	  tophad = top->momentum();
	}

	//	if  (!(isTopLep && isAntiTopLep)) 
	if (1 == 1)
	{
	
	  
	  //std::cout<<"Top = "<<isTopLep<<" Antitop = "<<isAntiTopLep<<std::endl;

	  FourMomentum ttbar = top->momentum() + antitop->momentum();
	  Vector3 z_versor(0,0,1);
	  Vector3 vtophad = tophad.vector3();
	  Vector3 vtoplep = toplep.vector3(); 

	  double ttbar_mass = ttbar.mass();
	  _h_ttbar_m_norm->fill(ttbar_mass, weight);
	  _h_ttbar_m->fill(ttbar_mass, weight);

	  double ttbar_pt = ttbar.pt();
	  _h_ttbar_pt_norm->fill(ttbar_pt, weight);
	  _h_ttbar_pt->fill(ttbar_pt, weight);

	  double ttbar_rap = ttbar.absrap();
	  _h_ttbar_rap_norm->fill(ttbar_rap, weight);
	  _h_ttbar_rap->fill(ttbar_rap, weight);
	  

	  double absPout = fabs(vtophad.dot((vtoplep.cross(z_versor))/(vtoplep.cross(z_versor).mod())));
	  _h_ttbar_ptout->fill(absPout);
	  _h_ttbar_ptout_norm->fill(absPout);
	  
	  double deltaPhi_ttbar = deltaPhi(toplep,tophad);
	  _h_ttbar_deltaphi->fill(deltaPhi_ttbar);
	  _h_ttbar_deltaphi_norm->fill(deltaPhi_ttbar);

	  double HT = tophad.pt() + toplep.pt();
	  _h_ttbar_ht->fill(HT);
	  _h_ttbar_ht_norm->fill(HT);

	  double rapboost = 0.5 * (fabs(tophad.rapidity() + toplep.rapidity()));
	  _h_ttbar_rapboost->fill(rapboost);
	  _h_ttbar_rapboost_norm->fill(rapboost);


	  double ystar = (tophad.pt() > toplep.pt()) ? 0.5 * (tophad.rap()-toplep.rap()) : 0.5*(toplep.rap()-tophad.rap()); 
	  double chi_ttbar = exp(2 * fabs(ystar));
	  _h_ttbar_chi->fill(chi_ttbar);
	  _h_ttbar_chi_norm->fill(chi_ttbar);
      

	  // fill only the hadronically decay top quark
	  _h_top_pt->fill(tophad.pt(), weight);
	  _h_top_pt_norm->fill(tophad.pt(), weight);
	  _h_top_rap->fill(tophad.absrap(), weight);
	  _h_top_rap_norm->fill(tophad.absrap(), weight);

	} else {
	  // std::cout<<"=================="<<std::endl;
	  // foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
	  //   if(p->status()==3) {
	  //     std::cout<<p->pdg_id()<<std::endl;
	  //   }
	  // }
	  // vetoEvent;
	}
      } else {
	MSG_ERROR("Did not find both top and anti-top!");
	vetoEvent;
      }
      
      if(top) delete top;
      if(antitop) delete antitop;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize to unity
      normalize(_h_top_pt_norm);
      normalize(_h_top_rap_norm);
      normalize(_h_ttbar_m_norm);
      normalize(_h_ttbar_pt_norm);
      normalize(_h_ttbar_rap_norm);
      normalize(_h_ttbar_ptout_norm);
      normalize(_h_ttbar_deltaphi_norm);
      normalize(_h_ttbar_ht_norm);
      normalize(_h_ttbar_rapboost_norm);
      normalize(_h_ttbar_chi_norm);


      // Normalize to cross-section
      const double norm = crossSection()/sumOfWeights()/picobarn;
      scale(_h_top_pt, norm); 
      scale(_h_top_rap, norm); 
      scale(_h_ttbar_m, norm);
      scale(_h_ttbar_pt, norm);
      scale(_h_ttbar_rap, norm);
      scale(_h_ttbar_ptout, norm);
      scale(_h_ttbar_deltaphi, norm);
      scale(_h_ttbar_ht, norm);
      scale(_h_ttbar_rapboost, norm);
      scale(_h_ttbar_chi, norm);
    }



  private:
    bool hasLeptonicDecay(Particle* p)
    {
      bool islep = false;

      // find out the decay modes
      vector<Particle> childs = p->children();
      for (unsigned int i=0; i<childs.size(); i++)
	if (childs[i].abspid() == 24) {
	  vector<Particle> Wchilds = childs[i].children();
	  while (Wchilds.size() == 1){
	    Wchilds = Wchilds[0].children();
	  }
	  if (Wchilds[0].abspid() > 10)
	    islep = true;
	}
      return islep;
    }
    

    Histo1DPtr _h_ttbar_m;
    Histo1DPtr _h_ttbar_m_norm;
    Histo1DPtr _h_ttbar_pt;
    Histo1DPtr _h_ttbar_pt_norm;
    Histo1DPtr _h_ttbar_rap;
    Histo1DPtr _h_ttbar_rap_norm;
    Histo1DPtr _h_top_pt;
    Histo1DPtr _h_top_pt_norm;
    Histo1DPtr _h_top_rap;
    Histo1DPtr _h_top_rap_norm;
    Histo1DPtr _h_ttbar_ptout;
    Histo1DPtr _h_ttbar_ptout_norm;
    Histo1DPtr _h_ttbar_deltaphi;
    Histo1DPtr _h_ttbar_deltaphi_norm;
    Histo1DPtr _h_ttbar_ht;
    Histo1DPtr _h_ttbar_ht_norm;
    Histo1DPtr _h_ttbar_rapboost;
    Histo1DPtr _h_ttbar_rapboost_norm;
    Histo1DPtr _h_ttbar_chi;
    Histo1DPtr _h_ttbar_chi_norm;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1404878_Parton);

}
