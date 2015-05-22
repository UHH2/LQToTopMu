#include "UHH2/LQToTopMu/include/LQGen.h"

using namespace std;
using namespace uhh2;

LQGen::LQGen(const vector<GenParticle> & genparticles, bool throw_on_failure)/*: m_type(e_notfound)*/ {    
    int n_LQ = 0, n_antiLQ = 0;
    for(unsigned int i=0; i<genparticles.size(); ++i) {
        const GenParticle & genp = genparticles[i];
        if (abs(genp.pdgId()) == 42){ // 42 = LQ's
            auto top = genp.daughter(&genparticles, 1);
            auto mu = genp.daughter(&genparticles, 2);
            if(!top || !mu){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has not ==2 daughters");
                return;
            }
            if(abs(top->pdgId()) != 6){
                std::swap(top, mu);
            }
            if(abs(top->pdgId()) != 6){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has no top daughter");
                return;
            }
            
            // NOTE: here, we could skip over intermediate W bosons. However,
            // this Pythia8-related problem is now fixed when creating ntuples already,
            // so this should not be necessary.
            
            if(abs(mu->pdgId()) != 13){
                if(throw_on_failure) throw runtime_error("LQGen: LQ has no muon daughter");
                return;
            }
            // now get W daughters:
            auto topd1 = top->daughter(&genparticles, 1);
            auto topd2 = top->daughter(&genparticles, 2);
            if(!topd1 || !topd2){
                if(throw_on_failure) throw runtime_error("LQGen: top has not ==2 daughters");
                return;
            }
            
            // now that we collected everything, fill the member variables. 
            // Use different member variables according to LQ charge.
            if(genp.pdgId() == 42){
                m_LQ = genp;
                m_TopLQ = *top;
                m_muLQ = *mu;
                m_Topdecay1 = *topd1;
                m_Topdecay2 = *topd2;
                ++n_LQ;
            }
            else{
                m_AntiLQ = genp;
                m_TopAntiLQ = *top;
                m_muAntiLQ = *mu;
                m_Antitopdecay1 = *topd1;
                m_Antitopdecay2 = *topd2;
                ++n_antiLQ;
            }
        }
    }
    if(n_LQ != 1 || n_antiLQ != 1){
        if(throw_on_failure)  throw runtime_error("LQGen: did not find exactly one LQ and one antiLQ in the event");
        return;
    }
    
    // calculate decay channel by counting the number of charged leptons
    // in the W daughters:
    /*       int n_e = 0, n_m = 0, n_t = 0;
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        int id = abs(wd.pdgId());
        if(id == 11) ++n_e;
        else if(id == 13) ++n_m;
        else if(id == 15) ++n_t;
    }

    // dilepton channels:
    if(n_e == 2){
        m_type = e_ee;
    }
    else if(n_e == 1 && n_m == 1){
        m_type = e_emu;
    }
    else if(n_e == 1 && n_t == 1){
        m_type = e_etau;
    }
    else if(n_m == 2){
        m_type = e_mumu;
    }
    else if(n_m == 1 && n_t == 1){
        m_type = e_mutau;
    }
    else if(n_t == 2){
        m_type = e_tautau;
    }
    // lepton+jet channels:
    else if(n_e == 1){
        m_type = e_ehad;
    }
    else if(n_m == 1){
        m_type = e_muhad;
    }
    else if(n_t == 1){
        m_type = e_tauhad;
    }
    // hadronic:
    else{
        m_type = e_had;
    }*/
}   


GenParticle LQGen::LQ() const{
    return m_LQ;
}

GenParticle LQGen::AntiLQ() const{
    return m_AntiLQ;
} 

GenParticle LQGen::TopLQ() const{
    return m_TopLQ;
}

GenParticle LQGen::TopAntiLQ() const{
    return m_TopAntiLQ;
}

GenParticle LQGen::muLQ() const{
    return m_muLQ;
}

GenParticle LQGen::muAntiLQ() const{
    return m_muAntiLQ;
} 

GenParticle LQGen::Topdecay1() const{
    return m_Topdecay1;
} 

GenParticle LQGen::Topdecay2() const{
    return m_Topdecay2;
} 

GenParticle LQGen::Antitopdecay1() const{
    return m_Antitopdecay1;
} 

GenParticle LQGen::Antitopdecay2() const{
    return m_Antitopdecay2;
} 

/*LQGen::E_DecayChannel LQGen::DecayChannel()  const{  
    return m_type;
    }*/

/*bool LQGen::IsTopHadronicDecay()  const{
    return abs(m_Wdecay1.pdgId()) <= 5;
    }*/

/*bool LQGen::IsAntiTopHadronicDecay()  const{
    return abs(m_WMinusdecay1.pdgId()) <= 5;
    }*/

namespace {
    
  /*bool is_charged_lepton(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 11 || id == 13 || id == 15;
}

bool is_neutrino(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 12 || id == 14 || id == 16;
    }*/

}

/*GenParticle LQGen::ChargedLepton() const{
    if (m_type != e_ehad &&  m_type != e_muhad  && m_type!= e_tauhad){
        throw runtime_error("LQGen::ChargedLepton called, but this is no l+jets ttbar event!");
    }
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        if(is_charged_lepton(wd)) return wd;
    }
    throw logic_error("logic error in LQGen::ChargedLepton");
    }*/

/*GenParticle LQGen::Neutrino() const{
    if (m_type != e_ehad &&  m_type != e_muhad  && m_type!= e_tauhad){
        throw runtime_error("LQGen::ChargedLepton called, but this is no l+jets ttbar event!");
    }
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        if(is_neutrino(wd)) return wd;
    }
    throw logic_error("logic error in LQGen::Neutrino");
    }*/

/*GenParticle LQGen::TopLep() const{
    if(ChargedLepton().charge()>0) return Top();
    else return Antitop();
}

GenParticle LQGen::TopHad() const{
    if(ChargedLepton().charge()<0) return Top();
    else return Antitop();
}

GenParticle LQGen::BLep() const{
    if(ChargedLepton().charge()>0) return bTop();
    else return bAntitop();
}

GenParticle LQGen::BHad() const{
    if(ChargedLepton().charge()<0) return bTop();
    else return bAntitop();
}

GenParticle LQGen::WLep() const{
    if(ChargedLepton().charge()>0) return WTop();
    else return WAntitop();
}

GenParticle LQGen::WHad() const{
    if(ChargedLepton().charge()<0) return WTop();
    else return WAntitop();
}

GenParticle LQGen::Q1() const{
    if(ChargedLepton().charge()>0) return WMinusdecay1();
    else return Wdecay1();
}

GenParticle LQGen::Q2() const{
    if(ChargedLepton().charge()>0) return WMinusdecay2();
    else return Wdecay2();
    }*/


LQGenProducer::LQGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
    h_LQLQbargen = ctx.get_handle<LQGen>(name);
}

bool LQGenProducer::process(Event & event){
    event.set(h_LQLQbargen, LQGen(*event.genparticles, throw_on_failure));
    return true;
}
