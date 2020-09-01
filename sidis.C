#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}

void sidis(){


  //some particles
  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.389,10.389);
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector g1(0,0,0,0);
  TLorentzVector g2(0,0,0,0);
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());

  //Create histogram
  TH1F hmass("pi0mass","Invariant Mass to 2#gamma",100,0,0.6);
  TH1F htime("DeltaTime","Time difference of 2#gamma",100,-10,10);
  auto* hmiss=new TH1F("missM","missM",200,-2,3);
  auto* hnphe = new TH1F("nphe", "nphi" , 100, 0,80);

  //create particles before looping to be more efficient
  TLorentzVector p4_gamma1;
  TLorentzVector p4_gamma2;
  TLorentzVector p4_pi1;
  TLorentzVector p4_pi2;
  TLorentzVector miss;
  TLorentzVector virtual_photon;   
  TLorentzVector hadronic_state;
  auto* hdipion_mass = new TH1F("h2pimass", "Invariant mass to 2 pi+/-", 25,0,2.0);
  auto* hdipion_missingmass = new TH1F("h2pimissingmass", "Missing mass for e pi- pi+ X", 25,0,1.5);
  auto* hQ2W = new TH2F("hQ2W", "Q2 vs W", 100, 0, 10.0, 100,0,5.0); 
  auto * hW = new TH1F("hW", "W [GeV]", 100, 0.0,10.0);
  auto * heoverp = new TH2F("heoverp", "heoverp", 100,0,10.0, 100,0,0.5);

  auto* hetheta = new TH1F("hetheta", "hetheta", 100,0.0,10.0);
  auto* h_betap = new TH2F("h_betap", "h_betap", 100,0,10,100,0.0,2.0);
  auto* h_vtp = new TH2F("h_vtp", "h_vtp", 100,0,10,100,0.0,2.0);

  gBenchmark->Start("timer");
  int counter=0;
   
  clas12root::HipoChain chain;

  //  TString runnumber = "011582"; //liquid deuterium
  //TString runnumber = "011592"; // liquid helium run
  TString runnumber = "011599"; //lead run
  chain.Add("/lustre19/expphy/volatile/clas12/rg-e/production/pass0/full/recon/"+runnumber+"/rec_clas_"+runnumber+".evio.000*.hipo"); 
  auto c12=chain.GetC12Reader();


  //c12->queryRcdb();

  while (chain.Next()){
    c12=chain.GetC12Reader();
    //std::cout << c12->event()->getHelicity() << " " <<  c12->event()->getHelicityRaw() << std::endl;
    // auto rcdbData=c12->getRcdbVals();//struct with all relevent rcdb values
  
    //The following run conditions can be returned directly by c12
    //cout<<"Event count: "<<rcdbData.event_count<<endl;
    //cout<<"Beam energy: "<<rcdbData.beam_energy<<endl;
    //cout<<"Beam current: "<<rcdbData.beam_current<<endl;

    auto gammas=c12->getByID(22);
    auto electrons=c12->getByID(11);
    auto protons=c12->getByID(2212);
    auto pips=c12->getByID(211);
    auto pims=c12->getByID(-211);

    if(electrons.size()<1) continue;
    //std::cout << "NPHE " << electrons[0]->par()->che(HTCC)->getNphe() <<std::endl;

    el.SetXYZM(electrons[0]->par()->getPx() , electrons[0]->par()->getPy() , electrons[0]->par()->getPz(), db->GetParticle(11)->Mass());    
    virtual_photon = beam-el;
    hadronic_state = beam + target-el;
   
    //std::cout << " Q2 " << -virtual_photon.M2() << std::endl; 
    //std::cout << " W " << hadronic_state.M() << std::endl;
    hQ2W->Fill(-virtual_photon.M2(), hadronic_state.M() );
    hW->Fill(hadronic_state.M());
    // std::cout << "theta" << electrons[0]->getTheta() << std::endl;
    hetheta->Fill(electrons[0]->getTheta());

    //if(electrons[0]->getRegion()==FD) std::cout << " Forward Detector " << std::endl;
    //if(electrons[0]->getRegion()==FT) std::cout << " Forward Tagger " << std::endl;

    if(electrons.size()==1 and gammas.size()>1){
      if(gammas[0]->getRegion()==FD && gammas[1]->getRegion()==FD){
        //Select photons velocity (to remove neutrons)
	//std::cout << "Velocity " << gammas[0]->par()->getBeta() <<  " " << gammas[1]->par()->getBeta() <<std::endl;
	//        if(gammas[0]->par()->getBeta()<0.9 or gammas[1]->par()->getBeta()<0.9) continue;
        //if(gammas[0]->par()->getBeta()>1.1 or gammas[1]->par()->getBeta()>1.1) continue; 

        p4_gamma1.SetXYZM(gammas[0]->par()->getPx(),gammas[0]->par()->getPy(),gammas[0]->par()->getPz(),0);
        p4_gamma2.SetXYZM(gammas[1]->par()->getPx(),gammas[1]->par()->getPy(),gammas[1]->par()->getPz(),0);
        el.SetXYZM(electrons[0]->par()->getPx() , electrons[0]->par()->getPy() , electrons[0]->par()->getPz(), db->GetParticle(11)->Mass());        
        auto pi0 = p4_gamma1 + p4_gamma2;
      
         hmass.Fill(pi0.M());
        //std::cout << pi0.M() << " " << gammas[0]->getTime() - gammas[1]->getTime() << "  " << std::endl;
	 //	std::cout <<"Electron momentum " <<  electrons[0]->par()->getP() << std::endl;
        htime.Fill(gammas[0]->getTime() - gammas[1]->getTime() );

	miss=beam+target-el-g1-g2;   
        hmiss->Fill(miss.M());
	//std::cout << " Missing mass " << miss.M() << std::endl;
        }
    }

    if(electrons.size()==1 && gammas.size()==2 && protons.size()==1 &&  pips.size()==1 &&pims.size() == 1){

      // set the particle momentum
      SetLorentzVector(el,electrons[0]);
      SetLorentzVector(pr,protons[0]);
      SetLorentzVector(g1,gammas[0]);
      SetLorentzVector(g2,gammas[1]);
      SetLorentzVector(pip,pips[0]);
      SetLorentzVector(pim,pims[0]);

      miss=beam+target-el-pr-g1-g2-pip-pim;
      //hmiss->Fill(miss.M2());
    }

    //std::cout << "pi- " << pims.size() << " pi+ " << pips.size() << " e" << electrons.size() << " gamma" << gammas.size() <<  std::endl;

    if(pips.size()>0 &&pims.size()>0){
      //  SetLorentzVector(el,electrons[0]);
      //SetLorentzVector(pip,pips[0]);
      p4_pi1.SetXYZM(pims[0]->par()->getPx(),pims[0]->par()->getPy(),pims[0]->par()->getPz(),db->GetParticle(211)->Mass());    
      p4_pi2.SetXYZM(pips[0]->par()->getPx(),pips[0]->par()->getPy(),pips[0]->par()->getPz(),db->GetParticle(211)->Mass());    

      // SetLorentzVector(pim,pims[0]);
      auto dipion = p4_pi1+p4_pi2;
      //std::cout << "positive pion" << pips[0]->par()->getPx() << " " << pips[0]->par()->getPy() << " " << pips[0]->par()->getPz() <<std::endl;
      //std::cout << "negative pion " << pims[0]->par()->getPx() << " " << pims[0]->par()->getPy() << " " << pims[0]->par()->getPz() <<std::endl;           
       
      hdipion_mass->Fill(dipion.M());

      //std::cout << pims.size() << " " << pips.size() << std::endl;
      //std::cout << " dipion invariant mass " << dipion.M() << std::endl;  
      miss = beam+target-el-pip-pim;
      //std::cout << " missing mass" << miss.M() << std::endl;
      hdipion_missingmass->Fill(dipion.M()); 
    }

    //           printf("pid = %8d time = %8.3f ec = %8.3f  pcal = %8.3f sf = %8.3f beta = %8.3f\n",
    //              pid,time-starttime,ecEnergy,pcalEnergy,sf, beta);
   
    // std::cout << " /////////////////////////////////////////////// Number of particles " << c12->getDetParticles().size() << std::endl;
    for(auto& p : c12->getDetParticles()){
      //  get predefined selected information
      //std::cout << " PID " << p->getPid() << std::endl;
      p->getTime();
      p->getDetEnergy();
      p->getDeltaEnergy();
      switch(p->getRegion()) {
      case FD :
	//std::cout << p->getCharge() << std::endl;
        h_betap->Fill(p->getP(), p->getBeta());
        //h_vtp->Fill(g->getP(), p->getvt());
	p->cal(PCAL)->getEnergy();
	p->cal(ECIN)->getEnergy();
	p->cal(ECOUT)->getEnergy();
	p->sci(FTOF1A)->getEnergy();
	p->sci(FTOF1B)->getEnergy();
	p->sci(FTOF2)->getEnergy();
	p->trk(DC)->getSector();
	p->che(HTCC)->getNphe();
	p->che(LTCC)->getNphe();
	//trajectories
	p->traj(LTCC)->getX();
	//std::cout << " p->che(HTCC)->getNphe() " << p->che(HTCC)->getNphe() << std::endl;
	//std::cout << "p->che(LTCC)->getNphe(); " << p->che(LTCC)->getNphe() <<std::endl;
	if(p->getPid()==11){
	  //std::cout << "//////////////// energy " << p->getDetEnergy() << std::endl;
	  // std::cout << " p->getP() <<" << p->getP() << std::endl;
	  // std::cout << "  p->cal(PCAL)->getEnergy(); " <<  p->cal(PCAL)->getEnergy() << std::endl;
	  //  std::cout << " p->cal(ECIN)->getEnergy(); " << p->cal(ECIN)->getEnergy() <<std::endl;
	  //  std::cout << " p->cal(ECOUT)->getEnergy(); " <<  p->cal(ECOUT)->getEnergy() <<std::endl;
	  //  std::cout << "all " << p->cal(PCAL)->getEnergy() + p->cal(ECIN)->getEnergy() +  p->cal(ECOUT)->getEnergy() << std::endl;
	  //std::cout <<  p->getP()  << " " << p->getP()/p->getDetEnergy() << std::endl;
            heoverp->Fill( p->getP() ,  1.0*p->getDetEnergy()/p->getP());
	}	
        //std::cout << "p->trk(DC)->getSector()" << p->trk(DC)->getSector() <<std::endl;
        if(p->getPid()==11) hnphe->Fill( p->che(HTCC)->getNphe());
	// p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	break;
      case FT :
	p->ft(FTCAL)->getEnergy();
	p->ft(FTHODO)->getEnergy();
	break;
      case CD:
	p->sci(CTOF)->getEnergy();
	p->sci(CND)->getEnergy();
	break;
      }
      //   covariance matrix (comment in to see!)
      // p->covmat()->print();
      p->cmat();
    }

  }

  TCanvas* can = new TCanvas("c","c",1200,900);
  can->Divide(2);
  can->cd(1);
  hmass.DrawNormalized();
  can->cd(2);
  htime.DrawNormalized();

  can->SaveAs("test_gammas"+runnumber+".pdf");
  can->Clear();
  can->Divide(2);
  can->cd(1);
  hdipion_mass->DrawNormalized();
  can->cd(2);
  hdipion_missingmass->DrawNormalized();
  can->SaveAs("test_chargedpion"+runnumber+".pdf");

  can->Clear();
  hnphe->Draw();
  can->SaveAs("nphe"+runnumber+".pdf");

  can->Clear();
  hQ2W->Draw("colz");
  gPad->SetLogz(1);
  can->SaveAs("Kinematics_Q2W"+runnumber+".pdf");  
  gPad->SetLogz(0);
  can->Clear();
  heoverp->Draw("colz");
  can->SaveAs("heoverp_"+runnumber+".pdf");
  can->Clear();
  hetheta->Draw();
  can->SaveAs("hetheta_"+runnumber+".pdf");
  can->Clear();
  h_betap->Draw();
  h_betap->SaveAs("h_betap_"+runnumber+".pdf");  
  TFile* fout  = new TFile("fout"+runnumber+".root","RECREATE");
  hmass.Write();
  htime.Write();
  hQ2W->Write();
  hnphe->Write();
  hW->Write();
  hmiss->Write();
  heoverp->Write();
  hdipion_mass->Write();
  hdipion_missingmass->Write();
  hetheta->Write();
  h_betap->Write();
  fout->Close();
  
}













