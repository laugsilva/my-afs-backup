#include <iostream>
#include <iomanip>
#include <fstream>
#include "myAnalysis.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TH2.h"
#include "TH3.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "TObjString.h"
#include "TH1F.h"
#include <string>
#include "TChain.h"
//FastJet stuff...
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
//#include <iostream>

#include "TTree.h"

#include "TMVA/Reader.h"

#include "Boost2011_simple.cc"

//ReadTruth stuff???
//#include "ReadTruth.h"


//==========================================================================================
//GoodRunList & TriggerReader
//==========================================================================================
#include "/afs/cern.ch/user/l/laugs/private/DataUtilities/GoodRunsLists-00-00-76/GoodRunsLists/TGoodRunsListReader.h"
#include "/afs/cern.ch/user/l/laugs/private/DataUtilities/GoodRunsLists-00-00-76/GoodRunsLists/TGoodRunsList.h"
#include "/afs/cern.ch/user/l/laugs/private/DataUtilities/GoodRunsLists-00-00-76/GoodRunsLists/DQHelperFunctions.h"

#include "TriggerReader.h"

//==========================================================================================

using namespace std;

ClassImp(myAnalysis) //Before the code body of first class member the ROOT's macro ClassImp() must be included (including class in ROOT?)


vector<int> bptot,bpini,bpfin,bmtot,bmini,bmfin; 
map<int,int> getind;

const int kPC = 11;
const int kFE = 12;
const int kDE = 14;
const int kFA = 15;
const int kGS = 13;

//Event  weight global vars
const float Lumi = 100.;
float Xsection = 99.; //45.6; 
float N0 = 1400000.;//99999.;

//Jet selection cuts
float jetEtaMax=2.1;
float jetPtMin=20.;
float lepPt = 15.0; 
float BtagCut = 5.85;
float DR_cut = 0.4;

//Quality cuts on tracks
double cutOnTrackPT = 1.; 
double cutOnChi2OverNDof = 3.;
int cutNBLayerHits=0; 
int cutNSCTHits=4;
int cutNPixelHits=2;
int cutNBLayerPlusPix = 0;
int cutNSiHits = 7;  
double cutd0PV = 2.0;
double cutz0PVsinTheta = 2.0;

int mc_channel_number = -1;

/* Constructor'n destructor */

myAnalysis::myAnalysis(TTree *tree,TTree* metaDataTree): AnalysisBase(tree) 
{
  cout << endl;
  cout << "------------------------------" << endl;
  cout << " myAnalysis " << endl;
  cout << " Number of events = " << tree->GetEntries() << endl;
  cout << "------------------------------" << endl;
  cout << endl;

 // Trigger info 
  m_confTree = metaDataTree;
  cout << "****  About to call TriggerReader ****** " << endl;
  m_triggerReader = TriggerReader(metaDataTree);
  cout << "****  Called TriggerReader ****** " << endl;
  //TriggerReader m_triggerReader(metaDataTree);


}

myAnalysis::~myAnalysis() {}




//=========================================================================================
// Fill Trees
//=========================================================================================

void myAnalysis::FillTrees(bool isDATA, bool isMC){
  

 
 


  //--------------- Selection cuts (pt/met in GeV!!) --------------
  //    float jetEtaMax=2.1;
  //    float jetPtMin=20;
  //    float lepPt = 15.0; //GeV
  //    float BtagCut = 4.50;//5.85;
  //    float DR_cut = 0.4;
  
  
   //    double cutOnTrackPT = 1.;//0.5;//
   //    double cutOnChi2OverNDof = 3.0; //10000
   //    int cutNBLayerHits=0; //
   //    int cutNSCTHits=4; //6;
   //    int cutNPixelHits=1;//2;
   //    int cutNBLayerPlusPix = 0;
   //    int cutNSiHits = 7;  
   //    double cutIPz=5.;//0.2;
   //    double cutD0 = 1.5;
   //    double cutZ0sinTheta = 1.5;
   //    double cutd0PV = 2.0;//1.5;
   //    double cutz0PVsinTheta = 2.0;
   //    double IPSignificanceCut = 1.;// <- This is a "less than". To test displaced trackas sigma>2.5; //2.; // 3.;    
   //    double IPzSignificanceCut = 5.;
   

  //---------------------- create a new root output file. -----------------

   char name[120];
   if(isMC) sprintf(name,"Trees_MC11_b_QCDjetjet_NO_ISO_r17.root");
   if(isDATA) sprintf(name,"Trees_DATA2011_QCDjetjet_NO_ISO_r17_PeriodGHIroot");

  TFile* file = TFile::Open(name,"RECREATE");

  TTree* mytree = new TTree("mytree","");
 
  mytree->SetDirectory(&*file); 
 
  float _pt, _drmax, _trkwidth,_width, _mass,_dr; // _jetTrackEcc, _jetCaloEcc;
  int _ntrk, _ncl, _nsubjets;
  double _tau1, _tau2, _tauratio, _drktaxis;
  float _weight;
  float _sv0weight;
  float _ip3dsv1weight;
  float _eta;//_phi;

  mytree->Branch("pt",&_pt,"pt/F");
  mytree->Branch("eta",&_eta,"eta/F");
  // mytree->Branch("phi",&_phi,"phi/F");
  //mytree->Branch("dr2closestjet",&_dr,"dr2closestjet/F");
  mytree->Branch("drmax",&_drmax,"drmax/F");
  mytree->Branch("width",&_width,"width/F");
  mytree->Branch("trkwidth",&_trkwidth,"trkwidth/F");
  mytree->Branch("mass",&_mass,"mass/F");
  //   mytree->Branch("jetTrackEcc",&_jetTrackEcc,"jetTrackEcc/F");
  //   mytree->Branch("jetCaloEcc",&_jetCaloEcc,"jetCaloEcc/F");
  mytree->Branch("ntrk",&_ntrk,"ntrk/I");
  mytree->Branch("ncl",&_ncl,"ncl/I");
  mytree->Branch("nsubjets",&_nsubjets,"nsubjets/I");
  mytree->Branch("weight",&_weight,"weight/F");
  mytree->Branch("sv0weight",&_sv0weight,"sv0weight/F");
  mytree->Branch("ip3dsv1weight",&_ip3dsv1weight,"ip3dsv1weight/F");
  mytree->Branch("tau1",&_tau1,"tau1/D");
  mytree->Branch("tau2",&_tau1,"tau2/D");
  mytree->Branch("tauratio",&_tauratio,"tauratio/D");
  mytree->Branch("drktaxis",&_drktaxis,"drktaxis/D");

  



  /**************************************************************/
  //Loop through entries
  /**************************************************************/

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
   nentries = 20000; //0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if(jentry%10000==0) cout << "Event #" << jentry << endl;

    bool isDEBUG = false;//true;//true;
    if(isDEBUG)cout<< __LINE__ <<endl;



  //--------------- Event weight ------------------
 
  int NowRunningOn = -1;
  float w = 1;
   if(isMC){

   NowRunningOn = mc_channel_number;
   bool isMC10b = false;
   bool isMC11 = false;
   bool isMC11b = true;

   if(isMC10b){
     if(NowRunningOn ==0){
      Xsection = 1000*9.8608e+06;
      N0 = 16388258;
    }else if(NowRunningOn==1){
      Xsection = 1000*6.7818e+05;
      N0 = 7377565;
    }else if(NowRunningOn==2){
      Xsection = 1000*4.0982e+04;
      N0 = 2796084;
    }else if(NowRunningOn==3){
      Xsection =1000*2.1929e+03;
      N0 = 2796879;
    }else if(NowRunningOn==4){
      Xsection = 1000*8.7701e+01;
      N0 = 2793179;
    }else if(NowRunningOn==5){
      Xsection = 1000*2.3501;
      N0 = 2790576;
    }else if(NowRunningOn==6){
      Xsection = 1000*3.361e-02;
      N0 = 2790601;
    }
   }else if(isMC11b){
    if(NowRunningOn ==0){
      Xsection = 1000*1.2032e+07 ;//9.8608e+06;
      N0 = 2798297;//91617;//3637;// 1363438;
    }else if(NowRunningOn==1){
      Xsection = 1000*8.0714e+05;//.7818e+05;
      N0 = 2798443;//7377565;//110951;//17111;//1395889;
    }else if(NowRunningOn==2){
      Xsection = 1000*4.8027e+04 ;//4.0982e+04;
      N0 = 2795891;//114381;//52253;//1396991;
    }else if(NowRunningOn==3){
      Xsection =1000*2.5364e+03 ;//2.1929e+03;
      N0 = 2797982;//213163;//99409;//1397590;
    }else if(NowRunningOn==4){
      Xsection = 1000*9.9605e+01;//8.7701e+01;
      N0 = 2797431;//154659;//148610;//1393487;
    }else if(NowRunningOn==5){
      Xsection = 1000*2.5947e+00;
      N0 = 2796405;//386425;//194467;//1392492l;
    }else if(NowRunningOn==6){
      Xsection = 1000*3.5457e-02 ;//3.361e-02;
      N0 = 2791826;//436390;//232837;//1391670;
    } 
    
    
   }

   float w = Xsection*Lumi/N0;
   }
   


   
    //--------------------- Vertex cut -----------------------------------
    unsigned int NEvtVtx =  vxp_n;
    int NGoodVtx = 0;//vxp_n; 
    for(int i=0; i<NEvtVtx; i++) {
      int vtx_ntrks = (*vxp_nTracks)[i];
      if(vtx_ntrks>4)NGoodVtx++;
    }
    

    if(NGoodVtx<1) continue;




    //------------------- MC10a Re-weighting ---------------------------------
    if(isDEBUG)cout << __LINE__ <<endl; 
    if(isMC){
    TLorentzVector lj;
    lj.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
    float Leadpt =lj.Pt(); 
    float mc_channel = mc_channel_number;
    if(mc_channel ==105009 && Leadpt >40) continue;
    else if(mc_channel==105010&& Leadpt >60)continue;
    else if(mc_channel==105011&& Leadpt >110)continue;
    else if(mc_channel==105012&& Leadpt >200)continue;
    else if(mc_channel==105013&& Leadpt >360)continue;
    else if(mc_channel==105014&& Leadpt >620)continue;
    else if(mc_channel==105015&& Leadpt >1200)continue;
    }



    //----------- Trigger Selection---------------------------------
    bool passes = true;
    if(isDATA){
      
      TLorentzVector leadJet;
      leadJet.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
      float LeadPt = leadJet.Pt();
     
      passes = PassedTrigger11 (RunNumber, 
				trig_DB_SMK, 
				trig_DB_L1PSK, 
				trig_DB_HLTPSK, 
				trig_L1_TAV,
				trig_EF_passedPhysics,
				LeadPt);
      }
    if (!passes) continue;

    





    //------------------ Loop through jets ------------------------------------

    int bjets[30] = {0}; int nj = 0; 
    int Nbjets = 0;
    int  jet_AntiKt4TopoEM_num = jet_AntiKt4TopoEM_pt->size();

    for(int i=0; i<jet_AntiKt4TopoEM_num; i++) {
      
     
      TLorentzVector j;
      j.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[i]/1000.,(*jet_AntiKt4TopoEM_eta)[i],(*jet_AntiKt4TopoEM_phi)[i],(*jet_AntiKt4TopoEM_E)[i]/1000.);
    
      //Cleaning cuts
      int goodness = 2;
      if(isDATA){
	goodness = (*jet_AntiKt4TopoEM_isBadLoose)[i];
	if(goodness != 0) continue;
      }

      float PT = j.Pt();
      float ETA = TMath::Abs(j.Eta());
      if(PT<20.|| ETA>jetEtaMax){
        j.Delete();
        continue;
      }
      float PHI = j.Phi();
      
      
     
      //-------------------- match to tracks --------------------------------

      int trk_num = trk_pt->size();
      TLorentzVector trackJet;
      TLorentzVector t1(0,0,0,0);
      TLorentzVector t2(0,0,0,0);
      float sumpt = 0.0; int ntrk = 0;
      float maxpt = 0.0;	
      vector<int> trklist; trklist.clear();
      //to compute jet eccentricity
      float sumEiEtai2 = 0;
      float sumEiPhii2 = 0;
      float sumEiEtaiPhii = 0;
      float jetTrackEccentricity = 99.;
      for(int trk=0; trk<trk_num; trk++) {
	
	if( ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk])>=cutNBLayerPlusPix && (*trk_nSCTHits)[trk]>= cutNSCTHits&&(*trk_nPixHits)[trk]>=cutNPixelHits && ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk]+(*trk_nSCTHits)[trk])>=cutNSiHits  &&(*trk_pt)[trk]/1000> cutOnTrackPT &&((*trk_chi2)[trk]/(*trk_ndof)[trk])<cutOnChi2OverNDof && TMath::Abs((*trk_d0_wrtPV)[trk])<cutd0PV && TMath::Abs((*trk_z0_wrtPV)[trk]*TMath::Sin((*trk_theta)[trk]))<cutz0PVsinTheta ){
	  
	  
	  TLorentzVector track;
	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trk]);
	  track.SetPtEtaPhiE((*trk_pt)[trk]/1000., (*trk_eta)[trk],(*trk_phi_wrtPV)[trk],track_e/1000);

	  float dr = j.DeltaR(track);
	  if(dr<0.4) {
	    if(ntrk==0) t1 = track;
	    if(ntrk==1) t2 = track;
	    sumpt += track.Pt();
	    sumEiEtai2 += track.Energy()*TMath::Power(track.Eta(),2);
	    sumEiPhii2 += track.Energy()*TMath::Power(track.Phi(),2);
	    sumEiEtaiPhii += track.Energy()*track.Phi()*track.Eta();
	    ntrk ++;
	    trackJet += track;
	    // to compute drmax/P1
	    if(track.Pt()>maxpt) maxpt = track.Pt();
	    trklist.push_back(trk);
	  }
	}
      }
      

      if(ntrk<1){
	cout<< " jet with less than 1 track, continue " <<endl;	
	continue;
      }

  
     
      
      
      //---------- Get calo/trk jet width ------------------------------
      float WIDTH = (*jet_AntiKt4TopoEM_WIDTH)[i];
      float trkWIDTH = 0;
      float trkWIDTH_num = 0;
      float trkWIDTH_den = 0;
      if(trklist.size()>1) {
        for(int k1=0; k1<trklist.size()-1; k1++) {
          TLorentzVector t1;
          float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
          t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
          float DR = t1.DeltaR(j);

          trkWIDTH_num+=t1.Pt()*DR;
          trkWIDTH_den+=t1.Pt(); 
          
        }
      }
      trkWIDTH = trkWIDTH_num/trkWIDTH_den;




      //---------- btag weight cut ----------------------------------------
      // bool isBjet = false;
      float Btag_w = (*jet_AntiKt4TopoEM_flavor_weight_SV0)[i]; 
      float Btag_combw = (*jet_AntiKt4TopoEM_flavor_weight_SV1)[i]+(*jet_AntiKt4TopoEM_flavor_weight_IP3D)[i]; 
      //if(Btag_w<BtagCut)continue;
      if(Btag_w>BtagCut){
	bjets[nj++]=i; Nbjets++;     
      }
   

      

     
      //---------- Isolation --------------------------------------------------------
      float closeby_unCalibPtCut = 7.;
      int Ncloseby7GeVJets = 0;
      int NclosebyBJets = 0;     

      TLorentzVector j2(0,0,0,0);
      float tmp_dr = 100.;
      for(int l=0; l<jet_AntiKt4TopoEM_n; l++) {
	
        float newJet_unCalibPt = (*jet_AntiKt4TopoEM_emscale_pt)[l]/1000.;
        TLorentzVector newJet;        
        newJet.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[l]/1000.,(*jet_AntiKt4TopoEM_eta)[l],(*jet_AntiKt4TopoEM_phi)[l],(*jet_AntiKt4TopoEM_E)[l]/1000.);    
	if(newJet==j)continue;
	
	if(newJet_unCalibPt>closeby_unCalibPtCut && j.DeltaR(newJet)<0.8)Ncloseby7GeVJets++;
	if(j.DeltaR(newJet)<0.8&&(*jet_AntiKt4TopoEM_flavor_weight_SV0)[l]<5.85 )NclosebyBJets++;

	float dr = j.DeltaR(newJet);
        if(dr<tmp_dr){
          tmp_dr = dr;
          j2=newJet;
        }
      }
      
      
      //   if(NclosebyBJets>0|| Ncloseby7GeVJets>0)continue; 
      float DRj12 = j.DeltaR(j2);
      

     

    
      //------------------- TEST FASTJETS ------------------------

      //std::cout << "Here we go..." << std::endl;
      fastjet::RecombinationScheme recombScheme = fastjet::pt2_scheme;
      fastjet::JetAlgorithm jetalgorithm = fastjet::kt_algorithm;
      double rParameter = 0.2;
      fastjet::Strategy ktstrategy = fastjet::Best;
      
      //std::cout << "Is this working?" << std::endl;
      fastjet::JetDefinition jetDef(jetalgorithm,rParameter, recombScheme, ktstrategy);
      vector<fastjet::PseudoJet> particles;

       //Loop to tracklist
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {
	  TLorentzVector t;
	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k]]);
	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],track_e/1000.);
	  particles.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
	}
      }
	
      //run the clustering, extract the jets
      fastjet::ClusterSequence cs(particles,jetDef);
      //      vector<PseudoJet> jets = cs.fastjet::inclusive_jets();
      vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
      
      unsigned int nsubjets = jets.size();
      


      
      //---------------- Get N-subjettiness variables -------------

      double Nsubjet1 = 99.;
      double Nsubjet2 = 99.;
      double Nsubjet3 = 99.;
      double Nsubjet12ratio = 99.;
      double DeltaRk2axes = 99.;


      vector<fastjet::PseudoJet> constituents;  // will use charged constituents
      //Loop to tracklist
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {
	  TLorentzVector t;
	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k]]);
	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],track_e/1000.);
	  constituents.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
	}
      }

      //      cout << " About to call Analysis() : tau1 = " << Nsubjet1 << " drktaxis " << DeltaRk2axes<< endl;

      Analysis(constituents, Nsubjet1, Nsubjet2, Nsubjet3,DeltaRk2axes);
      if(Nsubjet1!=99.&&Nsubjet2!=99.&&Nsubjet1!=0)Nsubjet12ratio = Nsubjet2/Nsubjet1;
      constituents.clear();

      //cout << " tau1 " << Nsubjet1 << " tau2 "<< Nsubjet2 << " tauratio " << Nsubjet12ratio << " DeltaRk2axes "<<DeltaRk2axes<<endl;




    
      //---------------- Get compute drmax ------------------------------

      float drmax = -1.0; 
      if(trklist.size()>1) {
	for(int k1=0; k1<trklist.size()-1; k1++) {
	  TLorentzVector t1;
	  float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
	  t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
	  for(int k2=k1+1; k2<trklist.size(); k2++) {
	    TLorentzVector t2;
	    float t2_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k2]]);
	    t2.SetPtEtaPhiE((*trk_pt)[trklist[k2]]/1000.,(*trk_eta)[trklist[k2]],(*trk_phi_wrtPV)[trklist[k2]],t2_e/1000.);
	    float DR = t1.DeltaR(t2);
	    if(DR>drmax) drmax = DR;
	  }
	}
      }
      

      //------------------- Get DR between two hardest tracks ---------------------

      TLorentzVector track1(0,0,0,0);
      TLorentzVector track2(0,0,0,0);
      int track1_index = 0;
      int track2_index = 0;
      float track1_pt = 0.;
      float track2_pt = 0.;
      float DeltaR_trk12 = 0.; 
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {    
	  if ((*trk_pt)[trklist[k]]/1000. > track1_pt) {
	    if (k>0) {
	      track2_pt  = track1_pt;
	      track2_index = track1_index;
	    }
	    track1_pt  = (*trk_pt)[trklist[k]]/1000.;
	    track1_index = trklist[k];
	  }
	  else if ( (*trk_pt)[trklist[k]]/1000.> track2_pt) {
	    track2_pt  = (*trk_pt)[trklist[k]]/1000.;
	    track2_index = trklist[k];
	  }
	}
      }
      
       float track1_e = TMath::Abs( 1./(*trk_qoverp)[track1_index]);       float track2_e = TMath::Abs( 1./(*trk_qoverp)[track2_index]);

      track1.SetPtEtaPhiE((*trk_pt)[track1_index]/1000.,(*trk_eta)[track1_index],(*trk_phi_wrtPV)[track1_index],track1_e/1000.);
      track2.SetPtEtaPhiE((*trk_pt)[track2_index]/1000.,(*trk_eta)[track2_index],(*trk_phi_wrtPV)[track2_index],track2_e/1000.);
      
      if(trklist.size()>1) DeltaR_trk12 = track1.DeltaR(track2); 
      
      
     
      //--------------- Get Jet mass ------------------------------------------

      float MASS = j.M();
     


      //--------------- Fill tree... ---------------------------------------------
      _drmax = drmax;
      _pt = PT;
      _eta = ETA;
      // _dr = DRj12;
      _width = WIDTH;
      _trkwidth = trkWIDTH;
      _mass = MASS;
      //       _jetTrackEcc = jetTrackEccentricity;
      //       _jetCaloEcc = jetCaloEccentricity;
      _ntrk = ntrk;
      int ncl = 99.;     
      _ncl = ncl;
      _nsubjets = nsubjets;
      _weight = w;
      _sv0weight = Btag_w;
      _ip3dsv1weight = Btag_combw;
      _tau1 = Nsubjet1;
      _tau2 = Nsubjet2;
      _tauratio = Nsubjet12ratio;
      _drktaxis =  DeltaRk2axes;
      
    
      mytree->Fill();
    
	
     
      
    }//jets
    
  }//entries
  
  //file.cd();
  file->cd();
  mytree->Write();
  //file.Close();
  file->Close();
  
 }//FillTree
	


//=========================================================================================
// Make histograms of input variables
//=========================================================================================
void myAnalysis::MakeHistos(bool isDATA, bool isMC){

  bool isdebug = false;
  if(isdebug)cout<< __LINE__ <<endl;

  //--------- Event  weight global vars ----------------------------
  
 //  const float Lumi = 100.;
//   float Xsection = 99.; //45.6; 
  
  
  //----------- Event Selection cuts (pt/met in GeV!!) --------------
 //  float jetEtaMax=2.1;
//   float jetPtMin=20.;
//   float lepPt = 15.0; //GeV
//   //float crackEtaMin = 1.37;
//   //float crackEtaMax = 1.52;

//   float BtagCut = 5.85;//5.72;//6.4;
//   float DR_cut = 0.4;
  
  //---------- Quality cuts on tracks -------------------------------

  // double cutOnTrackPT = 1.; 
//   double cutOnChi2OverNDof = 3.;
//   int cutNBLayerHits=0; 
//   int cutNSCTHits=4;
//   int cutNPixelHits=2;
//   int cutNBLayerPlusPix = 0;
//   int cutNSiHits = 7;  
//   double cutd0PV = 2.0;//1.5;
//   double cutz0PVsinTheta = 2.0;
//   //double IPSignificanceCut = 1.;// <- This is a "less than". To test displaced trackas sigma>2.5; //2.; // 3.;    
//   //double IPzSignificanceCut = 5.; //2.; // 3.;   

//   if(isdebug)cout<< __LINE__ <<endl;
  
  //------- PDGids for identifying B hadrons------------------------
  int NumBpdgids = 4;// 24;
  int B_pdgids[NumBpdgids];
  for(int i=0;i<NumBpdgids;i++)B_pdgids[i]=0;
  B_pdgids[0]=511;
  B_pdgids[1]=521;
  B_pdgids[2]=531;
  B_pdgids[3]=541;
  B_pdgids[4]=5122;
  B_pdgids[5]=5132;

  if(isdebug)cout<< __LINE__ <<endl;

  //------------------  Book histograms ---------------------------


  int nbins = 100; 
  TString separator = "_";


  // PT bins from Trigger Turn-ons

  const int N_PTbins = 8;

  float pt1[N_PTbins]={40.,60.,80.,110.,150.,200.,270.,360}; 
  float pt2[N_PTbins]={60.,80.,110.,150.,200.,270.,360,480}; 
  int Pt1[N_PTbins]={40.,60.,80.,110.,150.,200.,270.,360};


  //Tracking variables histos

  TH1F* d0PV_histo = new TH1F("d0PV_histo","",100,-5,5); 
  TH1F* z0PV_histo = new TH1F("z0PV_histo","",100,-25,25); 
  TH1F* z0PVsinTheta_histo = new TH1F("z0PVsinTheta_histo","",100,-5,5); 

  
  //Pt and eta distributions for all btagged jets 
  float minPt = 0.;
  float maxPt = 245.;
  float minEta = 0.;
  float maxEta = 2.0;

  TH1F* Pt_histos = new TH1F("Pt_histos","",300,minPt,600.); 
  Pt_histos->Sumw2();
  TH1F* Eta_histos = new TH1F("Eta_histos","",nbins,minEta,maxEta);
  Eta_histos->Sumw2();

  TH1F* Pt_L1J10 = new TH1F("Pt_L1J10","",nbins,minPt,200.); 
  Pt_L1J10->Sumw2();


  // DeltaR / DeltaPhi plots
  TH1F* DR_BB_all = new TH1F("DR_BB_all","",50,0.,4.);
  TH1F* DR_BB_allinEta = new TH1F("DR_BB_allinEta","",50,0.,4.);
  TH1F* DR_BB_allMatching = new TH1F("DR_BB_allMatching","",50,0.,4.);
  TH1F* DR_BB_merged = new TH1F("DR_BB_merged","",50,0.,4.);
  TH1F* DR_BB_single = new TH1F("DR_BB_single","",50,0.,4.);
  TH1F* DPhi_BB_all = new TH1F("DPhi_BB_all","",50,0.,4.);
  TH1F* DPhi_BB_allinEta = new TH1F("DPhi_BB_allinEta","",50,0.,4.);
  TH1F* DPhi_BB_allMatching = new TH1F("DPhi_BB_allMatching","",50,0.,4.);
  TH1F* DPhi_BB_merged = new TH1F("DPhi_BB_merged","",50,0.,4.);
  TH1F* DPhi_BB_single = new TH1F("DPhi_BB_single","",50,0.,4.);

  TH1F* DPhi_bjets = new TH1F("DPhi_bjets","",50,0.,4.);
  TH1F* DR_bjets = new TH1F("DR_bjets","",50,0.,4.);


  // Number of Bjets
  TH1F* NumBjets_histo = new TH1F("NumBjets_histo","",6,0,6);
  TH1F* NumIsolatedBjets_histo = new TH1F("NumIsolatedBjets_histo","",6,0,6);
  TH1F* NumNonIsolatedBjets_histo = new TH1F("NumNonIsolatedBjets_histo","",6,0,6);
  TH1F* NumIsolatedSingleBjets_histo = new TH1F("NumIsolatedSingleBjets_histo","",6,0,6);
  TH1F* NumIsolatedMergedBjets_histo = new TH1F("NumIsolatedMergedBjets_histo","",6,0,6);
  
  //counters
  int NumBjets = 0;
  int NumNonIsolatedBjets = 0;
  int NumIsolatedBjets = 0;
  int NumIsolatedSingleBjets =0;
  int NumIsolatedMergedBjets = 0;




  //------------- N-subjettiness histos ----------------------------

  TH1F* Tau1_histos[N_PTbins]; 
  TH1F* Tau2_histos[N_PTbins]; 
  TH1F* Tau12ratio_histos[N_PTbins]; 
  TH1F* kt2axisDeltaR_histos[N_PTbins];
  
  //correlations
  TH2F* Tau2Tau1Corr_histos[N_PTbins];
  TH2F* Tau1DRkT2axesCorr_histos[N_PTbins]; 
  TH2F* Tau2DRkT2axesCorr_histos[N_PTbins]; 
 for(int j=0;j<N_PTbins;j++){  
      Tau1_histos[j]=new TH1F("Tau1"+separator+"PT"+Form("%i",Pt1[j]),"Tau1"+separator+"PT"+Form("%i",Pt1[j]),50,-0.01,1.01);
      Tau1_histos[j]->Sumw2();
      Tau2_histos[j]=new TH1F("Tau2"+separator+"PT"+Form("%i",Pt1[j]),"Tau2"+separator+"PT"+Form("%i",Pt1[j]),50,-0.01,1.01);
      Tau2_histos[j]->Sumw2();

      Tau12ratio_histos[j]=new TH1F("Tau12ratio"+separator+"PT"+Form("%i",Pt1[j]),"Tau12ratio"+separator+"PT"+Form("%i",Pt1[j]),50,-0.01,2.01);
      Tau12ratio_histos[j]->Sumw2();
      
      kt2axisDeltaR_histos[j]= new TH1F("Kt2axisDeltaR"+separator+"PT"+Form("%i",Pt1[j]),"Kt2axisDeltaR"+separator+"PT"+Form("%i",Pt1[j]),50,0.,0.9);
      kt2axisDeltaR_histos[j]->Sumw2();

      Tau2Tau1Corr_histos[j]=new TH2F("Tau2Tau1Corr"+separator+"PT"+Form("%i",Pt1[j]),"Tau12ratio"+separator+"PT"+Form("%i",Pt1[j]),50,-0.01,1.01,50,-0.01,1.01);
      Tau1DRkT2axesCorr_histos[j]=new TH2F("Tau1DRkT2axesCorr"+separator+"PT"+Form("%i",Pt1[j]),"Tau1DRkT2axesCorr"+separator+"PT"+Form("%i",Pt1[j]),50,0.,0.9,50,-0.01,1.01);
      Tau2DRkT2axesCorr_histos[j]=new TH2F("Tau2DRkT2axesCorr"+separator+"PT"+Form("%i",Pt1[j]),"Tau2DRkT2axesCorr"+separator+"PT"+Form("%i",Pt1[j]),50,0.,0.9,50,-0.01,1.01);
 
 }

 //=====================================
 // pt histos for each bin of PT 
  TH1F* pt_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    pt_histos[j]=new TH1F("pt"+separator+"PT"+Form("%i",Pt1[j]),"pt"+separator+"PT"+Form("%i",Pt1[j]),200,minPt,400.);
    pt_histos[j]->Sumw2();
    pt_histos[j]->SetTitleSize(5.);
  }

  //===============================
  // Ntrk histos for each bin of PT 
  TH1F* Ntrk_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    Ntrk_histos[j]=new TH1F("Ntrk"+separator+"PT"+Form("%i",Pt1[j]),"Ntrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25);
    Ntrk_histos[j]->SetTitleSize(5.);
    Ntrk_histos[j]->Sumw2();
  }
  
  //====================================
  // DRmax histos for each bin of PT
  TH1F* DRmax_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRmax_histos[j]=new TH1F("DRmax"+separator+"PT"+Form("%i",Pt1[j]),"DRmax"+separator+"PT"+Form("%i",Pt1[j]),30,0.,1.);
    DRmax_histos[j]->Sumw2();
    DRmax_histos[j]->SetTitleSize(5.);
  }
  
  
  //====================================
  // DR_trk12 histos for each bin of PT
  TH1F* DRtrk12_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRtrk12_histos[j]=new TH1F("DRtrk12"+separator+"PT"+Form("%i",Pt1[j]),"DRtrk12"+separator+"PT"+Form("%i",Pt1[j]),50,0.,1.);
    DRtrk12_histos[j]->Sumw2();
    DRtrk12_histos[j]->SetTitleSize(5.);
  }

  
  //====================================
  // Width histos for each bin of PT
  TH1F* Width_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    Width_histos[j]=new TH1F("Width"+separator+"PT"+Form("%i",Pt1[j]),"Width"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4);
    Width_histos[j]->SetTitleSize(5.);
    Width_histos[j]->Sumw2();
  }


  
  //====================================
  // trkWidth histos for each bin of PT
  TH1F* trkWidth_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    trkWidth_histos[j]=new TH1F("trkWidth"+separator+"PT"+Form("%i",Pt1[j]),"Width"+separator+"PT"+Form("%i",Pt1[j]),100,0,0.8);//50,0,0.4);
    trkWidth_histos[j]->SetTitleSize(5.);
    trkWidth_histos[j]->Sumw2();
  }

  
  //====================================
  // NClus histos for each bin of PT
  TH1F* NClus_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    NClus_histos[j]=new TH1F("NClus"+separator+"PT"+Form("%i",Pt1[j]),"NClus"+separator+"PT"+Form("%i",Pt1[j]),40,0,40);
    NClus_histos[j]->SetTitleSize(5.);
    NClus_histos[j]->Sumw2();
  }
  
  //====================================
  // YFlip12 histos for each bin of PT
  TH1F* YFlip12_histos[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    YFlip12_histos[j]=new TH1F("YFlip12"+separator+"PT"+Form("%i",Pt1[j]),"YFlip12"+separator+"PT"+Form("%i",Pt1[j]),50,0.,2e-08); 
    YFlip12_histos[j]->SetTitleSize(5.);
    YFlip12_histos[j]->Sumw2();
  }
  
  // //====================================
  //   //Yscale histos for each bin of PT
  //   TH1F* Yscale_histos[N_PTbins];
  //   for(int k =0;k<n;k++){
  //     for(int j=0;j<N_PTbins;j++){  
  //       Yscale_histos[j]=new TH1F("Yscale"+separator+"PT"+Form("%i",Pt1[j]),"Yscale"+separator+"PT"+Form("%i",Pt1[j]),50,0.,1e-02);//-5e-8,4e-07); 
  //     Yscale_histos[j]->SetTitleSize(5.);
  //     }
  //   }
  
  //====================================
  //NumVtx histos for each bin of PT
  TH1F* NumVtx_histos[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    NumVtx_histos[j]=new TH1F("NumVtx"+separator+"PT"+Form("%i",Pt1[j]),"NumVtx"+separator+"PT"+Form("%i",Pt1[j]),15,0.,15);//-5e-8,4e-07); 
    NumVtx_histos[j]->SetTitleSize(5.);
    NumVtx_histos[j]->Sumw2();
  }

  
  //====================================
  //NumSubTrackJets histos for each bin of PT
  TH1F* NumSubTrackJets_histos[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    NumSubTrackJets_histos[j]=new TH1F("NumSubTrackJets"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJets"+separator+"PT"+Form("%i",Pt1[j]),15,0.,15);//-5e-8,4e-07); 
    NumSubTrackJets_histos[j]->SetTitleSize(5.);
    NumSubTrackJets_histos[j]->Sumw2();
  }

  //====================================
  //JetMass histos for each bin of PT
  TH1F* JetMass_histos[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    JetMass_histos[j]=new TH1F("JetMass"+separator+"PT"+Form("%i",Pt1[j]),"JetMass"+separator+"PT"+Form("%i",Pt1[j]),30,0.,20.);//-5e-8,4e-07); 
    JetMass_histos[j]->SetTitleSize(5.);
    JetMass_histos[j]->Sumw2();
  }

  //====================================
  //JetTrackEccentricity histos for each bin of PT
  TH1F* JetTrackEccentricity_histos[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    JetTrackEccentricity_histos[j]=new TH1F("JetTrackEccentricity"+separator+"PT"+Form("%i",Pt1[j]),"JetTrackEccentricity"+separator+"PT"+Form("%i",Pt1[j]),20,0.,1.);//-5e-8,4e-07); 
    JetTrackEccentricity_histos[j]->SetTitleSize(5.);
    JetTrackEccentricity_histos[j]->Sumw2();
  }

  //====================================
  //JetCaloEccentricity histos for each bin of PT
  TH1F* JetCaloEccentricity_histos[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    JetCaloEccentricity_histos[j]=new TH1F("JetCaloEccentricity"+separator+"PT"+Form("%i",Pt1[j]),"JetCaloEccentricity"+separator+"PT"+Form("%i",Pt1[j]),20,0.,1.);//-5e-8,4e-07); 
    JetCaloEccentricity_histos[j]->SetTitleSize(5.);
    JetCaloEccentricity_histos[j]->Sumw2();
  }


  
  //===================================
  // Correlation histos
  //==================================
  
  //===========
  // Width
  TH2F* WidthEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    WidthEta_Corr[j]=new TH2F("WidthEta"+separator+"PT"+Form("%i",Pt1[j]),"WidthEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,50,0,0.4);
    WidthEta_Corr[j]->SetTitleSize(5.);
  }
  
  TH2F* WidthNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    WidthNtrk_Corr[j]=new TH2F("WidthNtrk"+separator+"PT"+Form("%i",Pt1[j]),"WidthNtrk"+separator+"PT"+Form("%i",Pt1[j]),20,0,20,50,0,0.4);
    WidthNtrk_Corr[j]->SetTitleSize(5.);
  }
  

  //===========
  // trkWidth
  TH2F* trkWidthEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    trkWidthEta_Corr[j]=new TH2F("trkWidthEta"+separator+"PT"+Form("%i",Pt1[j]),"trkWidthEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,50,0,0.4);
    trkWidthEta_Corr[j]->SetTitleSize(5.);
  }
  
  TH2F* trkWidthNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    trkWidthNtrk_Corr[j]=new TH2F("trkWidthNtrk"+separator+"PT"+Form("%i",Pt1[j]),"trkWidthNtrk"+separator+"PT"+Form("%i",Pt1[j]),20,0,20,50,0,0.4);
    trkWidthNtrk_Corr[j]->SetTitleSize(5.);
  }
  
  
  //===========
  // DRmax
  TH2F* DRmaxEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRmaxEta_Corr[j]=new TH2F("DRmaxEta"+separator+"PT"+Form("%i",Pt1[j]),"DRmaxEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,30,0,1.);
    DRmaxEta_Corr[j]->SetTitleSize(5.);
  }
  
  
  TH2F* DRmaxNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRmaxNtrk_Corr[j]=new TH2F("DRmaxNtrk"+separator+"PT"+Form("%i",Pt1[j]),"DRmaxNtrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25,30,0,1.);
    DRmaxNtrk_Corr[j]->SetTitleSize(5.);
  }
  
  
  TH2F* DRmaxWidth_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRmaxWidth_Corr[j]=new TH2F("DRmaxWidth"+separator+"PT"+Form("%i",Pt1[j]),"DRmaxWidth"+separator+"PT"+Form("%i",Pt1[j]),50,0.,0.4,30,0,1.);
    DRmaxWidth_Corr[j]->SetTitleSize(5.);
  }

  
  //===========
  // DRtrk12
  TH2F* DRtrk12Eta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRtrk12Eta_Corr[j]=new TH2F("DRtrk12Eta"+separator+"PT"+Form("%i",Pt1[j]),"DRtrk12Eta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,50,0,1.);
    DRtrk12Eta_Corr[j]->SetTitleSize(5.);
  }
  
  
  TH2F* DRtrk12Ntrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRtrk12Ntrk_Corr[j]=new TH2F("DRtrk12Ntrk"+separator+"PT"+Form("%i",Pt1[j]),"DRtrk12Ntrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25,50,0,1.);
    DRtrk12Ntrk_Corr[j]->SetTitleSize(5.);
  }
  
  
  TH2F* DRtrk12Width_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    DRtrk12Width_Corr[j]=new TH2F("DRtrk12Width"+separator+"PT"+Form("%i",Pt1[j]),"DRtrk12Width"+separator+"PT"+Form("%i",Pt1[j]),50,0.,0.4,50,0,1.);
    DRtrk12Width_Corr[j]->SetTitleSize(5.);
  }

  

  //=========================
  // Jet cluster multiplicity 
  
  TH2F* NclusEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NclusEta_Corr[j]=new TH2F("NclusEta"+separator+"PT"+Form("%i",Pt1[j]),"NclusEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,40,0,40);
      NclusEta_Corr[j]->SetTitleSize(5.);
  }
  
  TH2F* NclusNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NclusNtrk_Corr[j]=new TH2F("NclusNtrk"+separator+"PT"+Form("%i",Pt1[j]),"NclusNtrk"+separator+"PT"+Form("%i",Pt1[j]),20,0,20,40,0,40);
      NclusNtrk_Corr[j]->SetTitleSize(5.);
  }

  TH2F* NclusWidth_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NclusWidth_Corr[j]=new TH2F("NclusWidth"+separator+"PT"+Form("%i",Pt1[j]),"NclusWidth"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4,40,0,40);
      NclusWidth_Corr[j]->SetTitleSize(5.);
  }
  

  TH2F* NclusDRmax_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    NclusDRmax_Corr[j]=new TH2F("NclusDRmax"+separator+"PT"+Form("%i",Pt1[j]),"NclusDRmax"+separator+"PT"+Form("%i",Pt1[j]),30,0.,1.,40,0,40);
    NclusDRmax_Corr[j]->SetTitleSize(5.);
  }
  
  
  TH2F* NclusDRtrk12_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    NclusDRtrk12_Corr[j]=new TH2F("NclusDRtrk12"+separator+"PT"+Form("%i",Pt1[j]),"NclusDRtrk12"+separator+"PT"+Form("%i",Pt1[j]),50,0.,1.,40,0,40);
    NclusDRtrk12_Corr[j]->SetTitleSize(5.);
  }

  

  TH2F* NclusYFlip12_Corr[N_PTbins];
  for(int j=0;j<N_PTbins;j++){  
    NclusYFlip12_Corr[j]=new TH2F("NclusYFlip12"+separator+"PT"+Form("%i",Pt1[j]),"NclusYFlip12"+separator+"PT"+Form("%i",Pt1[j]),50,0.,2e-09,40,0,40);//-5e-8,4e-07);
    NclusYFlip12_Corr[j]->SetTitleSize(5.);
  }

 //  TH2F* NclusYscale_Corr[N_PTbins];
//   for(int k =0;k<n;k++){
//     for(int j=0;j<N_PTbins;j++){  
//       NclusYscale_Corr[j]=new TH2F("NclusYscale"+separator+"PT"+Form("%i",Pt1[j]),"NclusYscale"+separator+"PT"+Form("%i",Pt1[j]),50,0.,1e-02,40,0,40);//-5e-8,4e-07);
//       NclusYscale_Corr[j]->SetTitleSize(5.);
//     }
//   }

  //=========================
  // Jet sub-track-jet multiplicity 
  
  TH2F* NumSubTrackJetsEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumSubTrackJetsEta_Corr[j]=new TH2F("NumSubTrackJetsEta"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJetsEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,15,0,15);
      NumSubTrackJetsEta_Corr[j]->SetTitleSize(5.);
  }
  TH2F* NumSubTrackJetsNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumSubTrackJetsNtrk_Corr[j]=new TH2F("NumSubTrackJetsNtrk"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJetsNtrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25,15,0,15);
      NumSubTrackJetsNtrk_Corr[j]->SetTitleSize(5.);
  }
  TH2F* NumSubTrackJetsWidth_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumSubTrackJetsWidth_Corr[j]=new TH2F("NumSubTrackJetsWidth"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJetsWidth"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4,15,0,15);
      NumSubTrackJetsWidth_Corr[j]->SetTitleSize(5.);
  }
  TH2F* NumSubTrackJetsDRmax_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumSubTrackJetsDRmax_Corr[j]=new TH2F("NumSubTrackJetsDRmax"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJetsDRmax"+separator+"PT"+Form("%i",Pt1[j]),30,0.,1.,15,0,15);
      NumSubTrackJetsDRmax_Corr[j]->SetTitleSize(5.);
  }
  TH2F* NumSubTrackJetsNclus_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumSubTrackJetsNclus_Corr[j]=new TH2F("NumSubTrackJetsNclus"+separator+"PT"+Form("%i",Pt1[j]),"NumSubTrackJetsNclus"+separator+"PT"+Form("%i",Pt1[j]),40,0,40,15,0,15);
      NumSubTrackJetsNclus_Corr[j]->SetTitleSize(5.);
  }


  //=========================
  // Jet Mass 
  
  TH2F* JetMassEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassEta_Corr[j]=new TH2F("JetMassEta"+separator+"PT"+Form("%i",Pt1[j]),"JetMassEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,30,0.,20.);
    JetMassEta_Corr[j]->SetTitleSize(5.);
  }
  TH2F* JetMassNtrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassNtrk_Corr[j]=new TH2F("JetMassNtrk"+separator+"PT"+Form("%i",Pt1[j]),"JetMassNtrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25,30,0.,20.);
    JetMassNtrk_Corr[j]->SetTitleSize(5.);
  }
  TH2F* JetMassWidth_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassWidth_Corr[j]=new TH2F("JetMassWidth"+separator+"PT"+Form("%i",Pt1[j]),"JetMassWidth"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4,30,0.,20.);
    JetMassWidth_Corr[j]->SetTitleSize(5.);
  }
  TH2F* JetMassDRmax_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassDRmax_Corr[j]=new TH2F("JetMassDRmax"+separator+"PT"+Form("%i",Pt1[j]),"JetMassDRmax"+separator+"PT"+Form("%i",Pt1[j]),30,0.,1.,30,0.,20.);
    JetMassDRmax_Corr[j]->SetTitleSize(5.);
  }
  TH2F* JetMassNclus_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassNclus_Corr[j]=new TH2F("JetMassNclus"+separator+"PT"+Form("%i",Pt1[j]),"JetMassNclus"+separator+"PT"+Form("%i",Pt1[j]),40,0,40,30,0.,20.);
    JetMassNclus_Corr[j]->SetTitleSize(5.);
  }
  TH2F* JetMassNumSubTrackJets_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    JetMassNumSubTrackJets_Corr[j]=new TH2F("JetMassNumSubTrackJets"+separator+"PT"+Form("%i",Pt1[j]),"JetMassNumSubTrackJets"+separator+"PT"+Form("%i",Pt1[j]),15,0,15,30,0.,20.);
    JetMassNumSubTrackJets_Corr[j]->SetTitleSize(5.);
  }
  

  //============================================
  // Jet YFlip12 
  TH2F* YFlip12Eta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    YFlip12Eta_Corr[j]=new TH2F("YFlip12Eta"+separator+"PT"+Form("%i",Pt1[j]),"YFlip12Eta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,50,0.,2e-09);
    YFlip12Eta_Corr[j]->SetTitleSize(5.);
  }
  
  TH2F* YFlip12Ntrk_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){
    YFlip12Ntrk_Corr[j]=new TH2F("YFlip12Ntrk"+separator+"PT"+Form("%i",Pt1[j]),"YFlip12Ntrk"+separator+"PT"+Form("%i",Pt1[j]),20,0,20,50,0.,2e-09);
    YFlip12Ntrk_Corr[j]->SetTitleSize(5.);
  }

  TH2F* YFlip12Width_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      YFlip12Width_Corr[j]=new TH2F("YFlip12Width"+separator+"PT"+Form("%i",Pt1[j]),"YFlip12Width"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4,50,0.,2e-09);
      YFlip12Width_Corr[j]->SetTitleSize(5.);
  }



  //YScale????!??!?!
  
  
  //=========================
  // vtx multiplicity 
  
  TH2F* NumVtxEta_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
      NumVtxEta_Corr[j]=new TH2F("NumVtxEta"+separator+"PT"+Form("%i",Pt1[j]),"NumVtxEta"+separator+"PT"+Form("%i",Pt1[j]),100,-0.1,2.1,15,0,15);
      NumVtxEta_Corr[j]->SetTitleSize(5.);
  }
  
  TH2F* NumVtxNtrk_Corr[N_PTbins]; 
    for(int j=0;j<N_PTbins;j++){  
      NumVtxNtrk_Corr[j]=new TH2F("NumVtxNtrk"+separator+"PT"+Form("%i",Pt1[j]),"NumVtxNtrk"+separator+"PT"+Form("%i",Pt1[j]),25,0,25,15,0,15);
      NumVtxNtrk_Corr[j]->SetTitleSize(5.);
    }
  

  TH2F* NumVtxWidth_Corr[N_PTbins]; 
  for(int j=0;j<N_PTbins;j++){  
    NumVtxWidth_Corr[j]=new TH2F("NumVtxWidth"+separator+"PT"+Form("%i",Pt1[j]),"NumVtxWidth"+separator+"PT"+Form("%i",Pt1[j]),50,0,0.4,15,0,15);
    NumVtxWidth_Corr[j]->SetTitleSize(5.);
  }
  


  //--------------------------------------------------------------------------------
  //--------------------------Loop through entries---------------------------------
  //--------------------------------------------------------------------------------
  if(isdebug)cout<< __LINE__ <<endl;
  
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("RunNumber",1);
   if(isMC)fChain->SetBranchStatus("mc_channel_number",1);
   fChain->SetBranchStatus("lbn",1);
   fChain->SetBranchStatus("trig_L1_TAV",1);
   fChain->SetBranchStatus("trig_L2_passedPhysics",1);
   fChain->SetBranchStatus("trig_EF_passedPhysics",1);
   fChain->SetBranchStatus("trig_DB_SMK",1);
   fChain->SetBranchStatus("trig_DB_L1PSK",1);
   fChain->SetBranchStatus("trig_DB_HLTPSK",1);
   fChain->SetBranchStatus("vxp_n",1);
   fChain->SetBranchStatus("vxp_nTracks",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_pt",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_emscale_pt",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_n",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_eta",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_phi",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_E",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_flavor_weight_SV0",1);
   fChain->SetBranchStatus("jet_AntiKt4TopoEM_isBadLoose",1);
   //if(isDATA)fChain->SetBranchStatus("jet_AntiKt4TopoEM_WIDTH",1);

  
   
   fChain->SetBranchStatus("trk_n",1);
   
   fChain->SetBranchStatus("trk_theta",1);
   fChain->SetBranchStatus("trk_qoverp",1);
   fChain->SetBranchStatus("trk_pt",1);
   fChain->SetBranchStatus("trk_eta",1);
   fChain->SetBranchStatus("trk_d0_wrtPV",1);
   fChain->SetBranchStatus("trk_z0_wrtPV",1);
   fChain->SetBranchStatus("trk_phi_wrtPV",1);
   //fChain->SetBranchStatus("trk_cov_d0_wrtPV",1);
   //fChain->SetBranchStatus("trk_cov_z0_wrtPV",1);
   fChain->SetBranchStatus("trk_chi2",1);
   fChain->SetBranchStatus("trk_ndof",1);
   fChain->SetBranchStatus("trk_nBLHits",1);
   fChain->SetBranchStatus("trk_nPixHits",1);
   fChain->SetBranchStatus("trk_nSCTHits",1);
  
   // link local variable to branch
   if(isdebug)cout<< __LINE__ <<endl;
   
   int  RunNumber = 0;
   if(isMC) int mc_channel_number =0;
   int  lbn = 0;
   vector<unsigned int> *trig_L1_TAV;
   vector<short>   *trig_EF_passedPhysics;
   UInt_t          trig_DB_SMK =0;
   UInt_t          trig_DB_L1PSK=0;
   UInt_t          trig_DB_HLTPSK=0;
   int      vxp_n = 0;
   vector<int>     *vxp_nTracks;
   int              jet_AntiKt4TopoEM_n = 0;
   vector<float>   *jet_AntiKt4TopoEM_E;
   vector<float>   *jet_AntiKt4TopoEM_pt;
   vector<float>   *jet_AntiKt4TopoEM_emscale_pt;
    vector<float>   *jet_AntiKt4TopoEM_eta;
   vector<float>   *jet_AntiKt4TopoEM_phi;
   vector<float>  *jet_AntiKt4TopoEM_flavor_weight_SV0;
   vector<int>     *jet_AntiKt4TopoEM_isBadLoose;
   int           trk_n =0;
   vector<float>   *trk_theta;
   vector<float>   *trk_qoverp;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_d0_wrtPV;
   vector<float>   *trk_z0_wrtPV;
   vector<float>   *trk_phi_wrtPV;
   // vector<float>   *trk_cov_d0_wrtPV;
   // vector<float>   *trk_cov_z0_wrtPV;
   vector<float>   *trk_chi2;
   vector<int>     *trk_ndof;
   vector<int>     *trk_nBLHits;
   vector<int>     *trk_nPixHits;
   vector<int>     *trk_nSCTHits;

   if(isdebug)cout<< __LINE__ <<endl;
   vxp_nTracks = 0;
   jet_AntiKt4TopoEM_E = 0;
   jet_AntiKt4TopoEM_pt = 0;
   jet_AntiKt4TopoEM_emscale_pt = 0;
   jet_AntiKt4TopoEM_eta = 0;
   jet_AntiKt4TopoEM_phi = 0;
   jet_AntiKt4TopoEM_flavor_weight_SV0 = 0;
   jet_AntiKt4TopoEM_isBadLoose=0;
   trk_theta = 0;
   trk_qoverp = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_d0_wrtPV = 0;
   trk_z0_wrtPV = 0;
   trk_phi_wrtPV = 0;
   //   trk_cov_d0_wrtPV = 0;
   //trk_cov_z0_wrtPV = 0;
   trk_chi2 = 0;
   trk_ndof = 0;
   trk_nBLHits = 0;
   trk_nPixHits = 0;
   trk_nSCTHits = 0;


   if(isdebug)cout<< __LINE__ <<endl;

   fChain->SetBranchAddress("RunNumber", &RunNumber);
   if(isMC)fChain->SetBranchAddress("mc_channel_number", &mc_channel_number);
   //fChain->SetBranchAddress("lbn", &lbn);
   fChain->SetBranchAddress("lbn", &lbn);
   fChain->SetBranchAddress("trig_L1_TAV",&trig_L1_TAV);
   fChain->SetBranchAddress("trig_L2_passedPhysics",&trig_L2_passedPhysics);
   fChain->SetBranchAddress("trig_EF_passedPhysics",&trig_EF_passedPhysics);
   fChain->SetBranchAddress("trig_DB_SMK",&trig_DB_SMK);
   fChain->SetBranchAddress("trig_DB_L1PSK",&trig_DB_L1PSK);
   fChain->SetBranchAddress("trig_DB_HLTPSK",&trig_DB_HLTPSK);
   fChain->SetBranchAddress("vxp_n",&vxp_n );
   fChain->SetBranchAddress("vxp_nTracks",&vxp_nTracks );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_pt",&jet_AntiKt4TopoEM_pt );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_pt",&jet_AntiKt4TopoEM_emscale_pt );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_n",&jet_AntiKt4TopoEM_n );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_eta",&jet_AntiKt4TopoEM_eta );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_phi",&jet_AntiKt4TopoEM_phi );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_E",&jet_AntiKt4TopoEM_E );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SV0",&jet_AntiKt4TopoEM_flavor_weight_SV0 );
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_isBadLoose",&jet_AntiKt4TopoEM_isBadLoose );

   fChain->SetBranchAddress("trk_theta",&trk_theta);
   fChain->SetBranchAddress("trk_qoverp",&trk_qoverp);
   fChain->SetBranchAddress("trk_pt",&trk_pt);
   fChain->SetBranchAddress("trk_eta",&trk_eta);
   fChain->SetBranchAddress("trk_d0_wrtPV",&trk_d0_wrtPV);
   fChain->SetBranchAddress("trk_z0_wrtPV",&trk_z0_wrtPV);
   fChain->SetBranchAddress("trk_phi_wrtPV",&trk_phi_wrtPV);
   //fChain->SetBranchAddress("trk_cov_d0_wrtPV",&trk_);
   //fChain->SetBranchAddress("trk_cov_z0_wrtPV",&trk_);
   fChain->SetBranchAddress("trk_chi2",&trk_chi2);
   fChain->SetBranchAddress("trk_ndof",&trk_ndof);
   fChain->SetBranchAddress("trk_nBLHits",&trk_nBLHits);
   fChain->SetBranchAddress("trk_nPixHits",&trk_nPixHits);
   fChain->SetBranchAddress("trk_nSCTHits",&trk_nSCTHits);
  
 

   if(isdebug)cout<< __LINE__ <<endl;

 if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  //  nentries = 200000;

  cout << " fChain entries " << nentries <<endl;

  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    //Long64_t ientry = LoadTree(jentry);
    //if (ientry < 0) break;
    //     nb = fChain->GetEntry(jentry);   nbytes += nb;

    //--------- Get branches ----------------------------
    if(isdebug)cout<< __LINE__ <<endl;

    int ientry = fChain->LoadTree(jentry);

    
    
    fChain->GetBranch("RunNumber")->GetEntry(ientry);
    if(isMC)fChain->GetBranch("mc_channel_number")->GetEntry(ientry);
    fChain->GetBranch("lbn")->GetEntry(ientry);
    fChain->GetBranch("trig_L1_TAV")->GetEntry(ientry);
    fChain->GetBranch("trig_L2_passedPhysics")->GetEntry(ientry);
    fChain->GetBranch("trig_EF_passedPhysics")->GetEntry(ientry);
    fChain->GetBranch("trig_DB_SMK")->GetEntry(ientry);
    fChain->GetBranch("trig_DB_L1PSK")->GetEntry(ientry);
    fChain->GetBranch("trig_DB_HLTPSK")->GetEntry(ientry);
    fChain->GetBranch("vxp_n")->GetEntry(ientry);
    fChain->GetBranch("vxp_nTracks")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_pt")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_n")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_eta")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_phi")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_E")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_flavor_weight_SV0")->GetEntry(ientry);
    fChain->GetBranch("jet_AntiKt4TopoEM_emscale_pt")->GetEntry(ientry);
    fChain->GetBranch("trk_n")->GetEntry(ientry);
    fChain->GetBranch("trk_theta")->GetEntry(ientry);
    fChain->GetBranch("trk_qoverp")->GetEntry(ientry);
    fChain->GetBranch("trk_pt")->GetEntry(ientry);
    fChain->GetBranch("trk_eta")->GetEntry(ientry);
    fChain->GetBranch("trk_d0_wrtPV")->GetEntry(ientry);
    fChain->GetBranch("trk_z0_wrtPV")->GetEntry(ientry);
    fChain->GetBranch("trk_phi_wrtPV")->GetEntry(ientry);
    //fChain->GetBranch("trk_cov_d0_wrtPV")->GetEntry(ientry);
    //fChain->GetBranch("trk_cov_z0_wrtPV")->GetEntry(ientry);
    fChain->GetBranch("trk_chi2")->GetEntry(ientry);
    fChain->GetBranch("trk_ndof")->GetEntry(ientry);
    fChain->GetBranch("trk_nBLHits")->GetEntry(ientry);
    fChain->GetBranch("trk_nPixHits")->GetEntry(ientry);
    fChain->GetBranch("trk_nSCTHits")->GetEntry(ientry);
    
    if(isdebug)cout<< __LINE__ <<endl;     

    if(jentry%10000==0) cout << " Event # " <<jentry<<endl;
    


    //------------------ Set Event Weight in JX (isMC) -----------------------
    
   
    int NowRunningOn = -1;
    float w = 1;
    
    if(isMC){

      NowRunningOn =mc_channel_number;
      bool isPythia = true;
      bool isPerugia = false;
      bool isHerwig = false;

      if(isPythia){
	
	if(NowRunningOn ==105009){
	  Xsection = 1000*1.2032e+07 ;//9.8608e+06;
	  N0 = 2798297;//91617;//3637;// 1363438;
	}else if(NowRunningOn==105010){
	  Xsection = 1000*8.0714e+05;//.7818e+05;
	  N0 = 2798443;//7377565;//110951;//17111;//1395889;
	}else if(NowRunningOn==105011){
	  Xsection = 1000*4.8027e+04 ;//4.0982e+04;
	  N0 = 2795891;//114381;//52253;//1396991;
	}else if(NowRunningOn==105012){
	  Xsection =1000*2.5364e+03 ;//2.1929e+03;
	  N0 = 2797982;//213163;//99409;//1397590;
	}else if(NowRunningOn==105013){
	  Xsection = 1000*9.9605e+01;//8.7701e+01;
	  N0 = 2797431;//154659;//148610;//1393487;
	}else if(NowRunningOn==105014){
	  Xsection = 1000*2.5947e+00;
	  N0 = 2796405;//386425;//194467;//1392492l;
	}else if(NowRunningOn==105015){
	  Xsection = 1000*3.5457e-02 ;//3.361e-02;
	  N0 = 2791826;//436390;//232837;//1391670;
	}
	
      }else if(isPerugia){
	
	// For Perugia!!
	if(NowRunningOn ==115849){
	  Xsection = 1000*7.7714E+06;
	  N0 = 399849;//16388258;
	}else if(NowRunningOn==115850){
	  Xsection =  1000*5.0385E+05;
	  N0 =1498393;// 7392565;
	}else if(NowRunningOn==115851){
	  Xsection =  1000*2.9358E+04;
	  N0 = 988144;//2796084;
	}else if(NowRunningOn==115852){
	  Xsection = 1000*1.5600E+03;
	  N0 = 394497;//2796879;
	}else if(NowRunningOn==115853){
	  Xsection =  1000*6.4393E+01;
	  N0 =399199 ;//2793179;
	}else if(NowRunningOn==115854){
	  Xsection =  1000*1.8764E+00;
	  N0 = 399046;//2790576;
	}else if(NowRunningOn==115855){
	  Xsection =  1000*3.0412E-02;
	  N0 =398900;// 2790601;
	}
	
      }else if(isHerwig){
	
	
	if(NowRunningOn ==113204){
	  Xsection = 1000*9.6139E+06 ;
	  N0 = 399799 ;
	}else if(NowRunningOn==113205){
	  Xsection = 1000*7.4366E+05 ;
	  N0 = 398897 ;
	}else if(NowRunningOn==113206){
	  Xsection = 1000*4.4307E+04 ;
	  N0 = 398498 ;
	}else if(NowRunningOn==113207){
	  Xsection =1000*2.3576E+03 ;
	  N0 = 399598 ;
	}else if(NowRunningOn==113208){
	  Xsection = 1000*9.4236E+01 ;
	  N0 = 399443 ;
	}else if(NowRunningOn==113209){
	  Xsection = 1000*2.5813E+00;
	  N0 = 399094 ;
	}else if(NowRunningOn==113210){
	  Xsection = 1000*3.9439E-02 ;
	  N0 = 398597 ;
	}
      }
      
      w = Xsection*Lumi/N0;
      
    }

   
   

    //-------------- Vertex cut --------------------------------------
    
    int NGoodVtx = 0; 
    int numvxt = vxp_nTracks->size();
    for(int i=0; i<numvxt; i++) {
      int vtx_ntrks = (*vxp_nTracks)[i];
      if(vtx_ntrks>4)NGoodVtx++;
    }
    

    if(NGoodVtx<1)continue;



    //------------- MC10a Re-weighting ---------------------------
   TLorentzVector lj;
    lj.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
    
    float LeadPt =lj.Pt(); 

  if(isMC){

   //  TLorentzVector lj;
//     lj.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
    
//     float LeadPt =lj.Pt(); 
    if(NowRunningOn ==105009 && LeadPt >40) continue;
    else if(NowRunningOn==105010 && LeadPt >60)continue;
    else if(NowRunningOn==105011 && LeadPt >110) continue;
    else if(NowRunningOn==105012 && LeadPt >200) continue;
    else if(NowRunningOn==105013 && LeadPt >360) continue;
    else if(NowRunningOn==105014 && LeadPt >620) continue;
    else if(NowRunningOn==105015 && LeadPt >1200) continue;
    }
    

    //----------- Trigger Selection---------------------------------
    bool passes = true;
    if(isDATA){
      
      TLorentzVector leadJet;
      leadJet.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
      float LeadPt = leadJet.Pt();
      
      passes = PassedTrigger11 (RunNumber, 
				trig_DB_SMK, 
				trig_DB_L1PSK, 
				trig_DB_HLTPSK, 
				trig_L1_TAV,
				trig_EF_passedPhysics,
				LeadPt);
    }
    if (!passes) continue;

    

    //--------------- Loop through jets ---------------------------
    int bjets[30] = {0}; int nj = 0; 
    
    int  jet_AntiKt4TopoEM_num = jet_AntiKt4TopoEM_pt->size();
    for(int i=0; i<jet_AntiKt4TopoEM_num; i++) {
     
      
      //Cleaning cuts
      int goodness = 2;
      if(isDATA){
	goodness= (*jet_AntiKt4TopoEM_isBadLoose)[i];
	if(goodness != 0) continue;
      }
      
      
      TLorentzVector j;
      j.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[i]/1000.,(*jet_AntiKt4TopoEM_eta)[i],(*jet_AntiKt4TopoEM_phi)[i],(*jet_AntiKt4TopoEM_E)[i]/1000.);
      
      float PT = j.Pt();
      float ETA = TMath::Abs(j.Eta()); 
      

      //Phase spaace cut
      bool isInGoodCaloRegion = RemoveDeadRegion(j.Eta(),j.Phi());
      if(PT<jetPtMin|| ETA>jetEtaMax || isInGoodCaloRegion==false ){
	j.Delete();
	continue;
      }
      float PHI = j.Phi();
     
      
      //------------- btag weight cut -----------------------------
     
     float Btag_w = (*jet_AntiKt4TopoEM_flavor_weight_SV0)[i]; 
     bool isBjet = false;
     if(Btag_w>BtagCut) isBjet = true; 
     else continue;

     
      //-------------- match to tracks -------------------------------
      
      int trk_num = trk_pt->size();
      TLorentzVector trackJet;
      TLorentzVector t1(0,0,0,0);
      TLorentzVector t2(0,0,0,0);
      float sumpt = 0.0; int ntrk = 0;
      float maxpt = 0.0;	
      vector<int> trklist;
      //to compute jet track eccentricity
      float sumEiEtai2 = 0;
      float sumEiPhii2 = 0;
      float sumEiEtaiPhii = 0;
      float jetTrackEccentricity = 99.;      
      for(int trk=0; trk<trk_num; trk++) {
	
	if( ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk])>=cutNBLayerPlusPix && (*trk_nSCTHits)[trk]>= cutNSCTHits&&(*trk_nPixHits)[trk]>=cutNPixelHits && ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk]+(*trk_nSCTHits)[trk])>=cutNSiHits  &&(*trk_pt)[trk]/1000> cutOnTrackPT &&((*trk_chi2)[trk]/(*trk_ndof)[trk])<cutOnChi2OverNDof && TMath::Abs((*trk_d0_wrtPV)[trk])<cutd0PV && TMath::Abs((*trk_z0_wrtPV)[trk]*TMath::Sin((*trk_theta)[trk]))<cutz0PVsinTheta ){ 
	  
	  TLorentzVector track;
	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trk]);
	  track.SetPtEtaPhiE((*trk_pt)[trk]/1000., (*trk_eta)[trk],(*trk_phi_wrtPV)[trk],track_e/1000);
	  
	  float dr = j.DeltaR(track);
	  if(dr<0.4) {
	    if(ntrk==0) t1 = track;
	    if(ntrk==1) t2 = track;
	    sumpt += track.Pt();
	    sumEiEtai2 += track.Energy()*TMath::Power(track.Eta(),2);
	    sumEiPhii2 += track.Energy()*TMath::Power(track.Phi(),2);
	    sumEiEtaiPhii += track.Energy()*track.Phi()*track.Eta();
	    //----------------
	    //Fill track histos
	    d0PV_histo->Fill((*trk_d0_wrtPV)[trk]);
	    z0PV_histo->Fill((*trk_z0_wrtPV)[trk]);
	    z0PVsinTheta_histo->Fill((*trk_z0_wrtPV)[trk]*TMath::Sin((*trk_theta)[trk]));
	    ntrk ++;
	    trackJet += track;
	    // to compute drmax/P1
	    if(track.Pt()>maxpt) maxpt = track.Pt();
	    trklist.push_back(trk);
	  }
	}
      }
      
      float drmax = -1.0; 
      if(ntrk<1){
	continue;
	cout<<" ntrk<1 - >continue"<<endl;
      }

 
     
    
      //Counting
      if(isBjet){
	bjets[nj++]=i;
	//	Nbjets++; //counting bjets in the event for further analysis
	//NumBjets++; //counting events for printing it out afterwards
      }
     
     
      //------------ Isolation --------------------
      float closeby_unCalibPtCut = 7.;
      int NclosebyJets = 0;
      int NclosebyBJets = 0;

      TLorentzVector j2(0,0,0,0);
      float tmp_dr = 100.;
      for(int l=0; l<jet_AntiKt4TopoEM_num; l++) {
	TLorentzVector newJet;        

	float newJet_unCalibPt = (*jet_AntiKt4TopoEM_emscale_pt)[l]/1000.;
	
	newJet.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[l]/1000.,(*jet_AntiKt4TopoEM_eta)[l],(*jet_AntiKt4TopoEM_phi)[l],(*jet_AntiKt4TopoEM_E)[l]/1000.);    
	
	
	if(newJet==j)continue;
	
	if(newJet_unCalibPt>closeby_unCalibPtCut && j.DeltaR(newJet)<0.8)NclosebyJets++;
        if(j.DeltaR(newJet)<0.8&&(*jet_AntiKt4TopoEM_flavor_weight_SV0)[l]>5.85 )NclosebyBJets++;
	
	
	float dr = j.DeltaR(newJet);
	if(dr<tmp_dr){
	  tmp_dr = dr;
	  j2=newJet;
	}
      }


      //------------------- Keep isolated jets ?? --------------

      if(NclosebyBJets>0|| NclosebyJets>0){
	//	continue;
      }     
      

      //-------------------------------------------------------------
      // GET VARIABLES
      //-------------------------------------------------------------


      //Jet TrackEccentricity
      if(ntrk>1){     
	float lambdaPlus = 0.5*(sumEiEtai2+sumEiPhii2+TMath::Sqrt(TMath::Power(sumEiEtai2-sumEiPhii2,2)+4*TMath::Power(sumEiEtaiPhii,2)));
	float lambdaMinus = 0.5*(sumEiEtai2+sumEiPhii2-TMath::Sqrt(TMath::Power(sumEiEtai2-sumEiPhii2,2)+4*TMath::Power(sumEiEtaiPhii,2)));
	float Ratio = lambdaMinus/lambdaPlus; 
	jetTrackEccentricity = TMath::Sqrt(1-TMath::Power(Ratio,2));
      }else if(ntrk==1){
	//	cout << " ntrk = 1 => can't compute eccentricity!!! *** JetEccentricity = 99.***"<<endl;
      }



      // get number of clusters in jet
      int ncl = 1;// jet_AntiKt4TopoEM_constituents_pt->at(i).size();
      

      //Jet calo-eccentricity
      float caloSumEiEtai2 = 0.;
      float caloSumEiPhii2 = 0.;
      float caloSumEiEtaiPhii = 0.;
      float jetCaloEccentricity = 99.;
      //cout<< "about to enter loop.."<<endl;
      //  for(int calo=0; calo<ncl; calo++){
      // 	//cout<< "in loop.."<<endl;
      // 	float Energy = (*jet_AntiKt4TopoEM_constituents_e)[i][calo];
      // 	//cout << " Energy " <<Energy << endl;
      // 	float Eta = (*jet_AntiKt4TopoEM_constituents_eta)[i][calo];
      // 	float Phi = (*jet_AntiKt4TopoEM_constituents_phi)[i][calo];
      // 	//cout << " Energy " <<Energy << " Eta " << Eta << " Phi " <<Phi <<endl;
      // 	caloSumEiEtai2 += Energy*TMath::Power(Eta,2);
      // 	caloSumEiPhii2 += Energy*TMath::Power(Phi,2);
      // 	caloSumEiEtaiPhii += Energy*Phi*Eta;
      //       }
      
      //       if(ncl>1){     
      // 	float lambdaPlus = 0.5*(caloSumEiEtai2+caloSumEiPhii2+TMath::Sqrt(TMath::Power(caloSumEiEtai2-caloSumEiPhii2,2)+4*TMath::Power(caloSumEiEtaiPhii,2)));
      // 	float lambdaMinus = 0.5*(caloSumEiEtai2+caloSumEiPhii2-TMath::Sqrt(TMath::Power(caloSumEiEtai2-caloSumEiPhii2,2)+4*TMath::Power(caloSumEiEtaiPhii,2)));
      // 	float Ratio = lambdaMinus/lambdaPlus; 
      // 	jetCaloEccentricity = TMath::Sqrt(1-TMath::Power(Ratio,2));
      //       }else if(ncl==1){
      // 	cout << " ncl = 1 => can't compute eccentricity!!! *** JetEccentricity = 99.***"<<endl;
      //       }



      //------------- TEST FASTJETS ------------------------------

      //std::cout << "FASTJET: Here we go..." << std::endl;

      fastjet::RecombinationScheme recombScheme = fastjet::pt2_scheme;
      fastjet::JetAlgorithm jetalgorithm = fastjet::kt_algorithm;
      double rParameter = 0.2;
      fastjet::Strategy ktstrategy = fastjet::Best;
      
      //std::cout << "Is this working?" << std::endl;
      fastjet::JetDefinition jetDef(jetalgorithm,rParameter, recombScheme, ktstrategy);
      vector<fastjet::PseudoJet> particles;
      
      //Loop to tracklist
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {
	  TLorentzVector t;
	  float track_e = TMath::Abs( (1./(*trk_qoverp)[trklist[k]])/1000.);
	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],track_e);
	  particles.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
	}
      }
      
      //run the clustering, extract the jets
      fastjet::ClusterSequence cs(particles,jetDef);
      //      vector<PseudoJet> jets = cs.fastjet::inclusive_jets();
      vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
      unsigned int nsubjets = jets.size();
      


      //--------- to compute N-subjettiness variables------------------
      

      double Nsubjet1 = 99.;
      double Nsubjet2 = 99.;
      double Nsubjet3 = 99.;
      double Nsubjet12ratio = 99.;
      double DeltaRk2axes = 99.;
      vector<fastjet::PseudoJet> constituents;  // will use charged constituents
      //Loop to tracklist
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {
	  float trk_energy = TMath::Abs( 1./(*trk_qoverp)[trklist[k]]);
	  TLorentzVector t;
	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],trk_energy/1000.);
	  constituents.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
	}
      }
      Analysis(constituents, Nsubjet1, Nsubjet2, Nsubjet3,DeltaRk2axes);
      if(Nsubjet1!=99.&&Nsubjet2!=99.&&Nsubjet1!=0)Nsubjet12ratio = Nsubjet2/Nsubjet1;
      constituents.clear();




      //to compute drmax
      
      if(trklist.size()>1) {
	for(int k1=0; k1<trklist.size()-1; k1++) {
	  TLorentzVector t1;
	  float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
	  t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
	  for(int k2=k1+1; k2<trklist.size(); k2++) {
	    TLorentzVector t2;
	    float t2_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k2]]);
	    t2.SetPtEtaPhiE((*trk_pt)[trklist[k2]]/1000.,(*trk_eta)[trklist[k2]],(*trk_phi_wrtPV)[trklist[k2]],t2_e/1000.);
	    float DR = t1.DeltaR(t2);
	    if(DR>drmax) drmax = DR;
	  }
	}
      }



      //compute DR between two hardest tracks
      
      TLorentzVector track1(0,0,0,0);
      TLorentzVector track2(0,0,0,0);
      int track1_index = 0;
      int track2_index = 0;
      float track1_pt = 0.;
      float track2_pt = 0.;
      float DeltaR_trk12 = 0.; 
      
      if(trklist.size()>1) {
	for(int k=0; k<trklist.size(); k++) {    
	  if ((*trk_pt)[trklist[k]]/1000. > track1_pt) {
	    if (k>0) {
	      track2_pt  = track1_pt;
	      track2_index = track1_index;
	    }
	    track1_pt  = (*trk_pt)[trklist[k]]/1000.;
	    track1_index = trklist[k];
	  }
	  else if ( (*trk_pt)[trklist[k]]/1000.> track2_pt) {
	    track2_pt  = (*trk_pt)[trklist[k]]/1000.;
	    track2_index = trklist[k];
	  }
	}
      }
      
      
      track1.SetPtEtaPhiE((*trk_pt)[track1_index]/1000.,(*trk_eta)[track1_index],(*trk_phi_wrtPV)[track1_index],(TMath::Abs( 1./(*trk_qoverp)[track1_index]))/1000.);
      track2.SetPtEtaPhiE((*trk_pt)[track2_index]/1000.,(*trk_eta)[track2_index],(*trk_phi_wrtPV)[track2_index],(TMath::Abs( 1./(*trk_qoverp)[track2_index]))/1000.);
      
      if(trklist.size()>1) DeltaR_trk12 = track1.DeltaR(track2); 
     

     


      //Get jet width 

      float WIDTH = 99;//(*jet_AntiKt4TopoEM_WIDTH)[i];
      float trkWIDTH = 0;
      float trkWIDTH_num = 0;
      float trkWIDTH_den = 0;
      if(trklist.size()>1) {
	for(int k1=0; k1<trklist.size()-1; k1++) {
	  TLorentzVector t1;
	  float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
	  t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
	  float DR = t1.DeltaR(j);

	  trkWIDTH_num+=t1.Pt()*DR;
	  trkWIDTH_den+=t1.Pt(); 
	  
	}
      }
      trkWIDTH = trkWIDTH_num/trkWIDTH_den;
      if(trkWIDTH>0.4)cout<<"********************************************* trkWIDTH = " <<endl;


      //Get jet YFlip12
      float YFlip12 = 99;//(*jet_AntiKt4TopoEM_YFlip12)[i];
      //float YFlip23 = (*jet_AntiKt4TopoEM_YFlip23)[i];

      

      //get number of vertices
      int num_vtx = 1;// jet_AntiKt4TopoEM_JetFitterComb_weight->size();
      

      //get Jet mass
      float MASS = j.M();
      // cout<<" jet mass = " << MASS <<endl;
             

	

      //---------- Fill histos -------------------------------
      
      for(int i_pt=0;i_pt<N_PTbins;i_pt++){
	
	if(isMC){
	if(PT>pt1[i_pt] && PT<=pt2[i_pt]){
	  // if(LeadPt>pt1[i_pt] && LeadPt<=pt2[i_pt] && PT>pt1[i_pt] && PT<=pt2[i_pt]){

	  Tau1_histos[i_pt]->Fill(Nsubjet1,w);
	  Tau2_histos[i_pt]->Fill(Nsubjet2,w);
	  kt2axisDeltaR_histos[i_pt]->Fill(DeltaRk2axes,w);	
	  Tau12ratio_histos[i_pt]->Fill(Nsubjet12ratio,w);

	  pt_histos[i_pt]->Fill(PT,w);     

	  Ntrk_histos[i_pt]->Fill(ntrk,w);     
	  NClus_histos[i_pt]->Fill(ncl,w);     
	  Width_histos[i_pt]->Fill(WIDTH,w);
	  trkWidth_histos[i_pt]->Fill(trkWIDTH,w);
	  YFlip12_histos[i_pt]->Fill(YFlip12,w);
	  DRmax_histos[i_pt]->Fill(drmax,w);
	  NumVtx_histos[i_pt]->Fill(num_vtx,w);
	  DRtrk12_histos[i_pt]->Fill(DeltaR_trk12,w);
	  NumSubTrackJets_histos[i_pt]->Fill(nsubjets,w);
	  JetMass_histos[i_pt]->Fill(MASS);
	  JetTrackEccentricity_histos[i_pt]->Fill(jetTrackEccentricity);
	  JetCaloEccentricity_histos[i_pt]->Fill(jetCaloEccentricity);

	  //Correlations
	  Tau2Tau1Corr_histos[i_pt]->Fill(Nsubjet1,Nsubjet2,w);
	  Tau1DRkT2axesCorr_histos[i_pt]->Fill(DeltaRk2axes,Nsubjet1,w);
	  Tau2DRkT2axesCorr_histos[i_pt]->Fill(DeltaRk2axes,Nsubjet2,w);

	  WidthEta_Corr[i_pt]->Fill(ETA,WIDTH,w);
	  WidthNtrk_Corr[i_pt]->Fill(ntrk,WIDTH,w);
	  trkWidthEta_Corr[i_pt]->Fill(ETA,trkWIDTH,w);
	  trkWidthNtrk_Corr[i_pt]->Fill(ntrk,trkWIDTH,w);
	  NclusEta_Corr[i_pt]->Fill(ETA,ncl,w);
	  NclusNtrk_Corr[i_pt]->Fill(ntrk,ncl,w);
	  NclusWidth_Corr[i_pt]->Fill(WIDTH,ncl,w);
	  NclusYFlip12_Corr[i_pt]->Fill(YFlip12,ncl,w);
	  NclusDRmax_Corr[i_pt]->Fill(drmax,ncl,w);
	  NclusDRtrk12_Corr[i_pt]->Fill(DeltaR_trk12,ncl,w);


	  NumSubTrackJetsEta_Corr[i_pt]->Fill(ETA,nsubjets,w);
	  NumSubTrackJetsNtrk_Corr[i_pt]->Fill(ntrk,nsubjets,w);
	  NumSubTrackJetsNclus_Corr[i_pt]->Fill(ncl,nsubjets,w);
	  NumSubTrackJetsWidth_Corr[i_pt]->Fill(WIDTH,nsubjets,w);
	  //NumSubTrackJetsYscale_Corr[i_pt]->Fill(Yscale,nsubjets,w);
	  //NumSubTrackJetsYFlip12_Corr[i_pt]->Fill(YFlip12,nsubjets,w);
	  NumSubTrackJetsDRmax_Corr[i_pt]->Fill(drmax,nsubjets,w);
	  //NumSubTrackJetsDRtrk12_Corr[i_pt]->Fill(DeltaR_trk12,nsubjets,w);

	  JetMassEta_Corr[i_pt]->Fill(ETA,MASS);
	  JetMassNtrk_Corr[i_pt]->Fill(ntrk,MASS);
	  JetMassNclus_Corr[i_pt]->Fill(ncl,MASS);
	  JetMassWidth_Corr[i_pt]->Fill(WIDTH,MASS);
	  JetMassDRmax_Corr[i_pt]->Fill(drmax,MASS);
	  JetMassNumSubTrackJets_Corr[i_pt]->Fill(nsubjets,MASS);

	  YFlip12Eta_Corr[i_pt]->Fill(ETA,YFlip12,w);
	  YFlip12Ntrk_Corr[i_pt]->Fill(ntrk,YFlip12,w);
	  YFlip12Width_Corr[i_pt]->Fill(WIDTH,YFlip12,w);
	  
	  DRmaxEta_Corr[i_pt]->Fill(ETA,drmax,w);
	  DRmaxNtrk_Corr[i_pt]->Fill(ntrk,drmax,w);
	  DRmaxWidth_Corr[i_pt]->Fill(WIDTH,drmax,w);

	  DRtrk12Eta_Corr[i_pt]->Fill(ETA,DeltaR_trk12,w);
	  DRtrk12Ntrk_Corr[i_pt]->Fill(ntrk,DeltaR_trk12,w);
	  DRtrk12Width_Corr[i_pt]->Fill(WIDTH,DeltaR_trk12,w);

	  NumVtxEta_Corr[i_pt]->Fill(ETA, num_vtx,w);
	  NumVtxNtrk_Corr[i_pt]->Fill(ntrk, num_vtx,w);
	  NumVtxWidth_Corr[i_pt]->Fill(WIDTH, num_vtx,w);

  
	}

	}else if(isDATA){


	  if(LeadPt>pt1[i_pt] && LeadPt<=pt2[i_pt] && PT>pt1[i_pt] && PT<=pt2[i_pt]){

	  Tau1_histos[i_pt]->Fill(Nsubjet1,w);
	  Tau2_histos[i_pt]->Fill(Nsubjet2,w);
	  kt2axisDeltaR_histos[i_pt]->Fill(DeltaRk2axes,w);	
	  Tau12ratio_histos[i_pt]->Fill(Nsubjet12ratio,w);

	  pt_histos[i_pt]->Fill(PT,w);     

	  Ntrk_histos[i_pt]->Fill(ntrk,w);     
	  NClus_histos[i_pt]->Fill(ncl,w);     
	  Width_histos[i_pt]->Fill(WIDTH,w);
	  trkWidth_histos[i_pt]->Fill(trkWIDTH,w);
	  YFlip12_histos[i_pt]->Fill(YFlip12,w);
	  DRmax_histos[i_pt]->Fill(drmax,w);
	  NumVtx_histos[i_pt]->Fill(num_vtx,w);
	  DRtrk12_histos[i_pt]->Fill(DeltaR_trk12,w);
	  NumSubTrackJets_histos[i_pt]->Fill(nsubjets,w);
	  JetMass_histos[i_pt]->Fill(MASS);
	  JetTrackEccentricity_histos[i_pt]->Fill(jetTrackEccentricity);
	  JetCaloEccentricity_histos[i_pt]->Fill(jetCaloEccentricity);

	  //Correlations
	  Tau2Tau1Corr_histos[i_pt]->Fill(Nsubjet1,Nsubjet2,w);
	  Tau1DRkT2axesCorr_histos[i_pt]->Fill(DeltaRk2axes,Nsubjet1,w);
	  Tau2DRkT2axesCorr_histos[i_pt]->Fill(DeltaRk2axes,Nsubjet2,w);

	  WidthEta_Corr[i_pt]->Fill(ETA,WIDTH,w);
	  WidthNtrk_Corr[i_pt]->Fill(ntrk,WIDTH,w);
	  trkWidthEta_Corr[i_pt]->Fill(ETA,trkWIDTH,w);
	  trkWidthNtrk_Corr[i_pt]->Fill(ntrk,trkWIDTH,w);
	  NclusEta_Corr[i_pt]->Fill(ETA,ncl,w);
	  NclusNtrk_Corr[i_pt]->Fill(ntrk,ncl,w);
	  NclusWidth_Corr[i_pt]->Fill(WIDTH,ncl,w);
	  NclusYFlip12_Corr[i_pt]->Fill(YFlip12,ncl,w);
	  NclusDRmax_Corr[i_pt]->Fill(drmax,ncl,w);
	  NclusDRtrk12_Corr[i_pt]->Fill(DeltaR_trk12,ncl,w);


	  NumSubTrackJetsEta_Corr[i_pt]->Fill(ETA,nsubjets,w);
	  NumSubTrackJetsNtrk_Corr[i_pt]->Fill(ntrk,nsubjets,w);
	  NumSubTrackJetsNclus_Corr[i_pt]->Fill(ncl,nsubjets,w);
	  NumSubTrackJetsWidth_Corr[i_pt]->Fill(WIDTH,nsubjets,w);
	  //NumSubTrackJetsYscale_Corr[i_pt]->Fill(Yscale,nsubjets,w);
	  //NumSubTrackJetsYFlip12_Corr[i_pt]->Fill(YFlip12,nsubjets,w);
	  NumSubTrackJetsDRmax_Corr[i_pt]->Fill(drmax,nsubjets,w);
	  //NumSubTrackJetsDRtrk12_Corr[i_pt]->Fill(DeltaR_trk12,nsubjets,w);

	  JetMassEta_Corr[i_pt]->Fill(ETA,MASS);
	  JetMassNtrk_Corr[i_pt]->Fill(ntrk,MASS);
	  JetMassNclus_Corr[i_pt]->Fill(ncl,MASS);
	  JetMassWidth_Corr[i_pt]->Fill(WIDTH,MASS);
	  JetMassDRmax_Corr[i_pt]->Fill(drmax,MASS);
	  JetMassNumSubTrackJets_Corr[i_pt]->Fill(nsubjets,MASS);

	  YFlip12Eta_Corr[i_pt]->Fill(ETA,YFlip12,w);
	  YFlip12Ntrk_Corr[i_pt]->Fill(ntrk,YFlip12,w);
	  YFlip12Width_Corr[i_pt]->Fill(WIDTH,YFlip12,w);
	  
	  DRmaxEta_Corr[i_pt]->Fill(ETA,drmax,w);
	  DRmaxNtrk_Corr[i_pt]->Fill(ntrk,drmax,w);
	  DRmaxWidth_Corr[i_pt]->Fill(WIDTH,drmax,w);

	  DRtrk12Eta_Corr[i_pt]->Fill(ETA,DeltaR_trk12,w);
	  DRtrk12Ntrk_Corr[i_pt]->Fill(ntrk,DeltaR_trk12,w);
	  DRtrk12Width_Corr[i_pt]->Fill(WIDTH,DeltaR_trk12,w);

	  NumVtxEta_Corr[i_pt]->Fill(ETA, num_vtx,w);
	  NumVtxNtrk_Corr[i_pt]->Fill(ntrk, num_vtx,w);
	  NumVtxWidth_Corr[i_pt]->Fill(WIDTH, num_vtx,w);

  
	}

	}



      }
      

      //Fill the other histos
      if(PT>40.)Pt_histos->Fill(PT,w);
      Eta_histos->Fill(ETA,w);
      
      
    }//jets

   
    
   }//entries

  
 
  
  
  //------------- Save histos-------------------------------------------------------
  char name[80]; 
  if(isDATA)sprintf(name,"AllHistos_Data2011_QCDjetjet_NO_ISO_r17_Jan132011_PeriodH.root");
  if(isMC)sprintf(name,"AllHistos_MC11b_QCDjetjet_NO_ISO_r17_Jan132011.root");
 
  TFile g(name,"recreate");

  Pt_L1J10->Write();

  Pt_histos->Write();

  Eta_histos->Write();
  NumBjets_histo->Write();  
  NumIsolatedBjets_histo->Write();//Fill(Evt_NIsolatedBjets,w);
  NumNonIsolatedBjets_histo->Write();
  NumIsolatedSingleBjets_histo->Write();//Fill(Evt_NIsolatedSingleBjets,w);
  NumIsolatedMergedBjets_histo->Write();
  
  d0PV_histo->Write();//((*trk_d0_wrtPV)[trk]);
  z0PV_histo->Write();//((*trk_z0_wrtPV)[trk]);
  z0PVsinTheta_histo->Write();

  for(int i_pt=0;i_pt<N_PTbins;i_pt++){
    
    Tau1_histos[i_pt]->Write();//(Nsubjet1,w);
    Tau2_histos[i_pt]->Write();//(Nsubjet2,w);
    kt2axisDeltaR_histos[i_pt]->Write();//(DeltaRk2axes,w);	
    Tau12ratio_histos[i_pt]->Write();
     
    pt_histos[i_pt]->Write();     
    Ntrk_histos[i_pt]->Write();     
    Width_histos[i_pt]->Write();
    trkWidth_histos[i_pt]->Write();
    NClus_histos[i_pt]->Write();     

    YFlip12_histos[i_pt]->Write();
    DRmax_histos[i_pt]->Write();//(drmax);
    NumVtx_histos[i_pt]->Write();//(num_vtx);
    DRtrk12_histos[i_pt]->Write();//(DeltaR_trk12); 
    NumSubTrackJets_histos[i_pt]->Write();
    JetMass_histos[i_pt]->Write();
    JetTrackEccentricity_histos[i_pt]->Write();//Fill(jetEccentricity);
    JetCaloEccentricity_histos[i_pt]->Write();


    Tau2Tau1Corr_histos[i_pt]->Write();//(Nsubjet1,Nsubjet2,w);
    Tau1DRkT2axesCorr_histos[i_pt]->Write();//(DeltaRk2axes,Nsubjet1,w);
    Tau2DRkT2axesCorr_histos[i_pt]->Write();//(DeltaRk2axes,Nsubjet2,w);

    WidthEta_Corr[i_pt]->Write();
    WidthNtrk_Corr[i_pt]->Write();
    trkWidthEta_Corr[i_pt]->Write();
    trkWidthNtrk_Corr[i_pt]->Write();

    NclusEta_Corr[i_pt]->Write();
    NclusNtrk_Corr[i_pt]->Write();
    NclusWidth_Corr[i_pt]->Write();
    NclusYFlip12_Corr[i_pt]->Write();//(YFlip12,ncl);
    NclusDRmax_Corr[i_pt]->Write();//(drmax,ncl);
    NclusDRtrk12_Corr[i_pt]->Write();//(DeltaR_trk12,ncl);

    NumSubTrackJetsEta_Corr[i_pt]->Write();
    NumSubTrackJetsNtrk_Corr[i_pt]->Write();
    NumSubTrackJetsNclus_Corr[i_pt]->Write();
    NumSubTrackJetsWidth_Corr[i_pt]->Write();
    // NumSubTrackJetsYFlip12_Corr[i_pt]->Write();//(YFlip12,ncl);
    NumSubTrackJetsDRmax_Corr[i_pt]->Write();//(drmax,ncl);
    //NumSubTrackJetsDRtrk12_Corr[i_pt]->Write();//(DeltaR_trk12,ncl);

    JetMassEta_Corr[i_pt]->Write();
    JetMassNtrk_Corr[i_pt]->Write();
    JetMassNclus_Corr[i_pt]->Write();
    JetMassWidth_Corr[i_pt]->Write();
    JetMassDRmax_Corr[i_pt]->Write();//(drmax,ncl);
    JetMassNumSubTrackJets_Corr[i_pt]->Write();//(drmax,ncl);

    YFlip12Eta_Corr[i_pt]->Write();
    YFlip12Ntrk_Corr[i_pt]->Write();
    YFlip12Width_Corr[i_pt]->Write();

    DRmaxNtrk_Corr[i_pt]->Write();//(ntrk,drmax);
    DRmaxWidth_Corr[i_pt]->Write();//(WIDTH,drmax);
    
    DRtrk12Eta_Corr[i_pt]->Write();//(ETA,DeltaR_trk12);
    DRtrk12Ntrk_Corr[i_pt]->Write();//(ntrk,DeltaR_trk12);
    DRtrk12Width_Corr[i_pt]->Write();//(WIDTH,DeltaR_trk12);
    NumVtxEta_Corr[i_pt]->Write();//(ETA, num_vtx);
    NumVtxNtrk_Corr[i_pt]->Write();//(ntrk, num_vtx);
    NumVtxWidth_Corr[i_pt]->Write();//(WIDTH, num_vtx);
    
  }
  
  //  DR_e->Write();
  //   DR_mu->Write();
  //   DR_particle->Write();
  //   DR_B1->Write();
  //   DR_B2->Write();
  //   NumBs->Write();
  
  DR_BB_all->Write();
  DR_BB_allinEta->Write();
  DR_BB_allMatching->Write();
  DR_BB_single->Write();
  DR_BB_merged->Write();
  DPhi_BB_all->Write();
  DPhi_BB_allinEta->Write();
  DPhi_BB_allMatching->Write();
  DPhi_BB_single->Write();
  DPhi_BB_merged->Write();
  DR_bjets->Write();
  DPhi_bjets->Write();

  g.Close();
  
  
  //====================================
  //delete histos to prevent memory leak  
  Pt_L1J10->Delete();

  Pt_histos->Delete();
  Eta_histos->Delete();
  NumBjets_histo->Delete();  
  NumIsolatedBjets_histo->Delete();//Fill(Evt_NIsolatedBjets,w);
  NumNonIsolatedBjets_histo->Delete();
  NumIsolatedSingleBjets_histo->Delete();//Fill(Evt_NIsolatedSingleBjets,w);
  NumIsolatedMergedBjets_histo->Delete();

  d0PV_histo->Delete();//((*trk_d0_wrtPV)[trk]);      
  z0PV_histo->Delete();//((*trk_z0_wrtPV)[trk]); 
  z0PVsinTheta_histo->Delete();


  for(int i_pt=0;i_pt<N_PTbins;i_pt++){
    
    Tau1_histos[i_pt]->Delete();//(Nsubjet1,w);
    Tau2_histos[i_pt]->Delete();//(Nsubjet2,w);
    kt2axisDeltaR_histos[i_pt]->Delete();
    Tau12ratio_histos[i_pt]->Delete();

    pt_histos[i_pt]->Delete();     
    Ntrk_histos[i_pt]->Delete();     
    Width_histos[i_pt]->Delete();
    trkWidth_histos[i_pt]->Delete();
    NClus_histos[i_pt]->Delete();     

    YFlip12_histos[i_pt]->Delete();
    DRmax_histos[i_pt]->Delete();//(drmax);
    NumVtx_histos[i_pt]->Delete();//(num_vtx);
    DRtrk12_histos[i_pt]->Delete();//(DeltaR_trk12); 
    NumSubTrackJets_histos[i_pt]->Delete();
    JetMass_histos[i_pt]->Delete();
    JetTrackEccentricity_histos[i_pt]->Delete();//Fill(jetEccentricity);
    JetCaloEccentricity_histos[i_pt]->Delete();

    Tau2Tau1Corr_histos[i_pt]->Delete();//(Nsubjet1,Nsubjet2,w);
    Tau1DRkT2axesCorr_histos[i_pt]->Delete();//(DeltaRk2axes,Nsubjet1,w);
    Tau2DRkT2axesCorr_histos[i_pt]->Delete();//(DeltaRk2axes,Nsubjet2,w);


    WidthEta_Corr[i_pt]->Delete();
    WidthNtrk_Corr[i_pt]->Delete();
    trkWidthEta_Corr[i_pt]->Delete();
    trkWidthNtrk_Corr[i_pt]->Delete();

    NclusEta_Corr[i_pt]->Delete();
    NclusNtrk_Corr[i_pt]->Delete();
    NclusWidth_Corr[i_pt]->Delete();
    NclusYFlip12_Corr[i_pt]->Delete();//(YFlip12,ncl);
    NclusDRmax_Corr[i_pt]->Delete();//(drmax,ncl);
    NclusDRtrk12_Corr[i_pt]->Delete();//(DeltaR_trk12,ncl);

    NumSubTrackJetsEta_Corr[i_pt]->Delete();
    NumSubTrackJetsNtrk_Corr[i_pt]->Delete();
    NumSubTrackJetsNclus_Corr[i_pt]->Delete();
    NumSubTrackJetsWidth_Corr[i_pt]->Delete();
    // NumSubTrackJetsYFlip12_Corr[i_pt]->Delete();//(YFlip12,ncl);
    NumSubTrackJetsDRmax_Corr[i_pt]->Delete();//(drmax,ncl);
    //NumSubTrackJetsDRtrk12_Corr[i_pt]->Delete();//(DeltaR_trk12,ncl);

    JetMassEta_Corr[i_pt]->Delete();
    JetMassNtrk_Corr[i_pt]->Delete();
    JetMassNclus_Corr[i_pt]->Delete();
    JetMassWidth_Corr[i_pt]->Delete();
    JetMassDRmax_Corr[i_pt]->Delete();//(drmax,ncl);
    JetMassNumSubTrackJets_Corr[i_pt]->Delete();//(drmax,ncl);

    YFlip12Eta_Corr[i_pt]->Delete();
    YFlip12Ntrk_Corr[i_pt]->Delete();
    YFlip12Width_Corr[i_pt]->Delete();

    DRmaxNtrk_Corr[i_pt]->Delete();//(ntrk,drmax);
    DRmaxWidth_Corr[i_pt]->Delete();//(WIDTH,drmax);
    
    DRtrk12Eta_Corr[i_pt]->Delete();//(ETA,DeltaR_trk12);
    DRtrk12Ntrk_Corr[i_pt]->Delete();//(ntrk,DeltaR_trk12);
    DRtrk12Width_Corr[i_pt]->Delete();//(WIDTH,DeltaR_trk12);
    NumVtxEta_Corr[i_pt]->Delete();//(ETA, num_vtx);
    NumVtxNtrk_Corr[i_pt]->Delete();//(ntrk, num_vtx);
    NumVtxWidth_Corr[i_pt]->Delete();//(WIDTH, num_vtx);
    
  }
  
  //  DR_e->Delete();
  //   DR_mu->Delete();
  //   DR_particle->Delete();
  //   DR_B1->Delete();
  //   DR_B2->Delete();
  //   NumBs->Delete();
  
  DR_BB_all->Delete();
  DR_BB_allinEta->Delete();
  DR_BB_allMatching->Delete();
  DR_BB_single->Delete();
  DR_BB_merged->Delete();
  DPhi_BB_all->Delete();
  DPhi_BB_allinEta->Delete();
  DPhi_BB_allMatching->Delete();
  DPhi_BB_single->Delete();
  DPhi_BB_merged->Delete();
  DR_bjets->Delete();
  DPhi_bjets->Delete();
  
 
  
}//MakeHistos



			
//=========================================================================================
// Make histograms of MVA output
//=========================================================================================						       
// void myAnalysis::EvaluateNN(){
 

//   cout << " In EvaluateNN() " <<endl;

//    //--------- Event global variables -------------------------
//    const float Lumi = 100.;
//    float Xsection = 99.; 
//    float N0 = 1400000.;//99999.; 

  

//    //--------- Selection cuts (pt/met in GeV!!) ----------------
//    float jetEtaMax=2.1;
//    float jetPtMin=20;//*GeV;
//    //float metMin=20.0;//*GeV;
//    float lepPt = 15.0; //GeV
//    //float crackEtaMin = 1.37;
//    //float crackEtaMax = 1.52;
//    float BtagCut = 5.85;//5.72;
//    float DR_cut = 0.4;
  
   
//    //--------- B hadrons pdgids ----------------------------

//    int NumBpdgids = 6;// 24;
//    int B_pdgids[NumBpdgids];
//    for(int i=0;i<NumBpdgids;i++)B_pdgids[i]=0;
  
//   B_pdgids[0]=511;
//   B_pdgids[1]=521;
//   B_pdgids[2]=531;
//   B_pdgids[3]=541;
//   B_pdgids[4]=5122;
//   B_pdgids[5]=5132;


//   //------- Quality cuts on tracks -----------------------
//   double cutOnTrackPT = 1.;//0.5;//
//   double cutOnChi2OverNDof = 3.0; //10000
//   int cutNBLayerHits=0; //
//   int cutNSCTHits=4; //6;
//   int cutNPixelHits=1;//2;
//   int cutNBLayerPlusPix = 0;
//   int cutNSiHits = 7;  
//   double cutd0PV = 2.0;//1.5;
//   double cutz0PVsinTheta = 2.0;
//   //double IPSignificanceCut = 1.;// <- This is a "less than". To test displaced trackas sigma>2.5; //2.; // 3.;    
//   //double IPzSignificanceCut = 5.; //2.; // 3.;    
  

//   //----------  Create a new root output file. ----------------
  
//   TFile* file = TFile::Open("Likelihood_MC11b_QCDjetjet_NOISO_drktaxis_p832_SV0.root","RECREATE");

  

//   //----------- Book histos ------------------------------------

//   TH1F* h_output = new TH1F("h_output","",100,-0.01,1.01);
//   h_output->SetFillColor(2);
//   h_output->SetLineWidth(2);
//   h_output->SetFillStyle(3001);
//   h_output->GetXaxis()->SetRangeUser(0.09,0.92);


//   TH1F* h_output_lowpt = new TH1F("h_output_lowpt","",100,-0.01,1.01);
//   TH1F* h_output_midpt = new TH1F("h_output_midpt","",100,-0.01,1.01);
//   TH1F* h_output_highpt = new TH1F("h_output_highpt","",100,-0.01,1.01);

//   TH1F *ptsig = new TH1F("ptsig", "", 50, 0,200);


//   const int N_PTbins = 8;
//   float pt1[N_PTbins]={40.,60.,80.,110.,150.,200.,270.,360}; 
//   float pt2[N_PTbins]={60.,80.,110.,150.,200.,270.,360,480}; 
//   int Pt1[N_PTbins]={40.,60.,80.,110.,150.,200.,270.,360};


//   TString separator = "_";

//   // NN output en bines de pt
//   TH1F* output_histo[N_PTbins]; 
//   for(int j=0;j<N_PTbins;j++){  
//     output_histo[j]=new TH1F("pt"+separator+"PT"+Form("%i",Pt1[j]),"pt"+separator+"PT"+Form("%i",Pt1[j]),50,-0.01,1.01);
//     output_histo[j]->Sumw2();
//   }


//   //------------- TMVA stuff -----------------------------------

//   vector<TString> inputVars;
//   // inputVars.push_back("ncl");
//   // inputVars.push_back("drmax");
//   //  inputVars.push_back("tauratio");
//   //inputVars.push_back("tau2");
//   inputVars.push_back("drktaxis");
//   inputVars.push_back("pt");
//   inputVars.push_back("ntrk");
//   inputVars.push_back("trkwidth");//nkt");
  

//   TMVA::Reader *reader = new TMVA::Reader(inputVars);  
//   cout << " YOU NEED TO HAVE weights.txt in your run directory to run this " <<endl;
 
//   // reader->BookMVA("kCFMlpANN", "./NNMC10b_weights/MVAnalysis_kCFMlpANN.weights.txt");
//   //reader->BookMVA("kMLP", "weights_MLP_newroot_fixOptionsMORE/MVAnalysis_MLP.weights.xml");
//   reader->BookMVA("kPDERS", "./PDERSMC10b_weights/MVAnalysis_PDERS.weights.txt");

  



//   /**************************************************************/
//   //Loop through entries
//   /**************************************************************/

//   cout << " About to load branches " <<endl;
//   bool isdebug =false;

 
//   if(isdebug) cout << __LINE__<<endl;

//   fChain->SetBranchStatus("*",0);
//  if(isdebug) cout << __LINE__<<endl;
//    fChain->SetBranchStatus("RunNumber",1);
//  if(isdebug) cout << __LINE__<<endl;
//    fChain->SetBranchStatus("mc_channel_number",1);
//    // fChain->SetBranchStatus("lbn",1);
//    fChain->SetBranchStatus("vxp_n",1);
//    fChain->SetBranchStatus("vxp_nTracks",1);
//  if(isdebug) cout << __LINE__<<endl;
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_pt",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_emscale_pt",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_n",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_eta",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_phi",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_E",1);
//    fChain->SetBranchStatus("jet_AntiKt4TopoEM_flavor_weight_SV0",1);
//    //  fChain->SetBranchStatus("jet_AntiKt4TopoEM_WIDTH",1);
//    //fChain->SetBranchStatus("jet_AntiKt4TopoEM_isBadLoose",1);
//  if(isdebug) cout << __LINE__<<endl;
//    fChain->SetBranchStatus("trk_n",1);
   
//    fChain->SetBranchStatus("trk_theta",1);
//    fChain->SetBranchStatus("trk_qoverp",1);
//    fChain->SetBranchStatus("trk_pt",1);
//    fChain->SetBranchStatus("trk_eta",1);
//    fChain->SetBranchStatus("trk_d0_wrtPV",1);
//    fChain->SetBranchStatus("trk_z0_wrtPV",1);
//    fChain->SetBranchStatus("trk_phi_wrtPV",1);
//    //fChain->SetBranchStatus("trk_cov_d0_wrtPV",1);
//    //fChain->SetBranchStatus("trk_cov_z0_wrtPV",1);
//    fChain->SetBranchStatus("trk_chi2",1);
//    fChain->SetBranchStatus("trk_ndof",1);
//    fChain->SetBranchStatus("trk_nBLHits",1);
//    fChain->SetBranchStatus("trk_nPixHits",1);
//    fChain->SetBranchStatus("trk_nSCTHits",1);
  
//    // link local variable to branch
//   if(isdebug) cout << __LINE__<<endl;

//    int  RunNumber = 0;
//    int mc_channel_number =0;
//    //int  lbn = 0;
//    int      vxp_n = 0;
//    vector<int>     *vxp_nTracks;
//    int              jet_AntiKt4TopoEM_n = 0;
//    vector<float>   *jet_AntiKt4TopoEM_E;
//    vector<float>   *jet_AntiKt4TopoEM_pt;
//    vector<float>   *jet_AntiKt4TopoEM_emscale_pt;
//     vector<float>   *jet_AntiKt4TopoEM_eta;
//    vector<float>   *jet_AntiKt4TopoEM_phi;
//    vector<float>  *jet_AntiKt4TopoEM_flavor_weight_SV0;
//    // vector<float>  *jet_AntiKt4TopoEM_WIDTH;
//    //vector<int>     *jet_AntiKt4TopoEM_isBadLoose;
//    int           trk_n =0;
//    vector<float>   *trk_theta;
//    vector<float>   *trk_qoverp;
//    vector<float>   *trk_pt;
//    vector<float>   *trk_eta;
//    vector<float>   *trk_d0_wrtPV;
//    vector<float>   *trk_z0_wrtPV;
//    vector<float>   *trk_phi_wrtPV;
//    // vector<float>   *trk_cov_d0_wrtPV;
//    // vector<float>   *trk_cov_z0_wrtPV;
//    vector<float>   *trk_chi2;
//    vector<int>     *trk_ndof;
//    vector<int>     *trk_nBLHits;
//    vector<int>     *trk_nPixHits;
//    vector<int>     *trk_nSCTHits;

  
//    vxp_nTracks = 0;
//    jet_AntiKt4TopoEM_E = 0;
//    jet_AntiKt4TopoEM_pt = 0;
//    jet_AntiKt4TopoEM_emscale_pt = 0;
//    jet_AntiKt4TopoEM_eta = 0;
//    jet_AntiKt4TopoEM_phi = 0;
//    jet_AntiKt4TopoEM_flavor_weight_SV0 = 0;
//    //  jet_AntiKt4TopoEM_WIDTH = 0;
//    //   jet_AntiKt4TopoEM_isBadLoose=0;
//    trk_theta = 0;
//    trk_qoverp = 0;
//    trk_pt = 0;
//    trk_eta = 0;
//    trk_d0_wrtPV = 0;
//    trk_z0_wrtPV = 0;
//    trk_phi_wrtPV = 0;
//    //   trk_cov_d0_wrtPV = 0;
//    //trk_cov_z0_wrtPV = 0;
//    trk_chi2 = 0;
//    trk_ndof = 0;
//    trk_nBLHits = 0;
//    trk_nPixHits = 0;
//    trk_nSCTHits = 0;

//   if(isdebug) cout << __LINE__<<endl;


//    fChain->SetBranchAddress("RunNumber", &RunNumber);
//   fChain->SetBranchAddress("mc_channel_number", &mc_channel_number);

//    //fChain->SetBranchAddress("lbn", &lbn);
//    fChain->SetBranchAddress("vxp_n",&vxp_n );
//    fChain->SetBranchAddress("vxp_nTracks",&vxp_nTracks );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_pt",&jet_AntiKt4TopoEM_pt );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_pt",&jet_AntiKt4TopoEM_emscale_pt );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_n",&jet_AntiKt4TopoEM_n );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_eta",&jet_AntiKt4TopoEM_eta );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_phi",&jet_AntiKt4TopoEM_phi );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_E",&jet_AntiKt4TopoEM_E );
//    fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SV0",&jet_AntiKt4TopoEM_flavor_weight_SV0 );
//    //  fChain->SetBranchAddress("jet_AntiKt4TopoEM_WIDTH",&jet_AntiKt4TopoEM_WIDTH );
//    // fChain->SetBranchAddress("jet_AntiKt4TopoEM_isBadLoose",&jet_AntiKt4TopoEM_isBadLoose );

//    fChain->SetBranchAddress("trk_theta",&trk_theta);
//    fChain->SetBranchAddress("trk_qoverp",&trk_qoverp);
//    fChain->SetBranchAddress("trk_pt",&trk_pt);
//    fChain->SetBranchAddress("trk_eta",&trk_eta);
//    fChain->SetBranchAddress("trk_d0_wrtPV",&trk_d0_wrtPV);
//    fChain->SetBranchAddress("trk_z0_wrtPV",&trk_z0_wrtPV);
//    fChain->SetBranchAddress("trk_phi_wrtPV",&trk_phi_wrtPV);
//    //fChain->SetBranchAddress("trk_cov_d0_wrtPV",&trk_);
//    //fChain->SetBranchAddress("trk_cov_z0_wrtPV",&trk_);
//    fChain->SetBranchAddress("trk_chi2",&trk_chi2);
//    fChain->SetBranchAddress("trk_ndof",&trk_ndof);
//    fChain->SetBranchAddress("trk_nBLHits",&trk_nBLHits);
//    fChain->SetBranchAddress("trk_nPixHits",&trk_nPixHits);
//    fChain->SetBranchAddress("trk_nSCTHits",&trk_nSCTHits);
  
//   if(isdebug) cout << __LINE__<<endl; 



//   if (fChain == 0) return;
//   Long64_t nentries = fChain->GetEntries();
//   //  nentries = 200; //0;
//   Long64_t nbytes = 0, nb = 0;
//   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   

//    //  Long64_t ientry = LoadTree(jentry);
// //     if (ientry < 0) break;
// //     nb = fChain->GetEntry(jentry);   nbytes += nb;
   

//     //--------- Get branches ----------------------------
//     int ientry = fChain->LoadTree(jentry);
//     if(isdebug) cout << __LINE__<<endl;
    
//     fChain->GetBranch("RunNumber")->GetEntry(ientry);
//     fChain->GetBranch("mc_channel_number")->GetEntry(ientry);
//     fChain->GetBranch("vxp_n")->GetEntry(ientry);
//     fChain->GetBranch("vxp_nTracks")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_pt")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_n")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_eta")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_phi")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_E")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_flavor_weight_SV0")->GetEntry(ientry);
//     //    fChain->GetBranch("jet_AntiKt4TopoEM_WIDTH")->GetEntry(ientry);
//     fChain->GetBranch("jet_AntiKt4TopoEM_emscale_pt")->GetEntry(ientry);
//     fChain->GetBranch("trk_n")->GetEntry(ientry);
//     fChain->GetBranch("trk_theta")->GetEntry(ientry);
//     fChain->GetBranch("trk_qoverp")->GetEntry(ientry);
//     fChain->GetBranch("trk_pt")->GetEntry(ientry);
//     fChain->GetBranch("trk_eta")->GetEntry(ientry);
//     fChain->GetBranch("trk_d0_wrtPV")->GetEntry(ientry);
//     fChain->GetBranch("trk_z0_wrtPV")->GetEntry(ientry);
//     fChain->GetBranch("trk_phi_wrtPV")->GetEntry(ientry);
//     //fChain->GetBranch("trk_cov_d0_wrtPV")->GetEntry(ientry);
//     //fChain->GetBranch("trk_cov_z0_wrtPV")->GetEntry(ientry);
//     fChain->GetBranch("trk_chi2")->GetEntry(ientry);
//     fChain->GetBranch("trk_ndof")->GetEntry(ientry);
//     fChain->GetBranch("trk_nBLHits")->GetEntry(ientry);
//     fChain->GetBranch("trk_nPixHits")->GetEntry(ientry);
//     fChain->GetBranch("trk_nSCTHits")->GetEntry(ientry);
    
//     if(isdebug) cout << __LINE__<<endl;

    
//     //    cout << " jentry " <<jentry<<endl;
//     if(jentry%10000==0){
//         cout << " jentry " <<jentry<<endl;
//       if(mc_channel_number == 105009) cout << "Event #" << jentry << " running on J0"<<endl;
//       if(mc_channel_number == 105010) cout << "Event #" << jentry << " running on J1"<<endl;
//       if(mc_channel_number == 105011) cout << "Event #" << jentry << " running on J2"<<endl;
//       if(mc_channel_number == 105012) cout << "Event #" << jentry << " running on J3"<<endl;
//       if(mc_channel_number == 105013) cout << "Event #" << jentry << " running on J4"<<endl;
//       if(mc_channel_number == 105014) cout << "Event #" << jentry << " running on J5"<<endl;
//       if(mc_channel_number == 105015) cout << "Event #" << jentry << " running on J6"<<endl;
//     }

    
  
//   if(isdebug) cout << __LINE__<<endl;

//     //------------- Set Event Weight in JX -------------------

//     int NowRunningOn = mc_channel_number;
//     float N0 = 1400000.;//99999.;
//     bool isPythia = true;
//     bool isPerugia = false;
//     bool isHerwig = false;
    
//     if(isPythia){
      

//     if(NowRunningOn ==105009){
//       Xsection = 1000*1.2032e+07 ;//9.8608e+06;
//       N0 = 2798297;//91617;//3637;// 1363438;
//     }else if(NowRunningOn==105010){
//       Xsection = 1000*8.0714e+05;//.7818e+05;
//       N0 = 2798443;//7377565;//110951;//17111;//1395889;
//     }else if(NowRunningOn==105011){
//       Xsection = 1000*4.8027e+04 ;//4.0982e+04;
//       N0 = 2795891;//114381;//52253;//1396991;
//     }else if(NowRunningOn==105012){
//       Xsection =1000*2.5364e+03 ;//2.1929e+03;
//       N0 = 2797982;//213163;//99409;//1397590;
//     }else if(NowRunningOn==105013){
//       Xsection = 1000*9.9605e+01;//8.7701e+01;
//       N0 = 2797431;//154659;//148610;//1393487;
//     }else if(NowRunningOn==105014){
//       Xsection = 1000*2.5947e+00;
//       N0 = 2796405;//386425;//194467;//1392492l;
//     }else if(NowRunningOn==105015){
//       Xsection = 1000*3.5457e-02 ;//3.361e-02;
//       N0 = 2791826;//436390;//232837;//1391670;
//     }
 
//     }else if(isPerugia){
      
//       // For Perugia!!
//       if(NowRunningOn ==115849){
// 	Xsection = 1000*7.7714E+06;
// 	N0 = 399849;//16388258;
//       }else if(NowRunningOn==115850){
// 	Xsection =  1000*5.0385E+05;
// 	N0 =1498393;// 7392565;
//       }else if(NowRunningOn==115851){
// 	Xsection =  1000*2.9358E+04;
// 	N0 = 988144;//2796084;
//       }else if(NowRunningOn==115852){
// 	Xsection = 1000*1.5600E+03;
// 	N0 = 394497;//2796879;
//       }else if(NowRunningOn==115853){
// 	Xsection =  1000*6.4393E+01;
// 	N0 =399199 ;//2793179;
//       }else if(NowRunningOn==115854){
// 	Xsection =  1000*1.8764E+00;
// 	N0 = 399046;//2790576;
//       }else if(NowRunningOn==115855){
// 	Xsection =  1000*3.0412E-02;
// 	N0 =398900;// 2790601;
//       }
      
//     }else if(isHerwig){
      
      
//       if(NowRunningOn ==113204){
// 	Xsection = 1000*9.6139E+06 ;
// 	N0 = 399799 ;
//       }else if(NowRunningOn==113205){
// 	Xsection = 1000*7.4366E+05 ;
// 	N0 = 398897 ;
//       }else if(NowRunningOn==113206){
// 	Xsection = 1000*4.4307E+04 ;
// 	N0 = 398498 ;
//       }else if(NowRunningOn==113207){
// 	Xsection =1000*2.3576E+03 ;
// 	N0 = 399598 ;
//       }else if(NowRunningOn==113208){
// 	Xsection = 1000*9.4236E+01 ;
// 	N0 = 399443 ;
//       }else if(NowRunningOn==113209){
// 	Xsection = 1000*2.5813E+00;
// 	N0 = 399094 ;
//       }else if(NowRunningOn==113210){
// 	Xsection = 1000*3.9439E-02 ;
// 	N0 = 398597 ;
//       }
//     }
    
    
        
//     float w = Xsection*Lumi/N0;

   



//     //----------------- Event Selection-------------------------

   
//     //------------------ Vertex cut -----------------------------

//     unsigned int NEvtVtx =  vxp_n;//vxp_x->size(); //Vertex_nTracks->size();
//        int NGoodVtx = 0;//vxp_n; 
//     for(int i=0; i<NEvtVtx; i++) {
//       int vtx_ntrks = (*vxp_nTracks)[i];
//       if(vtx_ntrks>4)NGoodVtx++;
//     }
    
//     //NPV==1
//     if(NGoodVtx<1) continue;
    
//     //      cout << __LINE__ <<endl;
    

//     //------------------ MC10a Re-weighting -------------------
//     TLorentzVector lj;
//     lj.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[0]/1000.,(*jet_AntiKt4TopoEM_eta)[0],(*jet_AntiKt4TopoEM_phi)[0],(*jet_AntiKt4TopoEM_E)[0]/1000.);
//     //float ljEMJES = (*jet_AntiKt4TopoEM_EMJES)[0];
//     //lj*=ljEMJES;
//     float Leadpt =lj.Pt(); 
//     if(NowRunningOn ==105009 && Leadpt >40) continue;
//     else if(NowRunningOn==105010&& Leadpt >60)continue;
//     else if(NowRunningOn==105011&& Leadpt >110)continue;
//     else if(NowRunningOn==105012&& Leadpt >200)continue;
//     else if(NowRunningOn==105013&& Leadpt >360)continue;
//     else if(NowRunningOn==105014&& Leadpt >620)continue;
//     else if(NowRunningOn==105015&& Leadpt >1200)continue;
//      // cout << __LINE__ <<endl;


 


//     //------------------ Loop through jets ------------------------------

//     //int  jet_AntiKt4TopoEM_num = jet_AntiKt4TopoEM_emscale_pt->size();
//      // cout << __LINE__ <<endl;
//     for(int i=0; i<jet_AntiKt4TopoEM_n; i++) {
     

//       //---------------------------------------------------
//       //Cleaning cuts
//       int goodness = 2;//(*jet_AntiKt4TopoEM_isGood)[i];
//       if(goodness != 2) continue;
      
//       TLorentzVector j;
//       j.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[i]/1000.,(*jet_AntiKt4TopoEM_eta)[i],(*jet_AntiKt4TopoEM_phi)[i],(*jet_AntiKt4TopoEM_E)[i]/1000.);


//       float PT = j.Pt();
//       float ETA = TMath::Abs(j.Eta());      


//       if(PT<jetPtMin || ETA>jetEtaMax){
// 	j.Delete();
// 	continue;
//       }


 
//       //---------- btag weight cut -----------------------------

//       float Btag_w = (*jet_AntiKt4TopoEM_flavor_weight_SV0)[i]; 
//       bool isBjet = false; 
//       if(Btag_w>BtagCut) isBjet = true; 
//       else continue;
      
//      // cout << __LINE__ <<endl;

//       //------------ match to tracks ----------------------------
//       int trk_num = trk_pt->size();
//       TLorentzVector trackJet;
//       TLorentzVector t1(0,0,0,0);
//       TLorentzVector t2(0,0,0,0);
//       float sumpt = 0.0; int ntrk = 0;
//       float maxpt = 0.0;	
//       vector<int> trklist;
//       //to compute jet eccentricity
//       float sumEiEtai2 = 0;
//       float sumEiPhii2 = 0;
//       float sumEiEtaiPhii = 0;
//       float jetTrackEccentricity = 99.; 
//       for(int trk=0; trk<trk_num; trk++) {
	
// 	if( ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk])>=cutNBLayerPlusPix && (*trk_nSCTHits)[trk]>= cutNSCTHits&&(*trk_nPixHits)[trk]>=cutNPixelHits && ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk]+(*trk_nSCTHits)[trk])>=cutNSiHits  &&(*trk_pt)[trk]/1000> cutOnTrackPT &&((*trk_chi2)[trk]/(*trk_ndof)[trk])<cutOnChi2OverNDof && TMath::Abs((*trk_d0_wrtPV)[trk])<cutd0PV &&TMath::Abs((*trk_z0_wrtPV)[trk]*TMath::Sin((*trk_theta)[trk]))<cutz0PVsinTheta){ 



// 	  TLorentzVector track;
// 	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trk]);
// 	   track.SetPtEtaPhiE((*trk_pt)[trk]/1000., (*trk_eta)[trk],(*trk_phi_wrtPV)[trk],track_e/1000);
	    
// 	  float dr = j.DeltaR(track);
// 	  if(dr<0.4) {
// 	    if(ntrk==0) t1 = track;
// 	    if(ntrk==1) t2 = track;
// 	    sumpt += track.Pt();
// 	    sumEiEtai2 += track.Energy()*TMath::Power(track.Eta(),2);
// 	    sumEiPhii2 += track.Energy()*TMath::Power(track.Phi(),2);
// 	    sumEiEtaiPhii += track.Energy()*track.Phi()*track.Eta();
	   
// 	    //-----------------
// 	    //Increase ntrk
// 	    ntrk ++;
// 	    trackJet += track;
// 	    // to compute drmax/P1
// 	    if(track.Pt()>maxpt) maxpt = track.Pt();
// 	    trklist.push_back(trk);
// 	  }
// 	}
//       }
      
	
//       // Tracks quantities
//       // float Pt1 = maxpt/sumpt; 
//       // float ftrk = sumpt/PT;
//       float drmax = -1.0; 
      
//       if(ntrk<1){
// 	continue;
// 	cout<<" ntrk<1 - >continue"<<endl;
//       }

//      // cout << __LINE__ <<endl;     

//       //---------------- Isolation ------------------------------
//       float closeby_unCalibPtCut = 4.;
//       int Ncloseby4GeVJets = 0;
//       int Ncloseby7GeVJets = 0;
//       int Ncloseby4GeVJetsNtrk = 0;
//       int Ncloseby4GeVBJets = 0;
//       int NclosebyBJets = 0;
//       TLorentzVector j2(0,0,0,0);
//       float tmp_dr = 100.;
//       for(int l=0; l<jet_AntiKt4TopoEM_n; l++) {

// 	TLorentzVector newJet;  

//       	float newJet_unCalibPt = (*jet_AntiKt4TopoEM_emscale_pt)[l]/1000.;
// 	newJet.SetPtEtaPhiE((*jet_AntiKt4TopoEM_pt)[l]/1000.,(*jet_AntiKt4TopoEM_eta)[l],(*jet_AntiKt4TopoEM_phi)[l],(*jet_AntiKt4TopoEM_E)[l]/1000.);


// 	if(newJet==j)continue;

// 	//track matching
// 	int newJet_ntrk = 0;
// 	for(int trk=0; trk<trk_num; trk++) {
	  
// 	  if( ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk])>=cutNBLayerPlusPix && (*trk_nSCTHits)[trk]>= cutNSCTHits&&(*trk_nPixHits)[trk]>=cutNPixelHits && ((*trk_nBLHits)[trk]+(*trk_nPixHits)[trk]+(*trk_nSCTHits)[trk])>=cutNSiHits  &&(*trk_pt)[trk]/1000> cutOnTrackPT &&((*trk_chi2)[trk]/(*trk_ndof)[trk])<cutOnChi2OverNDof && TMath::Abs((*trk_z0_wrtPV)[trk]*TMath::Sin((*trk_theta)[trk]))<cutz0PVsinTheta  ){ 

	    
// 	    TLorentzVector track;
// 	    float track_e = TMath::Abs( 1./(*trk_qoverp)[trk]);
// 	    track.SetPtEtaPhiE((*trk_pt)[trk]/1000., (*trk_eta)[trk],(*trk_phi_wrtPV)[trk],track_e/1000);
// 	    float dr = newJet.DeltaR(track);
// 	    if(dr<0.4) newJet_ntrk ++;
// 	  }
// 	}
	
// 	if(newJet_unCalibPt>closeby_unCalibPtCut && j.DeltaR(newJet)<0.8)Ncloseby4GeVJets++;
// 	if(newJet_unCalibPt>7. && j.DeltaR(newJet)<0.8)Ncloseby7GeVJets++;
//         if(newJet_unCalibPt>closeby_unCalibPtCut && j.DeltaR(newJet)<0.8 && newJet_ntrk>1 )Ncloseby4GeVJetsNtrk++;
//         if(newJet_unCalibPt>closeby_unCalibPtCut && j.DeltaR(newJet)<0.8&&(*jet_AntiKt4TopoEM_flavor_weight_SV0)[l]>5.85 )Ncloseby4GeVBJets++;
//         if(j.DeltaR(newJet)<0.8&&(*jet_AntiKt4TopoEM_flavor_weight_SV0)[l]>5.85 )NclosebyBJets++;


// 	float dr = j.DeltaR(newJet);
// 	if(dr<tmp_dr){
// 	  tmp_dr = dr;
// 	  j2=newJet;
// 	}
//       }
    
//       //------------ Keep isolated jets ??  --------------------------

//       // if(NclosebyBJets>0|| Ncloseby7GeVJets>0)continue;
   
//      // cout << __LINE__ <<endl;
      
//       //========================================================
//       //TEST FASTJETS
//       // std::cout << "Here we go..." << std::endl;
//       fastjet::RecombinationScheme recombScheme = fastjet::pt2_scheme;
//       fastjet::JetAlgorithm jetalgorithm = fastjet::kt_algorithm;
//       double rParameter = 0.2;
//       fastjet::Strategy ktstrategy = fastjet::Best;
      
//       //std::cout << "Is this working?" << std::endl;
//       fastjet::JetDefinition jetDef(jetalgorithm,rParameter, recombScheme, ktstrategy);
//       vector<fastjet::PseudoJet> particles;
      
//       //Loop to tracklist
//       if(trklist.size()>1) {
// 	for(int k=0; k<trklist.size(); k++) {
// 	  TLorentzVector t;
// 	  float track_e = TMath::Abs( (1./(*trk_qoverp)[trklist[k]])/1000.);
// 	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],track_e);
// 	  particles.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
// 	}
//       }
      
//       //run the clustering, extract the jets
//       fastjet::ClusterSequence cs(particles,jetDef);
//       //      vector<PseudoJet> jets = cs.fastjet::inclusive_jets();
//       vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
      
//       unsigned int nsubjets = jets.size();
//       //cout<< " jet " << i << " has " <<nsubjets << " sub-track jets" <<endl;
//       //std::cout << "All done!" << std::endl;
//       //=====================================================

//      // cout << __LINE__ <<endl;
//       //--------------------------------------------------
//       //to compute N-subjettiness variables
//       //--------------------------------------------------
//       double Nsubjet1 = 99.;
//       double Nsubjet2 = 99.;
//       double Nsubjet3 = 99.;
//       double Nsubjet12ratio = 99.;
//       double DeltaRk2axes = 99.;
//       vector<fastjet::PseudoJet> constituents;  // will use charged constituents
//       //Loop to tracklist
//       if(trklist.size()>1) {
// 	for(int k=0; k<trklist.size(); k++) {
// 	  TLorentzVector t;
// 	  float track_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k]]);
// 	  t.SetPtEtaPhiE((*trk_pt)[trklist[k]]/1000.,(*trk_eta)[trklist[k]],(*trk_phi_wrtPV)[trklist[k]],track_e/1000.);
// 	  constituents.push_back(fastjet::PseudoJet(t.Px(),t.Py(), t.Pz(),t.E()));
// 	}
//       }
//       Analysis(constituents, Nsubjet1, Nsubjet2, Nsubjet3,DeltaRk2axes);
//       if(Nsubjet1!=99.&&Nsubjet2!=99.&&Nsubjet1!=0)Nsubjet12ratio = Nsubjet2/Nsubjet1;
//       constituents.clear();


//  // cout << __LINE__ <<endl;

//       //=====================================================
//       //to compute drmax
//       if(trklist.size()>1) {
// 	for(int k1=0; k1<trklist.size()-1; k1++) {
// 	  TLorentzVector t1;
// 	  float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
// 	  t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
// 	  for(int k2=k1+1; k2<trklist.size(); k2++) {
// 	    TLorentzVector t2;
// 	    float t2_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k2]]);
// 	    t2.SetPtEtaPhiE((*trk_pt)[trklist[k2]]/1000.,(*trk_eta)[trklist[k2]],(*trk_phi_wrtPV)[trklist[k2]],t2_e/1000.);
// 	    float DR = t1.DeltaR(t2);
// 	    if(DR>drmax) drmax = DR;
// 	  }
// 	}
//       }

//       //=========================================================
//       //compute DR between two hardest tracks
//       TLorentzVector track1(0,0,0,0);
//       TLorentzVector track2(0,0,0,0);
//       int track1_index = 0;
//       int track2_index = 0;
//       float track1_pt = 0.;
//       float track2_pt = 0.;
//       float DeltaR_trk12 = 0.; 
      
//       if(trklist.size()>1) {
// 	for(int k=0; k<trklist.size(); k++) {    
// 	  if ((*trk_pt)[trklist[k]]/1000. > track1_pt) {
// 	    if (k>0) {
// 	      track2_pt  = track1_pt;
// 	      track2_index = track1_index;
// 	    }
// 	    track1_pt  = (*trk_pt)[trklist[k]]/1000.;
// 	    track1_index = trklist[k];
// 	  }
// 	  else if ( (*trk_pt)[trklist[k]]/1000.> track2_pt) {
// 	    track2_pt  = (*trk_pt)[trklist[k]]/1000.;
// 	    track2_index = trklist[k];
// 	  }
// 	}
//       }
      
      
//       track1.SetPtEtaPhiE((*trk_pt)[track1_index]/1000.,(*trk_eta)[track1_index],(*trk_phi_wrtPV)[track1_index],(TMath::Abs( 1./(*trk_qoverp)[track1_index]))/1000.);
//       track2.SetPtEtaPhiE((*trk_pt)[track2_index]/1000.,(*trk_eta)[track2_index],(*trk_phi_wrtPV)[track2_index],(TMath::Abs( 1./(*trk_qoverp)[track2_index]))/1000.);
      
//       if(trklist.size()>1) DeltaR_trk12 = track1.DeltaR(track2); 
//       //========================================================
//  // cout << __LINE__ <<endl;

//       //======================
//       //Get jet width 
//       float WIDTH = 1;//(*jet_AntiKt4TopoEM_WIDTH)[i];
//       float trkWIDTH = 0;
//       float trkWIDTH_num = 0;
//       float trkWIDTH_den = 0;
//       if(trklist.size()>1) {
// 	for(int k1=0; k1<trklist.size()-1; k1++) {
// 	  TLorentzVector t1;
// 	  float t1_e = TMath::Abs( 1./(*trk_qoverp)[trklist[k1]]);
// 	  t1.SetPtEtaPhiE((*trk_pt)[trklist[k1]]/1000.,(*trk_eta)[trklist[k1]],(*trk_phi_wrtPV)[trklist[k1]],t1_e/1000.);
// 	  float DR = t1.DeltaR(j);

// 	  trkWIDTH_num+=t1.Pt()*DR;
// 	  trkWIDTH_den+=t1.Pt(); 
	  
// 	}
//       }
//       trkWIDTH = trkWIDTH_num/trkWIDTH_den;
//       if(trkWIDTH!=trkWIDTH){
// 	//	cout<< " trkWIDTH NOT A NUMER NaN NaN NaN " << trkWIDTH <<endl;
//  	continue;
//       }


//      // cout << __LINE__ <<endl;

//       //======================
//       // get number of clusters in jet
//       int ncl = 1;//AntiKt4TopoNew_constituents_pt->at(i).size();

//       //=======================
//       //get Jet mass
//       float MASS = j.M();
//       //cout<<" jet mass = " << MASS <<endl;
      

//       //=======================
//       //Fill
      
      
//       vector<double> vars;
//       //vars.push_back(ncl);
//       //vars.push_back(drmax);
//       //vars.push_back(Nsubjet12ratio);
//       //vars.push_back(Nsubjet2);
//       vars.push_back(DeltaRk2axes);
//       vars.push_back(PT);
//       vars.push_back(ntrk);
//       vars.push_back(trkWIDTH);

//      // cout << __LINE__ <<endl;      

//       //float o = reader->EvaluateMVA(vars, "kCFMlpANN");//TMVA::Types::kCFMlpANN);
//       float o = reader->EvaluateMVA(vars, "kPDERS");
      
//       if(PT>40) h_output->Fill(o);//,w);
//       if(PT<40.){
// 	h_output_lowpt->Fill(o,w);
//       }
//       if(PT>40.&& PT<60.){
// 	h_output_midpt->Fill(o,w);
//       }
//       if(PT > 80.){ 
// 	h_output_highpt->Fill(o,w); 
//       }
      
//       ptsig->Fill(PT,w);       

//       //Fill new histos in bins of pt
//       for(int i_pt=0;i_pt<N_PTbins;i_pt++){
// 	if(PT>pt1[i_pt] && PT<=pt2[i_pt]){
// 	  output_histo[i_pt]->Fill(o,w);
// 	}
//       }
   
      
     
      
//       //  }
	
      
//     }//jets
    
//   }//entries
  
//   //file.cd();
//   file->cd();
 
//   h_output->Write();
 
//   h_output_lowpt->Write();
  
 
//   h_output_midpt->Write();
  
 
//   h_output_highpt->Write();

//   for(int i_pt=0;i_pt<N_PTbins;i_pt++){
//     output_histo[i_pt]->Write();
//   }

//   ptsig->Write();
//   //  ptbg->Write();

//   //file.Close();
//   file->Close();
  
// }//EvaluateNN







bool myAnalysis::PassedTrigger11 (UInt_t runNumber,					     
				  UInt_t trig_DB_SMK, 
				  UInt_t trig_DB_L1PSK, 
				  UInt_t trig_DB_HLTPSK, 
				  vector<unsigned int>* trig_L1_TAV,
				  vector<short>* trig_EF_passedPhysics,
				  double pt1) {
  
 //   cout << " In PassedTrigger11 : SMK " <<trig_DB_SMK<<endl;
//     cout << " In PassedTrigger11 : L1PSK " <<trig_DB_L1PSK<<endl;
//     cout << " In PassedTrigger11 : HLTPSK " <<trig_DB_HLTPSK<<endl;
  
  
  float avpt =  pt1;//(pt1 + pt2 )/2. ;
  
  const int ntrig = 11;
   double turnOnPtCuts[ntrig] = {20.,30.,40.,60.,80.,110.,150.,200.,270.,360.,480.}; // Always need to be checked using turn on curves
   // double turnOnPtCuts[ntrig] = {40.,60.,80.,110.,150.,200.,270.,360.,480.}; // Always need to be checked using turn on curves


  // For further info, take a look at:
  // https://twiki.cern.ch/twiki/bin/viewauth/Atlas/TrigJetMenu
  
  string trigWordEFincPeriodAB[ntrig] = {"EF_j10_a4_EFFS","EF_j15_a4_EFFS","EF_j20_a4_EFFS","EF_j30_a4_EFFS",
					 "EF_j40_a4_EFFS","EF_j55_a4_EFFS","EF_j75_a4_EFFS","EF_j100_a4_EFFS",
					 "EF_j135_a4_EFFS","EF_j180_a4_EFFS","EF_j240_a4_EFFS"};
  
  string trigWordEFincFromPeriodD[ntrig] = {"EF_j10_a4tc_EFFS","EF_j15_a4tc_EFFS","EF_j20_a4tc_EFFS","EF_j30_a4tc_EFFS",
					    "EF_j40_a4tc_EFFS","EF_j55_a4tc_EFFS","EF_j75_a4tc_EFFS","EF_j100_a4tc_EFFS",
					    "EF_j135_a4tc_EFFS","EF_j180_a4tc_EFFS","EF_j240_a4tc_EFFS"};
  
  
  string trigWordEFforwardPeriodAB[ntrig] = {"EF_j10_a4_EFFS","EF_j15_a4_EFFS","EF_j20_a4_EFFS","EF_fj30_a4_EFFS",
					     "EF_j40_a4_EFFS","EF_fj55_a4_EFFS","EF_fj75_a4_EFFS","EF_fj100_a4_EFFS",
					     "EF_j135_a4_EFFS,","EF_j180_a4_EFFS","EF_j240_a4_EFFS"};

  string trigWordEFforwardFromPeriodD[ntrig] = {"EF_j10_a4tc_EFFS","EF_j15_a4tc_EFFS","EF_j20_a4tc_EFFS","EF_fj30_a4tc_EFFS",
						"EF_j40_a4tc_EFFS","EF_fj55_a4tc_EFFS","EF_fj75_a4tc_EFFS","EF_fj100_a4tc_EFFS",
						"EF_j135_a4tc_EFFS,","EF_j180_a4tc_EFFS","EF_j240_a4tc_EFFS"};
  
  
  string trigWordEFinc[ntrig] = {""};
  string trigWordEFforward[ntrig] = {""};

  if (runNumber >= 177531 && runNumber <= 178109 ) {
    for (int i = 0;i<ntrig;i++) { trigWordEFinc[i] = trigWordEFincPeriodAB[i]; trigWordEFforward[i] = trigWordEFforwardPeriodAB[i]; }
  }
  else {
    for (int i = 0;i<ntrig;i++) { trigWordEFinc[i] = trigWordEFincFromPeriodD[i]; trigWordEFforward[i] = trigWordEFforwardFromPeriodD[i]; }
  }

  bool good = false; bool lockDebug = false;
  
  // Say hello: Trigger info
  if (lockDebug) {
    cout << " myAnalysis::PassedTrigger() called " << endl;
    cout << " run Number = " << runNumber << endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // Trigger strategy chosen for 2011 data (first term)
  //////////////////////////////////////////////////////////////////////////////
  
  // Trigger selection 
  if ( runNumber >= 177531 /* Periods A - B */ )
    {
      if ( avpt > turnOnPtCuts[0] && avpt <= turnOnPtCuts[1] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[0], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[0] << endl;
	}
      }
      else if ( avpt > turnOnPtCuts[1] && avpt <= turnOnPtCuts[2]  )  {
	if ( m_triggerReader.PassesEF(trigWordEFinc[1], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[0], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[1] << endl;
	}
      }
      else if ( avpt > turnOnPtCuts[2] && avpt <= turnOnPtCuts[3]  )  {
	if ( m_triggerReader.PassesEF(trigWordEFinc[2], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[1], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[2] << endl;
	} 
      }
      else if ( avpt > turnOnPtCuts[3] && avpt <= turnOnPtCuts[4]  ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[3], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[2], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[3], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[2], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[3] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[4] && avpt <= turnOnPtCuts[5] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[4], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[3], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[4], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[3], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[4] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[5] && avpt <= turnOnPtCuts[6] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[5], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[4], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[5], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[4], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[5] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[6] && avpt <= turnOnPtCuts[7]) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[6], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[5], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[6], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFforward[5], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  	  
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[6] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[7] && avpt <= turnOnPtCuts[8] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[7], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[6], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  	  	
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[7] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[8] && avpt <= turnOnPtCuts[9] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[8], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[7], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  	  	
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[7] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[9] && avpt <= turnOnPtCuts[10] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[9], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[8], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  	  	
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[7] << endl;	
	}
      }
      else if ( avpt > turnOnPtCuts[10] ) {
	if ( m_triggerReader.PassesEF(trigWordEFinc[10], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ||
	     m_triggerReader.PassesEF(trigWordEFinc[9], trig_DB_SMK, trig_DB_L1PSK, trig_DB_HLTPSK, trig_EF_passedPhysics) ) {
	  good = true;	  	  	  	  	  	
	  if (lockDebug && good) cout << "avpt = " << avpt << "\t" << " Trigger fired = " << trigWordEFinc[7] << endl;	
	}
      }
    }
  //else { if (lockDebug && good) cout << " Event has not passed any trigger item " << endl;  }
  else { if(good==false)cout << " Event has not passed any trigger item " << endl;  }
  
  
  
  return good;
}


  
bool myAnalysis::RemoveDeadRegion(float jeteta, float jetphi)
{

 //cout << "(eta,phi) = (" <<jeteta <<","<<jetphi<<")"<<endl;
 bool useMe = true;
 float step = TMath::TwoPi()/64.;   float bad[14] ={1,4,9,12,14,17,46,49,48,51,53,56,58,62};

 if (jeteta > 0.0 && jeteta < 0.8)
  {              
    if (jetphi > -TMath::Pi()+(bad[0]-1)*step && jetphi < -TMath::Pi()+(bad[1]-1)*step ) useMe = false;
    else if (jetphi > -TMath::Pi()+(bad[6]-1)*step && jetphi < -TMath::Pi()+(bad[7]-1)*step ) useMe = false;
    else if (jetphi > -TMath::Pi()+(bad[10]-1)*step && jetphi < -TMath::Pi()+(bad[11]-1)*step ) useMe = false;
    else {}
  }
 else if (jeteta > -0.8 && jeteta < 0.0)
  {
    if (jetphi > -TMath::Pi()+(bad[2]-1)*step && jetphi < -TMath::Pi()+(bad[3]-1)*step ) useMe = false;
    else if (jetphi > -TMath::Pi()+(bad[4]-1)*step && jetphi < -TMath::Pi()+(bad[5]-1)*step ) useMe = false;
    else if (jetphi > -TMath::Pi()+(bad[8]-1)*step && jetphi < -TMath::Pi()+(bad[9]-1)*step ) useMe = false;           else if (jetphi > -TMath::Pi()+(bad[12]-1)*step && jetphi < -TMath::Pi()+(bad[13]-1)*step ) useMe = false;           else {}
  }

 //if (!useMe) cout << "Bad region: " << "(eta,phi) = (" <<jeteta <<","<<jetphi<<")"<<endl;

 return useMe;

}
