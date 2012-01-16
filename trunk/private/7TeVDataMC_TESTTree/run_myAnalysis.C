{

  //==============================================================================
  // Run myAnalysis.cxx
  // For bjets - gluon splitting analysis
  //==============================================================================


  bool ispcatlas010 =false;//false;
  bool ispcuba001 = true;//true;

  bool isMakeHistos = false;
  bool isEvaluate = false;
  
  bool isDATA = true;
  bool isMC = false;

  
  //================================================================================
  //Fastjet stuff...
  //================================================================================
  if(ispcatlas010){

  gROOT->ProcessLine(".include /afs/cern.ch/sw/lcg/external/fastjet/2.4.2p3/i686-slc5-gcc43-opt/include/");
  gSystem->Load("/afs/cern.ch/sw/lcg/external/fastjet/2.4.2p3/i686-slc5-gcc43-opt/lib/libfastjet.so");


  }else if(ispcuba001){

    gROOT->ProcessLine(".include /home/laugs/fastjet-3.0.0/include/");
    gSystem->Load("/home/laugs/fastjet-3.0.0/lib/libfastjet.so");
   
  }
  //================================================================================
  //GoodRunListStuff


  if(isDATA){	
    if(ispcatlas010){
  
      if(isMakeHistos){
	std::string grlInclude = std::string("-I\"/afs/cern.ch/atlas/software/builds/AtlasEvent/16.6.3/DataQuality/GoodRunsLists/\" ");
	grlInclude = gSystem->GetIncludePath() + grlInclude;
	gSystem->SetIncludePath(grlInclude.c_str());
	gSystem->Load("/afs/cern.ch/atlas/software/builds/AtlasEvent/16.6.3/DataQuality/GoodRunsLists/i686-slc5-gcc43-opt/libGoodRunsListsLib.so");
	
      }else if(isEvaluate){
	std::string grlInclude = std::string("-I\"/afs/cern.ch/atlas/software/builds/AtlasEvent/15.6.12/DataQuality/GoodRunsLists/\" ");
	grlInclude = gSystem->GetIncludePath() + grlInclude;
	gSystem->SetIncludePath(grlInclude.c_str());
	gSystem->Load("/afs/cern.ch/atlas/software/builds/AtlasEvent/15.6.12/DataQuality/GoodRunsLists/i686-slc5-gcc43-opt/libGoodRunsListsLib.so");
      }
      
    }else if(ispcuba001){
      std::string grlInclude = std::string("-I\"/afs/cern.ch/user/l/laugs/private/DataUtilities/GoodRunsLists-00-00-76/\" ");
      grlInclude = gSystem->GetIncludePath() + grlInclude;
      gSystem->SetIncludePath(grlInclude.c_str());
      gSystem->Load("/afs/cern.ch/user/l/laugs/private/DataUtilities/Libs/goodRunsList.15.6.0-i686-slc4-gcc34/libGoodRunsLists.so");
    }
  }

  
  //================================================================================
  //Compile Macros 
  //================================================================================
  //gSystem->Load("libPhysics.so");
  gSystem->CompileMacro ("CollectionTree.C", "k");
  //gSystem->CompileMacro ("Helix.cxx", "k");
  gSystem->CompileMacro ("AnalysisBase.cxx", "k");
  gSystem->CompileMacro ("TriggerReader.cxx");
  gSystem->CompileMacro ("myAnalysis.cxx", "k");
 

  //================================================================================
  // Load TChains
  //===============================================================================

  if(isDATA==true){
    
  TChain* tree = new TChain("qcd");
  TChain* m_metaDataTree = new TChain("qcdMeta/TrigConfTree");


  //------ Load data in pcatlas010!! ----------------------------------

  //period B  
  
 tree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106141952/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178021.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142325/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178026.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142345/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178047.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142425/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178109.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142449/*.root*");
  
  //period D
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179710.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142512/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179739.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142552/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179771.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142615/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179804.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142634/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179938.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142651/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179939.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142719/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179940.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142742/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180122.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142803/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180139.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142856/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180144.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142916/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180149.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142937/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180153.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143011/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180225.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143120/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180241.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143137/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180309.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143215/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180400.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143342/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180481.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143437/*.root*");
 
  //period E
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180614.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143454/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180636.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143524/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180664.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143602_shadow/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180664.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143602_sub035297393/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180710.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143629/*.root*");
 tree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180776.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143650/*.root*");

  //period F
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182013.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143709/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182161.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143729/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182284.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143749/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182346.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143808/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182372.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143832/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182424.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143854/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182449.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144458/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182450.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144552/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182454.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144631/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182455.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144713/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182456.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144733/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182486.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144808/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182518.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145608/*.root*");
tree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182519.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145642/*.root*");

 


  //period G
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182726.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145706/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182747.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145730/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182766.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145755/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182787.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145820/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182796.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145842/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182879.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145907/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182886.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145935/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182997.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150000/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183021.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150104/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183038.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150133/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183045.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150201/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183054.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150230/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183078.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150301/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183079.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150335/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183112.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150432/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183127.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150500/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183130.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150643/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183216.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150714/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183272.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150753/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183286.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150820/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183391.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150920/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183407.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150947/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183412.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151014/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183426.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151656/*.root*");
tree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183462.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151739/*.root*");


  //period H
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183544.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151810/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183580.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151833/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183602.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151922/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183780.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151957/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183963.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152027/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184066.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152142/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184072.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152233/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184088.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152335/*.root*");
tree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184169.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152426/*.root*");

  //period I
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185353.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152513/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185518.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152538/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185536.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152608/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185644.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152637/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185649.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152702/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185731.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152727/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185761.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152830/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185976.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152940/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186049.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153036/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186156.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153109/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186169.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153142/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186179.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153250/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186180.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153318/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186216.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153420/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186275.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153531/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186361.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153610/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186396.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153637/*.root*");
tree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186456.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153735/*.root*");

  //period J

tree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186516.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153844/*.root*");
tree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186532.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153917/*.root*");
tree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186669.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154028/*.root*");
tree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186729.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154220/*.root*");
tree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186753.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154257/*.root*");

  //period K
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186873.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154359/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186877.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154430/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186878.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154505/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186965.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154740/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187014.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154824/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187196.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154905/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187453.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155038/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187543.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155238/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187552.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155321/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187763.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155412/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187811.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155446/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187812.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155520/*.root*");
tree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187815.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155547/*.root*");


  //period L
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00188921.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155613/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00188951.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106155708/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189011.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155734/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189027.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155804/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189079.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155937/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189090.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106160007/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189184.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106160036/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189205.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106160123/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189207.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160154/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189280.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160256/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189288.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160330/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189421.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160638/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189530.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106160937/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189561.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161104/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189602.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161220/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189639.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161359/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189692.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161604/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189693.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161638/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189751.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161738/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189774.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161806/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189813.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m991_p766.btag.v1.111106161923/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189836.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162030/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189845.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162055/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189875.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162120/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189963.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m997_p766.btag.v1.111106162143/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190046.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162244/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190116.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162311/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190120.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162413/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190127.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162447/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190157.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162533/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190162.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162609/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190236.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162641/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190256.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162723/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190295.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162752/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190297.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162824/*.root*");
tree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190300.physics_JetTauEtmiss.merge.NTUP_JETMET.f409_m1007_p766.btag.v1.111106162855/*.root*");


  //period M
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190504.physics_JetTauEtmiss.merge.NTUP_JETMET.f410_m1007_p766.btag.v1.111106163035/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190505.physics_JetTauEtmiss.merge.NTUP_JETMET.f410_m1007_p766.btag.v1.111106163102/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190618.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163253/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190643.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163323/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190644.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163358/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190689.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163458/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190878.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1014_p766.btag.v1.111106163619/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191138.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163918/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191150.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1019_p766.btag.v1.111106164129/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191190.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164206/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191235.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164423/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191239.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164451/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191343.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164633/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191373.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164758/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191376.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164824/*.root*");
tree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191425.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164911/*.root*");





  cout << "Number of events = " << tree->GetEntries() << endl;




  //==========================================================
  //meta data tree

  //m_metaDataM_MetaDataTree->Add(" /*.root*");
  

  //    m_metaDataTree->Add("/4C/laugs/Data2011/user.aia.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2276_p516_p523_p621.SKIM.V0.111014120538/*.root*");

  
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106141952/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178021.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142325/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178026.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142345/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178047.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142425/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00178109.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142449/*.root*");
  
  //period D
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179710.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142512/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179739.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142552/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179771.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142615/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179804.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142634/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179938.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142651/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179939.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142719/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00179940.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142742/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180122.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142803/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180139.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142856/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180144.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142916/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180149.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106142937/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180153.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143011/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180225.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143120/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180241.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143137/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180309.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143215/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180400.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143342/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-D/user.gotero.data11_7TeV.00180481.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143437/*.root*");

 
  //period E
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180614.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143454/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180636.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143524/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180664.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143602_shadow/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180664.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143602_sub035297393/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180710.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143629/*.root*");
  m_metaDataTree->Add("/4B/laugs/data11/r17/period-E/user.gotero.data11_7TeV.00180776.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143650/*.root*");

  //period F
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182013.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143709/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182161.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143729/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182284.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143749/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182346.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143808/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182372.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143832/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182424.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106143854/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182449.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144458/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182450.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144552/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182454.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144631/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182455.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144713/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182456.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144733/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182486.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106144808/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182518.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145608/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-F/user.gotero.data11_7TeV.00182519.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145642/*.root*");

 
  //period G
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182726.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145706/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182747.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145730/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182766.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145755/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182787.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145820/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182796.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145842/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182879.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145907/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182886.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106145935/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00182997.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150000/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183021.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150104/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183038.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150133/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183045.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150201/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183054.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150230/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183078.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150301/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183079.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150335/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183112.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150432/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183127.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150500/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183130.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150643/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183216.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150714/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183272.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150753/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183286.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150820/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183391.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150920/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183407.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106150947/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183412.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151014/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183426.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151656/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-G/user.gotero.data11_7TeV.00183462.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151739/*.root*");
 
   //period H
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183544.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151810/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183580.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151833/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183602.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151922/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183780.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106151957/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00183963.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152027/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184066.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152142/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184072.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152233/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184088.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152335/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-H/user.gotero.data11_7TeV.00184169.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106152426/*.root*");
 
  //period I
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185353.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152513/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185518.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152538/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185536.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152608/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185644.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152637/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185649.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152702/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185731.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152727/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185761.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152830/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00185976.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106152940/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186049.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153036/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186156.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153109/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186169.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153142/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186179.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153250/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186180.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153318/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186216.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153420/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186275.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153531/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186361.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153610/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186396.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153637/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-I/user.gotero.data11_7TeV.00186456.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153735/*.root*");
 
  //period J

 m_metaDataTree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186516.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153844/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186532.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106153917/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186669.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154028/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186729.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154220/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-J/user.gotero.data11_7TeV.00186753.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154257/*.root*");

  //period K
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186873.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154359/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186877.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154430/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186878.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154505/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00186965.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154740/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187014.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154824/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187196.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106154905/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187453.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155038/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187543.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155238/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187552.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155321/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187763.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155412/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187811.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155446/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187812.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155520/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-K/user.gotero.data11_7TeV.00187815.physics_JetTauEtmiss.merge.NTUP_JETMET.r2713_p705_p766.btag.v1.111106155547/*.root*");
 

  //period L
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00188921.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155613/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00188951.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106155708/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189011.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155734/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189027.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155804/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189079.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106155937/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189090.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m975_p766.btag.v1.111106160007/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189184.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106160036/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189205.physics_JetTauEtmiss.merge.NTUP_JETMET.f403_m980_p766.btag.v1.111106160123/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189207.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160154/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189280.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160256/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189288.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160330/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189421.physics_JetTauEtmiss.merge.NTUP_JETMET.f404_m980_p766.btag.v1.111106160638/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189530.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106160937/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189561.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161104/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189602.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161220/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189639.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161359/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189692.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161604/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189693.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161638/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189751.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161738/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189774.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m985_p766.btag.v1.111106161806/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189813.physics_JetTauEtmiss.merge.NTUP_JETMET.f405_m991_p766.btag.v1.111106161923/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189836.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162030/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189845.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162055/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189875.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m991_p766.btag.v1.111106162120/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00189963.physics_JetTauEtmiss.merge.NTUP_JETMET.f406_m997_p766.btag.v1.111106162143/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190046.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162244/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190116.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162311/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190120.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162413/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190127.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162447/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190157.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162533/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190162.physics_JetTauEtmiss.merge.NTUP_JETMET.f407_m997_p766.btag.v1.111106162609/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190236.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162641/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190256.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162723/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190295.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162752/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190297.physics_JetTauEtmiss.merge.NTUP_JETMET.f408_m1007_p766.btag.v1.111106162824/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-L/user.gotero.data11_7TeV.00190300.physics_JetTauEtmiss.merge.NTUP_JETMET.f409_m1007_p766.btag.v1.111106162855/*.root*");


  //period M
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190504.physics_JetTauEtmiss.merge.NTUP_JETMET.f410_m1007_p766.btag.v1.111106163035/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190505.physics_JetTauEtmiss.merge.NTUP_JETMET.f410_m1007_p766.btag.v1.111106163102/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190618.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163253/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190643.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163323/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190644.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163358/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190689.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163458/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00190878.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1014_p766.btag.v1.111106163619/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191138.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1007_p766.btag.v1.111106163918/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191150.physics_JetTauEtmiss.merge.NTUP_JETMET.f411_m1019_p766.btag.v1.111106164129/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191190.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164206/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191235.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164423/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191239.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164451/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191343.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164633/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191373.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164758/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191376.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164824/*.root*");
 m_metaDataTree->Add("/4B/laugs/data11/r17/period-M/user.gotero.data11_7TeV.00191425.physics_JetTauEtmiss.merge.NTUP_JETMET.f413_m1019_p766.btag.v1.111106164911/*.root*");
 


 // cout << "Number of events = " << m_metaDataTree->GetEntries() << endl;


  //tree = BranchSelector(tree);
  myAnalysis* y = new myAnalysis(tree,m_metaDataTree);
  y.FillTrees(isDATA,isMC);
  //  y.MakeHistos(isDATA,isMC);
  //  y.EvaluateNN();
  //y.EvaluateSamples();
  //y.EvaluateNN();




  }else if(isMC){




 TChain* tree = new TChain("qcd");
  tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105009.J0_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v1.111226181252/*.root*");
   tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105009.J0_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v1.111226181315/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105009.J0_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v1.111226181315_r1/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105010.J1_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226232832/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105010.J1_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226232905/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105011.J2_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226232951/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105011.J2_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233009/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105012.J3_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233036/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105012.J3_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233101/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105013.J4_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233144/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105013.J4_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233217/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105014.J5_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233311/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105015.J6_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233340/*.root*");
 tree->Add("/2B/laugs/MC11b/SKIM_p832/user.laugs.mc11_7TeV.105015.J6_pythia_jetjet.merge.NTUP_JETMET.e815_s1273_s1274_r2923_r2900_p832.slim.skim.v0.111226233355/*.root*");

  

   cout << "Number of events = " << tree->GetEntries() << endl;
   myAnalysis* y = new myAnalysis(tree);
   //  y.EvaluateNN();    
   //  y.MakeHistos();
   //y.MakeTemplates();
   //y.EvaluateNNmore();
   //y.JetSelectionStats(); 
 
  }

}
