//*****************************************************
// Class   : myAnalysis
// 
// Purpose : 
//           
//                      
// Author  :
//
// Date    : XX-XXX-2008 : Creation of class
//*****************************************************
#ifndef myAnalysis_h
#define myAnalysis_h

#include "AnalysisBase.h"
#include "TMatrixD.h"
  //#include "CollectionTree.h"
#include "TTree.h"
#include "TriggerReader.h"

  
class myAnalysis : public AnalysisBase
{
 public:
  myAnalysis(TTree *tree=0,TTree* metaDataTree=0);
  ~myAnalysis();
  
  
 
  //  void MakeHistos(int x);
  void MakeHistos(bool isdata,bool ismc);
  void FillTrees(bool isdata,bool ismc);
  // void EvaluateNN();
 
  bool RemoveDeadRegion(float jeteta, float jetphi);
  bool PassedTrigger11(UInt_t runNumber,
		       UInt_t SMK, 
		       UInt_t L1PSK, 
		       UInt_t HLTPSK, 
		       vector<unsigned int>* L1_TAV,
		       vector<short>* EF_passedPhysics,
		       double pt); // , double &prescale);
		     //		      double pt1, double pt2);
  private:

  TTree* m_confTree;
  TriggerReader  m_triggerReader; 

  ClassDef(myAnalysis,1)  //Adding class to ROOT?
};


#endif // #ifdef myAnalysis_cxx
