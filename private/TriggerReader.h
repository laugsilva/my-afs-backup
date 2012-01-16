#ifndef TRIGGER_READER_H
#define TRIGGER_READER_H

//////////////////////////////////////////////////////////
//Reader of D3PD trigger metadata
//based on example from Attila Krasznahorkay : 
//https://twiki.cern.ch/twiki/bin/view/Atlas/TriggerD3PDMaker
//and on the adaptation by Nils Ruthmann
//
// georgios.choudalakis@cern.ch , Oct. 2010
// peggiorato da V.Giangiobbe
////////////////////////////////////////////////////////////

#include <iostream>
#include <bitset>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"

class TriggerReader {
 public:
 
  typedef std::pair< UInt_t, std::pair< UInt_t, UInt_t > > DBKeys_t;
  
  TriggerReader(TTree* conftree=0);
  ~TriggerReader();
  
  void InitTrigMetaData();

  Bool_t ItemPassed( const std::vector< unsigned int >& ctp_bits, size_t bitpos );

  bool PassesEF(std::string triggerName, UInt_t SMK, UInt_t L1PSK, UInt_t HLTPSK, std::vector<short>* EF_passedPhysics);
  double GetPrescaleEF(std::string triggerName, UInt_t SMK, UInt_t L1PSK, UInt_t HLTPSK);
  
  bool PassesL1(std::string triggerName, UInt_t SMK, UInt_t L1PSK, UInt_t HLTPSK, std::vector<unsigned int>* L1_TAV);
  double GetPrescaleL1(std::string triggerName, UInt_t SMK, UInt_t L1PSK, UInt_t HLTPSK);
  
  void PrintNames(UInt_t SMK, UInt_t L1PSK, UInt_t HLTPSK); //print names and prescales

  void StoreTriggerTreeNoDuplicates(std::string outputFilename, std::string option="recreate");

  int currentMetaEntry;
  std::map< DBKeys_t, Int_t > configMap;
  TTree* MetaChain;

  Int_t SMK_conf;
  Int_t L1PSK_conf;
  Int_t HLTPSK_conf;
  std::map< std::string, int >*   LVL1NameMap;	 
  std::map< std::string, float >* LVL1PrescaleMap;
  std::map< std::string, int >*   HLTNameMap;	 
  std::map< std::string, float >* HLTPrescaleMap;
  
};

#endif // TRIGGER_READER_H


