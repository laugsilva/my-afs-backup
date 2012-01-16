//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 29 22:48:47 2011 by ROOT version 5.22/00i
// from TTree qcd/qcd
// found on file: /4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106141952/user.gotero.004075._00001.qcd.root
//////////////////////////////////////////////////////////

#ifndef CollectionTree_h
#define CollectionTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

class CollectionTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   UInt_t          bcid;
   Int_t           vxp_n;
   vector<float>   *vxp_z;
   vector<int>     *vxp_nTracks;
   Int_t           jet_AntiKt4TopoEM_n;
   vector<float>   *jet_AntiKt4TopoEM_E;
   vector<float>   *jet_AntiKt4TopoEM_pt;
   vector<float>   *jet_AntiKt4TopoEM_m;
   vector<float>   *jet_AntiKt4TopoEM_eta;
   vector<float>   *jet_AntiKt4TopoEM_phi;
   vector<float>   *jet_AntiKt4TopoEM_WIDTH;
   vector<float>   *jet_AntiKt4TopoEM_LArQuality;
   vector<float>   *jet_AntiKt4TopoEM_nTrk;
   vector<float>   *jet_AntiKt4TopoEM_sumPtTrk;
   vector<float>   *jet_AntiKt4TopoEM_NegativeE;
   vector<float>   *jet_AntiKt4TopoEM_BCH_CORR_CELL;
   vector<float>   *jet_AntiKt4TopoEM_BCH_CORR_DOTX;
   vector<float>   *jet_AntiKt4TopoEM_ENG_BAD_CELLS;
   vector<int>     *jet_AntiKt4TopoEM_SamplingMax;
   vector<int>     *jet_AntiKt4TopoEM_isUgly;
   vector<int>     *jet_AntiKt4TopoEM_isBadLoose;
   vector<int>     *jet_AntiKt4TopoEM_isBadMedium;
   vector<int>     *jet_AntiKt4TopoEM_isBadTight;
   vector<float>   *jet_AntiKt4TopoEM_emfrac;
   vector<float>   *jet_AntiKt4TopoEM_EMJES;
   vector<float>   *jet_AntiKt4TopoEM_EMJES_EtaCorr;
   vector<float>   *jet_AntiKt4TopoEM_EMJESnooffset;
   vector<float>   *jet_AntiKt4TopoEM_GCWJES;
   vector<float>   *jet_AntiKt4TopoEM_emscale_E;
   vector<float>   *jet_AntiKt4TopoEM_emscale_pt;
   vector<float>   *jet_AntiKt4TopoEM_emscale_eta;
   vector<float>   *jet_AntiKt4TopoEM_emscale_phi;
   vector<float>   *jet_AntiKt4TopoEM_jvtxf;
   vector<float>   *jet_AntiKt4TopoEM_GSCFactorF;
   vector<float>   *jet_AntiKt4TopoEM_WidthFraction;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_IP2D;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_IP3D;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_SV0;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_SV1;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_SV2;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_JetProb;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN;
   vector<float>   *jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_ip3d_pu;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_ip3d_pb;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv1_pu;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv1_pb;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv2_pu;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv2_pb;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_jfit_pu;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_jfit_pb;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_jfit_efrc;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_jfit_mass;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_jfit_sig3d;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_svp_ntrkv;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_svp_ntrkj;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_svp_n2t;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_svp_mass;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_svp_efrc;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_svp_ntrk;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_sv0p_n2t;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv0p_mass;
   vector<float>   *jet_AntiKt4TopoEM_flavor_component_sv0p_efrc;
   vector<int>     *jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk;
   vector<int>     *jet_AntiKt4TopoEM_L1_matched;
   vector<int>     *jet_AntiKt4TopoEM_L2_matched;
   vector<int>     *jet_AntiKt4TopoEM_EF_matched;
   Int_t           jet_AntiKt6TopoEM_n;
   vector<float>   *jet_AntiKt6TopoEM_E;
   vector<float>   *jet_AntiKt6TopoEM_pt;
   vector<float>   *jet_AntiKt6TopoEM_m;
   vector<float>   *jet_AntiKt6TopoEM_eta;
   vector<float>   *jet_AntiKt6TopoEM_phi;
   vector<float>   *jet_AntiKt6TopoEM_WIDTH;
   vector<float>   *jet_AntiKt6TopoEM_LArQuality;
   vector<float>   *jet_AntiKt6TopoEM_nTrk;
   vector<float>   *jet_AntiKt6TopoEM_sumPtTrk;
   vector<float>   *jet_AntiKt6TopoEM_NegativeE;
   vector<float>   *jet_AntiKt6TopoEM_BCH_CORR_CELL;
   vector<float>   *jet_AntiKt6TopoEM_BCH_CORR_DOTX;
   vector<float>   *jet_AntiKt6TopoEM_ENG_BAD_CELLS;
   vector<int>     *jet_AntiKt6TopoEM_SamplingMax;
   vector<int>     *jet_AntiKt6TopoEM_isUgly;
   vector<int>     *jet_AntiKt6TopoEM_isBadLoose;
   vector<int>     *jet_AntiKt6TopoEM_isBadMedium;
   vector<int>     *jet_AntiKt6TopoEM_isBadTight;
   vector<float>   *jet_AntiKt6TopoEM_emfrac;
   vector<float>   *jet_AntiKt6TopoEM_EMJES;
   vector<float>   *jet_AntiKt6TopoEM_EMJES_EtaCorr;
   vector<float>   *jet_AntiKt6TopoEM_EMJESnooffset;
   vector<float>   *jet_AntiKt6TopoEM_GCWJES;
   vector<float>   *jet_AntiKt6TopoEM_emscale_E;
   vector<float>   *jet_AntiKt6TopoEM_emscale_pt;
   vector<float>   *jet_AntiKt6TopoEM_emscale_eta;
   vector<float>   *jet_AntiKt6TopoEM_emscale_phi;
   vector<float>   *jet_AntiKt6TopoEM_jvtxf;
   vector<float>   *jet_AntiKt6TopoEM_GSCFactorF;
   vector<float>   *jet_AntiKt6TopoEM_WidthFraction;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_IP2D;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_IP3D;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_SV0;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_SV1;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_SV2;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_JetProb;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN;
   vector<float>   *jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_ip3d_pu;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_ip3d_pb;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv1_pu;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv1_pb;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv2_pu;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv2_pb;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_jfit_pu;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_jfit_pb;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_jfit_efrc;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_jfit_mass;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_jfit_sig3d;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_svp_ntrkv;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_svp_ntrkj;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_svp_n2t;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_svp_mass;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_svp_efrc;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_svp_ntrk;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_sv0p_n2t;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv0p_mass;
   vector<float>   *jet_AntiKt6TopoEM_flavor_component_sv0p_efrc;
   vector<int>     *jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk;
   vector<int>     *jet_AntiKt6TopoEM_L1_matched;
   vector<int>     *jet_AntiKt6TopoEM_L2_matched;
   vector<int>     *jet_AntiKt6TopoEM_EF_matched;
   Int_t           jet_AntiKt4LCTopo_n;
   vector<float>   *jet_AntiKt4LCTopo_E;
   vector<float>   *jet_AntiKt4LCTopo_pt;
   vector<float>   *jet_AntiKt4LCTopo_m;
   vector<float>   *jet_AntiKt4LCTopo_eta;
   vector<float>   *jet_AntiKt4LCTopo_phi;
   vector<float>   *jet_AntiKt4LCTopo_WIDTH;
   vector<float>   *jet_AntiKt4LCTopo_LArQuality;
   vector<float>   *jet_AntiKt4LCTopo_nTrk;
   vector<float>   *jet_AntiKt4LCTopo_sumPtTrk;
   vector<float>   *jet_AntiKt4LCTopo_NegativeE;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_CELL;
   vector<float>   *jet_AntiKt4LCTopo_BCH_CORR_DOTX;
   vector<float>   *jet_AntiKt4LCTopo_ENG_BAD_CELLS;
   vector<int>     *jet_AntiKt4LCTopo_SamplingMax;
   vector<int>     *jet_AntiKt4LCTopo_isUgly;
   vector<int>     *jet_AntiKt4LCTopo_isBadLoose;
   vector<int>     *jet_AntiKt4LCTopo_isBadMedium;
   vector<int>     *jet_AntiKt4LCTopo_isBadTight;
   vector<float>   *jet_AntiKt4LCTopo_emfrac;
   vector<float>   *jet_AntiKt4LCTopo_EMJES;
   vector<float>   *jet_AntiKt4LCTopo_EMJES_EtaCorr;
   vector<float>   *jet_AntiKt4LCTopo_EMJESnooffset;
   vector<float>   *jet_AntiKt4LCTopo_GCWJES;
   vector<float>   *jet_AntiKt4LCTopo_emscale_E;
   vector<float>   *jet_AntiKt4LCTopo_emscale_pt;
   vector<float>   *jet_AntiKt4LCTopo_emscale_eta;
   vector<float>   *jet_AntiKt4LCTopo_emscale_phi;
   vector<float>   *jet_AntiKt4LCTopo_jvtxf;
   vector<float>   *jet_AntiKt4LCTopo_GSCFactorF;
   vector<float>   *jet_AntiKt4LCTopo_WidthFraction;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_IP2D;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_IP3D;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV0;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV1;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SV2;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetProb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN;
   vector<float>   *jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip3d_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_ip3d_pb;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv1_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv1_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv2_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv2_pb;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_pu;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_pb;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_efrc;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_jfit_sig3d;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ntrkv;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ntrkj;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_n2t;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_svp_efrc;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_svp_ntrk;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_n2t;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_mass;
   vector<float>   *jet_AntiKt4LCTopo_flavor_component_sv0p_efrc;
   vector<int>     *jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk;
   vector<int>     *jet_AntiKt4LCTopo_L1_matched;
   vector<int>     *jet_AntiKt4LCTopo_L2_matched;
   vector<int>     *jet_AntiKt4LCTopo_EF_matched;
   Int_t           jet_AntiKt6LCTopo_n;
   vector<float>   *jet_AntiKt6LCTopo_E;
   vector<float>   *jet_AntiKt6LCTopo_pt;
   vector<float>   *jet_AntiKt6LCTopo_m;
   vector<float>   *jet_AntiKt6LCTopo_eta;
   vector<float>   *jet_AntiKt6LCTopo_phi;
   vector<float>   *jet_AntiKt6LCTopo_WIDTH;
   vector<float>   *jet_AntiKt6LCTopo_LArQuality;
   vector<float>   *jet_AntiKt6LCTopo_nTrk;
   vector<float>   *jet_AntiKt6LCTopo_sumPtTrk;
   vector<float>   *jet_AntiKt6LCTopo_NegativeE;
   vector<float>   *jet_AntiKt6LCTopo_BCH_CORR_CELL;
   vector<float>   *jet_AntiKt6LCTopo_BCH_CORR_DOTX;
   vector<float>   *jet_AntiKt6LCTopo_ENG_BAD_CELLS;
   vector<int>     *jet_AntiKt6LCTopo_SamplingMax;
   vector<int>     *jet_AntiKt6LCTopo_isUgly;
   vector<int>     *jet_AntiKt6LCTopo_isBadLoose;
   vector<int>     *jet_AntiKt6LCTopo_isBadMedium;
   vector<int>     *jet_AntiKt6LCTopo_isBadTight;
   vector<float>   *jet_AntiKt6LCTopo_emfrac;
   vector<float>   *jet_AntiKt6LCTopo_EMJES;
   vector<float>   *jet_AntiKt6LCTopo_EMJES_EtaCorr;
   vector<float>   *jet_AntiKt6LCTopo_EMJESnooffset;
   vector<float>   *jet_AntiKt6LCTopo_GCWJES;
   vector<float>   *jet_AntiKt6LCTopo_emscale_E;
   vector<float>   *jet_AntiKt6LCTopo_emscale_pt;
   vector<float>   *jet_AntiKt6LCTopo_emscale_eta;
   vector<float>   *jet_AntiKt6LCTopo_emscale_phi;
   vector<float>   *jet_AntiKt6LCTopo_jvtxf;
   vector<float>   *jet_AntiKt6LCTopo_GSCFactorF;
   vector<float>   *jet_AntiKt6LCTopo_WidthFraction;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_IP2D;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_IP3D;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_SV0;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_SV1;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_SV2;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_JetProb;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN;
   vector<float>   *jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_ip3d_pu;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_ip3d_pb;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv1_pu;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv1_pb;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv2_pu;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv2_pb;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_jfit_pu;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_jfit_pb;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_jfit_efrc;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_jfit_mass;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_jfit_sig3d;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_svp_ntrkv;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_svp_ntrkj;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_svp_n2t;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_svp_mass;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_svp_efrc;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_svp_ntrk;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_sv0p_n2t;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv0p_mass;
   vector<float>   *jet_AntiKt6LCTopo_flavor_component_sv0p_efrc;
   vector<int>     *jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk;
   vector<int>     *jet_AntiKt6LCTopo_L1_matched;
   vector<int>     *jet_AntiKt6LCTopo_L2_matched;
   vector<int>     *jet_AntiKt6LCTopo_EF_matched;
   Int_t           trk_n;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_phi;
   vector<float>   *trk_theta;
   vector<float>   *trk_qoverp;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_d0_wrtPV;
   vector<float>   *trk_z0_wrtPV;
   vector<float>   *trk_phi_wrtPV;
   vector<float>   *trk_cov_d0_wrtPV;
   vector<float>   *trk_cov_z0_wrtPV;
   vector<float>   *trk_cov_phi_wrtPV;
   vector<float>   *trk_cov_theta_wrtPV;
   vector<float>   *trk_cov_qoverp_wrtPV;
   vector<float>   *trk_d0_wrtBS;
   vector<float>   *trk_z0_wrtBS;
   vector<float>   *trk_phi_wrtBS;
   vector<float>   *trk_chi2;
   vector<int>     *trk_ndof;
   vector<int>     *trk_nBLHits;
   vector<int>     *trk_nPixHits;
   vector<int>     *trk_nSCTHits;
   vector<int>     *trk_nTRTHits;
   vector<int>     *trk_nMDTHits;
   vector<int>     *trk_nHits;
   Bool_t          L1_J15;
   vector<unsigned int> *trig_L1_TAV;
   vector<short>   *trig_L2_passedPhysics;
   vector<short>   *trig_EF_passedPhysics;
   vector<unsigned int> *trig_L1_TBP;
   vector<unsigned int> *trig_L1_TAP;
   vector<short>   *trig_L2_passedRaw;
   vector<short>   *trig_EF_passedRaw;
   vector<short>   *trig_L2_resurrected;
   vector<short>   *trig_EF_resurrected;
   vector<short>   *trig_L2_passedThrough;
   vector<short>   *trig_EF_passedThrough;
   UInt_t          trig_DB_SMK;
   UInt_t          trig_DB_L1PSK;
   UInt_t          trig_DB_HLTPSK;
   vector<int>     *trig_EF_jet_EF_j20_a4tc_EFFS;
   vector<int>     *trig_EF_jet_EF_j30_a4tc_EFFS;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_vxp_n;   //!
   TBranch        *b_vxp_z;   //!
   TBranch        *b_vxp_nTracks;   //!
   TBranch        *b_jet_AntiKt4TopoEM_n;   //!
   TBranch        *b_jet_AntiKt4TopoEM_E;   //!
   TBranch        *b_jet_AntiKt4TopoEM_pt;   //!
   TBranch        *b_jet_AntiKt4TopoEM_m;   //!
   TBranch        *b_jet_AntiKt4TopoEM_eta;   //!
   TBranch        *b_jet_AntiKt4TopoEM_phi;   //!
   TBranch        *b_jet_AntiKt4TopoEM_WIDTH;   //!
   TBranch        *b_jet_AntiKt4TopoEM_LArQuality;   //!
   TBranch        *b_jet_AntiKt4TopoEM_nTrk;   //!
   TBranch        *b_jet_AntiKt4TopoEM_sumPtTrk;   //!
   TBranch        *b_jet_AntiKt4TopoEM_NegativeE;   //!
   TBranch        *b_jet_AntiKt4TopoEM_BCH_CORR_CELL;   //!
   TBranch        *b_jet_AntiKt4TopoEM_BCH_CORR_DOTX;   //!
   TBranch        *b_jet_AntiKt4TopoEM_ENG_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt4TopoEM_SamplingMax;   //!
   TBranch        *b_jet_AntiKt4TopoEM_isUgly;   //!
   TBranch        *b_jet_AntiKt4TopoEM_isBadLoose;   //!
   TBranch        *b_jet_AntiKt4TopoEM_isBadMedium;   //!
   TBranch        *b_jet_AntiKt4TopoEM_isBadTight;   //!
   TBranch        *b_jet_AntiKt4TopoEM_emfrac;   //!
   TBranch        *b_jet_AntiKt4TopoEM_EMJES;   //!
   TBranch        *b_jet_AntiKt4TopoEM_EMJES_EtaCorr;   //!
   TBranch        *b_jet_AntiKt4TopoEM_EMJESnooffset;   //!
   TBranch        *b_jet_AntiKt4TopoEM_GCWJES;   //!
   TBranch        *b_jet_AntiKt4TopoEM_emscale_E;   //!
   TBranch        *b_jet_AntiKt4TopoEM_emscale_pt;   //!
   TBranch        *b_jet_AntiKt4TopoEM_emscale_eta;   //!
   TBranch        *b_jet_AntiKt4TopoEM_emscale_phi;   //!
   TBranch        *b_jet_AntiKt4TopoEM_jvtxf;   //!
   TBranch        *b_jet_AntiKt4TopoEM_GSCFactorF;   //!
   TBranch        *b_jet_AntiKt4TopoEM_WidthFraction;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_IP2D;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_IP3D;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_SV0;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_SV1;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_SV2;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_JetProb;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_ip3d_pu;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_ip3d_pb;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv1_pu;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv1_pb;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv2_pu;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv2_pb;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_pu;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_pb;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_efrc;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_mass;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_jfit_sig3d;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_ntrkv;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_ntrkj;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_n2t;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_mass;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_efrc;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_svp_ntrk;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_n2t;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_mass;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_efrc;   //!
   TBranch        *b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk;   //!
   TBranch        *b_jet_AntiKt4TopoEM_L1_matched;   //!
   TBranch        *b_jet_AntiKt4TopoEM_L2_matched;   //!
   TBranch        *b_jet_AntiKt4TopoEM_EF_matched;   //!
   TBranch        *b_jet_AntiKt6TopoEM_n;   //!
   TBranch        *b_jet_AntiKt6TopoEM_E;   //!
   TBranch        *b_jet_AntiKt6TopoEM_pt;   //!
   TBranch        *b_jet_AntiKt6TopoEM_m;   //!
   TBranch        *b_jet_AntiKt6TopoEM_eta;   //!
   TBranch        *b_jet_AntiKt6TopoEM_phi;   //!
   TBranch        *b_jet_AntiKt6TopoEM_WIDTH;   //!
   TBranch        *b_jet_AntiKt6TopoEM_LArQuality;   //!
   TBranch        *b_jet_AntiKt6TopoEM_nTrk;   //!
   TBranch        *b_jet_AntiKt6TopoEM_sumPtTrk;   //!
   TBranch        *b_jet_AntiKt6TopoEM_NegativeE;   //!
   TBranch        *b_jet_AntiKt6TopoEM_BCH_CORR_CELL;   //!
   TBranch        *b_jet_AntiKt6TopoEM_BCH_CORR_DOTX;   //!
   TBranch        *b_jet_AntiKt6TopoEM_ENG_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt6TopoEM_SamplingMax;   //!
   TBranch        *b_jet_AntiKt6TopoEM_isUgly;   //!
   TBranch        *b_jet_AntiKt6TopoEM_isBadLoose;   //!
   TBranch        *b_jet_AntiKt6TopoEM_isBadMedium;   //!
   TBranch        *b_jet_AntiKt6TopoEM_isBadTight;   //!
   TBranch        *b_jet_AntiKt6TopoEM_emfrac;   //!
   TBranch        *b_jet_AntiKt6TopoEM_EMJES;   //!
   TBranch        *b_jet_AntiKt6TopoEM_EMJES_EtaCorr;   //!
   TBranch        *b_jet_AntiKt6TopoEM_EMJESnooffset;   //!
   TBranch        *b_jet_AntiKt6TopoEM_GCWJES;   //!
   TBranch        *b_jet_AntiKt6TopoEM_emscale_E;   //!
   TBranch        *b_jet_AntiKt6TopoEM_emscale_pt;   //!
   TBranch        *b_jet_AntiKt6TopoEM_emscale_eta;   //!
   TBranch        *b_jet_AntiKt6TopoEM_emscale_phi;   //!
   TBranch        *b_jet_AntiKt6TopoEM_jvtxf;   //!
   TBranch        *b_jet_AntiKt6TopoEM_GSCFactorF;   //!
   TBranch        *b_jet_AntiKt6TopoEM_WidthFraction;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_IP2D;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_IP3D;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_SV0;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_SV1;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_SV2;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_JetProb;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_ip3d_pu;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_ip3d_pb;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv1_pu;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv1_pb;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv2_pu;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv2_pb;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_pu;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_pb;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_efrc;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_mass;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_jfit_sig3d;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_ntrkv;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_ntrkj;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_n2t;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_mass;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_efrc;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_svp_ntrk;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_n2t;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_mass;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_efrc;   //!
   TBranch        *b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk;   //!
   TBranch        *b_jet_AntiKt6TopoEM_L1_matched;   //!
   TBranch        *b_jet_AntiKt6TopoEM_L2_matched;   //!
   TBranch        *b_jet_AntiKt6TopoEM_EF_matched;   //!
   TBranch        *b_jet_AntiKt4LCTopo_n;   //!
   TBranch        *b_jet_AntiKt4LCTopo_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_pt;   //!
   TBranch        *b_jet_AntiKt4LCTopo_m;   //!
   TBranch        *b_jet_AntiKt4LCTopo_eta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_WIDTH;   //!
   TBranch        *b_jet_AntiKt4LCTopo_LArQuality;   //!
   TBranch        *b_jet_AntiKt4LCTopo_nTrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_sumPtTrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_NegativeE;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_CELL;   //!
   TBranch        *b_jet_AntiKt4LCTopo_BCH_CORR_DOTX;   //!
   TBranch        *b_jet_AntiKt4LCTopo_ENG_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt4LCTopo_SamplingMax;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isUgly;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadLoose;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadMedium;   //!
   TBranch        *b_jet_AntiKt4LCTopo_isBadTight;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emfrac;   //!
   TBranch        *b_jet_AntiKt4LCTopo_EMJES;   //!
   TBranch        *b_jet_AntiKt4LCTopo_EMJES_EtaCorr;   //!
   TBranch        *b_jet_AntiKt4LCTopo_EMJESnooffset;   //!
   TBranch        *b_jet_AntiKt4LCTopo_GCWJES;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_E;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_pt;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_eta;   //!
   TBranch        *b_jet_AntiKt4LCTopo_emscale_phi;   //!
   TBranch        *b_jet_AntiKt4LCTopo_jvtxf;   //!
   TBranch        *b_jet_AntiKt4LCTopo_GSCFactorF;   //!
   TBranch        *b_jet_AntiKt4LCTopo_WidthFraction;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_IP2D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_IP3D;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV0;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV1;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SV2;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetProb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv1_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv1_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv2_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv2_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_pu;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_pb;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_jfit_sig3d;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkv;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkj;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_n2t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_svp_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_n2t;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_mass;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_efrc;   //!
   TBranch        *b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk;   //!
   TBranch        *b_jet_AntiKt4LCTopo_L1_matched;   //!
   TBranch        *b_jet_AntiKt4LCTopo_L2_matched;   //!
   TBranch        *b_jet_AntiKt4LCTopo_EF_matched;   //!
   TBranch        *b_jet_AntiKt6LCTopo_n;   //!
   TBranch        *b_jet_AntiKt6LCTopo_E;   //!
   TBranch        *b_jet_AntiKt6LCTopo_pt;   //!
   TBranch        *b_jet_AntiKt6LCTopo_m;   //!
   TBranch        *b_jet_AntiKt6LCTopo_eta;   //!
   TBranch        *b_jet_AntiKt6LCTopo_phi;   //!
   TBranch        *b_jet_AntiKt6LCTopo_WIDTH;   //!
   TBranch        *b_jet_AntiKt6LCTopo_LArQuality;   //!
   TBranch        *b_jet_AntiKt6LCTopo_nTrk;   //!
   TBranch        *b_jet_AntiKt6LCTopo_sumPtTrk;   //!
   TBranch        *b_jet_AntiKt6LCTopo_NegativeE;   //!
   TBranch        *b_jet_AntiKt6LCTopo_BCH_CORR_CELL;   //!
   TBranch        *b_jet_AntiKt6LCTopo_BCH_CORR_DOTX;   //!
   TBranch        *b_jet_AntiKt6LCTopo_ENG_BAD_CELLS;   //!
   TBranch        *b_jet_AntiKt6LCTopo_SamplingMax;   //!
   TBranch        *b_jet_AntiKt6LCTopo_isUgly;   //!
   TBranch        *b_jet_AntiKt6LCTopo_isBadLoose;   //!
   TBranch        *b_jet_AntiKt6LCTopo_isBadMedium;   //!
   TBranch        *b_jet_AntiKt6LCTopo_isBadTight;   //!
   TBranch        *b_jet_AntiKt6LCTopo_emfrac;   //!
   TBranch        *b_jet_AntiKt6LCTopo_EMJES;   //!
   TBranch        *b_jet_AntiKt6LCTopo_EMJES_EtaCorr;   //!
   TBranch        *b_jet_AntiKt6LCTopo_EMJESnooffset;   //!
   TBranch        *b_jet_AntiKt6LCTopo_GCWJES;   //!
   TBranch        *b_jet_AntiKt6LCTopo_emscale_E;   //!
   TBranch        *b_jet_AntiKt6LCTopo_emscale_pt;   //!
   TBranch        *b_jet_AntiKt6LCTopo_emscale_eta;   //!
   TBranch        *b_jet_AntiKt6LCTopo_emscale_phi;   //!
   TBranch        *b_jet_AntiKt6LCTopo_jvtxf;   //!
   TBranch        *b_jet_AntiKt6LCTopo_GSCFactorF;   //!
   TBranch        *b_jet_AntiKt6LCTopo_WidthFraction;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_IP2D;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_IP3D;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_SV0;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_SV1;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_SV2;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_JetProb;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_ip3d_pu;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_ip3d_pb;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv1_pu;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv1_pb;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv2_pu;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv2_pb;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_pu;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_pb;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_efrc;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_mass;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_jfit_sig3d;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_ntrkv;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_ntrkj;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_n2t;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_mass;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_efrc;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_svp_ntrk;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_n2t;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_mass;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_efrc;   //!
   TBranch        *b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk;   //!
   TBranch        *b_jet_AntiKt6LCTopo_L1_matched;   //!
   TBranch        *b_jet_AntiKt6LCTopo_L2_matched;   //!
   TBranch        *b_jet_AntiKt6LCTopo_EF_matched;   //!
   TBranch        *b_trk_n;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_theta;   //!
   TBranch        *b_trk_qoverp;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_d0_wrtPV;   //!
   TBranch        *b_trk_z0_wrtPV;   //!
   TBranch        *b_trk_phi_wrtPV;   //!
   TBranch        *b_trk_cov_d0_wrtPV;   //!
   TBranch        *b_trk_cov_z0_wrtPV;   //!
   TBranch        *b_trk_cov_phi_wrtPV;   //!
   TBranch        *b_trk_cov_theta_wrtPV;   //!
   TBranch        *b_trk_cov_qoverp_wrtPV;   //!
   TBranch        *b_trk_d0_wrtBS;   //!
   TBranch        *b_trk_z0_wrtBS;   //!
   TBranch        *b_trk_phi_wrtBS;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_ndof;   //!
   TBranch        *b_trk_nBLHits;   //!
   TBranch        *b_trk_nPixHits;   //!
   TBranch        *b_trk_nSCTHits;   //!
   TBranch        *b_trk_nTRTHits;   //!
   TBranch        *b_trk_nMDTHits;   //!
   TBranch        *b_trk_nHits;   //!
   TBranch        *b_L1_J15;   //!
   TBranch        *b_trig_L1_TAV;   //!
   TBranch        *b_trig_L2_passedPhysics;   //!
   TBranch        *b_trig_EF_passedPhysics;   //!
   TBranch        *b_trig_L1_TBP;   //!
   TBranch        *b_trig_L1_TAP;   //!
   TBranch        *b_trig_L2_passedRaw;   //!
   TBranch        *b_trig_EF_passedRaw;   //!
   TBranch        *b_trig_L2_resurrected;   //!
   TBranch        *b_trig_EF_resurrected;   //!
   TBranch        *b_trig_L2_passedThrough;   //!
   TBranch        *b_trig_EF_passedThrough;   //!
   TBranch        *b_trig_DB_SMK;   //!
   TBranch        *b_trig_DB_L1PSK;   //!
   TBranch        *b_trig_DB_HLTPSK;   //!
   TBranch        *b_trig_EF_jet_EF_j20_a4tc_EFFS;   //!
   TBranch        *b_trig_EF_jet_EF_j30_a4tc_EFFS;   //!

   CollectionTree(TTree *tree=0);
   virtual ~CollectionTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CollectionTree_cxx
CollectionTree::CollectionTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106141952/user.gotero.004075._00001.qcd.root");
      if (!f) {
         f = new TFile("/4B/laugs/data11/r17/period-B/user.gotero.data11_7TeV.00177986.physics_JetTauEtmiss.merge.NTUP_JETMET.r2603_p659_p766.btag.v1.111106141952/user.gotero.004075._00001.qcd.root");
      }
      tree = (TTree*)gDirectory->Get("qcd");

   }
   Init(tree);
}

CollectionTree::~CollectionTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CollectionTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CollectionTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CollectionTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vxp_z = 0;
   vxp_nTracks = 0;
   jet_AntiKt4TopoEM_E = 0;
   jet_AntiKt4TopoEM_pt = 0;
   jet_AntiKt4TopoEM_m = 0;
   jet_AntiKt4TopoEM_eta = 0;
   jet_AntiKt4TopoEM_phi = 0;
   jet_AntiKt4TopoEM_WIDTH = 0;
   jet_AntiKt4TopoEM_LArQuality = 0;
   jet_AntiKt4TopoEM_nTrk = 0;
   jet_AntiKt4TopoEM_sumPtTrk = 0;
   jet_AntiKt4TopoEM_NegativeE = 0;
   jet_AntiKt4TopoEM_BCH_CORR_CELL = 0;
   jet_AntiKt4TopoEM_BCH_CORR_DOTX = 0;
   jet_AntiKt4TopoEM_ENG_BAD_CELLS = 0;
   jet_AntiKt4TopoEM_SamplingMax = 0;
   jet_AntiKt4TopoEM_isUgly = 0;
   jet_AntiKt4TopoEM_isBadLoose = 0;
   jet_AntiKt4TopoEM_isBadMedium = 0;
   jet_AntiKt4TopoEM_isBadTight = 0;
   jet_AntiKt4TopoEM_emfrac = 0;
   jet_AntiKt4TopoEM_EMJES = 0;
   jet_AntiKt4TopoEM_EMJES_EtaCorr = 0;
   jet_AntiKt4TopoEM_EMJESnooffset = 0;
   jet_AntiKt4TopoEM_GCWJES = 0;
   jet_AntiKt4TopoEM_emscale_E = 0;
   jet_AntiKt4TopoEM_emscale_pt = 0;
   jet_AntiKt4TopoEM_emscale_eta = 0;
   jet_AntiKt4TopoEM_emscale_phi = 0;
   jet_AntiKt4TopoEM_jvtxf = 0;
   jet_AntiKt4TopoEM_GSCFactorF = 0;
   jet_AntiKt4TopoEM_WidthFraction = 0;
   jet_AntiKt4TopoEM_flavor_weight_IP2D = 0;
   jet_AntiKt4TopoEM_flavor_weight_IP3D = 0;
   jet_AntiKt4TopoEM_flavor_weight_SV0 = 0;
   jet_AntiKt4TopoEM_flavor_weight_SV1 = 0;
   jet_AntiKt4TopoEM_flavor_weight_SV2 = 0;
   jet_AntiKt4TopoEM_flavor_weight_JetProb = 0;
   jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag = 0;
   jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN = 0;
   jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt4TopoEM_flavor_component_ip3d_pu = 0;
   jet_AntiKt4TopoEM_flavor_component_ip3d_pb = 0;
   jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk = 0;
   jet_AntiKt4TopoEM_flavor_component_sv1_pu = 0;
   jet_AntiKt4TopoEM_flavor_component_sv1_pb = 0;
   jet_AntiKt4TopoEM_flavor_component_sv2_pu = 0;
   jet_AntiKt4TopoEM_flavor_component_sv2_pb = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_pu = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_pb = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_efrc = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_mass = 0;
   jet_AntiKt4TopoEM_flavor_component_jfit_sig3d = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_ntrkv = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_ntrkj = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_n2t = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_mass = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_efrc = 0;
   jet_AntiKt4TopoEM_flavor_component_svp_ntrk = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_n2t = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_mass = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_efrc = 0;
   jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk = 0;
   jet_AntiKt4TopoEM_L1_matched = 0;
   jet_AntiKt4TopoEM_L2_matched = 0;
   jet_AntiKt4TopoEM_EF_matched = 0;
   jet_AntiKt6TopoEM_E = 0;
   jet_AntiKt6TopoEM_pt = 0;
   jet_AntiKt6TopoEM_m = 0;
   jet_AntiKt6TopoEM_eta = 0;
   jet_AntiKt6TopoEM_phi = 0;
   jet_AntiKt6TopoEM_WIDTH = 0;
   jet_AntiKt6TopoEM_LArQuality = 0;
   jet_AntiKt6TopoEM_nTrk = 0;
   jet_AntiKt6TopoEM_sumPtTrk = 0;
   jet_AntiKt6TopoEM_NegativeE = 0;
   jet_AntiKt6TopoEM_BCH_CORR_CELL = 0;
   jet_AntiKt6TopoEM_BCH_CORR_DOTX = 0;
   jet_AntiKt6TopoEM_ENG_BAD_CELLS = 0;
   jet_AntiKt6TopoEM_SamplingMax = 0;
   jet_AntiKt6TopoEM_isUgly = 0;
   jet_AntiKt6TopoEM_isBadLoose = 0;
   jet_AntiKt6TopoEM_isBadMedium = 0;
   jet_AntiKt6TopoEM_isBadTight = 0;
   jet_AntiKt6TopoEM_emfrac = 0;
   jet_AntiKt6TopoEM_EMJES = 0;
   jet_AntiKt6TopoEM_EMJES_EtaCorr = 0;
   jet_AntiKt6TopoEM_EMJESnooffset = 0;
   jet_AntiKt6TopoEM_GCWJES = 0;
   jet_AntiKt6TopoEM_emscale_E = 0;
   jet_AntiKt6TopoEM_emscale_pt = 0;
   jet_AntiKt6TopoEM_emscale_eta = 0;
   jet_AntiKt6TopoEM_emscale_phi = 0;
   jet_AntiKt6TopoEM_jvtxf = 0;
   jet_AntiKt6TopoEM_GSCFactorF = 0;
   jet_AntiKt6TopoEM_WidthFraction = 0;
   jet_AntiKt6TopoEM_flavor_weight_IP2D = 0;
   jet_AntiKt6TopoEM_flavor_weight_IP3D = 0;
   jet_AntiKt6TopoEM_flavor_weight_SV0 = 0;
   jet_AntiKt6TopoEM_flavor_weight_SV1 = 0;
   jet_AntiKt6TopoEM_flavor_weight_SV2 = 0;
   jet_AntiKt6TopoEM_flavor_weight_JetProb = 0;
   jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag = 0;
   jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN = 0;
   jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt6TopoEM_flavor_component_ip3d_pu = 0;
   jet_AntiKt6TopoEM_flavor_component_ip3d_pb = 0;
   jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk = 0;
   jet_AntiKt6TopoEM_flavor_component_sv1_pu = 0;
   jet_AntiKt6TopoEM_flavor_component_sv1_pb = 0;
   jet_AntiKt6TopoEM_flavor_component_sv2_pu = 0;
   jet_AntiKt6TopoEM_flavor_component_sv2_pb = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_pu = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_pb = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_efrc = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_mass = 0;
   jet_AntiKt6TopoEM_flavor_component_jfit_sig3d = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_ntrkv = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_ntrkj = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_n2t = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_mass = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_efrc = 0;
   jet_AntiKt6TopoEM_flavor_component_svp_ntrk = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_n2t = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_mass = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_efrc = 0;
   jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk = 0;
   jet_AntiKt6TopoEM_L1_matched = 0;
   jet_AntiKt6TopoEM_L2_matched = 0;
   jet_AntiKt6TopoEM_EF_matched = 0;
   jet_AntiKt4LCTopo_E = 0;
   jet_AntiKt4LCTopo_pt = 0;
   jet_AntiKt4LCTopo_m = 0;
   jet_AntiKt4LCTopo_eta = 0;
   jet_AntiKt4LCTopo_phi = 0;
   jet_AntiKt4LCTopo_WIDTH = 0;
   jet_AntiKt4LCTopo_LArQuality = 0;
   jet_AntiKt4LCTopo_nTrk = 0;
   jet_AntiKt4LCTopo_sumPtTrk = 0;
   jet_AntiKt4LCTopo_NegativeE = 0;
   jet_AntiKt4LCTopo_BCH_CORR_CELL = 0;
   jet_AntiKt4LCTopo_BCH_CORR_DOTX = 0;
   jet_AntiKt4LCTopo_ENG_BAD_CELLS = 0;
   jet_AntiKt4LCTopo_SamplingMax = 0;
   jet_AntiKt4LCTopo_isUgly = 0;
   jet_AntiKt4LCTopo_isBadLoose = 0;
   jet_AntiKt4LCTopo_isBadMedium = 0;
   jet_AntiKt4LCTopo_isBadTight = 0;
   jet_AntiKt4LCTopo_emfrac = 0;
   jet_AntiKt4LCTopo_EMJES = 0;
   jet_AntiKt4LCTopo_EMJES_EtaCorr = 0;
   jet_AntiKt4LCTopo_EMJESnooffset = 0;
   jet_AntiKt4LCTopo_GCWJES = 0;
   jet_AntiKt4LCTopo_emscale_E = 0;
   jet_AntiKt4LCTopo_emscale_pt = 0;
   jet_AntiKt4LCTopo_emscale_eta = 0;
   jet_AntiKt4LCTopo_emscale_phi = 0;
   jet_AntiKt4LCTopo_jvtxf = 0;
   jet_AntiKt4LCTopo_GSCFactorF = 0;
   jet_AntiKt4LCTopo_WidthFraction = 0;
   jet_AntiKt4LCTopo_flavor_weight_IP2D = 0;
   jet_AntiKt4LCTopo_flavor_weight_IP3D = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV0 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV1 = 0;
   jet_AntiKt4LCTopo_flavor_weight_SV2 = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetProb = 0;
   jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN = 0;
   jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_sv1_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_sv1_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_sv2_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_sv2_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_pu = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_pb = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_jfit_sig3d = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ntrkv = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ntrkj = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_n2t = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_svp_ntrk = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_n2t = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_mass = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_efrc = 0;
   jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk = 0;
   jet_AntiKt4LCTopo_L1_matched = 0;
   jet_AntiKt4LCTopo_L2_matched = 0;
   jet_AntiKt4LCTopo_EF_matched = 0;
   jet_AntiKt6LCTopo_E = 0;
   jet_AntiKt6LCTopo_pt = 0;
   jet_AntiKt6LCTopo_m = 0;
   jet_AntiKt6LCTopo_eta = 0;
   jet_AntiKt6LCTopo_phi = 0;
   jet_AntiKt6LCTopo_WIDTH = 0;
   jet_AntiKt6LCTopo_LArQuality = 0;
   jet_AntiKt6LCTopo_nTrk = 0;
   jet_AntiKt6LCTopo_sumPtTrk = 0;
   jet_AntiKt6LCTopo_NegativeE = 0;
   jet_AntiKt6LCTopo_BCH_CORR_CELL = 0;
   jet_AntiKt6LCTopo_BCH_CORR_DOTX = 0;
   jet_AntiKt6LCTopo_ENG_BAD_CELLS = 0;
   jet_AntiKt6LCTopo_SamplingMax = 0;
   jet_AntiKt6LCTopo_isUgly = 0;
   jet_AntiKt6LCTopo_isBadLoose = 0;
   jet_AntiKt6LCTopo_isBadMedium = 0;
   jet_AntiKt6LCTopo_isBadTight = 0;
   jet_AntiKt6LCTopo_emfrac = 0;
   jet_AntiKt6LCTopo_EMJES = 0;
   jet_AntiKt6LCTopo_EMJES_EtaCorr = 0;
   jet_AntiKt6LCTopo_EMJESnooffset = 0;
   jet_AntiKt6LCTopo_GCWJES = 0;
   jet_AntiKt6LCTopo_emscale_E = 0;
   jet_AntiKt6LCTopo_emscale_pt = 0;
   jet_AntiKt6LCTopo_emscale_eta = 0;
   jet_AntiKt6LCTopo_emscale_phi = 0;
   jet_AntiKt6LCTopo_jvtxf = 0;
   jet_AntiKt6LCTopo_GSCFactorF = 0;
   jet_AntiKt6LCTopo_WidthFraction = 0;
   jet_AntiKt6LCTopo_flavor_weight_IP2D = 0;
   jet_AntiKt6LCTopo_flavor_weight_IP3D = 0;
   jet_AntiKt6LCTopo_flavor_weight_SV0 = 0;
   jet_AntiKt6LCTopo_flavor_weight_SV1 = 0;
   jet_AntiKt6LCTopo_flavor_weight_SV2 = 0;
   jet_AntiKt6LCTopo_flavor_weight_JetProb = 0;
   jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag = 0;
   jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN = 0;
   jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN = 0;
   jet_AntiKt6LCTopo_flavor_component_ip3d_pu = 0;
   jet_AntiKt6LCTopo_flavor_component_ip3d_pb = 0;
   jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk = 0;
   jet_AntiKt6LCTopo_flavor_component_sv1_pu = 0;
   jet_AntiKt6LCTopo_flavor_component_sv1_pb = 0;
   jet_AntiKt6LCTopo_flavor_component_sv2_pu = 0;
   jet_AntiKt6LCTopo_flavor_component_sv2_pb = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_pu = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_pb = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_efrc = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_mass = 0;
   jet_AntiKt6LCTopo_flavor_component_jfit_sig3d = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_ntrkv = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_ntrkj = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_n2t = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_mass = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_efrc = 0;
   jet_AntiKt6LCTopo_flavor_component_svp_ntrk = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_n2t = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_mass = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_efrc = 0;
   jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk = 0;
   jet_AntiKt6LCTopo_L1_matched = 0;
   jet_AntiKt6LCTopo_L2_matched = 0;
   jet_AntiKt6LCTopo_EF_matched = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_phi = 0;
   trk_theta = 0;
   trk_qoverp = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_d0_wrtPV = 0;
   trk_z0_wrtPV = 0;
   trk_phi_wrtPV = 0;
   trk_cov_d0_wrtPV = 0;
   trk_cov_z0_wrtPV = 0;
   trk_cov_phi_wrtPV = 0;
   trk_cov_theta_wrtPV = 0;
   trk_cov_qoverp_wrtPV = 0;
   trk_d0_wrtBS = 0;
   trk_z0_wrtBS = 0;
   trk_phi_wrtBS = 0;
   trk_chi2 = 0;
   trk_ndof = 0;
   trk_nBLHits = 0;
   trk_nPixHits = 0;
   trk_nSCTHits = 0;
   trk_nTRTHits = 0;
   trk_nMDTHits = 0;
   trk_nHits = 0;
   trig_L1_TAV = 0;
   trig_L2_passedPhysics = 0;
   trig_EF_passedPhysics = 0;
   trig_L1_TBP = 0;
   trig_L1_TAP = 0;
   trig_L2_passedRaw = 0;
   trig_EF_passedRaw = 0;
   trig_L2_resurrected = 0;
   trig_EF_resurrected = 0;
   trig_L2_passedThrough = 0;
   trig_EF_passedThrough = 0;
   trig_EF_jet_EF_j20_a4tc_EFFS = 0;
   trig_EF_jet_EF_j30_a4tc_EFFS = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("vxp_n", &vxp_n, &b_vxp_n);
   fChain->SetBranchAddress("vxp_z", &vxp_z, &b_vxp_z);
   fChain->SetBranchAddress("vxp_nTracks", &vxp_nTracks, &b_vxp_nTracks);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_n", &jet_AntiKt4TopoEM_n, &b_jet_AntiKt4TopoEM_n);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_E", &jet_AntiKt4TopoEM_E, &b_jet_AntiKt4TopoEM_E);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_pt", &jet_AntiKt4TopoEM_pt, &b_jet_AntiKt4TopoEM_pt);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_m", &jet_AntiKt4TopoEM_m, &b_jet_AntiKt4TopoEM_m);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_eta", &jet_AntiKt4TopoEM_eta, &b_jet_AntiKt4TopoEM_eta);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_phi", &jet_AntiKt4TopoEM_phi, &b_jet_AntiKt4TopoEM_phi);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_WIDTH", &jet_AntiKt4TopoEM_WIDTH, &b_jet_AntiKt4TopoEM_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_LArQuality", &jet_AntiKt4TopoEM_LArQuality, &b_jet_AntiKt4TopoEM_LArQuality);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_nTrk", &jet_AntiKt4TopoEM_nTrk, &b_jet_AntiKt4TopoEM_nTrk);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_sumPtTrk", &jet_AntiKt4TopoEM_sumPtTrk, &b_jet_AntiKt4TopoEM_sumPtTrk);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_NegativeE", &jet_AntiKt4TopoEM_NegativeE, &b_jet_AntiKt4TopoEM_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_BCH_CORR_CELL", &jet_AntiKt4TopoEM_BCH_CORR_CELL, &b_jet_AntiKt4TopoEM_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_BCH_CORR_DOTX", &jet_AntiKt4TopoEM_BCH_CORR_DOTX, &b_jet_AntiKt4TopoEM_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_ENG_BAD_CELLS", &jet_AntiKt4TopoEM_ENG_BAD_CELLS, &b_jet_AntiKt4TopoEM_ENG_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_SamplingMax", &jet_AntiKt4TopoEM_SamplingMax, &b_jet_AntiKt4TopoEM_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_isUgly", &jet_AntiKt4TopoEM_isUgly, &b_jet_AntiKt4TopoEM_isUgly);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_isBadLoose", &jet_AntiKt4TopoEM_isBadLoose, &b_jet_AntiKt4TopoEM_isBadLoose);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_isBadMedium", &jet_AntiKt4TopoEM_isBadMedium, &b_jet_AntiKt4TopoEM_isBadMedium);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_isBadTight", &jet_AntiKt4TopoEM_isBadTight, &b_jet_AntiKt4TopoEM_isBadTight);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emfrac", &jet_AntiKt4TopoEM_emfrac, &b_jet_AntiKt4TopoEM_emfrac);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_EMJES", &jet_AntiKt4TopoEM_EMJES, &b_jet_AntiKt4TopoEM_EMJES);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_EMJES_EtaCorr", &jet_AntiKt4TopoEM_EMJES_EtaCorr, &b_jet_AntiKt4TopoEM_EMJES_EtaCorr);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_EMJESnooffset", &jet_AntiKt4TopoEM_EMJESnooffset, &b_jet_AntiKt4TopoEM_EMJESnooffset);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_GCWJES", &jet_AntiKt4TopoEM_GCWJES, &b_jet_AntiKt4TopoEM_GCWJES);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_E", &jet_AntiKt4TopoEM_emscale_E, &b_jet_AntiKt4TopoEM_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_pt", &jet_AntiKt4TopoEM_emscale_pt, &b_jet_AntiKt4TopoEM_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_eta", &jet_AntiKt4TopoEM_emscale_eta, &b_jet_AntiKt4TopoEM_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_emscale_phi", &jet_AntiKt4TopoEM_emscale_phi, &b_jet_AntiKt4TopoEM_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_jvtxf", &jet_AntiKt4TopoEM_jvtxf, &b_jet_AntiKt4TopoEM_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_GSCFactorF", &jet_AntiKt4TopoEM_GSCFactorF, &b_jet_AntiKt4TopoEM_GSCFactorF);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_WidthFraction", &jet_AntiKt4TopoEM_WidthFraction, &b_jet_AntiKt4TopoEM_WidthFraction);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_IP2D", &jet_AntiKt4TopoEM_flavor_weight_IP2D, &b_jet_AntiKt4TopoEM_flavor_weight_IP2D);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_IP3D", &jet_AntiKt4TopoEM_flavor_weight_IP3D, &b_jet_AntiKt4TopoEM_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SV0", &jet_AntiKt4TopoEM_flavor_weight_SV0, &b_jet_AntiKt4TopoEM_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SV1", &jet_AntiKt4TopoEM_flavor_weight_SV1, &b_jet_AntiKt4TopoEM_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SV2", &jet_AntiKt4TopoEM_flavor_weight_SV2, &b_jet_AntiKt4TopoEM_flavor_weight_SV2);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_JetProb", &jet_AntiKt4TopoEM_flavor_weight_JetProb, &b_jet_AntiKt4TopoEM_flavor_weight_JetProb);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag", &jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag, &b_jet_AntiKt4TopoEM_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN", &jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN, &b_jet_AntiKt4TopoEM_flavor_weight_JetFitterTagNN);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN", &jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt4TopoEM_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_ip3d_pu", &jet_AntiKt4TopoEM_flavor_component_ip3d_pu, &b_jet_AntiKt4TopoEM_flavor_component_ip3d_pu);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_ip3d_pb", &jet_AntiKt4TopoEM_flavor_component_ip3d_pb, &b_jet_AntiKt4TopoEM_flavor_component_ip3d_pb);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk", &jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk, &b_jet_AntiKt4TopoEM_flavor_component_ip3d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv1_pu", &jet_AntiKt4TopoEM_flavor_component_sv1_pu, &b_jet_AntiKt4TopoEM_flavor_component_sv1_pu);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv1_pb", &jet_AntiKt4TopoEM_flavor_component_sv1_pb, &b_jet_AntiKt4TopoEM_flavor_component_sv1_pb);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv2_pu", &jet_AntiKt4TopoEM_flavor_component_sv2_pu, &b_jet_AntiKt4TopoEM_flavor_component_sv2_pu);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv2_pb", &jet_AntiKt4TopoEM_flavor_component_sv2_pb, &b_jet_AntiKt4TopoEM_flavor_component_sv2_pb);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_pu", &jet_AntiKt4TopoEM_flavor_component_jfit_pu, &b_jet_AntiKt4TopoEM_flavor_component_jfit_pu);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_pb", &jet_AntiKt4TopoEM_flavor_component_jfit_pb, &b_jet_AntiKt4TopoEM_flavor_component_jfit_pb);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t", &jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t, &b_jet_AntiKt4TopoEM_flavor_component_jfit_nvtx1t);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_efrc", &jet_AntiKt4TopoEM_flavor_component_jfit_efrc, &b_jet_AntiKt4TopoEM_flavor_component_jfit_efrc);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_mass", &jet_AntiKt4TopoEM_flavor_component_jfit_mass, &b_jet_AntiKt4TopoEM_flavor_component_jfit_mass);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_jfit_sig3d", &jet_AntiKt4TopoEM_flavor_component_jfit_sig3d, &b_jet_AntiKt4TopoEM_flavor_component_jfit_sig3d);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_ntrkv", &jet_AntiKt4TopoEM_flavor_component_svp_ntrkv, &b_jet_AntiKt4TopoEM_flavor_component_svp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_ntrkj", &jet_AntiKt4TopoEM_flavor_component_svp_ntrkj, &b_jet_AntiKt4TopoEM_flavor_component_svp_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_n2t", &jet_AntiKt4TopoEM_flavor_component_svp_n2t, &b_jet_AntiKt4TopoEM_flavor_component_svp_n2t);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_mass", &jet_AntiKt4TopoEM_flavor_component_svp_mass, &b_jet_AntiKt4TopoEM_flavor_component_svp_mass);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_efrc", &jet_AntiKt4TopoEM_flavor_component_svp_efrc, &b_jet_AntiKt4TopoEM_flavor_component_svp_efrc);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_svp_ntrk", &jet_AntiKt4TopoEM_flavor_component_svp_ntrk, &b_jet_AntiKt4TopoEM_flavor_component_svp_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv", &jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj", &jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_n2t", &jet_AntiKt4TopoEM_flavor_component_sv0p_n2t, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_n2t);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_mass", &jet_AntiKt4TopoEM_flavor_component_sv0p_mass, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_mass);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_efrc", &jet_AntiKt4TopoEM_flavor_component_sv0p_efrc, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_efrc);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk", &jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk, &b_jet_AntiKt4TopoEM_flavor_component_sv0p_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_L1_matched", &jet_AntiKt4TopoEM_L1_matched, &b_jet_AntiKt4TopoEM_L1_matched);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_L2_matched", &jet_AntiKt4TopoEM_L2_matched, &b_jet_AntiKt4TopoEM_L2_matched);
   fChain->SetBranchAddress("jet_AntiKt4TopoEM_EF_matched", &jet_AntiKt4TopoEM_EF_matched, &b_jet_AntiKt4TopoEM_EF_matched);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_n", &jet_AntiKt6TopoEM_n, &b_jet_AntiKt6TopoEM_n);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_E", &jet_AntiKt6TopoEM_E, &b_jet_AntiKt6TopoEM_E);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_pt", &jet_AntiKt6TopoEM_pt, &b_jet_AntiKt6TopoEM_pt);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_m", &jet_AntiKt6TopoEM_m, &b_jet_AntiKt6TopoEM_m);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_eta", &jet_AntiKt6TopoEM_eta, &b_jet_AntiKt6TopoEM_eta);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_phi", &jet_AntiKt6TopoEM_phi, &b_jet_AntiKt6TopoEM_phi);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_WIDTH", &jet_AntiKt6TopoEM_WIDTH, &b_jet_AntiKt6TopoEM_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_LArQuality", &jet_AntiKt6TopoEM_LArQuality, &b_jet_AntiKt6TopoEM_LArQuality);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_nTrk", &jet_AntiKt6TopoEM_nTrk, &b_jet_AntiKt6TopoEM_nTrk);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_sumPtTrk", &jet_AntiKt6TopoEM_sumPtTrk, &b_jet_AntiKt6TopoEM_sumPtTrk);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_NegativeE", &jet_AntiKt6TopoEM_NegativeE, &b_jet_AntiKt6TopoEM_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_BCH_CORR_CELL", &jet_AntiKt6TopoEM_BCH_CORR_CELL, &b_jet_AntiKt6TopoEM_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_BCH_CORR_DOTX", &jet_AntiKt6TopoEM_BCH_CORR_DOTX, &b_jet_AntiKt6TopoEM_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_ENG_BAD_CELLS", &jet_AntiKt6TopoEM_ENG_BAD_CELLS, &b_jet_AntiKt6TopoEM_ENG_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_SamplingMax", &jet_AntiKt6TopoEM_SamplingMax, &b_jet_AntiKt6TopoEM_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_isUgly", &jet_AntiKt6TopoEM_isUgly, &b_jet_AntiKt6TopoEM_isUgly);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_isBadLoose", &jet_AntiKt6TopoEM_isBadLoose, &b_jet_AntiKt6TopoEM_isBadLoose);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_isBadMedium", &jet_AntiKt6TopoEM_isBadMedium, &b_jet_AntiKt6TopoEM_isBadMedium);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_isBadTight", &jet_AntiKt6TopoEM_isBadTight, &b_jet_AntiKt6TopoEM_isBadTight);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_emfrac", &jet_AntiKt6TopoEM_emfrac, &b_jet_AntiKt6TopoEM_emfrac);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_EMJES", &jet_AntiKt6TopoEM_EMJES, &b_jet_AntiKt6TopoEM_EMJES);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_EMJES_EtaCorr", &jet_AntiKt6TopoEM_EMJES_EtaCorr, &b_jet_AntiKt6TopoEM_EMJES_EtaCorr);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_EMJESnooffset", &jet_AntiKt6TopoEM_EMJESnooffset, &b_jet_AntiKt6TopoEM_EMJESnooffset);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_GCWJES", &jet_AntiKt6TopoEM_GCWJES, &b_jet_AntiKt6TopoEM_GCWJES);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_emscale_E", &jet_AntiKt6TopoEM_emscale_E, &b_jet_AntiKt6TopoEM_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_emscale_pt", &jet_AntiKt6TopoEM_emscale_pt, &b_jet_AntiKt6TopoEM_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_emscale_eta", &jet_AntiKt6TopoEM_emscale_eta, &b_jet_AntiKt6TopoEM_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_emscale_phi", &jet_AntiKt6TopoEM_emscale_phi, &b_jet_AntiKt6TopoEM_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_jvtxf", &jet_AntiKt6TopoEM_jvtxf, &b_jet_AntiKt6TopoEM_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_GSCFactorF", &jet_AntiKt6TopoEM_GSCFactorF, &b_jet_AntiKt6TopoEM_GSCFactorF);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_WidthFraction", &jet_AntiKt6TopoEM_WidthFraction, &b_jet_AntiKt6TopoEM_WidthFraction);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_IP2D", &jet_AntiKt6TopoEM_flavor_weight_IP2D, &b_jet_AntiKt6TopoEM_flavor_weight_IP2D);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_IP3D", &jet_AntiKt6TopoEM_flavor_weight_IP3D, &b_jet_AntiKt6TopoEM_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_SV0", &jet_AntiKt6TopoEM_flavor_weight_SV0, &b_jet_AntiKt6TopoEM_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_SV1", &jet_AntiKt6TopoEM_flavor_weight_SV1, &b_jet_AntiKt6TopoEM_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_SV2", &jet_AntiKt6TopoEM_flavor_weight_SV2, &b_jet_AntiKt6TopoEM_flavor_weight_SV2);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_JetProb", &jet_AntiKt6TopoEM_flavor_weight_JetProb, &b_jet_AntiKt6TopoEM_flavor_weight_JetProb);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag", &jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag, &b_jet_AntiKt6TopoEM_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN", &jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN, &b_jet_AntiKt6TopoEM_flavor_weight_JetFitterTagNN);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN", &jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt6TopoEM_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_ip3d_pu", &jet_AntiKt6TopoEM_flavor_component_ip3d_pu, &b_jet_AntiKt6TopoEM_flavor_component_ip3d_pu);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_ip3d_pb", &jet_AntiKt6TopoEM_flavor_component_ip3d_pb, &b_jet_AntiKt6TopoEM_flavor_component_ip3d_pb);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk", &jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk, &b_jet_AntiKt6TopoEM_flavor_component_ip3d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv1_pu", &jet_AntiKt6TopoEM_flavor_component_sv1_pu, &b_jet_AntiKt6TopoEM_flavor_component_sv1_pu);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv1_pb", &jet_AntiKt6TopoEM_flavor_component_sv1_pb, &b_jet_AntiKt6TopoEM_flavor_component_sv1_pb);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv2_pu", &jet_AntiKt6TopoEM_flavor_component_sv2_pu, &b_jet_AntiKt6TopoEM_flavor_component_sv2_pu);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv2_pb", &jet_AntiKt6TopoEM_flavor_component_sv2_pb, &b_jet_AntiKt6TopoEM_flavor_component_sv2_pb);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_pu", &jet_AntiKt6TopoEM_flavor_component_jfit_pu, &b_jet_AntiKt6TopoEM_flavor_component_jfit_pu);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_pb", &jet_AntiKt6TopoEM_flavor_component_jfit_pb, &b_jet_AntiKt6TopoEM_flavor_component_jfit_pb);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t", &jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t, &b_jet_AntiKt6TopoEM_flavor_component_jfit_nvtx1t);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_efrc", &jet_AntiKt6TopoEM_flavor_component_jfit_efrc, &b_jet_AntiKt6TopoEM_flavor_component_jfit_efrc);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_mass", &jet_AntiKt6TopoEM_flavor_component_jfit_mass, &b_jet_AntiKt6TopoEM_flavor_component_jfit_mass);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_jfit_sig3d", &jet_AntiKt6TopoEM_flavor_component_jfit_sig3d, &b_jet_AntiKt6TopoEM_flavor_component_jfit_sig3d);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_ntrkv", &jet_AntiKt6TopoEM_flavor_component_svp_ntrkv, &b_jet_AntiKt6TopoEM_flavor_component_svp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_ntrkj", &jet_AntiKt6TopoEM_flavor_component_svp_ntrkj, &b_jet_AntiKt6TopoEM_flavor_component_svp_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_n2t", &jet_AntiKt6TopoEM_flavor_component_svp_n2t, &b_jet_AntiKt6TopoEM_flavor_component_svp_n2t);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_mass", &jet_AntiKt6TopoEM_flavor_component_svp_mass, &b_jet_AntiKt6TopoEM_flavor_component_svp_mass);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_efrc", &jet_AntiKt6TopoEM_flavor_component_svp_efrc, &b_jet_AntiKt6TopoEM_flavor_component_svp_efrc);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_svp_ntrk", &jet_AntiKt6TopoEM_flavor_component_svp_ntrk, &b_jet_AntiKt6TopoEM_flavor_component_svp_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv", &jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj", &jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_n2t", &jet_AntiKt6TopoEM_flavor_component_sv0p_n2t, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_n2t);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_mass", &jet_AntiKt6TopoEM_flavor_component_sv0p_mass, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_mass);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_efrc", &jet_AntiKt6TopoEM_flavor_component_sv0p_efrc, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_efrc);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk", &jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk, &b_jet_AntiKt6TopoEM_flavor_component_sv0p_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_L1_matched", &jet_AntiKt6TopoEM_L1_matched, &b_jet_AntiKt6TopoEM_L1_matched);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_L2_matched", &jet_AntiKt6TopoEM_L2_matched, &b_jet_AntiKt6TopoEM_L2_matched);
   fChain->SetBranchAddress("jet_AntiKt6TopoEM_EF_matched", &jet_AntiKt6TopoEM_EF_matched, &b_jet_AntiKt6TopoEM_EF_matched);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_n", &jet_AntiKt4LCTopo_n, &b_jet_AntiKt4LCTopo_n);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_E", &jet_AntiKt4LCTopo_E, &b_jet_AntiKt4LCTopo_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_pt", &jet_AntiKt4LCTopo_pt, &b_jet_AntiKt4LCTopo_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_m", &jet_AntiKt4LCTopo_m, &b_jet_AntiKt4LCTopo_m);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_eta", &jet_AntiKt4LCTopo_eta, &b_jet_AntiKt4LCTopo_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_phi", &jet_AntiKt4LCTopo_phi, &b_jet_AntiKt4LCTopo_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_WIDTH", &jet_AntiKt4LCTopo_WIDTH, &b_jet_AntiKt4LCTopo_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_LArQuality", &jet_AntiKt4LCTopo_LArQuality, &b_jet_AntiKt4LCTopo_LArQuality);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_nTrk", &jet_AntiKt4LCTopo_nTrk, &b_jet_AntiKt4LCTopo_nTrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_sumPtTrk", &jet_AntiKt4LCTopo_sumPtTrk, &b_jet_AntiKt4LCTopo_sumPtTrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_NegativeE", &jet_AntiKt4LCTopo_NegativeE, &b_jet_AntiKt4LCTopo_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_CELL", &jet_AntiKt4LCTopo_BCH_CORR_CELL, &b_jet_AntiKt4LCTopo_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_BCH_CORR_DOTX", &jet_AntiKt4LCTopo_BCH_CORR_DOTX, &b_jet_AntiKt4LCTopo_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_ENG_BAD_CELLS", &jet_AntiKt4LCTopo_ENG_BAD_CELLS, &b_jet_AntiKt4LCTopo_ENG_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_SamplingMax", &jet_AntiKt4LCTopo_SamplingMax, &b_jet_AntiKt4LCTopo_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isUgly", &jet_AntiKt4LCTopo_isUgly, &b_jet_AntiKt4LCTopo_isUgly);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadLoose", &jet_AntiKt4LCTopo_isBadLoose, &b_jet_AntiKt4LCTopo_isBadLoose);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadMedium", &jet_AntiKt4LCTopo_isBadMedium, &b_jet_AntiKt4LCTopo_isBadMedium);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_isBadTight", &jet_AntiKt4LCTopo_isBadTight, &b_jet_AntiKt4LCTopo_isBadTight);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emfrac", &jet_AntiKt4LCTopo_emfrac, &b_jet_AntiKt4LCTopo_emfrac);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EMJES", &jet_AntiKt4LCTopo_EMJES, &b_jet_AntiKt4LCTopo_EMJES);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EMJES_EtaCorr", &jet_AntiKt4LCTopo_EMJES_EtaCorr, &b_jet_AntiKt4LCTopo_EMJES_EtaCorr);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EMJESnooffset", &jet_AntiKt4LCTopo_EMJESnooffset, &b_jet_AntiKt4LCTopo_EMJESnooffset);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_GCWJES", &jet_AntiKt4LCTopo_GCWJES, &b_jet_AntiKt4LCTopo_GCWJES);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_E", &jet_AntiKt4LCTopo_emscale_E, &b_jet_AntiKt4LCTopo_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_pt", &jet_AntiKt4LCTopo_emscale_pt, &b_jet_AntiKt4LCTopo_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_eta", &jet_AntiKt4LCTopo_emscale_eta, &b_jet_AntiKt4LCTopo_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_emscale_phi", &jet_AntiKt4LCTopo_emscale_phi, &b_jet_AntiKt4LCTopo_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_jvtxf", &jet_AntiKt4LCTopo_jvtxf, &b_jet_AntiKt4LCTopo_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_GSCFactorF", &jet_AntiKt4LCTopo_GSCFactorF, &b_jet_AntiKt4LCTopo_GSCFactorF);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_WidthFraction", &jet_AntiKt4LCTopo_WidthFraction, &b_jet_AntiKt4LCTopo_WidthFraction);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_IP2D", &jet_AntiKt4LCTopo_flavor_weight_IP2D, &b_jet_AntiKt4LCTopo_flavor_weight_IP2D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_IP3D", &jet_AntiKt4LCTopo_flavor_weight_IP3D, &b_jet_AntiKt4LCTopo_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV0", &jet_AntiKt4LCTopo_flavor_weight_SV0, &b_jet_AntiKt4LCTopo_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV1", &jet_AntiKt4LCTopo_flavor_weight_SV1, &b_jet_AntiKt4LCTopo_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SV2", &jet_AntiKt4LCTopo_flavor_weight_SV2, &b_jet_AntiKt4LCTopo_flavor_weight_SV2);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetProb", &jet_AntiKt4LCTopo_flavor_weight_JetProb, &b_jet_AntiKt4LCTopo_flavor_weight_JetProb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag", &jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag, &b_jet_AntiKt4LCTopo_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN", &jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterTagNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN", &jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt4LCTopo_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_pu", &jet_AntiKt4LCTopo_flavor_component_ip3d_pu, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_pb", &jet_AntiKt4LCTopo_flavor_component_ip3d_pb, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk", &jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_ip3d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv1_pu", &jet_AntiKt4LCTopo_flavor_component_sv1_pu, &b_jet_AntiKt4LCTopo_flavor_component_sv1_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv1_pb", &jet_AntiKt4LCTopo_flavor_component_sv1_pb, &b_jet_AntiKt4LCTopo_flavor_component_sv1_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv2_pu", &jet_AntiKt4LCTopo_flavor_component_sv2_pu, &b_jet_AntiKt4LCTopo_flavor_component_sv2_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv2_pb", &jet_AntiKt4LCTopo_flavor_component_sv2_pb, &b_jet_AntiKt4LCTopo_flavor_component_sv2_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_pu", &jet_AntiKt4LCTopo_flavor_component_jfit_pu, &b_jet_AntiKt4LCTopo_flavor_component_jfit_pu);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_pb", &jet_AntiKt4LCTopo_flavor_component_jfit_pb, &b_jet_AntiKt4LCTopo_flavor_component_jfit_pb);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t", &jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t, &b_jet_AntiKt4LCTopo_flavor_component_jfit_nvtx1t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_efrc", &jet_AntiKt4LCTopo_flavor_component_jfit_efrc, &b_jet_AntiKt4LCTopo_flavor_component_jfit_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_mass", &jet_AntiKt4LCTopo_flavor_component_jfit_mass, &b_jet_AntiKt4LCTopo_flavor_component_jfit_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_jfit_sig3d", &jet_AntiKt4LCTopo_flavor_component_jfit_sig3d, &b_jet_AntiKt4LCTopo_flavor_component_jfit_sig3d);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ntrkv", &jet_AntiKt4LCTopo_flavor_component_svp_ntrkv, &b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ntrkj", &jet_AntiKt4LCTopo_flavor_component_svp_ntrkj, &b_jet_AntiKt4LCTopo_flavor_component_svp_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_n2t", &jet_AntiKt4LCTopo_flavor_component_svp_n2t, &b_jet_AntiKt4LCTopo_flavor_component_svp_n2t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_mass", &jet_AntiKt4LCTopo_flavor_component_svp_mass, &b_jet_AntiKt4LCTopo_flavor_component_svp_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_efrc", &jet_AntiKt4LCTopo_flavor_component_svp_efrc, &b_jet_AntiKt4LCTopo_flavor_component_svp_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_svp_ntrk", &jet_AntiKt4LCTopo_flavor_component_svp_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_svp_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv", &jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj", &jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_n2t", &jet_AntiKt4LCTopo_flavor_component_sv0p_n2t, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_n2t);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_mass", &jet_AntiKt4LCTopo_flavor_component_sv0p_mass, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_mass);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_efrc", &jet_AntiKt4LCTopo_flavor_component_sv0p_efrc, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_efrc);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk", &jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk, &b_jet_AntiKt4LCTopo_flavor_component_sv0p_ntrk);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_L1_matched", &jet_AntiKt4LCTopo_L1_matched, &b_jet_AntiKt4LCTopo_L1_matched);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_L2_matched", &jet_AntiKt4LCTopo_L2_matched, &b_jet_AntiKt4LCTopo_L2_matched);
   fChain->SetBranchAddress("jet_AntiKt4LCTopo_EF_matched", &jet_AntiKt4LCTopo_EF_matched, &b_jet_AntiKt4LCTopo_EF_matched);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_n", &jet_AntiKt6LCTopo_n, &b_jet_AntiKt6LCTopo_n);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_E", &jet_AntiKt6LCTopo_E, &b_jet_AntiKt6LCTopo_E);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_pt", &jet_AntiKt6LCTopo_pt, &b_jet_AntiKt6LCTopo_pt);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_m", &jet_AntiKt6LCTopo_m, &b_jet_AntiKt6LCTopo_m);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_eta", &jet_AntiKt6LCTopo_eta, &b_jet_AntiKt6LCTopo_eta);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_phi", &jet_AntiKt6LCTopo_phi, &b_jet_AntiKt6LCTopo_phi);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_WIDTH", &jet_AntiKt6LCTopo_WIDTH, &b_jet_AntiKt6LCTopo_WIDTH);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_LArQuality", &jet_AntiKt6LCTopo_LArQuality, &b_jet_AntiKt6LCTopo_LArQuality);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_nTrk", &jet_AntiKt6LCTopo_nTrk, &b_jet_AntiKt6LCTopo_nTrk);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_sumPtTrk", &jet_AntiKt6LCTopo_sumPtTrk, &b_jet_AntiKt6LCTopo_sumPtTrk);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_NegativeE", &jet_AntiKt6LCTopo_NegativeE, &b_jet_AntiKt6LCTopo_NegativeE);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_BCH_CORR_CELL", &jet_AntiKt6LCTopo_BCH_CORR_CELL, &b_jet_AntiKt6LCTopo_BCH_CORR_CELL);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_BCH_CORR_DOTX", &jet_AntiKt6LCTopo_BCH_CORR_DOTX, &b_jet_AntiKt6LCTopo_BCH_CORR_DOTX);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_ENG_BAD_CELLS", &jet_AntiKt6LCTopo_ENG_BAD_CELLS, &b_jet_AntiKt6LCTopo_ENG_BAD_CELLS);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_SamplingMax", &jet_AntiKt6LCTopo_SamplingMax, &b_jet_AntiKt6LCTopo_SamplingMax);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_isUgly", &jet_AntiKt6LCTopo_isUgly, &b_jet_AntiKt6LCTopo_isUgly);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_isBadLoose", &jet_AntiKt6LCTopo_isBadLoose, &b_jet_AntiKt6LCTopo_isBadLoose);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_isBadMedium", &jet_AntiKt6LCTopo_isBadMedium, &b_jet_AntiKt6LCTopo_isBadMedium);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_isBadTight", &jet_AntiKt6LCTopo_isBadTight, &b_jet_AntiKt6LCTopo_isBadTight);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_emfrac", &jet_AntiKt6LCTopo_emfrac, &b_jet_AntiKt6LCTopo_emfrac);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_EMJES", &jet_AntiKt6LCTopo_EMJES, &b_jet_AntiKt6LCTopo_EMJES);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_EMJES_EtaCorr", &jet_AntiKt6LCTopo_EMJES_EtaCorr, &b_jet_AntiKt6LCTopo_EMJES_EtaCorr);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_EMJESnooffset", &jet_AntiKt6LCTopo_EMJESnooffset, &b_jet_AntiKt6LCTopo_EMJESnooffset);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_GCWJES", &jet_AntiKt6LCTopo_GCWJES, &b_jet_AntiKt6LCTopo_GCWJES);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_emscale_E", &jet_AntiKt6LCTopo_emscale_E, &b_jet_AntiKt6LCTopo_emscale_E);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_emscale_pt", &jet_AntiKt6LCTopo_emscale_pt, &b_jet_AntiKt6LCTopo_emscale_pt);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_emscale_eta", &jet_AntiKt6LCTopo_emscale_eta, &b_jet_AntiKt6LCTopo_emscale_eta);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_emscale_phi", &jet_AntiKt6LCTopo_emscale_phi, &b_jet_AntiKt6LCTopo_emscale_phi);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_jvtxf", &jet_AntiKt6LCTopo_jvtxf, &b_jet_AntiKt6LCTopo_jvtxf);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_GSCFactorF", &jet_AntiKt6LCTopo_GSCFactorF, &b_jet_AntiKt6LCTopo_GSCFactorF);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_WidthFraction", &jet_AntiKt6LCTopo_WidthFraction, &b_jet_AntiKt6LCTopo_WidthFraction);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_IP2D", &jet_AntiKt6LCTopo_flavor_weight_IP2D, &b_jet_AntiKt6LCTopo_flavor_weight_IP2D);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_IP3D", &jet_AntiKt6LCTopo_flavor_weight_IP3D, &b_jet_AntiKt6LCTopo_flavor_weight_IP3D);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_SV0", &jet_AntiKt6LCTopo_flavor_weight_SV0, &b_jet_AntiKt6LCTopo_flavor_weight_SV0);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_SV1", &jet_AntiKt6LCTopo_flavor_weight_SV1, &b_jet_AntiKt6LCTopo_flavor_weight_SV1);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_SV2", &jet_AntiKt6LCTopo_flavor_weight_SV2, &b_jet_AntiKt6LCTopo_flavor_weight_SV2);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_JetProb", &jet_AntiKt6LCTopo_flavor_weight_JetProb, &b_jet_AntiKt6LCTopo_flavor_weight_JetProb);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag", &jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag, &b_jet_AntiKt6LCTopo_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN", &jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN, &b_jet_AntiKt6LCTopo_flavor_weight_JetFitterTagNN);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN", &jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN, &b_jet_AntiKt6LCTopo_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_ip3d_pu", &jet_AntiKt6LCTopo_flavor_component_ip3d_pu, &b_jet_AntiKt6LCTopo_flavor_component_ip3d_pu);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_ip3d_pb", &jet_AntiKt6LCTopo_flavor_component_ip3d_pb, &b_jet_AntiKt6LCTopo_flavor_component_ip3d_pb);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk", &jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk, &b_jet_AntiKt6LCTopo_flavor_component_ip3d_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv1_pu", &jet_AntiKt6LCTopo_flavor_component_sv1_pu, &b_jet_AntiKt6LCTopo_flavor_component_sv1_pu);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv1_pb", &jet_AntiKt6LCTopo_flavor_component_sv1_pb, &b_jet_AntiKt6LCTopo_flavor_component_sv1_pb);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv2_pu", &jet_AntiKt6LCTopo_flavor_component_sv2_pu, &b_jet_AntiKt6LCTopo_flavor_component_sv2_pu);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv2_pb", &jet_AntiKt6LCTopo_flavor_component_sv2_pb, &b_jet_AntiKt6LCTopo_flavor_component_sv2_pb);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_pu", &jet_AntiKt6LCTopo_flavor_component_jfit_pu, &b_jet_AntiKt6LCTopo_flavor_component_jfit_pu);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_pb", &jet_AntiKt6LCTopo_flavor_component_jfit_pb, &b_jet_AntiKt6LCTopo_flavor_component_jfit_pb);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t", &jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t, &b_jet_AntiKt6LCTopo_flavor_component_jfit_nvtx1t);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_efrc", &jet_AntiKt6LCTopo_flavor_component_jfit_efrc, &b_jet_AntiKt6LCTopo_flavor_component_jfit_efrc);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_mass", &jet_AntiKt6LCTopo_flavor_component_jfit_mass, &b_jet_AntiKt6LCTopo_flavor_component_jfit_mass);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_jfit_sig3d", &jet_AntiKt6LCTopo_flavor_component_jfit_sig3d, &b_jet_AntiKt6LCTopo_flavor_component_jfit_sig3d);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_ntrkv", &jet_AntiKt6LCTopo_flavor_component_svp_ntrkv, &b_jet_AntiKt6LCTopo_flavor_component_svp_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_ntrkj", &jet_AntiKt6LCTopo_flavor_component_svp_ntrkj, &b_jet_AntiKt6LCTopo_flavor_component_svp_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_n2t", &jet_AntiKt6LCTopo_flavor_component_svp_n2t, &b_jet_AntiKt6LCTopo_flavor_component_svp_n2t);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_mass", &jet_AntiKt6LCTopo_flavor_component_svp_mass, &b_jet_AntiKt6LCTopo_flavor_component_svp_mass);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_efrc", &jet_AntiKt6LCTopo_flavor_component_svp_efrc, &b_jet_AntiKt6LCTopo_flavor_component_svp_efrc);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_svp_ntrk", &jet_AntiKt6LCTopo_flavor_component_svp_ntrk, &b_jet_AntiKt6LCTopo_flavor_component_svp_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv", &jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkv);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj", &jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrkj);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_n2t", &jet_AntiKt6LCTopo_flavor_component_sv0p_n2t, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_n2t);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_mass", &jet_AntiKt6LCTopo_flavor_component_sv0p_mass, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_mass);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_efrc", &jet_AntiKt6LCTopo_flavor_component_sv0p_efrc, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_efrc);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk", &jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk, &b_jet_AntiKt6LCTopo_flavor_component_sv0p_ntrk);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_L1_matched", &jet_AntiKt6LCTopo_L1_matched, &b_jet_AntiKt6LCTopo_L1_matched);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_L2_matched", &jet_AntiKt6LCTopo_L2_matched, &b_jet_AntiKt6LCTopo_L2_matched);
   fChain->SetBranchAddress("jet_AntiKt6LCTopo_EF_matched", &jet_AntiKt6LCTopo_EF_matched, &b_jet_AntiKt6LCTopo_EF_matched);
   fChain->SetBranchAddress("trk_n", &trk_n, &b_trk_n);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_theta", &trk_theta, &b_trk_theta);
   fChain->SetBranchAddress("trk_qoverp", &trk_qoverp, &b_trk_qoverp);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_d0_wrtPV", &trk_d0_wrtPV, &b_trk_d0_wrtPV);
   fChain->SetBranchAddress("trk_z0_wrtPV", &trk_z0_wrtPV, &b_trk_z0_wrtPV);
   fChain->SetBranchAddress("trk_phi_wrtPV", &trk_phi_wrtPV, &b_trk_phi_wrtPV);
   fChain->SetBranchAddress("trk_cov_d0_wrtPV", &trk_cov_d0_wrtPV, &b_trk_cov_d0_wrtPV);
   fChain->SetBranchAddress("trk_cov_z0_wrtPV", &trk_cov_z0_wrtPV, &b_trk_cov_z0_wrtPV);
   fChain->SetBranchAddress("trk_cov_phi_wrtPV", &trk_cov_phi_wrtPV, &b_trk_cov_phi_wrtPV);
   fChain->SetBranchAddress("trk_cov_theta_wrtPV", &trk_cov_theta_wrtPV, &b_trk_cov_theta_wrtPV);
   fChain->SetBranchAddress("trk_cov_qoverp_wrtPV", &trk_cov_qoverp_wrtPV, &b_trk_cov_qoverp_wrtPV);
   fChain->SetBranchAddress("trk_d0_wrtBS", &trk_d0_wrtBS, &b_trk_d0_wrtBS);
   fChain->SetBranchAddress("trk_z0_wrtBS", &trk_z0_wrtBS, &b_trk_z0_wrtBS);
   fChain->SetBranchAddress("trk_phi_wrtBS", &trk_phi_wrtBS, &b_trk_phi_wrtBS);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_ndof", &trk_ndof, &b_trk_ndof);
   fChain->SetBranchAddress("trk_nBLHits", &trk_nBLHits, &b_trk_nBLHits);
   fChain->SetBranchAddress("trk_nPixHits", &trk_nPixHits, &b_trk_nPixHits);
   fChain->SetBranchAddress("trk_nSCTHits", &trk_nSCTHits, &b_trk_nSCTHits);
   fChain->SetBranchAddress("trk_nTRTHits", &trk_nTRTHits, &b_trk_nTRTHits);
   fChain->SetBranchAddress("trk_nMDTHits", &trk_nMDTHits, &b_trk_nMDTHits);
   fChain->SetBranchAddress("trk_nHits", &trk_nHits, &b_trk_nHits);
   fChain->SetBranchAddress("L1_J15", &L1_J15, &b_L1_J15);
   fChain->SetBranchAddress("trig_L1_TAV", &trig_L1_TAV, &b_trig_L1_TAV);
   fChain->SetBranchAddress("trig_L2_passedPhysics", &trig_L2_passedPhysics, &b_trig_L2_passedPhysics);
   fChain->SetBranchAddress("trig_EF_passedPhysics", &trig_EF_passedPhysics, &b_trig_EF_passedPhysics);
   fChain->SetBranchAddress("trig_L1_TBP", &trig_L1_TBP, &b_trig_L1_TBP);
   fChain->SetBranchAddress("trig_L1_TAP", &trig_L1_TAP, &b_trig_L1_TAP);
   fChain->SetBranchAddress("trig_L2_passedRaw", &trig_L2_passedRaw, &b_trig_L2_passedRaw);
   fChain->SetBranchAddress("trig_EF_passedRaw", &trig_EF_passedRaw, &b_trig_EF_passedRaw);
   fChain->SetBranchAddress("trig_L2_resurrected", &trig_L2_resurrected, &b_trig_L2_resurrected);
   fChain->SetBranchAddress("trig_EF_resurrected", &trig_EF_resurrected, &b_trig_EF_resurrected);
   fChain->SetBranchAddress("trig_L2_passedThrough", &trig_L2_passedThrough, &b_trig_L2_passedThrough);
   fChain->SetBranchAddress("trig_EF_passedThrough", &trig_EF_passedThrough, &b_trig_EF_passedThrough);
   fChain->SetBranchAddress("trig_DB_SMK", &trig_DB_SMK, &b_trig_DB_SMK);
   fChain->SetBranchAddress("trig_DB_L1PSK", &trig_DB_L1PSK, &b_trig_DB_L1PSK);
   fChain->SetBranchAddress("trig_DB_HLTPSK", &trig_DB_HLTPSK, &b_trig_DB_HLTPSK);
   fChain->SetBranchAddress("trig_EF_jet_EF_j20_a4tc_EFFS", &trig_EF_jet_EF_j20_a4tc_EFFS, &b_trig_EF_jet_EF_j20_a4tc_EFFS);
   fChain->SetBranchAddress("trig_EF_jet_EF_j30_a4tc_EFFS", &trig_EF_jet_EF_j30_a4tc_EFFS, &b_trig_EF_jet_EF_j30_a4tc_EFFS);
   Notify();
}

Bool_t CollectionTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CollectionTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CollectionTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CollectionTree_cxx
