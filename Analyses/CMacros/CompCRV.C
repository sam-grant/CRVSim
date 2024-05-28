//
//  Compare track extrapolation with CRV hits in extracted mode
//  Original author: Dave Brown 
// 

#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "THStack.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TChain.h"
#include "FillChain.C"
#include "CRVCuts.C"

using namespace std;

void CompCRV(TTree* ta, const char* savesuffix="") {

    TH2D* tcrvpg = new TH2D("tcrvpg","KKInter TCRV Layer 1 Position, Has CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    TH2D* tcrvpb = new TH2D("tcrvpb","KKInter TCRV Layer 1 Position,  No CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    
    tcrvpg->SetStats(0);
    tcrvpb->SetStats(0);

    TH1D* tcrvts = new TH1D("tcrvts","KKInter Time - CRV Time, Layer 1 #Delta T;T_{KKInter} - T_{CRV} (ns)",250,-30,30);
    TH1D* tcrvy = new TH1D("tcrvy","KKInter TCRV Layer 1 #Delta Y;Y_{KKInter} - Y_{CRV} (mm)",100,-50,50);

    TH2D* dtvx = new TH2D("dtvx","KKInter TCRV Layer 1 #Delta T vs X;KKInter X (mm);T_{KKInter}-T_{CRV} (ns)",50,-3400,3400,100,-15,15);
    TH2D* dtvz = new TH2D("dtvz","KKInter TCRV Layer 1 #Delta T vs Z;KKInter Z (mm);T_{KKInter}-T_{CRV} (ns)",50,-2500,1500,100,-15,15);
    TH2D* xvx = new TH2D("xvx","CRV Layer1 X vs KKInter X;KKInter X (mm);CRV X (mm)",50,-3400,3400,50,-3400,3400);
    TH2D* zvz = new TH2D("zvz","CRV Layer1 Z vs KKInter Z;KKInter Z (mm);CRV Z (mm)",50,-2500,1500,50,-2500,1500);

    dtvx->SetStats(0);
    dtvz->SetStats(0);
    xvx->SetStats(0);
    zvz->SetStats(0);

    tcrvts->SetLineColor(kBlue);

    // I had no idea such a thing existed
    // Is it worth translating this into Python? 

    ta->Project("tcrvpg","klfit.pos.X():klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("tcrvpb","klfit.pos.X():klfit.pos.Z()",goodtrk+KLCRV1+bestfit+noCRV);
    ta->Project("tcrvy","klfit.pos.Y()-crvcoincs.pos.Y()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("tcrvts","klfit.time-crvcoincs.time",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("dtvz","klfit.time-crvcoincs.time:klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("dtvx","klfit.time-crvcoincs.time:klfit.pos.X()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("xvx","crvcoincs.pos.X():klfit.pos.X()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("zvz","crvcoincs.pos.Z():klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);

    TCanvas* tcan = new TCanvas("tcan","tcan",1200,1200);
    tcan->Divide(2,2);
    tcan->cd(1);
    tcrvts->Fit("gaus");
    tcan->cd(3);
    // gPad->SetLogz();
    zvz->Draw("colorz");
    tcan->cd(2);
    // gPad->SetLogz();
    dtvx->Draw("colorz");
    tcan->cd(4);
    // gPad->SetLogz();
    xvx->Draw("colorz");

    TCanvas* pcan = new TCanvas("pcan","pcan",1200,800);
    pcan->Divide(2,1);
    pcan->cd(1);
    // gPad->SetLogz();
    tcrvpg->Draw("colorz");
    pcan->cd(2);
    // gPad->SetLogz();
    tcrvpb->Draw("colorz");

    string ssuf(savesuffix);
    string term(".png");
    if(true) { // !ssuf.empty()){
        string tcanfile = string("tcan_") + ssuf + term;
        tcan->Draw();
        tcan->SaveAs(tcanfile.c_str());
        string pcanfile = string("pcan_") + ssuf + term;
        pcan->Draw();
        pcan->SaveAs(pcanfile.c_str());
    }
    return;
}

void CompCRVFile(const char* file,const char* ssuf="") {
    TFile* tf = new TFile(file);
    TTree* ta = (TTree*)tf->Get("TrkAnaExt/trkana");
    cout<<"\n---> Opened file "<<file<<", "<<tf<<endl;
    cout<<"---> Opened tree "<<ta<<endl;
    CompCRV(ta,ssuf);
    return;
}

void CompCRVChain(const char* files,const char* cname="TAKK/trkana",const char* ssuf="") {
    TChain* ta = new TChain(cname);
    FillChain(ta,files);
    CompCRV(ta,ssuf);
    return;
}

int main() {

    // Dataset: nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.root 
    // First file: /pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/ac/68/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000457.root
    const char* file = "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/ac/68/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000457.root";
    CompCRVFile(file);

    return 0;

}