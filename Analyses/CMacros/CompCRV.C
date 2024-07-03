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

void TCutCheck(TTree* ta) {

    // Debugging/sanity checks! 
    Int_t klstatus[3] = {0, 0, 0}; 
    Int_t crvhits[3] = {0, 0, 0}; 

    ta->SetBranchAddress("kl.status", &klstatus);
    ta->SetBranchAddress("crvcoincs.nHits", &crvhits);
    
    vector<int> hits; 

    // Loop over all entries and count the number of values in each event
    for (Long64_t i(0); ta->GetEntries(); i++) {
        ta->GetEntry(i);
        for (int i(0); i<3; i++) {
            // Evaluate the cut condition manually
            if (klstatus[i] > 0) {
                hits.push_back(crvhits[i]);
            }
        }

    }

    // Print the number of entries that pass the cut
    std::cout << "crv hits passing: " << hits.size() << std::endl;
}

void SanityPlots(TTree* ta, bool cut=false) { 

    cout<<"---> Running SanityPlots with cut (on/off) = "<<cut<<endl;

    string cutStr = "";
    if (cut) cutStr += "_cut";


    TH1D* klstatus = new TH1D("klstatus","All;Track status;Counts",2,0,2);
    if (!cut) ta->Project("klstatus","kl.status");
    else ta->Project("klstatus","kl.status", goodtrk);
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    klstatus->Draw();
    c1->SaveAs(("klstatus"+cutStr+".png").c_str());
    delete klstatus;
    delete c1;

    TH1D* crvhits = new TH1D("crvhits","All;Hits / CRV coincidence;Counts;",100,0,100);
    if (!cut) ta->Project("crvhits","crvcoincs.nHits");
    else ta->Project("crvhits","crvcoincs.nHits", goodtrk); // goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    // else ta->Project("crvhits","crvcoincs.nHits", "goodtrk && KLCRV1 && CRV1 && bestfit && goodCRV");
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    crvhits->Draw();
    c2->SaveAs(("crvhits"+cutStr+".png").c_str());
    delete crvhits;
    delete c2;

    return;

    // TH1D* klhits = new TH1D("klhits","klhits;nhits;",100,0,100);
    // if (!cut) ta->Project("klhits","kl.nhits");
    // else ta->Project("klhits","kl.nhits", goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    // TCanvas* c2 = new TCanvas("c2","c2",800,600);
    // klhits->Draw();
    // c2->SaveAs(("klhits"+cutStr+".png").c_str());
    // delete klhits;
    // delete c2;

    TH1D* klfittime = new TH1D("klfittime","klfittime;time;",100,0,5e4);
    if (!cut) ta->Project("klfittime","klfit.time");
    else ta->Project("klfittime","klfit.time", goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    klfittime->Draw();
    c3->SaveAs(("klfittime"+cutStr+".png").c_str());
    delete klfittime;
    delete c3;

    TH1D* klklz0err = new TH1D("klklz0err","klklz0err;z0err;",100,0,1);
    ta->Project("klklz0err","klkl.z0err");
    if (!cut) ta->Project("klklz0err","klkl.z0err");
    else ta->Project("klklz0err","klkl.z0err", goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    TCanvas* c4 = new TCanvas("c4","c4",800,600);
    klklz0err->Draw();
    c4->SaveAs(("klklz0err"+cutStr+".png").c_str());
    delete klklz0err;
    delete c4;


    return;
}

void CompCRV(TTree* ta, const char* savesuffix="") {

    TCutCheck(ta);

    return;

    // // Debugging/sanity checks! 
    // Int_t crvhits[3] = {0, 0, 0}; 
    // ta->SetBranchAddress("crvcoincs.nHits", &crvhits);
    // // Counters
    // int n_crvhits = 0;
    // // Loop over all entries and count the number of values in each event
    // for (Long64_t i(0); i<10; i++) { //ta->GetEntries(); i++) {
    //     ta->GetEntry(i);
    //     cout<<"["<<crvhits[0]<<", "<<crvhits[1]<<", "<<crvhits[2]<<"]"<<endl;
    //     for (int j(0); j<3; j++) {
    //         if (crvhits[j] > 0){
    //             n_crvhits++;
    //         }
    //         cout<<n_crvhits<<endl;
    //     }
    //     cout<<endl;
    // }

    // Debugging/sanity checks! 
    // // TBits klstatus; 
    // Int_t klstatus[3] = {0, 0, 0}; 
    // ta->SetBranchAddress("kl.status", &klstatus);
    // // Loop over all entries and count the number of values in each event
    // for (Long64_t i(0); i<5; i++) { //ta->GetEntries(); i++) {
    //     ta->GetEntry(i);
    //     // cout<<klstatus<<endl;
    //     cout<<"["<<klstatus[0]<<", "<<klstatus[1]<<", "<<klstatus[2]<<"]"<<endl;
    //     // for (int j(0); j<3; j++) {
    //     //     if (crvhits[j] > 0){
    //     //         n_crvhits++;
    //     //     }
    //     //     cout<<n_crvhits<<endl;
    //     // }
    //     cout<<endl;
    // }

    // return;
    // cout<<"n_crvhits = "<<n_crvhits<<endl;

    // delete ta;

    // return;

    // sanity checks! 

    // TH1D* crvhits = new TH1D("crvhits","crvhits;nhits;",100,0,100);
    // ta->Project("crvhits","crvcoincs.nHits");
    // TCanvas* c1 = new TCanvas("c1","c1",800,600);
    // crvhits->Draw();
    // c1->SaveAs("crvhits.png");

    // TH1D* klhits = new TH1D("klhits","klhits;nhits;",100,0,100);
    // ta->Project("klhits","kl.nhits");
    // TCanvas* c2 = new TCanvas("c2","c2",800,600);
    // klhits->Draw();
    // c2->SaveAs("klhits.png");

    // TH1D* klfittime = new TH1D("klfittime","klfittime;time;",100,0,5e4);
    // ta->Project("klfittime","klfit.time");
    // TCanvas* c3 = new TCanvas("c3","c3",800,600);
    // klfittime->Draw();
    // c3->SaveAs("klfittime.png");

    // TH1D* klklz0err = new TH1D("klklz0err","klklz0err;z0err;",100,0,1);
    // ta->Project("klklz0err","klkl.z0err");
    // TCanvas* c4 = new TCanvas("c4","c4",800,600);
    // klklz0err->Draw();
    // c4->SaveAs("klklz0err.png");

    // // Now with cuts
    // ta->Project("tcrvy","klfit.pos.Y()-crvcoincs.pos.Y()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);

    SanityPlots(ta, false);
    // SanityPlots(ta, true);

    return;

    //TH2D* tcrvpb = new TH2D("tcrvpb","KKInter TCRV Layer 1 Position,  No CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    

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

    ta->Project("tcrvpg","klfit.pos.X():klfit.pos.Z()");//,goodCRV); // +KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("tcrvpb","klfit.pos.X():klfit.pos.Z()");//,noCRV); //+KLCRV1+bestfit+noCRV);

    std::cout<<"\n---> tcrvpg->GetEntries()\t"<<tcrvpg->GetEntries()<<std::endl;
    std::cout<<"---> tcrvpb->GetEntries()\t"<<tcrvpb->GetEntries()<<std::endl;

    return;

    ta->Project("tcrvy","klfit.pos.Y()-crvcoincs.pos.Y()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("tcrvts","klfit.time-crvcoincs.time",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("dtvz","klfit.time-crvcoincs.time:klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("dtvx","klfit.time-crvcoincs.time:klfit.pos.X()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("xvx","crvcoincs.pos.X():klfit.pos.X()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    ta->Project("zvz","crvcoincs.pos.Z():klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);

    

    std::cout<<"\n---> tcrvpg->GetEntries()\t"<<tcrvpg->GetEntries()<<std::endl;
    std::cout<<"---> tcrvpb->GetEntries()\t"<<tcrvpb->GetEntries()<<std::endl;

    std::cout<<"\n---> tcrvts->GetEntries()\t"<<tcrvts->GetEntries()<<std::endl;
    std::cout<<"---> dtvx->GetEntries()\t"<<dtvx->GetEntries()<<std::endl;

    std::cout<<"\n---> zvz->GetEntries()\t"<<zvz->GetEntries()<<std::endl;
    std::cout<<"---> xvx->GetEntries()\t"<<xvx->GetEntries()<<std::endl;
    
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
    tf->Close();
    return;
}

void CompCRVChain(const char* files,const char* cname="TAKK/trkana",const char* ssuf="") {
// void CompCRVChain(const char* files,const char* cname="TrkAnaExt/trkana",const char* ssuf="") {
    TChain* ta = new TChain(cname);
    FillChain(ta,files);
    CompCRV(ta,ssuf);
    delete ta;
    return;
}

int main() {

    // Dataset: nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.root 
    // First file: /pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/ac/68/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000457.root
    const char* file = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00033/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000080.root"; //  #  "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/ac/68/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000457.root";
    CompCRVFile(file);
    // delete file;

    return 0;

}
