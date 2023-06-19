#include "FancyDraw.h"

/*void DrawTwoTH1(TH1* h1, TH1* h2, std::string title, std::string fname) {

	TCanvas *c = new TCanvas("c","c",800,600);

	h1->SetTitle(title.c_str());

	//h1->SetStats(0);
	gStyle->SetOptStat(0);//2210);
			
	h1->GetXaxis()->SetTitleSize(.04);
	h1->GetYaxis()->SetTitleSize(.04);
	h1->GetXaxis()->SetTitleOffset(1.1);
	h1->GetYaxis()->SetTitleOffset(1.1);
	h1->GetXaxis()->CenterTitle(1);
	h1->GetYaxis()->CenterTitle(1);
	h1->GetYaxis()->SetMaxDigits(4);
	h1->SetLineWidth(3); h2->SetLineWidth(3);
	h1->SetLineColor(kBlack);
	h1->SetLineColor(kRed);

	//c->SetRightMargin(0.13);
	h1->Draw("HIST");
	h2->Draw("HIST SAME");

	// TLegend *l = new TLegend(0.70,0.65,0.89,0.85);
	TLegend *l = new TLegend(0.69,0.69,0.89,0.89);
	h1->SetName(h1->GetName());//name1.c_str());
	h2->SetName(h2->GetName());//name2.c_str());
	gPad->Update();
	l->SetBorderSize(0);
	l->SetTextSize(26);
	l->SetTextFont(44);
	l->AddEntry(h1, h1->GetName());
	l->AddEntry(h2, h2->GetName());
	l->Draw("SAME");

	c->SaveAs((fname+".pdf").c_str());
	c->SaveAs((fname+".png").c_str());

	delete c;

	return;
}
*/
void TrkCRVPlots() {

	TString finName = "../../../plots/trkana-reco.test.root";	
	TFile *fin = TFile::Open(finName);

	cout<<"---> Got input file "<<finName<<", "<<fin<<endl;

	// TH1Fs
	TH1F *h1_trkTime = (TH1F*)fin->Get("trkCRVPlots/h1_trkTime"); 
	TH1F *h1_crvTime = (TH1F*)fin->Get("trkCRVPlots/h1_crvTime"); 
  	TH1F *h1_trkX = (TH1F*)fin->Get("trkCRVPlots/h1_trkX");
  	TH1F *h1_crvX = (TH1F*)fin->Get("trkCRVPlots/h1_crvX");
   	TH1F *h1_trkY = (TH1F*)fin->Get("trkCRVPlots/h1_trkY");
  	TH1F *h1_crvY = (TH1F*)fin->Get("trkCRVPlots/h1_crvY");
  	TH1F *h1_trkZ = (TH1F*)fin->Get("trkCRVPlots/h1_trkZ");
  	TH1F *h1_crvZ = (TH1F*)fin->Get("trkCRVPlots/h1_crvZ");
  	TH1F *h1_deltaTime = (TH1F*)fin->Get("trkCRVPlots/h1_deltaTime");
   	TH1F *h1_deltaX = (TH1F*)fin->Get("trkCRVPlots/h1_deltaX"); 
    TH1F *h1_deltaY = (TH1F*)fin->Get("trkCRVPlots/h1_deltaY");
  	TH1F *h1_deltaZ = (TH1F*)fin->Get("trkCRVPlots/h1_deltaZ"); 

/*	DrawTH1(h1_trkTime, "", "../../../images/h1_trkTime");
	DrawTH1(h1_crvTime, "", "../../../images/h1_crvTime");
	DrawTH1(h1_trkX, "", "../../../images/h1_trkX");
	DrawTH1(h1_crvX, "", "../../../images/h1_crvX");
	DrawTH1(h1_trkY, "", "../../../images/h1_trkY");
	DrawTH1(h1_crvY, "", "../../../images/h1_crvY");
	DrawTH1(h1_trkZ, "", "../../../images/h1_trkZ");
	DrawTH1(h1_crvZ, "", "../../../images/h1_crvZ");
	DrawTH1(h1_deltaTime, "", "../../../images/h1_deltaTime");
	DrawTH1(h1_deltaX, "", "../../../images/h1_deltaX");
	DrawTH1(h1_deltaY, "", "../../../images/h1_deltaY");
	DrawTH1(h1_deltaZ, "", "../../../images/h1_deltaZ");*/

  	// TH2Fs

 	TH2F *h2_trkYX = (TH2F*)fin->Get("trkCRVPlots/h2_trkYX");
  	TH2F *h2_crvYX = (TH2F*)fin->Get("trkCRVPlots/h2_crvYX");
	TH2F *h2_deltaYX = (TH2F*)fin->Get("trkCRVPlots/h2_deltaYX");

	TH2F *h2_trkZX = (TH2F*)fin->Get("trkCRVPlots/h2_trkZX");
  	TH2F *h2_crvZX = (TH2F*)fin->Get("trkCRVPlots/h2_crvZX");
	TH2F *h2_deltaZX = (TH2F*)fin->Get("trkCRVPlots/h2_deltaZX");

/*	DrawTH2(h2_trkYX, "", "../../../images/h2_trkYX");
	DrawTH2(h2_crvYX, "", "../../../images/h2_crvYX");
	DrawTH2(h2_deltaYX, "", "../../../images/h2_deltaYX");

	DrawTH2(h2_trkZX, "", "../../../images/h2_trkZX");
	DrawTH2(h2_crvZX, "", "../../../images/h2_crvZX");
	DrawTH2(h2_deltaZX, "", "../../../images/h2_deltaZX");*/

	// Overlays

/*	DrawTwoTH1(h1_trkX, h1_crvX, "", "../../../images/h2_trkCRVX_overlay");
	DrawTwoTH1(h1_trkY, h1_crvY, "", "../../../images/h2_trkCRVY_overlay");
	DrawTwoTH1(h1_trkZ, h1_crvZ, "", "../../../images/h2_trkCRVZ_overlay");
*/
	fin->Close();

	return;
}