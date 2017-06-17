#include <memory>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <TH2.h>
#include <THStack.h>
#include <TCanvas.h>
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TText.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <utility> 
#include <TLorentzVector.h>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLine.h"
#include "TEfficiency.h"
#include "TBinomialEfficiencyFitter.h"
#include "TList.h"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"
#include "TGraph2D.h"
#include "TMatrixD.h"

#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>

Double_t myfunc(double * x , double * par){
	double a = (1.-x[0])*(1.-x[0]);
	double b = (1.-x[0]*x[0]);
	double c = (1.+x[0])*(1.+x[0]);  
	return par[2]*(3./8.*a*par[0]+3./4.*b*par[1]+3./8.*c*(1.-par[0]-par[1])); 
}
Double_t linearFit(double * x , double * par){
	return (par[0]*x[0]+par[1]);
}

void makeCard(double N, double S, double dS, double B, double dB, string sOut) {

    // === DATA CARD ===
    ofstream fOut(sOut.c_str());
    fOut.precision(3);
    fOut << "imax 1  number of channels" << std::endl;
    fOut << "jmax 1  number of backgrounds" << std::endl;
    fOut << "kmax 2  number of nuisance parameters (sources of systematic uncertainties)" << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin b1" << std::endl;
    fOut << "observation " << N << std::endl;
    fOut << "---" << std::endl;
    fOut << "bin              b1     b1" << std::endl;
    fOut << "process         SMS    All" << std::endl;
    fOut << "process          0     1  " << std::endl;
    fOut << "rate           " << S << "\t" << B << std::endl;
    fOut << "---" << std::endl;
    fOut << "dS  lnN    " << 1 + dS << "\t-" << std::endl;
    fOut << "dB  lnN    - \t " << 1 + dB << std::endl;
    fOut.close();

}

TGraphAsymmErrors* plotSig(TH1 *hSgn, TH1 *hBkg, TString xtitle, TString cutType, int type, double sys) {

    int nbins = hSgn->GetXaxis()->GetNbins();
    float *x = new float[nbins + 1];
    float *ex = new float[nbins + 1];
    float *y = new float[nbins + 1];
    float *ey = new float[nbins + 1];
    float *eyp = new float[nbins + 1];
    float *eym = new float[nbins + 1];

    for (int i = 1; i <= nbins+1; i++) {

        x[i - 1] = hSgn->GetBinLowEdge(i);
        ex[i - 1] = hSgn->GetBinWidth(i)/2.0;
	
        double s = (cutType == "Lower_Cut") ? hSgn->Integral(i, nbins + 1) : hSgn->Integral(0, i-1);
        double ds = sqrt(s) + s * sys;
        double b = (cutType == "Lower_Cut") ? hBkg->Integral(i, nbins + 1) : hBkg->Integral(0, i-1);
        double db = sqrt(b) + b * sys;
	double S_sum = hSgn->Integral();

	if (b < 0) b = 0; // added this to avoid nan, on Sep 21 

        if (b == 0 || s == 0) {
            y[i - 1] = .0;
            ey[i - 1] = .0;
        } else {
            if (type == 0) {
                y[i - 1] = s / sqrt(b);
                ey[i - 1] = y[i - 1] * fabs(ds / s - db / (2 * b)); // changed the + with -
            }
            if (type == 1) {
                y[i - 1] = s / sqrt(s + b);
                ey[i - 1] = y[i - 1] * fabs(ds / s - (db + ds) / (2 * (b + s))); // changed the + with -
            }
            /*if (type == 2) {
                y[i - 1] = s / b;
                ey[i - 1] = y[i - 1] * (ds / s + db / b);
            }*/
            /*if (type == 2) {
		double Eff = s / S_sum;
		double FoM = 2.1 * sqrt(b + 2.04) / Eff;

                y[i - 1] = FoM;
		double dEff = sqrt(Eff * (1 - Eff) / S_sum);
                ey[i - 1] = y[i - 1] * (dEff / Eff + db / (2 * b + 4.08));
	    }*/
            if (type == 2) {
                y[i - 1] = sqrt(2 * ( (s + b) * TMath::Log(1 + s / b) - s ));
                ey[i - 1] = fabs((ds + db) * TMath::Log(1 + s / b) / y[i - 1] + (ds - db * s / b) / y[i - 1] - ds / y[i - 1]);
            }
	    if (type == 3) {
                cout<<cutType<<" bin "<<i<<endl;
                makeCard(b, s, sys, b, sys, "datacard");
                if (!(std::ifstream("datacard")).good()) continue;
                system("combine -M Asymptotic datacard");
                TTree* tree;

                TFile * flimit = new TFile("higgsCombineTest.Asymptotic.mH120.root");

                flimit->GetObject("limit", tree);

                Double_t limit;
                TBranch *b_limit; //!
                tree->SetBranchAddress("limit", &limit, &b_limit);

                Float_t quantileExpected;
                TBranch *b_quantileExpected; //!
                tree->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);

                std::vector<double> vLimit;
                Long64_t nEntrs = tree->GetEntriesFast();
                for (Long64_t iEntr = 0; iEntr < nEntrs; iEntr++) {
                    tree->GetEntry(iEntr);
                    cout << ">> quantileExpected: " << quantileExpected << "\tlimit: " << limit << endl;
                    vLimit.push_back(limit);
                }

                double SgmP2(vLimit[0]), SgmP1(vLimit[1]), Mdn(vLimit[2]), SgmM1(vLimit[3]), SgmM2(vLimit[4]), Obs(vLimit[5]);

                y[i - 1] = Mdn;
                eyp[i - 1] = SgmM1 - y[i - 1];
                eym[i - 1] = y[i - 1] - SgmP1;

                system("rm -f higgsCombineTest.Asymptotic.mH120.root");
                system("rm -f datacard");
                system("rm -f roostat*");
	  }
        }
    }

    TGraphAsymmErrors *sig = new TGraphAsymmErrors(nbins +1 , x, y, ex, ex, ey, ey);
    if (type == 3) sig = new TGraphAsymmErrors(nbins+1 , x, y, ex, ex, eym, eyp);

    sig->SetTitle("");
    sig->GetXaxis()->SetTitle(xtitle + "_" + cutType);
    sig->SetMarkerStyle(20);
    sig->SetFillColor(kBlue-7); 
    sig->SetFillStyle(3005);
////
    sig->GetXaxis()->SetLabelFont(42);
    sig->GetXaxis()->SetLabelOffset(0.007);
    sig->GetXaxis()->SetLabelSize(0.045);
    sig->GetXaxis()->SetTitleSize(0.06);
    sig->GetXaxis()->SetTitleOffset(0.95);
    sig->GetXaxis()->SetTitleFont(42);

    sig->GetYaxis()->SetLabelFont(42);
    sig->GetYaxis()->SetLabelSize(0.045);
    sig->GetYaxis()->SetTitleSize(0.06);
    sig->GetYaxis()->SetTitleOffset(1.45);
    sig->GetYaxis()->SetTitleFont(42);
///
    if (type == 0) sig->GetYaxis()->SetTitle("S/#sqrt{B}");
    if (type == 1) sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
    //if (type == 2) sig->GetYaxis()->SetTitle("S/B");
    //if (type == 2) sig->GetYaxis()->SetTitle("FoM");
    if (type == 2) sig->GetYaxis()->SetTitle("AMS");
    if (type == 3) sig->GetYaxis()->SetTitle("signal strength (r)");

    return sig;
}

TCanvas * canvasCreator(TString name, TGraphAsymmErrors * gr1, TGraphAsymmErrors * gr2, TGraphAsymmErrors * gr3, TGraphAsymmErrors * gr4) {
   TCanvas *c1 = new TCanvas("c1", "MT2Optimization_Upper_Cut",53,26,900,700);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: c1_1
   TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.51,0.49,0.99);
   c1_1->Draw();
   c1_1->cd();
   c1_1->Range(-50.0852,-0.01934188,152.869,0.1442097);
   c1_1->SetFillColor(0);
   c1_1->SetBorderMode(0);
   c1_1->SetBorderSize(2);
   c1_1->SetLeftMargin(0.1605283);
   c1_1->SetBottomMargin(0.1213526);
   c1_1->SetFrameBorderMode(0);
   c1_1->SetGridx();
   c1_1->SetGridy();
   
   gr1->Draw("alp");
   c1_1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_2
   TPad *c1_2 = new TPad("c1_2", "c1_2",0.51,0.51,0.99,0.99);
   c1_2->Draw();
   c1_2->cd();
   c1_2->Range(-50.0852,-0.01934188,152.869,0.1442097);
   c1_2->SetFillColor(0);
   c1_2->SetBorderMode(0);
   c1_2->SetBorderSize(2);
   c1_2->SetLeftMargin(0.1605283);
   c1_2->SetBottomMargin(0.1213526);
   c1_2->SetGridx();
   c1_2->SetGridy();
   
   gr2->Draw("alp");
   c1_2->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_3
   TPad *c1_3 = new TPad("c1_3", "c1_3",0.01,0.01,0.49,0.49);
   c1_3->Draw();
   c1_3->cd();
   //c1_3->Range(-48.92522,-0.0002109698,152.6758,0.001527516);
   c1_3->Range(-50.0852,-0.01934188,152.869,0.1442097);
   c1_3->SetFillColor(0);
   c1_3->SetBorderMode(0);
   c1_3->SetBorderSize(2);
   c1_3->SetLeftMargin(0.155878);
   c1_3->SetBottomMargin(0.1213526);
   c1_3->SetFrameBorderMode(0);
   c1_3->SetGridx();
   c1_3->SetGridy();
   
   gr3->Draw("alp");
   c1_3->Modified();
   c1->cd();
  
// ------------>Primitives in pad: c1_4
   TPad *c1_4 = new TPad("c1_4", "c1_4",0.51,0.01,0.99,0.49);
   c1_4->Draw();
   c1_4->cd();
   //c1_4->Range(-50.06627,-0.0002109698,152.8027,0.001527516);
   c1_4->Range(-50.0852,-0.01934188,152.869,0.1442097);
   c1_4->SetFillColor(0);
   c1_4->SetBorderMode(0);
   c1_4->SetBorderSize(2);
   c1_4->SetLeftMargin(0.1605283);
   c1_4->SetBottomMargin(0.1213526);
   c1_4->SetFrameBorderMode(0);
   c1_4->SetGridx();
   c1_4->SetGridy();
   c1_4->SetLogy();
   
   gr4->Draw("alp");
   c1_4->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);

   c1->SaveAs(name+".png");
   return c1;
}

TCanvas * canvasCreator(TString name, TGraphAsymmErrors * graphs[3][4]) {
   TCanvas *c1 = new TCanvas("c1", "optimization",53,26,900,700);

   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
   TMultiGraph *mg[4];
   TString names[3] = {"SMS_DM_400","SMS_DM_250","SMS_DM_200"};	

   TPad * c1_[4];
   c1_[0] = new TPad("c1_1", "c1_1",0.01,0.51,0.49,0.99);
   c1_[1] = new TPad("c1_2", "c1_2",0.51,0.51,0.99,0.99);
   c1_[2] = new TPad("c1_3", "c1_3",0.01,0.01,0.49,0.49);
   c1_[3] = new TPad("c1_4", "c1_4",0.51,0.01,0.99,0.49);

	for (unsigned int i = 0; i < 4; i++) {

		c1_[i]->Draw();
                c1_[i]->cd();
                c1_[i]->Range(-50.0852,-0.01934188,152.869,0.1442097);
                c1_[i]->SetFillColor(0);
                c1_[i]->SetBorderMode(0);
                c1_[i]->SetBorderSize(2);
                c1_[i]->SetLeftMargin(0.1605283);
                c1_[i]->SetBottomMargin(0.1213526);
                c1_[i]->SetFrameBorderMode(0);
                c1_[i]->SetGridx();
                c1_[i]->SetGridy();

		mg[i] = new TMultiGraph();

		for (unsigned int s = 0; s < 3; s++) {
			graphs[s][i]->SetFillStyle(3005);
			graphs[s][i]->SetMarkerStyle(s+20);
			graphs[s][i]->SetMarkerColor(s+2);
			graphs[s][i]->SetLineColor(s+2);
			graphs[s][i]->SetTitle(names[s]);
			mg[i]->Add(graphs[s][i]);
	
		}
		c1_[i]->cd();
		mg[i]->Draw("apl");
		c1_[i]->Update();
		mg[i]->GetXaxis()->SetTitle(graphs[0][i]->GetXaxis()->GetTitle());
		mg[i]->GetYaxis()->SetTitle(graphs[0][i]->GetYaxis()->GetTitle());

                mg[i]->GetXaxis()->SetLabelFont(42);
                mg[i]->GetXaxis()->SetLabelOffset(0.007);
                mg[i]->GetXaxis()->SetLabelSize(0.045);
                mg[i]->GetXaxis()->SetTitleSize(0.06);
                mg[i]->GetXaxis()->SetTitleOffset(0.95);
                mg[i]->GetXaxis()->SetTitleFont(42);

                mg[i]->GetYaxis()->SetLabelFont(42);
                mg[i]->GetYaxis()->SetLabelSize(0.045);
                mg[i]->GetYaxis()->SetTitleSize(0.06);
                mg[i]->GetYaxis()->SetTitleOffset(1.45);
                mg[i]->GetYaxis()->SetTitleFont(42);

		c1_[i]->Modified();
		if (i == 3) c1_[i]->BuildLegend();
		c1_[i]->Modified();
		c1->cd();

	}
   
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);

   c1->SaveAs(name+".png");
   return c1;
}

void Macro() {

	TH1::SetDefaultSumw2();

	TFile * f =  new TFile("pfmet_sumMT.root");

	TH1F * Bkg = (TH1F *)f->Get("h_Bkg");
	TH1F * SMS_DM_400 = (TH1F *)f->Get("h_SMS_TChipmStauSnuMStau400");
	TH1F * SMS_DM_250 = (TH1F *)f->Get("h_SMS_TChipmStauSnuMStau250");
	TH1F * SMS_DM_200 = (TH1F *)f->Get("h_SMS_TChipmStauSnuMStau200");

	std::vector<TH1F*> SMS_histos;
	SMS_histos.push_back(SMS_DM_400);
	SMS_histos.push_back(SMS_DM_250);
	SMS_histos.push_back(SMS_DM_200);

/*
	Bkg->Rebin(5);
	SMS_DM_400->Rebin(5);
	SMS_DM_250->Rebin(5);
	SMS_DM_200->Rebin(5);

	TH1F * BKg = (TH1F *)f->Get("_h_BKg_80_");
	TH1F * SMS_DM_ = (TH1F *)f->Get("_h_SMS_DM__80_");

	// when running on minDphi, it was initially 32 bin, so no sence if we want to rebin the histogram via the following commands

	BKg->Rebin(10);
	SMS_DM_->Rebin(10);

	cout<<"number of bins of all SM: "<<BKg->GetNbinsX()<<endl;
	cout<<"number of bins of SMS_DM_nal: "<<SMS_DM_->GetNbinsX()<<endl;
*/

	TGraphAsymmErrors* myGraphs[3][4];

	//TString cutTypes[2] = {"Lower_Cut","Upper_Cut"}; 
	TString cutType = "Lower_Cut"; 

	//for (unsigned int j = 0; j < 2; j++) {

		//if (cutTypes[j] == "Upper_Cut") continue;

		TString var = "sumMT";
		TString name = "Optimization";
		//s.Contains("SMS_DM_400")
		//name = var + name + "_" + cutTypes[j];
		name = var + name + "_" + cutType;
		

		for (int s = 0; s < 3; s++) {
			for (int i = 0; i < 4; i++) {
					myGraphs[s][i] = plotSig(SMS_histos[s], Bkg, var, cutType, i, 0.1);}
		}
		//canvasCreator(name, myGraphs[0], myGraphs[1], myGraphs[2], myGraphs[3]);
		canvasCreator(name, myGraphs);

	//}
}
