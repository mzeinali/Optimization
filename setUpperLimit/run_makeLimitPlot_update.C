#include <TString.h>
#include <TGraph.h>
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
#include "treeBase.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <TLorentzVector.h>
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLine.h"
#include <stdlib.h>
#include <stdio.h>


TLegend* legend() {
    TLegend *leg = new TLegend(0.2, 0.5, 0.4, 0.8, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(102);
    return leg;
}

TH2D* getXsecUp(TH2D* hRatio, TH2D* hXsec) {
    TH2D* hXsecUp = (TH2D*) hRatio->Clone("hXsecUp");
    
    hXsecUp->Multiply(hRatio, hXsec);

    //    hXsecUp->GetXaxis()->SetRangeUser(0, 600);
    //    hXsecUp->GetYaxis()->SetRangeUser(0, 600);
    hXsecUp->GetXaxis()->SetRangeUser(100, 700);
    hXsecUp->GetYaxis()->SetRangeUser(0, 700);
    hXsecUp->GetZaxis()->SetRangeUser(0.001, 100);

    hXsecUp->GetXaxis()->SetTitle("m_{#tilde{#chi}} [GeV/c]");
    hXsecUp->GetYaxis()->SetTitle("m_{LSP} [GeV/c]");
    hXsecUp->GetZaxis()->SetTitle("#sigma upper limit [pb]");

    return hXsecUp;
}

void fixTH2D(TH2D* h, int lw, int lc, int ls) {
    h->SetLineWidth(lw);
    h->SetLineColor(lc);
    h->SetLineStyle(ls);

    int nbnx = h->GetNbinsX();
    int nbny = h->GetNbinsY();
    for (int i = 1; i <= nbnx; i++) {
        for (int j = 1; j <= nbny; j++) {
	//  int glob = h->FindBin(i,j);
	  
	//  double z = h->GetBinContent(glob);
	//  if (z == 0) h->SetBinContent(glob, 10);
            double z = h->GetBinContent(i, j);
           if (z == 0) h->SetBinContent(i, j, 10);
        }
    }
            h->SetMinimum(0.01);
            h->SetMaximum(100);
            h->SetEntries(1980);
            h->SetContour(1);
            h->SetContourLevel(0, 1);

}

TGraph* accRange() {

    // ===
    TGraph *g00 = new TGraph(4);
    g00->SetFillColor(kWhite);
    g00->SetLineColor(kBlack);
    g00->SetLineWidth(4);
    g00->SetPoint(0, 200, 0);
    g00->SetPoint(1, 120, 80);
    g00->SetPoint(2, 120, 100);
    g00->SetPoint(3, 700, 680);

    return g00;
}

TGraph* atlasGraph() {

    // ===
    TGraph *g00 = new TGraph(7);
    g00->SetFillColor(kWhite);
    g00->SetLineColor(kRed);
    g00->SetLineWidth(3);
    g00->SetPoint(0, 345, 0);
    g00->SetPoint(1, 350, 25);
    g00->SetPoint(2, 335, 50);
    g00->SetPoint(3, 270, 70);
    g00->SetPoint(4, 225, 60);
    g00->SetPoint(5, 150, 25);
    g00->SetPoint(6, 120, 0);

    return g00;
}

void text(float x, float y, float a, float w, TString s) {
    TLatex Tl;
//    Tl.SetTextAlign(0);
    Tl.SetTextSize(w);
    Tl.SetTextAngle(a);
    Tl.DrawLatex(x, y, s);
}


struct vGraph {
    TGraph** g;
    Int_t n;

    void draw() {
        for (int i = 0; i < n; i++)
            g[i]->Draw("same");
            //g[i]->Draw("");
    }
};

void Contour2Graph(TH2D* h, vGraph &vg) {
    h->Draw("CONT Z LIST");
    gPad->Update();
    TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

    Int_t TotalConts = 0;
    if (!(conts == NULL)) TotalConts = conts->GetSize();

    for (int i = 0; i < TotalConts; i++) {
        TList* contLevel = (TList*) conts->At(i);
        vg.n = contLevel->GetSize();

        vg.g = new TGraph*[vg.n];
        for (int j = 0; j < vg.n; j++) {

            TGraph* curv = (TGraph*) contLevel->At(j);
            Int_t np = curv->GetN();
            vg.g[j] = new TGraph(np);
            Double_t x0, y0;
            for (int k = 0; k < np; k++) {
                curv->GetPoint(k, x0, y0);
                vg.g[j]->SetPoint(k, x0 - 10, y0 - 10);
            }// k over np
            vg.g[j]->SetLineColor(h->GetLineColor());
            vg.g[j]->SetLineStyle(h->GetLineStyle());
            vg.g[j]->SetLineWidth(h->GetLineWidth());
        }// j over ngc
    }// i over TotalConts
}

void run_makeLimitPlot_update() {
    
    gStyle->SetOptStat(0);	

    std::vector< TString > sin;
    int argc = gApplication->Argc();
    for (int i = 0; i < argc; i++) {
        TString argvi = gApplication->Argv(i);
        if (argvi == "run_makeLimitPlot_update.C")
            for (int j = i + 1; j < argc; j++) {
                TString argvj = gApplication->Argv(j);
                sin.push_back(argvj);
                i++;
            }
    }


//     TH2D* hXsecUp = (TH2D*) TFile::Open("/dataLOCAL/MT2Tau/share/muTauNewOptimizationsPreSelectionMinMetJDPhi40/Bin1_3_MT2g90_tauMTg200.root")->Get("h_PN_MLSP_MChi");
//     hXsecUp->Add( (TH2D*) TFile::Open("/dataLOCAL/MT2Tau/share/eleTauNewOptimizationsPreSelectionMinMetJDPhi40/Bin1_3.root")->Get("h_PN_MLSP_MChi") );
//     hXsecUp->Add( (TH2D*) TFile::Open("/dataLOCAL/MT2Tau/share/TauTau/counting_TBTTalk_23Jan/counting_firstBin/maxMT_gt_200/countingForExclusion_TBD_Histos.root")->Get("h_PN_MLSP_MChi") );
//     hXsecUp->Add( (TH2D*) TFile::Open("/dataLOCAL/MT2Tau/share/TauTau/counting_TBTTalk_23Jan/counting_secondBin/countingForExclusion_TBD_Histos.root")->Get("h_PN_MLSP_MChi") );


//     for (int iy = 0; iy < hXsecUp->GetNbinsY(); iy++) {
//       for (int ix = 0; ix < hXsecUp->GetNbinsX(); ix++) {
// 	double content = hXsecUp->GetBinContent(ix, iy);
// 	double error   = hXsecUp->GetBinError(ix, iy);
	
// 	if(content != 0)
// 	  hXsecUp->SetBinContent(ix, iy, error/content);
//       }
//     }
    
 //    hXsecUp->GetXaxis()->SetRangeUser(100, 500);
//     hXsecUp->GetYaxis()->SetRangeUser(0, 400);
//     hXsecUp->GetZaxis()->SetRangeUser(0.1, 0.6);


    TString sgnStrFile = sin[0];
    TString sgnStrName = sin[1];

    TH2D* hXsec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("C1C1_13TeV_NLONLL_LSP");

//     TH2D* hXsec = (TH2D*) TFile::Open("referenceXSecs.root")->Get("StSt_8TeV_NLONLL_LSP"); hXsec->Rebin2D(2,2); hXsec->Scale(0.25);

    TH2D* h00 = (TH2D*) TFile::Open(sgnStrFile)->Get("hMdn");
    TH2D* hOb = (TH2D*) TFile::Open(sgnStrFile)->Get("hMdn");
    TH2D* h01 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmP1");
    TH2D* h02 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmP2");
    TH2D* h_1 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmM1");
    TH2D* h_2 = (TH2D*) TFile::Open(sgnStrFile)->Get("hSgmM2");

    
    //For HybridNew hSgmP2 is the median
 //     h00 = h02 ;
//      h02 = h01 ;
  
    //For HybridNew hSgmP2 is the hSgmP1
//     h01 = h02 ;
//     h02 = h00 ;
    float a = 1.0;//0.2;
    h00->Scale(a);
    hOb->Scale(a);
    h01->Scale(a);
    h02->Scale(a);
    h_1->Scale(a);
    h_2->Scale(a);

    
    TH2D* hXsecUp = getXsecUp(h00, hXsec);

    h00->Rebin2D(25, 25);	
    hOb->Rebin2D(25, 25);	
    h01->Rebin2D(25, 25);	
    h02->Rebin2D(25, 25);	
    h_1->Rebin2D(25, 25);	
    h_2->Rebin2D(25, 25);	

    fixTH2D(h00, 3, kBlue, 1);
    fixTH2D(hOb, 2, kBlack, 9);
    //    h00->SetFillStyle(3007);
    //    h00->SetFillColor(kBlue - 7);
    fixTH2D(h01, 2, kBlue, 2);
    fixTH2D(h02, 1, kBlue, 7);
    fixTH2D(h_1, 2, kBlue, 2);
    fixTH2D(h_2, 1, kBlue, 7);

    vGraph vg00, vgOb, vg01, vg02, vg_1, vg_2;
    Contour2Graph(h00, vg00);
    Contour2Graph(hOb, vgOb);
    Contour2Graph(h01, vg01);
    Contour2Graph(h02, vg02);
    Contour2Graph(h_1, vg_1);
    Contour2Graph(h_2, vg_2);
    
    TGraph* g00 = accRange();
    TGraph* g01 = atlasGraph();
    
    // ===  ===
    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    c1->SetGrid();
    c1->cd();

    gPad->SetLogz();
    gPad->SetRightMargin(0.2);

    hXsecUp->Rebin2D(25, 25);	
    hXsecUp->Draw("colz");

//    h00->Draw("cont3same");
//    h01->Draw("cont3same");
//    h02->Draw("cont3same");
    vg00.draw();
    vgOb.draw();
    vg01.draw();
    vg02.draw();
    vg_1.draw();
    vg_2.draw();

////    g00->Draw();
////    g01->Draw("C");
    //text(110,70,-43, 0.02, "m(#tilde{#tau}) < 95 GeV"); // commented on 29 May

    TLegend* leg = legend();
    leg->AddEntry(h00, sgnStrName, "lpf");
    leg->AddEntry(h01, (sgnStrName+" #pm 1#sigma"), "lpf");
//      leg->AddEntry(h02, (sgnStrName+" + 2#sigma"), "lpf");
    //leg->AddEntry(hOb, "Observed", "lpf"); // commented on 29 May
    //leg->AddEntry(g01, ("ATLAS Expected"), "lpf");
    //leg->AddEntry(g00, "Accepted Range", "lpf");
    leg->Draw();

}
