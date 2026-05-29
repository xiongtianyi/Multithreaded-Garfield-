//------------------------------------------------------------------//
// Author: QXHusky
// Final editable version (paper-style Z distribution plot)
//------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

#include <TChain.h>
#include <TLegend.h>
#include <THStack.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

// THGEM 尺寸 [cm]
const double pitch  = 0.1;
const double metal  = 0.0012;
const double thick  = 0.08;
const double induce = 0.2;
const double driftz = 0.5;

//------------------------------------------------------------
// Global plot style: cleaner, more paper-like
//------------------------------------------------------------
void SetPaperStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTextFont(42);

    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.055, "XYZ");

    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.20, "Y");

    gStyle->SetNdivisions(506, "XY");
    gStyle->SetLegendBorderSize(0);
}

//------------------------------------------------------------
// Main analysis
//------------------------------------------------------------
void ExtractData(TString inputFile, const char* outputFile, int N) {

    std::vector<float> z1, z2, x2, y2;
    std::vector<double> data;

    //--------------------------------------------------------
    // Read existing output file (preserve original logic)
    //--------------------------------------------------------
    std::ifstream inputFileStream(outputFile);
    if (inputFileStream.is_open()) {
        string header;
        string line;

        // 读取第一行
        getline(inputFileStream, header);

        // 读取第二行数据
        if (std::getline(inputFileStream, line)) {
            stringstream ss(line);
            double value;
            while (ss >> value) {
                data.push_back(value);
            }
        }
        inputFileStream.close();
    } else {
        // 如果文件不存在，则 data 为空，后续按 0 处理
        std::cerr << "Warning: " << outputFile
                  << " not found. Previous accumulated data will be treated as zero."
                  << std::endl;
    }

    //--------------------------------------------------------
    // Open output file for overwrite
    //--------------------------------------------------------
    std::ofstream outputFileStream(outputFile);
    if (!outputFileStream.is_open()) {
        std::cerr << "Error: Failed to open output file " << outputFile << std::endl;
        return;
    }

    //--------------------------------------------------------
    // Open ROOT file and tree
    //--------------------------------------------------------
    TString RootFile = inputFile + to_string(N) + TString(".root");

    TFile* f = TFile::Open(RootFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Failed to open input file " << RootFile << std::endl;
        outputFileStream.close();
        return;
    }

    TTree* t = (TTree*)f->Get("Tree");
    if (!t) {
        std::cerr << "Error: Failed to find TTree in input file " << RootFile << std::endl;
        f->Close();
        outputFileStream.close();
        return;
    }

    int ne = 0, ni = 0;
    vector<double>* e1hit = nullptr;
    vector<double>* e2hit = nullptr;
    vector<double>* i1hit = nullptr;
    vector<double>* i2hit = nullptr;

    TBranch* be1hit = nullptr;
    TBranch* be2hit = nullptr;
    TBranch* bi1hit = nullptr;
    TBranch* bi2hit = nullptr;

    t->SetBranchAddress("e1hit", &e1hit, &be1hit);
    t->SetBranchAddress("e2hit", &e2hit, &be2hit);
    t->SetBranchAddress("i1hit", &i1hit, &bi1hit);
    t->SetBranchAddress("i2hit", &i2hit, &bi2hit);
    t->SetBranchAddress("ne", &ne);
    t->SetBranchAddress("ni", &ni);

    double eps = 1e-6 * thick;

    //--------------------------------------------------------
    // Extract electron / ion z positions
    //--------------------------------------------------------
    for (Int_t i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);

        // Electron hits
        for (UInt_t j = 2; j < e2hit->size(); j += 5) {
            double val = e2hit->at(j);
            if (val >= induce + metal - eps && val <= induce + metal + thick) {
                z1.push_back(val);
            }
        }

        // Ion hits
        for (UInt_t j = 2; j < i2hit->size(); j += 4) {
            double val = i2hit->at(j);
            if (val >= induce + metal - eps && val <= induce + metal + thick) {
                x2.push_back(i2hit->at(j - 2));
                y2.push_back(i2hit->at(j - 1));
                z2.push_back(val);
            }
        }
    }

    t->ResetBranchAddresses();
    f->Close();

    //--------------------------------------------------------
    // Figure output path
    //--------------------------------------------------------
    TString figPdf = TString("./fig/ZDistribution_") + to_string(N) + TString(".pdf");
    TString figSvg = TString("./fig/ZDistribution_") + to_string(N) + TString(".svg");

    //--------------------------------------------------------
    // Histogram binning for output accumulation (original logic)
    // These histograms remain in cm
    //--------------------------------------------------------
    const int nbins = 22;
    double edges[nbins + 1];

    // bin #1: [induce+metal-eps, induce+metal]
    edges[0] = induce + metal - eps;
    edges[1] = induce + metal;

    // intermediate bins
    for (int i = 2; i < nbins; i++) {
        double fraction = double(i - 1) / 20.0;
        edges[i] = induce + metal + eps + fraction * (thick - eps);
    }

    // last edge
    edges[nbins] = induce + metal + thick;

    TH1D* h1 = new TH1D("h1", "", nbins, edges); // electron counts for output
    TH1D* h2 = new TH1D("h2", "", nbins, edges); // ion counts for output

    //--------------------------------------------------------
    // Histograms for plotting only
    // Convert x-axis to mm here, while keeping original logic above
    //--------------------------------------------------------
    TH1D* h3 = new TH1D("h3", "", 20,
                        10.0 * (induce + metal),
                        10.0 * (induce + metal + thick));

    TH1D* h4 = new TH1D("h4", "", 20,
                        10.0 * (induce + metal),
                        10.0 * (induce + metal + thick));

    //--------------------------------------------------------
    // Fill histograms
    //--------------------------------------------------------
    int en = static_cast<int>(z1.size());
    for (int i = 0; i < en; ++i) {
        h1->Fill(z1[i]);            // original logic in cm
        h3->Fill(10.0 * z1[i]);     // plotting in mm
    }

    int in = static_cast<int>(z2.size());
    for (int i = 0; i < in; ++i) {
        h2->Fill(z2[i]);            // original logic in cm
        h4->Fill(10.0 * z2[i]);     // plotting in mm
    }

    //--------------------------------------------------------
    // Plot
    //--------------------------------------------------------
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("fig", kTRUE);
    SetPaperStyle();

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 650);
    c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.13);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetTicks(1, 1);
    c1->SetLogy(1);

    // Electrons: blue hatched
    h3->SetLineColor(kBlue + 1);
    h3->SetLineWidth(2);
    h3->SetFillColor(kBlue + 1);
    h3->SetFillStyle(3354);

    // Ions: red hatched
    h4->SetLineColor(kRed + 1);
    h4->SetLineWidth(2);
    h4->SetFillColor(kRed + 1);
    h4->SetFillStyle(3345);

    h3->GetXaxis()->SetTitle("z [mm]");
    h3->GetYaxis()->SetTitle("Counts");

    h3->GetXaxis()->SetTitleFont(42);
    h3->GetYaxis()->SetTitleFont(42);
    h3->GetXaxis()->SetLabelFont(42);
    h3->GetYaxis()->SetLabelFont(42);

    h3->GetXaxis()->SetTitleSize(0.060);
    h3->GetYaxis()->SetTitleSize(0.055);
    h3->GetXaxis()->SetLabelSize(0.045);
    h3->GetYaxis()->SetLabelSize(0.045);

    // y-axis starts from 10
    double ymax = std::max(h3->GetMaximum(), h4->GetMaximum());
    h3->SetMinimum(500.0);
    h3->SetMaximum(std::max(20.0, ymax * 5.0));

    h3->Draw("HIST");
    h4->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.62, 0.78, 0.88, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->AddEntry(h3, "Electrons", "f");
    leg->AddEntry(h4, "Ions", "f");
    leg->Draw();

    c1->SaveAs(figPdf);
    c1->SaveAs(figSvg);

    //--------------------------------------------------------
    // Write updated output data (preserve original logic)
    //--------------------------------------------------------
    int numBins = h1->GetNbinsX();

    for (int i = 1; i <= numBins; ++i) {
        outputFileStream << "n" << i << (i < numBins ? "," : "\n");
    }

    for (int i = 1; i <= numBins; ++i) {
        double e = h1->GetBinContent(i);
        double n = h2->GetBinContent(i);
        double dataValue = (i - 1 < static_cast<int>(data.size())) ? data[i - 1] : 0.0;
        double count = (i == 1 || i == numBins) ? 0.0 : ((n - e) * 30.0 + dataValue);
        outputFileStream << count << (i < numBins ? " " : "\n");
    }

    outputFileStream.close();

    //--------------------------------------------------------
    // Cleanup
    //--------------------------------------------------------
    delete leg;
    delete c1;
    delete h1;
    delete h2;
    delete h3;
    delete h4;
}

//------------------------------------------------------------
// main
//------------------------------------------------------------
int main(int argc, char* argv[]) {
    char Path[100] = "./result";

    if (argc > 1) {
        int N = std::stoi(argv[1]);
        TString inputFile = TString(Path) + TString("/result");
        const char* outputFile = "output.txt";
        ExtractData(inputFile, outputFile, N);
    } else {
        std::cerr << "Usage: " << argv[0] << " N" << std::endl;
        return 1;
    }

    return 0;
}