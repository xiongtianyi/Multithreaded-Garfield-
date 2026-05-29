//------------------------------------------------------------------//
// Draw result1: gain evolution vs iteration number
// Fit: G(n) = G_inf + A1 exp(-n/tau1) + A2 exp(-n/tau2)
// Error bar: 1% relative error for each gain point
// Author: Tianyi Xiong, revised for JINST response figures
//------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TMath.h>
#include <TString.h>

using namespace std;

// ------------------------------------------------------------------
// Global style close to the reference charging-up figure
// ------------------------------------------------------------------
void SetResult1Style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetFrameLineWidth(1);
    gStyle->SetLineWidth(2);
    gStyle->SetEndErrorSize(4);

    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTextFont(42);

    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.05, "Y");
    gStyle->SetNdivisions(510, "XY");
}

// ------------------------------------------------------------------
// Mean value of the tail region, used as the initial estimate of G_inf
// ------------------------------------------------------------------
double TailMean(const vector<double>& y, int nTail) {
    if (y.empty()) return 0.0;
    const int n = static_cast<int>(y.size());
    const int i0 = std::max(0, n - nTail);
    double sum = 0.0;
    int cnt = 0;
    for (int i = i0; i < n; ++i) {
        sum += y[i];
        ++cnt;
    }
    return (cnt > 0) ? sum / cnt : y.back();
}

// ------------------------------------------------------------------
// Main function
// inputPrefix example: "./result/result"
// It reads ./result/result1.root, ./result/result2.root, ...
// ------------------------------------------------------------------
void ExtractData(TString inputPrefix = "./result/result",
                 int N = 160,
                 int totalEvents = 10000) {

    SetResult1Style();
    gSystem->mkdir("./fig", kTRUE);

    vector<double> iter;
    vector<double> gain;
    vector<double> gainErr;

    iter.reserve(N);
    gain.reserve(N);
    gainErr.reserve(N);

    for (int i = 1; i <= N; ++i) {
        TString rootFile = inputPrefix + to_string(i) + TString(".root");
        TFile* f = TFile::Open(rootFile, "READ");
        if (!f || f->IsZombie()) {
            cerr << "Error: failed to open " << rootFile << endl;
            if (f) f->Close();
            break;
        }

        TTree* t = dynamic_cast<TTree*>(f->Get("Tree"));
        if (!t) {
            cerr << "Error: failed to find TTree named Tree in " << rootFile << endl;
            f->Close();
            break;
        }

        vector<double>* e2hit = nullptr;
        TBranch* be2hit = nullptr;
        t->SetBranchAddress("e2hit", &e2hit, &be2hit);

        long long totalElectronEndpoints = 0;
        const Long64_t nEntries = t->GetEntries();
        for (Long64_t a = 0; a < nEntries; ++a) {
            t->GetEntry(a);
            if (!e2hit) continue;

            // e2hit stores endpoint information in groups of five:
            // x, y, z, t, energy. One endpoint therefore corresponds to one stride.
            totalElectronEndpoints += static_cast<long long>(e2hit->size() / 5);
        }

        t->ResetBranchAddresses();
        f->Close();

        const double g = static_cast<double>(totalElectronEndpoints) /
                         static_cast<double>(totalEvents);

        // Use iteration number starting from 0, consistent with charging-up plots.
        // If you prefer the first point to be 1, change i - 1 to i.
        iter.push_back(static_cast<double>(i - 1));
        gain.push_back(g);

        // 1% relative error requested by the user.
        gainErr.push_back(0.01 * g);
    }

    const int nPoints = static_cast<int>(gain.size());
    if (nPoints < 6) {
        cerr << "Error: too few valid points for a double-exponential fit." << endl;
        return;
    }

    vector<double> xErr(nPoints, 0.0);

    // Determine y-axis range and leave enough blank space for the parameter text.
    double yMin = gain[0] - gainErr[0];
    double yMax = gain[0] + gainErr[0];
    for (int i = 0; i < nPoints; ++i) {
        yMin = std::min(yMin, gain[i] - gainErr[i]);
        yMax = std::max(yMax, gain[i] + gainErr[i]);
    }
    const double ySpan = std::max(yMax - yMin, 1.0e-9);
    const double yLow  = std::max(0.0, yMin - 0.08 * ySpan);
    const double yHigh = yMax + 0.45 * ySpan;

    // ------------------------------------------------------------------
    // Build graph with 1% error bars
    // ------------------------------------------------------------------
    TGraphErrors* gr = new TGraphErrors(nPoints,
                                        iter.data(), gain.data(),
                                        xErr.data(), gainErr.data());
    gr->SetTitle("");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.9);
    gr->SetMarkerColor(kBlack);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(1);

    gr->GetXaxis()->SetTitle("Iteration number");
    gr->GetYaxis()->SetTitle("Gain");
    gr->GetXaxis()->SetLimits(0.0, iter.back() + 8.0);
    gr->GetYaxis()->SetRangeUser(yLow, yHigh);
    gr->GetXaxis()->SetTitleOffset(1.08);
    gr->GetYaxis()->SetTitleOffset(0.95);
    gr->GetXaxis()->SetNdivisions(510);
    gr->GetYaxis()->SetNdivisions(510);

    // ------------------------------------------------------------------
    // Double-exponential fit
    // G(n) = G_inf + A1 exp(-n/tau1) + A2 exp(-n/tau2)
    // ------------------------------------------------------------------
    const double xFitMin = iter.front();
    const double xFitMax = iter.back();

    TF1* fit = new TF1("fit_result1_biexp",
                       "[0] + [1]*exp(-x/[2]) + [3]*exp(-x/[4])",
                       xFitMin, xFitMax);
    fit->SetParNames("Ginf", "A1", "tau1", "A2", "tau2");

    const double g0     = gain.front();
    const double gTail  = TailMean(gain, std::max(8, nPoints / 8));
    const double amp    = std::max(g0 - gTail, 1.0e-6);

    fit->SetParameters(gTail,
                       0.60 * amp, 2.0,
                       0.40 * amp, 15.0);

    // The following limits keep the two exponential components stable.
    // If your data require a much slower component, increase the upper limit of tau2.
    fit->SetParLimits(0, 0.0, std::max(10.0 * g0, 1.0));
    fit->SetParLimits(1, 0.0, std::max(10.0 * amp, 1.0));
    fit->SetParLimits(2, 0.05, 20.0);      // fast component
    fit->SetParLimits(3, 0.0, std::max(10.0 * amp, 1.0));
    fit->SetParLimits(4, 2.0, 5000.0);     // slow component

    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    fit->SetLineStyle(1);

    // ------------------------------------------------------------------
    // Draw figure
    // ------------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1", "result1", 760, 650);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.04);
    c1->SetTopMargin(0.04);
    c1->SetTicks(1, 1);

    gr->Draw("AP E1");
    gr->Fit(fit, "RQ");
    fit->Draw("SAME");
    gr->Draw("P E1 SAME");

    // ------------------------------------------------------------------
    // Text block inside the plot, arranged in the same style as the
    // reference THGEM charging-up figure. Modify these setup strings when
    // drawing a different condition.
    // ------------------------------------------------------------------
    const TString detectorGasText = "THGEM; Ne/CH_{4}(5%); #it{pen}=0.4";
    const TString geomText        = "#it{t}=0.8; #it{a}=1; #it{d}=0.5; noRim [mm]";
    const TString fieldText       = "#it{E}_{d}=0.5 kV/cm, #it{E}_{i}=0.5 kV/cm";
    const TString voltageText     = "#DeltaV=600; n_{AV}=360";
    const TString rateText        = "Rate=10[Hz]; E=8[keV]";

    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    const double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
    const double g0fit = fit->Eval(0.0);
    const double gInffit = fit->GetParameter(0);

    TString stepText;
    if (totalEvents >= 1000 && totalEvents % 1000 == 0) {
        stepText = Form("step=%dk", totalEvents / 1000*30);
    } else {
        stepText = Form("step=%d", totalEvents);
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.037);
    latex.SetTextAlign(22);

    const double xText = 0.58;
    double yText = 0.90;
    const double dy = 0.049;

    latex.DrawLatex(xText, yText, detectorGasText);
    yText -= dy;
    latex.DrawLatex(xText, yText, geomText);
    yText -= dy;
    latex.DrawLatex(xText, yText, fieldText);
    yText -= dy;
    latex.DrawLatex(xText, yText, voltageText);
    yText -= dy;
    latex.DrawLatex(xText, yText, rateText);
    yText -= dy;
    latex.DrawLatex(xText, yText, Form("G_{0}^{FIT}=%.2f; G_{#infty}^{FIT}=%.2f", g0fit, gInffit));
    yText -= dy;
    latex.DrawLatex(xText, yText, Form("Fit #chi^{2}/n.d.f.=%.2f", chi2ndf));
    yText -= dy;
    latex.DrawLatex(xText, yText,
                    Form("#tau=%.2f+%.2f [iteration]; %s",
                         fit->GetParameter(2), fit->GetParameter(4), stepText.Data()));

    c1->SaveAs("./fig/result1.pdf");
    c1->SaveAs("./fig/result1.svg");
    c1->SaveAs("./fig/result1.png");

    // ------------------------------------------------------------------
    // Print and save fit summary
    // ------------------------------------------------------------------
    cout << "\n===== result1 double-exponential fit =====" << endl;
    cout << "Function: G(n) = Ginf + A1 exp(-n/tau1) + A2 exp(-n/tau2)" << endl;
    cout << "Points: " << nPoints << endl;
    cout << "Error: 1% relative error for each point" << endl;
    cout << "chi2/ndf = " << chi2 << " / " << ndf << " = " << chi2ndf << endl;
    for (int ip = 0; ip < fit->GetNpar(); ++ip) {
        cout << fit->GetParName(ip) << " = "
             << fit->GetParameter(ip) << " +/- "
             << fit->GetParError(ip) << endl;
    }

    ofstream fout("./fig/result1_fit_summary.txt");
    fout << "result1 double-exponential fit\n";
    fout << "Function: G(n) = Ginf + A1 exp(-n/tau1) + A2 exp(-n/tau2)\n";
    fout << "Error: 1% relative error for each gain point\n";
    fout << "Number of points: " << nPoints << "\n";
    fout << "chi2/ndf = " << chi2 << " / " << ndf << " = " << chi2ndf << "\n";
    for (int ip = 0; ip < fit->GetNpar(); ++ip) {
        fout << fit->GetParName(ip) << " = "
             << fit->GetParameter(ip) << " +/- "
             << fit->GetParError(ip) << "\n";
    }
    fout.close();
}

// ------------------------------------------------------------------
// Standalone usage:
//   g++ result1_biexp_reference_legend.C $(root-config --cflags --libs) -o draw_result1
//   ./draw_result1 160 ./result/result 10000
//
// ROOT macro usage:
//   root -l -q 'result1_biexp_reference_legend.C("./result/result",160,10000)'
// ------------------------------------------------------------------
#ifndef __CLING__
int main(int argc, char* argv[]) {
    int N = 160;
    TString inputPrefix = "./result/result";
    int totalEvents = 10000;

    if (argc > 1) N = std::stoi(argv[1]);
    if (argc > 2) inputPrefix = TString(argv[2]);
    if (argc > 3) totalEvents = std::stoi(argv[3]);

    ExtractData(inputPrefix, N, totalEvents);
    return 0;
}
#endif
