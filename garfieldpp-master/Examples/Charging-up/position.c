#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>

// ---- 与 analysis.c 一致的几何（单位：cm）
static const double pitch  = 0.1;
static const double metal  = 0.0012;
static const double thick  = 0.04;
static const double induce = 0.2;

static inline double cm2um(double v_cm) { return v_cm * 1e4; }

void PlotEndPoints3D(const TString& filePrefix, int N) {
  // 打开 ROOT 文件
  TString rootFile = filePrefix + TString::Format("%d.root", N);
  TFile* f = TFile::Open(rootFile, "READ");
  if (!f || f->IsZombie()) { std::cerr << "Open failed: " << rootFile << "\n"; return; }

  TTree* t = (TTree*)f->Get("Tree");
  if (!t) { std::cerr << "Tree not found.\n"; f->Close(); return; }

  // 绑定分支
  std::vector<double>* e2hit = nullptr; // 每个电子端点: x,y,z,t,E  共5个值
  std::vector<double>* i2hit = nullptr; // 每个离子端点:   x,y,z,t    共4个值
  int ne = 0, ni = 0;
  t->SetBranchAddress("e2hit", &e2hit);
  t->SetBranchAddress("i2hit", &i2hit);
  t->SetBranchAddress("ne", &ne);
  t->SetBranchAddress("ni", &ni);

  // 收集端点（只要落在电介质板厚范围内）
  const double zmin = 0;
  const double zmax = induce + metal + thick + 1e-6 * thick;
  const double zc   = induce + metal + 0.5 * thick;   // 中面 (cm)

  std::vector<double> ex, ey, ez; // µm, 以中面为0
  std::vector<double> ix, iy, iz; // µm, 以中面为0

  Long64_t nent = t->GetEntries();
  ex.reserve(200000); ey.reserve(200000); ez.reserve(200000);
  ix.reserve(200000); iy.reserve(200000); iz.reserve(200000);

  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);

    // 电子最终端点 e2hit: [x,y,z,t,E] 步长 5
    for (size_t j = 0; j + 2 < e2hit->size(); j += 5) {
      double z = e2hit->at(j + 2);
      if (z >= zmin && z <= zmax) {
        ex.push_back(cm2um(e2hit->at(j)));
        ey.push_back(cm2um(e2hit->at(j + 1)));
        ez.push_back(cm2um(z - zc));
      }
    }
    // 离子最终端点 i2hit: [x,y,z,t]   步长 4
    for (size_t j = 0; j + 2 < i2hit->size(); j += 4) {
      double z = i2hit->at(j + 2);
      if (z >= zmin && z <= zmax) {
        ix.push_back(cm2um(i2hit->at(j)));
        iy.push_back(cm2um(i2hit->at(j + 1)));
        iz.push_back(cm2um(z - zc));
      }
    }
  }
  f->Close();

  // 准备 3D 点云
  TGraph2D* gEle = new TGraph2D((int)ex.size());
  for (int i = 0; i < (int)ex.size(); ++i) gEle->SetPoint(i, ex[i], ey[i], ez[i]);

  TGraph2D* gIon = new TGraph2D((int)ix.size());
  for (int i = 0; i < (int)ix.size(); ++i) gIon->SetPoint(i, ix[i], iy[i], iz[i]);

  // 画布：左右两个子图 (a)(b)
  TCanvas* c = new TCanvas("c", "Final endpoints", 1200, 600);
  c->Divide(2, 1);
  gStyle->SetCanvasPreferGL(true);     // OpenGL 渲染更顺滑
  gStyle->SetOptStat(0);

  // ------- (a) 电子 -------
  c->cd(1);
  gEle->SetTitle(";X (#mum);Y (#mum);Z (#mum)");
  gEle->SetMarkerStyle(20);
  gEle->SetMarkerSize(0.5);
  // 半透明紫色（ROOT>=6 可用 SetMarkerColorAlpha；若版本旧用 TColor::GetColorTransparent）
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,00)
  gEle->SetMarkerColorAlpha(kViolet+6, 0.55);
#else
  gEle->SetMarkerColor(TColor::GetColorTransparent(kViolet+6, 0.55));
#endif
  gEle->Draw("p0");  // p0 = 只画点

  // ------- (b) 离子 -------
  c->cd(2);
  gIon->SetTitle(";X (#mum);Y (#mum);Z (#mum)");
  gIon->SetMarkerStyle(20);
  gIon->SetMarkerSize(0.5);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,00)
  gIon->SetMarkerColorAlpha(kAzure-4, 0.45);
#else
  gIon->SetMarkerColor(TColor::GetColorTransparent(kAzure-4, 0.45));
#endif
  gIon->Draw("p0");

  // （可选）加两行 “(a) (b)” 标注
  c->cd(1); TLatex* la = new TLatex(); la->SetNDC(); la->SetTextSize(0.05); la->DrawLatex(0.12,0.88,"(a) Electron");
  c->cd(2); TLatex* lb = new TLatex(); lb->SetNDC(); lb->SetTextSize(0.05); lb->DrawLatex(0.12,0.88,"(b) Ion");

  // 输出
  gSystem->mkdir("fig", true);
  TString out = TString::Format("fig/EndPoints3D-%d.pdf", N);
  c->SaveAs(out);
}

int main(int argc, char** argv) {
  int N = (argc > 1) ? atoi(argv[1]) : 1;
  // 你的结果在 ./result/resultN.root
  PlotEndPoints3D("./result/result", N);
  return 0;
}
