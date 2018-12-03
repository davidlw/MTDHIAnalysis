void tof()
{
  TFile* f = new TFile("../test/mtdhi_30ps.root");
  TH2D* h = (TH2D*)f->Get("mtdhiana/BetavsP");
  
  TCanvas* c = new TCanvas("c","c",600,600);
  h->Draw("colz");

  double ptbins[] = {0.7,0.8,0.9,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,5.0,6.0};

  TH1D* h1D[100];

  for(int i=0;i<12;i++)
  {
    int xmin = h->GetXaxis()->FindBin(ptbins[i]);
    int xmax = h->GetXaxis()->FindBin(ptbins[i+1]);

    h1D[i] = (TH1D*)h->ProjectionY(Form("h1D_%d",i),xmin,xmax,"e");
  }
  
  TCanvas* cc = new TCanvas("cc","cc",1200,900); 
  cc->Divide(4,3);
  for(int i=0;i<12;i++)
  {
    cc->cd(i+1);
    h1D[i]->SetMarkerStyle(20);
    h1D[i]->Draw("PE");
  }

}
