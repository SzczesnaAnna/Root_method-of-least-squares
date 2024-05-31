#include "TH1D.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TGraph.h"

double dopasuj(Int_t st,Int_t N,Double_t *X,Double_t *Y,Double_t *sigmaj,Double_t *w,Double_t *uw)
{
    TMatrixD A(N,st);
    TMatrixD H(N,N);
    TMatrixD G(N,N);
    TMatrixD y(N,1);
    TMatrixD c(N,1);
    for(int i = 0; i < N; i++)
    {
        G(i,i) = 1/TMath::Power(sigmaj[i], 2);
        H(i,i) = 1/sigmaj[i];
        y(i, 0) = Y[i];
        c(i, 0) = Y[i];
        for(int j = 0; j < st; j++)
        {
            A(i,j) = TMath::Power(X[i], j);
        }
    }
   /* cout<< "A:" << endl;
    A.Print();
    cout<< "H:" << endl;
    H.Print();
    cout<< "G:" << endl;
    G.Print();
    cout<< "Y:" << endl;
    y.Print();*/

    TMatrixD Apr(H,TMatrixD::kMult,A);
    TMatrixD cpr(H,TMatrixD::kMult,c);
    TMatrixD AprTr(TMatrixD::kTransposed,Apr);
    TMatrixD AprTrApr(AprTr,TMatrixD::kMult,Apr);
    AprTrApr.Invert();
    TMatrixD AprTrcpr(AprTr,TMatrixD::kMult,cpr);
    TMatrixD xfl(AprTrApr,TMatrixD::kMult,AprTrcpr);
     for(int j = 0; j < st; j++)
        {
            w[j] = xfl(j,0);
        }

     for(int j = 0; j < st; j++)
        {
            uw[j]=TMath::Sqrt(AprTrApr(j,j));
        }

    TMatrixD nfl(A,TMatrixD::kMult,xfl);
    cout<< "x~:" << endl;
   // xfl.Print();
   // nfl.Print();

    TMatrixD y_nfl = TMatrixD(y-nfl).T();
    auto M = y_nfl*G*(y-nfl);
    M.Print();
double x=M(0,0);
    return x;
}
void least_squares_method(){
    double x[] = { -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
    double y[] = { 86, 50, 26, 24, 20, 55, 113, 186, 339, 601};
    double ux[] = {0,0,0,0,0,0,0,0,0,0,0};
    double uy[] = {10, 14, 13, 9, 13, 12, 12, 13, 11, 14};
    const int N = 10;
   
    double W[N];
    double UW[N];

    auto c = new TCanvas();      
    c->cd(1);
    auto g = new TGraphErrors(N,x,y,ux,uy);
    g->Draw();
     
    for(int i=1; i<7; i++)
    {
        dopasuj(i,N,x,y,uy,W,UW);
        auto dop = new TF1(Form("pol%i",i),Form("pol%i",i), -1,1);
        dop->SetParameters(W);
        dop->SetLineColor(i+1);
        dop->Draw("SAME"); 
        double ndf=N-1-i;
        
    }


   // c->Print("least_squares_method.pdf");



}