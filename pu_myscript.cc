#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include <TProfile.h>
#include <TF1.h>
#include <TRandom.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include <iostream>
#include <algorithm>
#include <string>

using namespace std;

#include "Pulse.h"

#include "IniEventFile.C"
#include "FillWaveform.C"
#include "extractsigma.C"
#include "CreateNSamples.C"



int main (int argc, char**argv)
{

  TString inputfile=argv[1];
  //  TString shapingtime=argv[2];
  TString addendum="ns.root";
  //  inputfile+=shapingtime;
  inputfile+=addendum;
  cout<<inputfile<<endl;
  string method=argv[1];
  transform(method.begin(), method.end(),method.begin(),::toupper);
  TString methodupper=method;
  cout<<"methodupper"<<methodupper<<endl;

  TString inieventfile=IniEventFile(inputfile,methodupper);//,argv[2]);
  //  TString inieventfile=IniEventFile("crrc43ns.root");
  cout<<inieventfile<<endl;

  
  double ampli_uncerts [20];
  double tim_uncerts [20];
  double ampli_uncerts_err [20];
  double tim_uncerts_err [20];
  double pileup [20];
  double *pileup_err=0;
  cout <<"1"<<endl;
  cout<<argv[2]<<endl;
  cout<<argv[3]<<endl;
  cout<<argv[4]<<endl;
  cout<<argv[5]<<endl;

  int nsmpl=atoi(argv[2]);
  cout<<nsmpl<<endl;
  int nfreq=atoi(argv[3]);
  int num_points=0;


 
  char canvasname1[120];
  sprintf(canvasname1,"fit_samples_4 events_event");
  TCanvas* canv1=new TCanvas(canvasname1,canvasname1,1000,1000);
  canv1->Divide (2,2);

  int canvaspad=1;
  int peakphase=atoi(argv[4]);
  int risestart=atoi(argv[5]);
  
  cout <<"2"<<endl;

  for (int npu=0;npu<20;npu+=1)
    {
      num_points+=1;
      TString waveforms =FillWaveform(inieventfile,nsmpl,nfreq,npu*10);
      TString samples_file =CreateNSamples(waveforms,inieventfile,peakphase,nsmpl,nfreq);
      double* unc;
      unc=extractsigma(samples_file,inieventfile,npu,canv1,canvaspad,nsmpl,nfreq,risestart);
      ampli_uncerts[npu]=unc[0];
      tim_uncerts[npu]=unc[1];
      ampli_uncerts_err[npu]=unc[2];
      tim_uncerts_err[npu]=unc[3];
      pileup[npu]=npu*10;
    }

  cout <<"3"<<endl;

  TGraphErrors* ampU_v_pileup=new TGraphErrors(num_points,pileup,ampli_uncerts,pileup_err,ampli_uncerts_err);
  TGraphErrors* timU_v_pileup=new TGraphErrors(num_points,pileup,tim_uncerts,pileup_err,tim_uncerts_err);
  
  char canvasname[120];
  sprintf(canvasname,"%sns_Uncertainity_vs_Pileup_%i_samples_%i_frequency",argv[1],nsmpl,nfreq);
  TCanvas* canv=new TCanvas(canvasname,canvasname,1000,500);

  cout<<"average_amplitude_uncertainity = "<<ampU_v_pileup->GetMean(2)<<endl;
  cout<<"average_timing_uncertainity = "<<timU_v_pileup->GetMean(2)<<endl;

  canv->Divide(2,1);
  canv->cd();
  canv->cd(1);
  ampU_v_pileup->Draw("ap");
  ampU_v_pileup->GetXaxis()->SetTitle("pileup");
  ampU_v_pileup->GetYaxis()->SetTitle("amplitute_uncertainity");
  ampU_v_pileup->SetTitle("ampl_uncertainity vs pileup");

  canv->cd(2);
  timU_v_pileup->Draw("ap");
  timU_v_pileup->GetXaxis()->SetTitle("pileup");
  timU_v_pileup->GetYaxis()->SetTitle("timing_uncertainity");
  timU_v_pileup->SetTitle("tim_uncertainity vs pileup");
  
  char savename[120];
  char savename1[120];

  cout<<"jadhjas"<<endl;
  sprintf(savename,"all_plots/pileup/Uncertainity/%sns_Uncertainity_vs_Pileup_%i_samples_%i_frequency_efS.png",argv[1],nsmpl,nfreq);
	  

  sprintf(savename1,"all_plots/pileup/Fits/%sns_Fits_%i_samples_%i_frequency_efS.png",argv[1],nsmpl,nfreq); 

  cout<<"asd"<<endl;
  canv->SaveAs(savename);

  cout<<"chumprasta"<<endl;
  canv1->SaveAs(savename1);

  return 2;
}

  
  

