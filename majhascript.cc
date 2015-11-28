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
#include "extractsigma2.C"
#include "CreateNSamples.C"



int main (int argc, char**argv)
{

  TString inputfile=argv[1];
  TString shapingtime=argv[2];
  TString addendum="ns.root";
  inputfile+=shapingtime;
  inputfile+=addendum;
  cout<<inputfile<<endl;
  string method=argv[1];
  transform(method.begin(), method.end(),method.begin(),::toupper);
  TString methodupper=method;
  cout<<"methodupper"<<methodupper<<endl;

  TString inieventfile=IniEventFile(inputfile,methodupper,argv[2]);
  //  TString inieventfile=IniEventFile("crrc43ns.root");
  cout<<inieventfile<<endl;

  double ampli_uncerts [101];
  double tim_uncerts [101];
  double ampli_uncerts_err [101];
  double tim_uncerts_err [101];
  double phase [101];
  double *phase_err=0;
  int nsmpl=atoi(argv[3]);
  int nfreq=atoi(argv[4]);
  int num_points=0;

  TString waveforms =FillWaveform(inieventfile,nsmpl,nfreq);
  cout<<waveforms<<endl;
  


  char canvasname1[120];
  sprintf(canvasname1,"fit_samples_4 events_event");
  TCanvas* canv1=new TCanvas(canvasname1,canvasname1,1000,1000);
  canv1->Divide (2,2);

  int canvaspad=1;


  for (int phaseshift=0;phaseshift<=101;phaseshift++)
    {
      num_points+=1;
      TString samples_file =CreateNSamples(waveforms,inieventfile,0.5*phaseshift,nsmpl,nfreq);
      double* unc;
      cout<<endl;
      cout<<"performing fit for phaseshift:"<<phaseshift<<endl;
      cout<<endl;
      unc=extractsigma2(samples_file,inieventfile,phaseshift,canv1,canvaspad,nsmpl,nfreq);
      //cout<<"guju"<<endl;
      ampli_uncerts[phaseshift]=unc[0];
      tim_uncerts[phaseshift]=unc[1];
      ampli_uncerts_err[phaseshift]=unc[2];
      tim_uncerts_err[phaseshift]=unc[3];
      phase[phaseshift]=IDSTART-0.5*phaseshift;
    }
  TGraphErrors* ampU_v_phase=new TGraphErrors(num_points,phase,ampli_uncerts,phase_err,ampli_uncerts_err);
  TGraphErrors* timU_v_phase=new TGraphErrors(num_points,phase,tim_uncerts,phase_err,tim_uncerts_err);
  
  char canvasname[120];
  sprintf(canvasname,"%s%sns_Uncertainity_vs_Phase_shift_%i_samples_%i_frequency",argv[1],argv[2],nsmpl,nfreq);
  TCanvas* canv=new TCanvas(canvasname,canvasname,1000,500);

  cout<<"average_amplitude_uncertainity = "<<ampU_v_phase->GetMean(2)<<endl;
  cout<<"average_timing_uncertainity = "<<timU_v_phase->GetMean(2)<<endl;

  canv->Divide(2,1);
  canv->cd();
  canv->cd(1);
  ampU_v_phase->Draw("ap");
  ampU_v_phase->GetXaxis()->SetTitle("phase");
  ampU_v_phase->GetYaxis()->SetTitle("amplitute_uncertainity");
  ampU_v_phase->SetTitle("ampl_uncertainity vs phase");

  canv->cd(2);
  timU_v_phase->Draw("ap");
  timU_v_phase->GetXaxis()->SetTitle("phase");
  timU_v_phase->GetYaxis()->SetTitle("timing_uncertainity");
  timU_v_phase->SetTitle("tim_uncertainity vs phase");
  
  char savename[120];
  char savename1[120];

  sprintf(savename,"all_plots/Uncertainity/%s%sns_Uncertainity_vs_Phase_shift_%i_samples_%i_frequency_efS.png",argv[1],argv[2],nsmpl,nfreq);
	  
  sprintf(savename1,"all_plots/Fits/%s%sns_Fits_%i_samples_%i_frequency_efS.png",argv[1],argv[2],nsmpl,nfreq); 


  char canvtitle[120];
  char canvtitle1[120];

  sprintf(canvtitle,"all_plots/Uncertainity/Uncertainity_vs_Phase_shift");
	  
  sprintf(canvtitle1,"all_plots/Fits/Fits"); 
  

  canv->SaveAs(savename);
  canv1->SaveAs(savename1);
  

  canv->SetTitle(canvtitle);
  canv1->SetTitle(canvtitle1);
  

  return 2;
}

  
  
