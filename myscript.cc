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

using namespace std;

#include "Pulse.h"

#include "IniEventFile.C"
#include "FillWaveform.C"
#include "FitPulseShape.C"
#include "CreateNSamples.C"



int main (int argc, char**argv)
{
  TString inieventfile=IniEventFile("crrc43ns.root");
  cout<<inieventfile<<endl;
  TString waveforms =FillWaveform(inieventfile);
  cout<<waveforms<<endl;
  
  double ampli_uncerts [101];
  double tim_uncerts [101];
  double ampli_uncerts_err [101];
  double tim_uncerts_err [101];
  double phase [101];
  double *phase_err=0;
  int nsmpl=atoi(argv[1]);
  int nfreq=atoi(argv[2]);
  int num_points=0;

  char canvasname1[120];
  sprintf(canvasname1,"fit_samples_4 events_event");
  TCanvas* canv1=new TCanvas(canvasname1,canvasname1,1000,1000);
  char canvasname2[120];
  char canvasname3[120];
  sprintf(canvasname2,"amplitude_histos");
  sprintf(canvasname3,"timing_histos");
  TCanvas* canv2=new TCanvas(canvasname2,canvasname2,1000,1000);
  TCanvas* canv3=new TCanvas(canvasname3,canvasname3,1000,1000);
  canv2->Divide(2,2);
  canv3->Divide(2,2);
  canv1->Divide (2,2);

  int canvaspad=1;


  for (int phaseshift=0;phaseshift<=100;phaseshift++)
    {
      num_points+=1;
      TString samples_file =CreateNSamples(waveforms,0.5*phaseshift,nsmpl,nfreq);
      double* unc;
      cout<<endl;
      cout<<"performing fit for phaseshift:"<<phaseshift<<endl;
      cout<<endl;
      unc=FitPulseShape(samples_file,phaseshift,canv1,canv2,canv3,canvaspad,nsmpl,nfreq);
      ampli_uncerts[phaseshift]=unc[0];
      tim_uncerts[phaseshift]=unc[1];
      ampli_uncerts_err[phaseshift]=unc[2];
      tim_uncerts_err[phaseshift]=unc[3];
      phase[phaseshift]=IDSTART-0.5*phaseshift;
    }
  TGraphErrors* ampU_v_phase=new TGraphErrors(num_points,phase,ampli_uncerts,phase_err,ampli_uncerts_err);
  TGraphErrors* timU_v_phase=new TGraphErrors(num_points,phase,tim_uncerts,phase_err,tim_uncerts_err);
  
  char canvasname[120];
  sprintf(canvasname,"Uncertainity_vs_Phase_shift_%i_samples_%i_frequency", nsmpl,nfreq);
  TCanvas* canv=new TCanvas(canvasname,canvasname,1000,500);

    

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
  char savename2[120];
  char savename3[120];
  sprintf(savename,"all_plots/Uncertainity/Uncertainity_vs_Phase_shift_%i_samples_%i_frequency.png",nsmpl,nfreq);
	  
	  sprintf(savename1,"all_plots/Fits/Fits_%i_samples_%i_frequency.png",nsmpl,nfreq); 
	  sprintf(savename2,"all_plots/Ampli_histos/Amplitude_histos%i_samples_%i_frequency.png",nsmpl,nfreq); 
	  sprintf(savename3,"all_plots/Timing_histos/Timing_histos%i_samples_%i_frequency.png",nsmpl,nfreq);

  char canvtitle[120];
  char canvtitle1[120];
  char canvtitle2[120];
  char canvtitle3[120];
  sprintf(canvtitle,"all_plots/Uncertainity/Uncertainity_vs_Phase_shift");
	  
	  sprintf(canvtitle1,"all_plots/Fits/Fits"); 
	  sprintf(canvtitle2,"all_plots/Ampli_histos/Amplitude_histos"); 
	  sprintf(canvtitle3,"all_plots/Timing_histos/Timing_histos");

  canv->SaveAs(savename);
  canv1->SaveAs(savename1);
  canv2->SaveAs(savename2);
  canv3->SaveAs(savename3);

  canv->SetTitle(canvtitle);
  canv1->SetTitle(canvtitle1);
  canv2->SetTitle(canvtitle2);
  canv3->SetTitle(canvtitle3);

  return 2;
}
  
  

