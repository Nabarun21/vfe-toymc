//
// Fit one event with a pulse shape
// To run:
// > root -l Example05.C+
//



Pulse pSh;

double uncertainities[4];


double funcPulseShape( double *x, double *par )
{
  double t = x[0] - par[0];
  return par[1] + par[2] * pSh.fShape(t);
}



double* FitPulseShape(TString samples_file,int iter,TCanvas* canv1,TCanvas* canv2,TCanvas* canv3,int& canvaspad,int nsmpl=NSAMPLES,int nfrq=NFREQ) 
{
  // Fit function
  
  TF1 *fPulseShape = new TF1("fPulseShape", funcPulseShape, -500, 500, 3);
  fPulseShape->SetNpx(1000);
  
  

  TString fileinput=samples_file;

  
  // Get one event with samples

  TFile *file2 = new TFile(fileinput);
  //  TFile *file2 = new TFile("data/samples_signal_10GeV_eta_0.0_pu_140.root");

  double samples[nsmpl];
  double amplitudeTruth;
  TTree *tree = (TTree*)file2->Get("Samples");
  tree->SetBranchAddress("amplitudeTruth",      &amplitudeTruth);
  tree->SetBranchAddress("samples",             samples);

  

  if (iter%12==0 && iter!=0)
    { 
      tree->GetEntry(iter);
      
      // Create TGraphErrors with the pulse to fit
      
      TGraphErrors *gr = new TGraphErrors();
      for(int i=0; i<nsmpl; i++){
	double x = i * nfrq;  
	gr->SetPoint(i, x, samples[i]);
	gr->SetPointError(i, 0., 0.044);  // 44 MeV for all samples
      }
      // Fit 
      canv1->cd(canvaspad);
      fPulseShape->SetParameters(70., 0., amplitudeTruth);
      fPulseShape->SetLineColor(2);
      
      char grname[120];
      sprintf(grname,"fit_%ith_iter_%ith_event",iter,iter);    
      gr->Fit("fPulseShape","E");
      gr->Draw("AP");
      gr->SetTitle(grname);
      
    }
  
  // Errors in reconstructed timing and amplitude

  TH1D * h1 = new TH1D("h1", "amplitude", 10000/nsmpl, -5., 5.);
  TH1D * h2 = new TH1D("h2", "timing", 5000, 30.0-iter, 100.0);
  int nentries = tree->GetEntries();
  for(int ievt=0; ievt<nentries; ievt++){
    tree->GetEntry(ievt);
    TGraphErrors *gr = new TGraphErrors();
    for(int i=0; i<nsmpl; i++){
      double x = i * nfrq;  
      gr->SetPoint(i, x, samples[i]);
      gr->SetPointError(i, 0., 0.044);  // 44 MeV for all samples
    }
    
    fPulseShape->SetParameters(70., 0., amplitudeTruth);
    gr->Fit("fPulseShape","QE");
    h1->Fill(fPulseShape->GetParameter(2)-amplitudeTruth);
    h2->Fill(fPulseShape->GetParameter(0));
    //    if (ievt%10==0)cout<<fPulseShape->GetParameter(0)<<endl;
  }
  cout << " amplitude uncertainty   = " << h1->GetRMS() << " GeV" << endl;
  cout << " timing uncertainty      = " << h2->GetRMS() << " ns" << endl;


  if (iter%12==0 && iter!=0)
    { 
      /*      canv2->cd(canvaspad);
      h1->Draw();
      char histname1[120];
      sprintf(histname1,"amplitude_%ith_iter",iter);
      h1->SetTitle(histname1);
      */
      canv3->cd(canvaspad);
      h2->Draw();
      char histname2[120];
      sprintf(histname2,"timing_%ith_iter",iter);
      h2->SetTitle(histname2);
      canvaspad+=1;
    }
  
    
  uncertainities[0]=h1->GetRMS();
  uncertainities[1]=h2->GetRMS();
  uncertainities[2]=h1->GetRMSError();
  uncertainities[3]=h2->GetRMSError();
  return uncertainities;
}
