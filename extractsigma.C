//
// Fit one event with a pulse shape
// To run:
// > root -l Example07.C+
//


Pulse pSh;

double uncertainities[4];

pair<double,double> effectiveSigma(vector<double> a)
{
  std::sort(a.begin(),a.end());
  unsigned int nInInterval = int(0.68*a.size());
  double minInterval = 1e+9;
  double amin = 0;
  double amax = 0;
  for (unsigned int i=0; i<a.size(); i++) {
    double interval_size=minInterval;
    if ((i+nInInterval)<a.size()) interval_size = (a[i+nInInterval]-a[i]);
    if (interval_size < minInterval){
      minInterval = interval_size;
      amin = a[i];
      amax = amin + minInterval;
    }
  }
  pair<double,double> result;
  result.first  = amin;
  result.second = amax;
  return result;     
}



double funcPulseShape( double *x, double *par )
{
  double t = x[0] - par[0];
  return par[1] + par[2] * pSh.fShape(t);
}



double * extractsigma(TString samples_file,int iter,TCanvas* canv1,int& canvaspad,int nsmpl=NSAMPLES,int nfrq=NFREQ) 
{
  // Fit function
  // cout<<"1"<<endl;
  
  TF1 *fPulseShape = new TF1("fPulseShape", funcPulseShape, -500, 500, 3);
  fPulseShape->SetNpx(1000);

  TString fileinput=samples_file;

  TFile *file2 = new TFile(fileinput);

  // cout<<"2"<<endl;
  // Get one event with samples

  //  TFile *file2 = new TFile("data/samples_qie25_signal_10GeV_pu_0.root");
  //  TFile *file2 = new TFile("data/samples_signal_10GeV_pu_0.root");
  //  TFile *file2 = new TFile("data/samples_signal_10GeV_eta_0.0_pu_140.root");

  double samples[NSAMPLES];
  double amplitudeTruth;
  TTree *tree = (TTree*)file2->Get("Samples");
  tree->SetBranchAddress("amplitudeTruth",      &amplitudeTruth);
  tree->SetBranchAddress("samples",             samples);
  //cout<<"3"<<endl;

  if (iter%24==0 && iter!=0)
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
      canvaspad+=1;
      
      //      delete gr;
    }
  
  // cout<<"4"<<endl;

  // Reconstructed amplitude and timing
  vector<double> resultsA;
  vector<double> resultsT;

  int nentries = tree->GetEntries();
  for(int ievt=0; ievt<nentries; ievt++){
    //cout<<"kije"<<endl;
    tree->GetEntry(ievt);
    TGraphErrors *gr = new TGraphErrors();
    //cout<<"5"<<endl;
    for(int i=0; i<nsmpl; i++){
      double x = i * nfrq;  
      //cout<<"i"<<i<<endl;
      gr->SetPoint(i, x, samples[i]);
      //cout<<"6"<<endl;
      gr->SetPointError(i, 0., 0.044);  // 44 MeV for all samples
      //cout<<"7"<<endl;
    }
    //cout<<"kije2"<<endl;
    fPulseShape->SetParameters(70., 0., amplitudeTruth);
    //cout<<"kije3"<<endl;    
    gr->Fit("fPulseShape","QE");
    //cout<<"kije3"<<endl;    
    //cout<<"param "<<fPulseShape->GetParameter(2)<<endl;
    //cout<<"amptr "<<amplitudeTruth<<endl;
    double inp=fPulseShape->GetParameter(2)-amplitudeTruth;
    double inp2=fPulseShape->GetParameter(0);

    //cout<<"inp "<<inp<<endl;
    resultsA.push_back(inp);
    //cout<<"kije4"<<endl;    
     resultsT.push_back(inp2);
    //cout<<"kije3"<<endl;     
 
}


   if(resultsA.size()>0){
    pair<double,double> effA = effectiveSigma(resultsA);
    double avgValue_A = 0.5 * (effA.first + effA.second);
    double sigValue_A = 0.5 * fabs(effA.second - effA.first) ;
    double avgError_A = sigValue_A / sqrt(double(resultsA.size()));
    double sigError_A = sigValue_A / sqrt(2.0 * double(resultsA.size()));
    cout << " Mean of amplitude      " << avgValue_A << " +- " << avgError_A << " GeV" << endl;
    cout << " Sigma_eff of amplitude " << sigValue_A << " +- " << sigError_A << " GeV" << endl;
    uncertainities[0]=sigValue_A;
    uncertainities[2]=sigError_A;
      
  }

  if(resultsT.size()>0){
    pair<double,double> effT = effectiveSigma(resultsT);
    double avgValue_T = 0.5 * (effT.first + effT.second);
    double sigValue_T = 0.5 * fabs(effT.second - effT.first) ;
    double avgError_T = sigValue_T / sqrt(double(resultsT.size()));
    double sigError_T = sigValue_T / sqrt(2.0 * double(resultsT.size()));
    cout << " Mean of timing      " << avgValue_T << " +- " << avgError_T << " ns" << endl;
    cout << " Sigma_eff of timing " << sigValue_T << " +- " << sigError_T << " ns" << endl;
    uncertainities[1]=sigValue_T;
    uncertainities[3]=sigError_T;

}
  
  return uncertainities;

  
}
