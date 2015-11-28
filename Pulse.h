#include <TRandom.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>

TRandom rnd;

// total number of bunches in "LHC" bunch train
const int NBXTOTAL = 2800;

// length of a waveform in 1ns steps
const int WFLENGTH  = 500;

// number of samples per hit
const int NSAMPLES   = 10;

// distance between samples in 1ns steps
const int NFREQ      = 25;

// position of a 1st sample inside waveform
const int IDSTART    = 180;

// CRRC shaping time in ns. It is used to calculate noise
// correlations. For QIE, set it to 1e-1
const double TAU = 43.0;

// filename with pulse shape
const TString FNAMESHAPE_DEF       = "data/EmptyFileCRRC43.root";



class Pulse{

  double *weights_;
  double *mC_;
  double *mL_;
  TGraph *grPS_;
  float tMin_;
  float fPar0_;
  float fPar1_;
  TFile *filePS_;
  int nsamples_;
  int nfreq_;

 public:

  Pulse(TString fnameshape=FNAMESHAPE_DEF,int nsmpl=NSAMPLES,int nfreq=NFREQ);
  ~Pulse();
  TGraph* grPS() {return grPS_; };
  float tMin() const { return tMin_; };
  float fPar0() const { return fPar0_; };
  float fPar1() const { return fPar1_; };
  double weight(int i) const { return weights_[i]; };
  double corr(int i) const { return mC_[i]; };
  double cholesky(int i, int j) const; 
  void NoiseInit();
  double fShape(double);  
};


Pulse::Pulse(TString fnameshape, int nsmpl, int nfreq)
{
  nsamples_=nsmpl;
  nfreq_=nfreq;

  mC_ = new double [nsamples_];
  mL_ = new double [nsamples_*nsamples_];  
  weights_=new double [nsamples_];
  filePS_ = new TFile(fnameshape.Data());
  grPS_ = (TGraph*)filePS_->Get("PulseShape/grPulseShape");
  TTree *trPS = (TTree*)filePS_->Get("PulseShape/Tail");
  trPS->SetBranchAddress("timeMin",      &tMin_);
  trPS->SetBranchAddress("expAmplitude", &fPar0_);
  trPS->SetBranchAddress("expTime",      &fPar1_);
  trPS->GetEntry(0);

  // In-time sample is i=5
  
  for(int i=0; i<nsamples_; i++){
    double x = double( IDSTART + nfreq_ * i - WFLENGTH / 2);
    weights_[i] = fShape(x);
  }
  NoiseInit();
}


Pulse::~Pulse()
{
  free(weights_); weights_ = NULL;
  free(mC_); mC_ = NULL;
  free(mL_); mL_ = NULL;

}


double Pulse::cholesky(int i, int j) const 

{ 
  return mL_[(i*nsamples_+j)]; 
}


double Pulse::fShape(double x)
{
  if( grPS_ !=0 && x > 0.){
    if(x<800.){
      return grPS_->Eval(x);
    }else{
      return fPar0_ * exp( -x * fPar1_ );
    }
  }else{
    return 0.;
  }
}



void Pulse::NoiseInit()
{

  //cout<<"jugu"<<endl;
  for(int i=0; i<nsamples_; i++){
    double y = 1. - exp( -double(nfreq_ * i) / (sqrt(2.) * TAU));
    mC_[i] = 1. - y * y;
    //cout<<mC_[i]<<endl;
  }
  //cout<<"juga"<<endl;
  // initialize
  //cout<<"nsamples "<<nsamples_<<endl;

  //cout<<"pugu"<<endl;
  for(int i=0; i<nsamples_; ++i){
    //cout<<"i "<<i<<endl;
    for(int j=0; j<nsamples_; ++j){
      //cout<<"j "<<j<<endl;
      mL_[(i*nsamples_+j)]=0;
    }
  }
  //cout<<"jugi"<<endl;
  // decomposition
 
  mL_[0] = sqrt(mC_[0]);
  
  for( int col=1; col<nsamples_; col++){
    mL_[col]=0;
  } 

  //cout<<"pugi"<<endl;
  
  for( int row=1; row<nsamples_; row++)
    {
      for( int col=0; col<row; col++ )
	{
	  double sum1 = 0;
	  int m=abs(row-col);
	
	  for( int k=0; k<col; ++k) 
	    {
	      sum1 += mL_[row*nsamples_+k]*mL_[col*nsamples_+k];
	    }
	
	  mL_[row*nsamples_+col] = (mC_[m] - sum1)/mL_[col*nsamples_+col];
	}
     
      //cout<<"pugu"<<endl;
      double sum2 = 0;
    
      for( int k=0; k<row; ++k)
	{
	  sum2 += mL_[row*nsamples_+k]*mL_[row*nsamples_+k];
	}
    
      mL_[row*nsamples_+row] = sqrt( mC_[0] - sum2 );
    
      for( int col=row+1; col<nsamples_; col++ )
	{
	  mL_[row*nsamples_+col] = 0;
	}
    }
}


