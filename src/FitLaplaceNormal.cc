#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "src/LaplaceGauss.h"
#include "src/Replicates.h"


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> repCmmd("-r","replicate config file");
  commandArg<int> firstCmmd("-f","first column with data (0-based)",1);
  commandArg<int> lastCmmd("-l","last column with data (0-based)",-1);
  commandLineParser P(argc,argv);
  P.SetDescription("Computes a distribution based on FPKM values/replicates.");
  P.registerArg(fileCmmd);
  P.registerArg(repCmmd);
  P.registerArg(firstCmmd);
  P.registerArg(lastCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string repName = P.GetStringValueFor(repCmmd);
  int first =  P.GetIntValueFor(firstCmmd);
  int last =  P.GetIntValueFor(lastCmmd);

  Replicates rep;
  rep.ReadConfig(repName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;
  int k = 0;
  
  parser.ParseLine();
 
  for (i=0; i<parser.GetItemCount(); i++)
    rep.Configure(parser.AsString(i), i);

  svec<double> hist;
  hist.resize(300);
  double range = 3.;
  double scale = range/(double)(hist.isize()/2);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    svec<double> val;

     
    for (i=0; i<rep.Count(); i++) {
      val.clear();
      const svec<int> & c = rep.Get(i);
      //const svec<int> & c = rep.Get(i);
      //cout << "Rep " << i << endl;
      for (j=0; j<c.isize(); j++) {
	if (c[j] < 0) {
	  cout << "ERROR: Unassigned replicate!!! " << rep.Label(i, j) << endl;
	  return -1;
	}
	//cout << "Push " << c[j] << endl;
	val.push_back(parser.AsFloat(c[j]));
      }
      AddFold(hist, val, scale);
      val.clear();
    }
  }


  double sum = 0.;
  for (i=0; i<hist.isize(); i++) {
    sum += hist[i];
  }
  for (i=0; i<hist.isize(); i++) {
    hist[i] /= sum;
  }

  LaplaceGauss lg;
  double min = 999999999999.;
  double sigma = 0.;
  double epsilon = 0.;

  double start_sigma = 0.001;
  double end_sigma = 1.;
  double start_epsilon = 0.001;
  double end_epsilon = 1.;
  double inc = 0.1;
  double inc_init = 0.1;
  double n = 1.;
  while (inc > 0.0005) {
    
    for (double x=start_sigma; x<=end_sigma; x+= inc) {
      for (double y=start_epsilon; y<=end_epsilon; y+= inc) {     
	lg.SetSigma(x);
	lg.SetEpsilon(y);
	lg.Compute(range, hist.isize());
	double s = lg.Dist(hist);
	if (s < min) {
	  min = s;
	  sigma = x;
	  epsilon = y;
	}
      }
    }
    start_sigma = sigma - 2*inc;
    start_epsilon = epsilon - 2*inc;

    if (start_sigma < 0)
      start_sigma = 0.;
    if (start_epsilon < 0)
      start_epsilon = 0.;

    end_sigma = start_sigma + 4*inc;
    end_epsilon = start_epsilon + 4*inc;


    n += 1.;
    inc = inc_init/n;
    cerr << "sigma: " << sigma << " epsilon: " << epsilon << " inc: " << inc;
    cerr << endl; 
  }

  lg.SetSigma(sigma);
  lg.SetEpsilon(epsilon);
  lg.Compute(range, hist.isize());

  cout << "Sigma (Laplace): " << sigma << " Epsilon (Gauss): " << epsilon << endl;
  for (i=0; i<hist.isize(); i++) {
    cout << i-hist.isize()/2 << " " << hist[i] << " " << lg.F(i) << endl;
  }

  return 0;
}
