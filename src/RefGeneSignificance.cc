#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>
#include "src/LaplaceGauss.h"
#include "src/Replicates.h"

class CDFLaplaceGauss
{
public:
  CDFLaplaceGauss() {}

  void SetParams(double sigma, double epsilon) {
    LaplaceGauss lg;
    lg.SetSigma(sigma);
    lg.SetEpsilon(epsilon);
    
    int size = 10000;
    double range = 20.;
    svec<double> cdf;
    m_cdf.resize(size/2, 0.);
    lg.Compute(range, size);
    
    
    int i;
    double sum = 0.;
    for (i=size-1; i>=size/2; i--) {
      //double fold = range*(i-size/2);
      double f = lg.F(i);
      sum += f;
      m_cdf[i-size/2] = sum;
    }
    m_fold.resize(m_cdf.isize(), 0.);
    for (i=0; i<size/2; i++) {
      double fold = range*i/(double)(size/2);
      m_fold[i] = fold;
      m_cdf[i] = 0.5*m_cdf[i]/sum;    
      //cout << m_fold[i] << " " << m_cdf[i] << endl;
    }
  }

  double PVal(double fold)  {
    //int index = BinSearch(m_fold, fold);
    //cout << "Size: " << m_fold.isize() << " " << m_cdf.isize() << endl;

    int i;
    for (i=0; i<m_fold.isize(); i++) {
      //cout << "Search " << fold << " " << m_fold[i] << endl;
      if (fold > m_fold[i])
	continue;
      break;
    }
    //cout << "Fold=" << fold << " Index= " << i << endl;
    //cout << "cdf=" << m_cdf[i] << endl;
    return m_cdf[i];
  }


private:
  svec<double> m_fold;
  svec<double> m_cdf;
};

bool Good(const svec<double> & a)
{
  int i;
  bool b;
  for (i=0; i<a.isize(); i++) {
    if (a[i] > 2.)
      b = true;
  }
  return b;
}

double EValueUp(CDFLaplaceGauss & cdf, 
		const svec<double> & a, 
		const svec<double> & b,
		double corr) 
{
  if (!Good(a) && !Good(b))
    return 1.;

  int i, j;
  double n = 0.;
  double p = 1.;
  double off = 1.;

  for (i=0; i<a.isize(); i++) {
    for (j=0; j<b.isize(); j++) {
      double fold = log((b[j]+off)*corr) - log(a[i]+off);
      n += 1.;
      if (fold < 0)
	continue;
      double pval = 2*cdf.PVal(fold);
      //cout << "Fold: " << fold << " p=" << pval << endl;
      p *= pval;
    }
  }
  return p * n;
}

double Limit(double a)
{
  if (a > 1.)
    return 1.;
  return a;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> repCmmd("-r","replicates file");
  commandArg<string> distCmmd("-d","distribution file", "");
  commandArg<double> sigmaCmmd("-s","sigma (Laplace)", 0.);
  commandArg<double> epsilonCmmd("-e","epsilon (Gauss)", 0.);
  commandArg<int> readsCmmd("-a","average read counts per sample", 0);
  commandArg<int> lenCmmd("-l","column that specifies the transcript length", -1);
  commandLineParser P(argc,argv);
  P.SetDescription("Computes sigificance for differential expression based on a Laplace/Normal distribution.");
  P.registerArg(fileCmmd);
  P.registerArg(repCmmd);
  P.registerArg(sigmaCmmd);
  P.registerArg(epsilonCmmd);
  P.registerArg(distCmmd);
  P.registerArg(readsCmmd);
  P.registerArg(lenCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string distName = P.GetStringValueFor(distCmmd);
  string repName = P.GetStringValueFor(repCmmd);
  double sigma = P.GetDoubleValueFor(sigmaCmmd);
  double epsilon = P.GetDoubleValueFor(epsilonCmmd);
  double total = P.GetIntValueFor(readsCmmd);
  int lencol = P.GetIntValueFor(lenCmmd);

  if (total == 0)
    total = 5000000;

  Replicates rep;
  rep.ReadConfig(repName);

  if (distName == "" && sigma == 0. && epsilon == 0.) {
    cout << "ERROR: you must specify sigma & epsilon OR a distribution file!" << endl;
    return -1;
  }
  if (distName != "") {
    FlatFileParser parser2;  
    parser2.Open(distName);
    while (parser2.ParseLine()) {
      if (parser2.AsString(0) == "Sigma") {
	sigma = parser2.AsFloat(2);
	epsilon = parser2.AsFloat(5);
	break;
      }
    }
  }

  //cout << "sigma " << sigma << " epsilon " << epsilon << endl;
 
  // Count genes
  FlatFileParser parser1;  
  parser1.Open(fileName);
  parser1.ParseLine();
  double n = 0.;
  while (parser1.ParseLine()) {
    if (parser1.GetItemCount() == 0)
      continue;
    n += 1.;
  }

 
  CDFLaplaceGauss cdf;

  cdf.SetParams(sigma, epsilon);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;
  int k = 0;
  
  parser.ParseLine();
  //cout << parser.Line() << endl;

  for (i=0; i<parser.GetItemCount(); i++)
    rep.Configure(parser.AsString(i), i);

  n *= rep.Count();

  svec<double> hist;
  hist.resize(300);
  double range = 3.;
  double scale = range/(double)(hist.isize()/2);


  //cout << "Gene\t0vs30+\t0vs30-\t0vs90+\t0vs90-\t30vs90+\t30vs90-" << endl;
  cout << "Gene";
  for (i=0; i<rep.Count(); i++) {
    for (j=i+1; j<rep.Count(); j++) {
      cout << "\t" << rep.Name(i) << "_vs_" << rep.Name(j) << "+";
      cout << "\t" << rep.Name(i) << "_vs_" << rep.Name(j) << "-";
    }
  }
  cout << endl;


  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    cout << parser.AsString(0);
 

    for (i=0; i<rep.Count(); i++) {
      for (j=i+1; j<rep.Count(); j++) {
    

	const svec<int> & one = rep.Get(i);
	const svec<int> & two = rep.Get(j);

	svec<double> val1, val2;
	int x;
	
	double avg = 0.;
	double div = 0.;
	for (x=0; x<one.isize(); x++) {
	  val1.push_back(parser.AsFloat(one[x]));
	  avg += parser.AsFloat(one[x]);
	  div += 1.;
	}
	for (x=0; x<two.isize(); x++) {
	  val2.push_back(parser.AsFloat(two[x]));
	  avg += parser.AsFloat(two[x]);
	  div += 1.;
	}

	avg /= div;

	double corr = 1.;
	if (total > 0.) {
	  double len = 1;
	  if (lencol >0)
	    len = parser.AsFloat(lencol)/1000;
	  double reads = avg*len*total/1000000;
	  
	  corr = 1.-sqrt(reads)/reads;
	  //cout << endl << "reads " << reads << " avg " << avg << " len " << len << " corr: " << corr << endl;
	}
	//cout << endl << "Testing " << val1.isize() << " " << val2.isize() << " " << n << " " << corr << endl;
	cout << "\t" << Limit(EValueUp(cdf, val1, val2, corr)*n);
	cout << "\t" << Limit(EValueUp(cdf, val2, val1, corr)*n);
       }
    }
    cout << endl;
  }


  return 0;
}
