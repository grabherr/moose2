#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include <math.h>

double Norm(double d) 
{
  if (d < 0)
    return -d;
  return d;
}


class NormVal
{
public:
  NormVal() {
    m_norm = 0.;
    m_val = 0.;
  }
  double Norm() const {return m_norm;}
  double Value() const {return m_val;}
  void SetNorm(double d) {m_norm = d;}
  void SetVal(double d) {m_val = d;}

  
  bool operator < (const NormVal & v) const {
    return m_val < v.m_val;
  }

private:
  double m_norm;
  double m_val;
};

class RefGene
{
public:
  RefGene() {
  }

  void Add(double val, double norm) {
    NormVal tmp;
    tmp.SetNorm(norm);
    tmp.SetVal(val);
    m_val.push_back(tmp);
  }
  

  void Done() {
    Sort(m_val);
    //for (int i=0; i<m_val.isize(); i++) {
    //cout << "Print " << i << "\t" << m_val[i].Value() << " " << m_val[i].Norm() << endl;
    //}
  }

  double Normalize(double r) const {
    int i;
    if (r < m_val[0].Value())
      return r/m_val[0].Norm();
    if (r >= m_val[m_val.isize()-1].Value())
      return r/m_val[m_val.isize()-1].Norm();

    for (i=0; i<m_val.isize()-1; i++) {
      if (log(r+1) >= log(m_val[i].Value()+1.) && log(r+1) <= log(m_val[i+1].Value()+1))
	break;
    }
       

    //double d1 = r-m_val[i].Value();
    //double d2 = m_val[i+1].Value()-r;
    //double div = m_val[i+1].Value()-m_val[i].Value();
    
    double d1 = log(r+1)-log(m_val[i].Value()+1);
    double d2 = log(m_val[i+1].Value()+1)-log(r+1);
    double div = log(m_val[i+1].Value()+1)-log(m_val[i].Value()+1);
    
    //cout << endl;
    //cout << "NORM(1): " << m_val[i].Value() << " " << r << " " << m_val[i+1].Value() << endl;

    d1 /= div;
    d2 /= div;
    double f = m_val[i].Norm()*(1.-d1) + m_val[i+1].Norm()*(1.-d2);
    //double w = 1. - 1./div;
    
    //cout << endl << "Norm: " << r << " " << m_val[i].Value() << " - " << m_val[i+1].Value() << " " << f << endl;
    //cout << "NORM(2): " << d1 << " " << d2 << " " << f << " " << r << endl;
    double cf = 1.0;
    return cf * r/f + (1.-cf)*r;
  } 

private:
  svec<NormVal> m_val;
 
};


class NormSignal
{
public:
  NormSignal() {
  }

  void Setup(int n, double max, const RefGene & ref) {
    m_table.resize(n, 0.);
    m_max = log(max);
    int i;
    double step = m_max/(double)m_table.isize();
    double x = 0.;
    for (i=0; i<m_table.isize(); i++) {
      double y = exp(x);
      m_table[i] = ref.Normalize(y)/y;
      x += step;
    }
    //Smooth();
  }

  double NormFactor(double v) const {
    if (v <= 1.)
      return m_table[0];
    double step = m_max/(double)m_table.isize();
    int i = (int)(log(v)/step + 0.5);
    if (i >= m_table.isize())
      return m_table[m_table.isize()-1];
    return m_table[i];
  }
private:
  void Smooth() {
    svec<double> tmp;
    tmp.resize(m_table.isize(), 0.);
    int i, j;
    int win = 50;
    for (j=0; j<50; j++) {
      tmp = m_table;
       
      for (i=win; i<m_table.isize()-win; i++) {
	double n = 0.;
	double d = 0.;
	for (int k=i-win; k<=i+win; k++) {
	  d += m_table[k];
	  n += 1.;
	}
	tmp[i] = d/n;
	// tmp[i] = 0.2*m_table[i] + 0.8*(m_table[i-1] + m_table[i+1])/2.;
	//tmp[i] = (m_table[i] + m_table[i-1] + m_table[i+1] + m_table[i-2] + m_table[i+2])/5.;
      }
      m_table = tmp;
    }
  }


  svec<double> m_table;
  double m_max;
};


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> refCmmd("-r","reference file");
  commandLineParser P(argc,argv);
  P.SetDescription("Normalizes FPKM/RPKM values by reference genes.");
  P.registerArg(fileCmmd);
  P.registerArg(refCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string refName = P.GetStringValueFor(refCmmd);
  

  //comment. ???
  FlatFileParser parserR;
  
  parserR.Open(refName);
  int i;
  svec<RefGene> ref;

  double hi = 0.;
  while (parserR.ParseLine()) {
    if (parserR.GetItemCount() == 0)
      continue;
    if (ref.isize() == 0)
      ref.resize(parserR.GetItemCount()-2);
    svec<double> val, norm;
    double sum = 0.;
    double nn = 0.;
    for (i=1; i<parserR.GetItemCount()-1; i++) {
      val.push_back(parserR.AsFloat(i));
      norm.push_back(parserR.AsFloat(i));
      nn += 1.;
      sum += parserR.AsFloat(i);
    }
    for (i=0; i<norm.isize(); i++) {
      norm[i] = nn*norm[i]/sum;
      ref[i].Add(val[i], norm[i]);
      if (val[i] > hi)
	hi = val[i];
    }
    
  }

  svec<NormSignal> norm;
  norm.resize(ref.isize());

  for (i=0; i<ref.isize(); i++) {
    ref[i].Done();
    norm[i].Setup(5000, hi, ref[i]);
  }

  /*
  double max = 0.1 + log(hi);
  double step = max/1000.;
  for (double x=0; x<max; x += step) {
    cout << x << " ";
    for (i=0; i<norm.isize(); i++) {
      //cout << exp(x)/ref[i].Normalize(exp(x)) << " ";
      cout << 1./norm[i].NormFactor(exp(x)) << " ";
    }
    cout << endl;
  }
  return 0;
  */
  
  FlatFileParser parser;
  parser.Open(fileName);
  parser.ParseLine();
  cout << parser.Line() << endl;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    cout << parser.AsString(0);
    double sum = 0.;
    double nnn = 0.;
    for (i=1; i<parser.GetItemCount()-2; i++) {
      double orig = parser.AsFloat(i);
      double dd = parser.AsFloat(i)*norm[i-1].NormFactor(parser.AsFloat(i));
      //dd = (dd + orig)/2.;
      cout << "\t" << dd;
      sum += dd;
      nnn += 1.;
      //cout << "\t" << ref[i-1].Normalize(parser.AsFloat(i)); // + parser.AsFloat(i))/2.;
    }

    cout << "\t" << parser.AsString(parser.GetItemCount()-2) << "\t" << sum/nnn << endl;
  }
  return 0;
}
