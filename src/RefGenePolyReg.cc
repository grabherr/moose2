#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>

class Equation
{
public:
  Equation() {
    m_data.resize(3);
    for (int i=0; i<m_data.isize(); i++)
      m_data[i].resize(4);
    m_a = m_b = m_c = 0.;
  }

  void Set(int column, int row, double d) {
    (m_data[row])[column] = d;
    //cout << "Set: " << column << " " << row << " " << d << endl;
  }
  double Get(int column, int row) const {
    return (m_data[row])[column];
  }

  void Solve();
  
  double a() const {return m_a;}
  double b() const {return m_b;}
  double c() const {return m_c;}

private:
  svec < svec < double > > m_data; 
  double m_a;
  double m_b;
  double m_c;

 

};




void Equation::Solve() 
{
  double a0, a2, a3;

  double n = Get(1, 0);

  //cout << "n=" << n << endl;
  int i;

  // Make a 1
  Set(0, 0, Get(0, 0)/n);
  Set(1, 0, Get(1, 0)/n);
  Set(2, 0, Get(2, 0)/n);
  Set(3, 0, Get(3, 0)/n);

  n = Get(1, 1);
  //cout << "n=" << n << endl;
  Set(0, 1, Get(0, 1)/n);
  Set(1, 1, Get(1, 1)/n);
  Set(2, 1, Get(2, 1)/n);
  Set(3, 1, Get(3, 1)/n);

  n = Get(1, 2);
  //cout << "n=" << n << endl;
  Set(0, 2, Get(0, 2)/n);
  Set(1, 2, Get(1, 2)/n);
  Set(2, 2, Get(2, 2)/n);
  Set(3, 2, Get(3, 2)/n);

 
  a0 = Get(0, 0);
  a2 = Get(2, 0);
  a3 = Get(3, 0);

  for (i=0; i<4; i++)  
    Set(i, 0, Get(i, 0)-Get(i, 1));

  for (i=0; i<4; i++)  
    Set(i, 2, Get(i, 2)-Get(i, 1));
  

  n = Get(2, 0);
  Set(0, 0, Get(0, 0)/n);
  Set(1, 0, Get(1, 0)/n);
  Set(2, 0, Get(2, 0)/n);
  Set(3, 0, Get(3, 0)/n);
  
  n = Get(2, 2);
  Set(0, 2, Get(0, 2)/n);
  Set(1, 2, Get(1, 2)/n);
  Set(2, 2, Get(2, 2)/n);
  Set(3, 2, Get(3, 2)/n);

  double b0, b3;
  b0 = Get(0, 0);
  b3 = Get(3, 0);

  for (i=0; i<4; i++)  
    Set(i, 0, Get(i, 0)-Get(i, 2));

  m_c = -Get(0, 0)/Get(3, 0);

  m_b = -b0 -m_c*b3;
  m_a = -a0 - m_b*a2 - m_c* a3;
  
  
  cerr << "a=" << m_a << " b=" << m_b << " c=" << m_c << endl;

 
}




class Linear
{
public:
  Linear() {}

  double Value(const svec<double> & a, const svec<double> & avg) {
    int i;
    double sum = 0.;
    double div = 0.001;
    for (i=0; i<a.isize(); i++) {
      if (a[i] < 0.001 || avg[i] < 0.001)
	continue;
      sum += a[i]/avg[i];
      div += 1.;
    }
    double v = sum/div;
    double lim = 2.1; 
    if (v > lim)
      v = lim;
    if (v < 1/lim)
      v = 1/lim;
    return v;
  }

};





class Waypoints
{
public:
  Waypoints(int n) {
    m_data.resize(n);
  }

  void Add(const string & name, const string & tail) {
    m_names.push_back(name);
    m_tails.push_back(tail);
    svec<double> d;   
  }
  void AddData(int i, double val) {
    m_data[i].push_back(val);
  }
  
  void Print() const {
  }

  const svec<double> & Get(int i) const {return m_data[i];}

private:
  svec< svec < double > > m_data;
  svec<string> m_names;
  svec<string> m_tails;
};

double Sum(const svec<double> & x) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i];
  return d;
}
double Sum2(const svec<double> & x) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i]*x[i];
  return d;
}
double Sum3(const svec<double> & x) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i]*x[i]*x[i];
  return d;
}
double Sum4(const svec<double> & x) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i]*x[i]*x[i]*x[i];
  return d;
}
double Sum(const svec<double> & x, const svec<double> & y) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i]*y[i];
  return d;
}
double SumY2(const svec<double> & x, const svec<double> & y) {
  double d = 0;
  for (int i=0; i<x.isize(); i++)
    d += x[i]*y[i]*y[i];
  return d;
}

bool PolyReg(double & a, double & b, double & c, const svec<double> & x, const svec<double> & y) 
{
  double n = x.isize();

  //cout << "SUM " << Sum(x) << " " << Sum(y) << endl;
  
  Equation eq;
  eq.Set(0, 0, -Sum(y));
  eq.Set(1, 0, n);
  eq.Set(2, 0, Sum(x));
  eq.Set(3, 0, Sum2(x));

  eq.Set(0, 1, -Sum(x, y));
  eq.Set(1, 1, Sum(x));
  eq.Set(2, 1, Sum2(x));
  eq.Set(3, 1, Sum3(x));

  eq.Set(0, 2, -SumY2(y, x));
  eq.Set(1, 2, Sum2(x));
  eq.Set(2, 2, Sum3(x));
  eq.Set(3, 2, Sum4(x));

  eq.Solve();
  a = eq.a();
  b = eq.b();
  c = eq.c();

  //cout << "RESULT: " << a << " " << b << " " << c << endl;
  
  bool bSuccess = true;

  if (b < 0 || b > 3.) {
    a = 0.;
    c = 0.;
    Linear backoff;
    b = backoff.Value(x, y);
    cerr << "Fallback to a=" << a << " b=" << b << " c=" << c << endl;
    bSuccess = false;
  }

  //a = (Sum(y)*Sum2(x)-Sum(x)*Sum(x,y))/(n*Sum2(x)-Sum(x)*Sum(x));
  //b = (n*Sum(x,y)-Sum(x)*Sum(y))/(n*Sum2(x)-Sum(x)*Sum(x));
  
  //cout << "a=" << a << " b=" << b << endl;
  /*  int i;
  for (i=0; i<x.isize(); i++) {
    cout << x[i] << " " << y[i] << " " << a+b*x[i]+c*x[i]*x[i];
    double a1 = a-y[i];
    double inv = (-b + sqrt(b*b-4*a1*c))/(2.*a1);
    cout << " " << inv << endl;
 
    }*/

  return bSuccess;
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file", "");
  commandArg<string> wayCmmd("-w","waypoint gene file");
  commandArg<int> firstCmmd("-f","first column with data (0-based)",1);
  commandArg<int> lastCmmd("-l","last column with data (0-based)",-1);
  commandArg<bool> linCmmd("-linear","use linear regression",false);
  commandArg<bool> fCmmd("-force","force polynomial fit",false);
  commandArg<double> tCmmd("-t","threshold",1.);
  commandLineParser P(argc,argv);
  P.SetDescription("Normalizes a list of genes.");
  P.registerArg(fileCmmd);
  P.registerArg(wayCmmd);
  P.registerArg(firstCmmd);
  P.registerArg(lastCmmd);
  P.registerArg(linCmmd);
  P.registerArg(fCmmd);
  P.registerArg(tCmmd);
 
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string wayName = P.GetStringValueFor(wayCmmd);
  int first =  P.GetIntValueFor(firstCmmd);
  int last =  P.GetIntValueFor(lastCmmd);
  bool bLin =  P.GetBoolValueFor(linCmmd);
  bool bForce =  P.GetBoolValueFor(fCmmd);
  double thresh =  P.GetDoubleValueFor(tCmmd);

  Waypoints w(last-first+1);

  FlatFileParser parser1;
  
  parser1.Open(wayName);

  int i, j;
  
  svec<double> avg;
  
  while (parser1.ParseLine()) {
    if (parser1.GetItemCount() == 0)
      continue;
    if (last <= 0)
      last = parser1.GetItemCount()-2;

    bool bSkip = false;
    for (i=first; i<=last; i++) {
      if (parser1.AsFloat(i) < thresh)
	bSkip = true;
    }
    if (bSkip)
      continue;

    string trail;
    for (i=last+1; i<parser1.GetItemCount(); i++) {
      trail += parser1.AsString(i);
      trail += " ";
    }


    w.Add(parser1.AsString(0), trail);
    double m = 0.;
    double div = 0.;
    for (i=first; i<=last; i++) {
      w.AddData(i-first, log(1.+parser1.AsFloat(i)));
      m += log(1.+parser1.AsFloat(i));
      div += 1.;
    }
    //cout << m << " " << div << endl;
    avg.push_back(m/div);
  }

  int ref = 6;
  svec<double> a, b, c, l;
  a.resize(last-first+1);
  b.resize(last-first+1);
  c.resize(last-first+1);
  l.resize(last-first+1);

  FILE * p = fopen("factors.out", "w");

  Linear lin;
  for (i=0; i<a.isize(); i++) {     
    l[i] = lin.Value(w.Get(i), avg);
    if (!bLin) {
      bool bYes = PolyReg(a[i], b[i], c[i], w.Get(i), avg);
      if (!bYes && !bForce) {
	cout << "ERROR: Sample " << i << " failed! Use the -linear or -force option to proceed anyway." << endl;
	exit(-1);
      }
      //FILE * p = fopen("factors.out", "a");
      fprintf(p, "%f %f %f\n", a[i], b[i], c[i]);
      //fclose(p);
      //cerr << "OUCH!" << endl;
    } else {
      a[i] = 0.;
      if (l[i] > 0.)
	b[i] = 1./l[i];
      else
	b[i] = 1;
      c[i] = 0.;
      cerr << "Using linear model a=" << a[i] << " b=" << b[i] << " c=" << c[i] << endl;
      //FILE * p = fopen("factors.out", "a");
      fprintf(p, "%f %f %f\n", a[i], b[i], c[i]);
      //fclose(p);

    }
  }
  
  fclose(p);
 
  FlatFileParser parser;
  
  parser.Open(fileName);
  parser.ParseLine();
  cout << parser.Line() << endl;
  double bound = 1.;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    cout << parser.AsString(0);
    for (i=first; i<=last; i++) {
      double d = parser.AsFloat(i)+1.;
   
      double val;
      if (d >= bound) {
	double x, y;
	y = log(d);
 

	x = a[i-first]+b[i-first]*y+c[i-first]*y*y;
	val = exp(x)-1;
	
      } else {
	// Backoff to linear
	double bb = log(bound);
	double vv = a[i-first]+b[i-first]*bb+c[i-first]*bb*bb;
	vv = exp(vv)-1;
	val = d / bound * vv;
	//cerr << endl << "Before: " << d << " after " << val << " vv: " << vv << endl;
      }


      cout << "\t" << val;
    }

    for (i=last+1; i<parser.GetItemCount(); i++)
      cout << " " << parser.AsString(i);
 
    cout << endl;

  }
  
  return 0;
}
