#ifndef LAPLACEGAUSS_H
#define LAPLACEGAUSS_H

#include "base/SVector.h"
#include <string>
#include <math.h>


#define PI 3.1415926535897932384626433832795

double Fold3(double a, double b) { 
  double f = log(a) - log(b);
  return f;
}

void AddFold(svec<double> & hist, svec<double> & val, double scale) {
  int i, j;
  for (i=0; i<val.isize(); i++) {
    for (j=0; j<val.isize(); j++) {
      if (i == j)
	continue;
      if (val[i] < 2. || val[j] < 2.)
	continue;
      double f = Fold3(val[i], val[j])/scale;
      int index = (int)(f+hist.isize()/2+0.5);
      if (index < 0 || index >= hist.isize()) {
	//cout << "OUT OF RANGE: " << val[i] << " " << val[j] << " " << f << endl;
      } else {
	hist[index] += 1.;
      }
 
    }
  }
}


double Gauss(double center, double sigma, double x)
{
  double f = (x-center);
  if (f < 0.)
    f = -f;

  
  double f1 = x-center+1.;
  if (f1 < 0.)
    f1 = -f1;
  double f2 = x-center-1.;
  if (f2 < 0.)
    f2 = -f2;
  

  double y = exp(-f/sigma);
  double y1 = exp(-f1/sigma);
  double y2 = exp(-f2/sigma);
 //double y = sigma/((x-center)*(x-center) + sigma);
  //cout << x << " " << center << " " << sigma << " " << y << endl;

  return (y + y1 + y2)/3.;
}


class LaplaceGauss
{
public:
  LaplaceGauss() {
    m_sigma = m_epsilon = 1.;
    m_range = 1.;
  }
  
  void SetSigma(double sigma) {
    m_sigma = sigma;
  }
  void SetEpsilon(double epsilon) {
    m_epsilon = epsilon;
  }

  // Range is plus/minus
  void Compute(double range, int size) {
    m_data.clear();
    m_data.resize(size, 0.);
    m_range = range;

    int i, j;
    for (i=0; i<m_data.isize(); i++) {
      int ii = i-m_data.isize()/2;
      //if (ii != 0)
      //continue;
      //for (i=m_data.isize()/2; i<=m_data.isize()/2; i++) {
      double x = m_range*(i-m_data.isize()/2)/(double)(m_data.isize()/2);
      double f = 0.;
      double n = 0.;
      //cout << "i=" << i <<  " x=" << x << endl;
      for (j=0; j<m_data.isize(); j++) {
	int jj = j-m_data.isize()/2;
	//if (jj != 0)
	//continue;

	double y = x-m_range*(j-m_data.isize()/2)/(double)(m_data.isize()/2);
	//cout << "  y=" << y << endl;
	int index = j;
	//if (index >= 0 && index <= m_data.isize()) {
	m_data[index] += Laplace(y)*Gauss(x);
	  //cout << "  apply=" << y << " index=" << index << endl;
	//cout << index << " x=" << x << " y=" << y << " " << Laplace(x)*Gauss(y) << endl;
	  n++;
	  //}
      }
      
      //m_data[i] = f/n;
    }
    
    double sum = 0.;
    for (i=0; i<m_data.isize(); i++)
      sum += m_data[i];
    for (i=0; i<m_data.isize(); i++)
      m_data[i] /= sum;    

    
  }

  double F(double x) const {
    int index = (int)(m_data.isize()/2 + x*m_data.isize()/2/m_range + 0.5);
    return m_data[index];
  }

  double F(int i) const {
    return m_data[i];
  }
  
  double ErrfNum(double x) const {
    int index = (int)(m_data.isize()/2 + x*m_data.isize()/2/m_range + 0.5);
    double p = 0.;
    for (int i=index; i<m_data.isize(); i++) {
      p += m_data[i];
    }
    return p;
  }

  double Dist(const svec<double> & val) {
    double sum = 0.;
    for (int i=0; i<m_data.isize(); i++) {
      sum += (m_data[i] - val[i])*(m_data[i] - val[i]);
    }
    return sqrt(sum);
  }

private:
  double Gauss(double x) {
    //if (x < 0.1 && x > -0.1)
    //  return 1.;
    //return 0;

    double f = exp(-0.5*x*x/(m_epsilon*m_epsilon))/(m_epsilon*sqrt(2*PI));
    return f;
  }

  double Laplace(double x) {
    if (x < 0.)
      x = -x;
    double f = exp(-x/m_sigma)/(2*m_sigma);
    return f;
    
  }

  double m_sigma;
  double m_epsilon;
  double m_range;
  svec<double> m_data;
};





#endif //LAPLACEGAUSS_H

