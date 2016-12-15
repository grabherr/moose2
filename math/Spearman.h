#ifndef SPEARMAN_H
#define SPEARMAN_H

#include "base/SVector.h"
#include <math.h>



class DoublePair
{
 public:
  DoublePair() {
    m_x = m_y = m_z = 0.;
  }
  DoublePair(double x, double y, double z) {
    m_x = x;
    m_y = y;
    m_z = z;
  }

  void Set(double x, double y, double z) {
    m_x = x;
    m_y = y;
    m_z = z;
  }

  double X() const {return m_x;}
  double Y() const {return m_y;}
  double Z() const {return m_z;}

  void SetZ(double z) {
    m_z = z;
  } 
  void SetX(double x) {
    m_x = x;
  } 
  void SetY(double y) {
    m_y = y;
  } 

  bool operator < (const DoublePair & d) const {
    return m_x < d.m_x;
  }


 private:
  double m_x;
  double m_y; 
  double m_z;
};


class SpearmansRho
{
 public:
  SpearmansRho() {
    m_sig = 0.;
  }

  double Compute(const svec<double> & x, const svec<double> &y) {
    int i, j;
    svec<DoublePair> data;
    data.resize(x.isize());
    svec<DoublePair> dataInv;
    dataInv.resize(x.isize());
    
    svec<double> rank;
    svec<double> d;
    rank.resize(data.size(), 0.);
    d.resize(data.isize(), 0.);

    for (i=0; i<data.isize(); i++) {
      data[i].Set(x[i], y[i], i);
      dataInv[i].Set(y[i], x[i], i);
    }

    Sort(dataInv);

    for (i=0; i<dataInv.isize(); i++) {
      int k = i;
      int n = i;
      while (i+1<dataInv.isize() && dataInv[i+1].X() == dataInv[i].X()) {
	i++;
	n+=i;
      }

      //cout << "i=" << i << " n=" << n << " k=" << k << endl;
      for (j=k; j<=i; j++) {
	double val = (double)(n)/(double)(i-k+1);
	//cout << "Set j=" << j << "  " << val << endl;
	dataInv[j].SetY(val);
      }
    }

    
    for (i=0; i<dataInv.isize(); i++) {
      int index = (int)(dataInv[i].Z()+0.5);
      data[index].SetZ(dataInv[i].Y());
    }


    Sort(data);

    for (i=0; i<data.isize(); i++) {
      //cout << data[i].X() << "\t" << data[i].Y() << "\t" << rank[i];
      //cout << "\t" << data[i].Z() << endl;

      rank[i] = i;
      int k = i;
      int n = i;
      while (i+1<data.isize() && data[i+1].X() == data[i].X()) {
	i++;
	n+=i;
      }
      for (j=k; j<=i; j++) 
	rank[j] = (double)(n)/(double)(i-k+1);
    }

    double sum = 0.;
    for (i=0; i<d.isize(); i++) {
      //cout << data[i].X() << "\t" << data[i].Y() << "\t" << rank[i];
      //cout << "\t" << data[i].Z() << "\t";
      double diff = rank[i] - data[i].Z();
      double sq = diff*diff;
      sum += sq;
      //cout << diff << "\t" << sq << endl;

    }
    double n = d.isize();
    double rho = 1. - 6*sum/(n*(n*n-1));

    double r = rho;
    if (r < 0.)
      r = -rho;

    double z = 0.5*log((1+rho)/(1-rho))*sqrt((n-3.)/1.06);
    if (z < 0.)
      z = -z;
    
    //cout << "z=" << z << " " << z/sqrt(2.) << endl;
    
    m_sig = erfc(z/sqrt(2.));
    
    //cout << "Sig=" << m_sig << endl;

    //m_sig = rho*sqrt((n-2.)/(1-rho*rho));

    return rho;
  }

  double Significance() const {return m_sig;}

 private:
  double m_sig;
};

//=================================================
class PearsonsRho
{
 public:
  PearsonsRho() {
    m_sig = 0.;
  }

  double Compute(const svec<double> & x, const svec<double> &y) {
    int i, j;
    double n = x.isize();

    double xm = 0.;
    double ym = 0.;

    for (i=0; i<x.isize(); i++) {
      xm += x[i];
      ym += y[i];
    }
    xm /= n;
    ym /= n;

    double a = 0.;
    double b = 0.;
    double c = 0.;

    for (i=0; i<x.isize(); i++) {
      a += (x[i]-xm)*(y[i]-ym);
      b += (x[i]-xm)*(x[i]-xm);
      c += (y[i]-ym)*(y[i]-ym);
    }
    double rho = a / (sqrt(b)*sqrt(c));
  
    double z = 0.5*log((1+rho)/(1-rho))*sqrt((n-3.)/1.06);
    if (z < 0.)
      z = -z;
    
    m_sig = erfc(z/sqrt(2.));
    
    return rho;
  }

  double Significance() const {return m_sig;}

 private:
  double m_sig;
};


#endif


