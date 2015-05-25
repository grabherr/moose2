#ifndef REPLICATES_H
#define REPLICATES_H

#include <string>
#include "base/FileParser.h"

class Replicates
{
public:
  Replicates() {}

  void ReadConfig(const string & fileName) {
    FlatFileParser parser;
  
    parser.Open(fileName);
    int i;
    while (parser.ParseLine()) {
      if (parser.GetItemCount() == 0)
	continue;
      svec<string> cond;
      svec<int> col;
      m_name.push_back(parser.AsString(0));
      for (i=1; i<parser.GetItemCount(); i++) {
	cond.push_back(parser.AsString(i));
	col.push_back(-1);
      }
      m_labels.push_back(cond);
      m_rep.push_back(col);
    }
    
  }
  void Configure(const string & label, int column) {
    int i, j;
    for (i=0; i<m_labels.isize(); i++) {
      const svec<string> & l = m_labels[i];
      svec<int> & c = m_rep[i];
      for (j=0; j<l.isize(); j++) {
	if (l[j] == label) {
	  c[j] = column;
	}
      }
    }
  }

  int Count() const {return m_rep.isize();}
  const svec<int> & Get(int i) {return m_rep[i];}
  const string & Name(int i) const {return m_name[i];}
  const string & Label(int i, int j) const {
    return (m_labels[i])[j];
  }
private:
  svec < svec < int > > m_rep;
  svec < svec < string > > m_labels;
  svec<string> m_name;
};

#endif 
