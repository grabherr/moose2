#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include <math.h>

#define CRAZY_COST 99999999999.
#define REWARD 4.
#define VIOLATION_PENALTY 5.0

//#define REWARD 10.
//#define VIOLATION_PENALTY 3.2

class Node
{
public:
  Node() {
    m_score = CRAZY_COST;
    m_back = -1;
    m_reward = REWARD;
    m_dist = 0.;
    m_bYes = false;
    m_flat = 1.;

    m_violation = VIOLATION_PENALTY;
  }

  void SetViolation(double v) {
    m_violation = v;
  }

  void SetName(const string & name) {
    m_name = name;
  }
  const string & Name() const {return m_name;}

  void SetWaypoint() {    
    m_reward = REWARD * 200;
  }
  void SetTrail(const string & name) {
    m_trail = name;
  }
  const string & Trail() const {return m_trail;}

  void SetHead(const string & name) {
    m_head = name;
  }
  const string & Head() const {return m_head;}

  int isize() const {return m_pos.isize();}
  double & operator [] (int i) {return m_pos[i];}
  const double & operator [] (int i) const {return m_pos[i];}
  void resize(int n) {m_pos.resize(n, 0.);}
  void ComputeDist() {
    m_norm.resize(m_pos.isize(), 0.);
    m_dist = 0.;
    double sum = 0.;
    int i;
    for (i=0; i<m_pos.isize(); i++) {
      m_dist += m_pos[i]*m_pos[i];
      sum += m_pos[i];      
    }
    m_dist = sqrt(m_dist/(double)m_pos.isize());
    m_flat = 0.;
    double avg = 1. / (double)m_pos.isize();
    sum /= (double)m_norm.isize();
    for (i=0; i<m_norm.isize(); i++) {
      m_norm[i] = m_pos[i]/sum;
      m_flat += (m_norm[i]-avg)*(m_norm[i]-avg);
    }
  }
  double Norm(int i) const {return m_norm[i];}

  bool operator < (const Node & n) const {
    return m_dist < n.m_dist;
  }


  double Score() const {return m_score;}
  int Back() const {return m_back;}

  void SetScore(double d) {m_score = d;
    m_bYes = true;
  }
  
  bool IsGood() const {return m_bYes;}
  void Reset() {m_bYes = false;}

  bool Minimize(const Node & n, int index) {
    int i;
    double dist = 0.;
    /* for (i=0; i<isize(); i++) {
      //if (n[i] >= m_pos[i])
      //return CRAZY_COST;
      double a = log(n[i]+1.);
      double b = log(m_pos[i]+1.);
      double dd = (a-b)*(a-b);
      dist += dd;
      }
    // Cost
    

    */

    for (i=0; i<isize(); i++) {
      dist += (m_norm[i] - n.Norm(i))*(m_norm[i] - n.Norm(i));
     }
    //dist = 15.*sqrt(dist/(double)isize());
    dist = 15.*dist/(double)isize();

    for (i=0; i<isize(); i++) {
      if (m_pos[i] < n.m_pos[i])
	dist += m_violation;
    }

    //if (m_reward > 100)
    //cout << "Considering waypoint." << endl;

    double score = n.Score() + dist - m_reward /*+ m_flat * 0.1*/;
    //score -= log(m_dist+1.)/25.;
    //score -= log(m_dist+1.)/25.;
    if (score < m_score) {
      m_score = score;
      m_back = index;
      m_bYes = true;
      return true;
    }

    return false;
  }

  void Print() const {
    cout << m_name << "\t";
    if (m_head != "")
      cout << m_head << "\t";

    for (int i=0; i<m_pos.isize(); i++)
      cout << m_pos[i] << "\t";
    cout << m_trail;

    //cout << " d=" << m_dist;
    cout << endl;
  }
  
private:
  double m_dist;
  double m_score;
  int m_back;
  svec<double> m_pos;
  svec<double> m_norm;
  double m_reward;
  string m_name;
  string m_trail;
  string m_head;
  bool m_bYes;
  double m_flat;
  double m_violation;
};


void LoadWaypoint(svec<string> & way, const string & fileName)
{
  FlatFileParser parser;
  parser.Open(fileName);
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    way.push_back(parser.AsString(0));
  }
}


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> wayCmmd("-w","waypoint gene file","");
  commandArg<int> firstCmmd("-f","first column with data (0-based)",1);
  commandArg<int> lastCmmd("-l","last column with data (0-based)",-1);
  commandArg<double> penCmmd("-p","violation penalty",VIOLATION_PENALTY);
  commandLineParser P(argc,argv);
  P.SetDescription("Prepares a list of normalized reference genes for RPKM normalization.");
  P.registerArg(fileCmmd);
  P.registerArg(wayCmmd);
  P.registerArg(firstCmmd);
  P.registerArg(lastCmmd);
  P.registerArg(penCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string wayName = P.GetStringValueFor(wayCmmd);

  int first =  P.GetIntValueFor(firstCmmd);
  int last =  P.GetIntValueFor(lastCmmd);
  double pen =  P.GetDoubleValueFor(penCmmd);


  svec<string> waypoints;
  if (wayName != "")
    LoadWaypoint(waypoints, wayName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i, j;
  svec<Node> nodes;

  parser.ParseLine();
  int lo = -1;
  int hi = -1;
  double lowscore = 99999999.;
  int hiscore = 0.;
    

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (last <= 0)
      last = parser.GetItemCount()-2;

    Node tmp;
    tmp.SetName(parser.AsString(0));

    string head;
    for (i=1; i<first; i++) {
      head += parser.AsString(i);
      head += " ";
    }
    string trail;
    for (i=last+1; i<parser.GetItemCount(); i++) {
      trail += parser.AsString(i);
      trail += " ";
    }

    bool bWay = false;
    for (i=0; i<waypoints.isize(); i++) {
      if (tmp.Name() == waypoints[i]) {
	bWay = true;
      }
    }
    if (bWay)
      tmp.SetWaypoint();
    
 
    tmp.resize(last-first+1);

    tmp.SetHead(head);
    tmp.SetTrail(trail);

    double s = 0.;
    bool bSkip = false;
    for (i=first; i<=last; i++) {
      tmp[i-first] = parser.AsFloat(i);
      s += tmp[i-first];
      if (tmp[i-first] < 1.)
	bSkip = true;
    }

    /*
    if (bWay) {
      if (bSkip) 
	cout << "Skip waypoint (1) " << tmp.Name() << endl;
      if (s < tmp.isize()) 
	cout << "Skip waypoint (2) " << tmp.Name() << endl;
	}*/

    if (s < tmp.isize() || bSkip)
      continue;

    if (s < lowscore) {
      lowscore = s;
      lo = nodes.isize();
    }
    if (s > hiscore) {
      hiscore = s;
      hi = nodes.isize();
    }
    tmp.SetViolation(pen);
    tmp.ComputeDist();
    nodes.push_back(tmp);     
  }

  
  nodes[lo].SetScore(0);
  hi = nodes.isize()-1;
  lo = 0;

  Sort(nodes);
  //cout << "Lo: " << lo << " " << nodes[lo].Score() << endl;
  //nodes[lo].Print();
  int k;
  bool bCont = true;
  for (j=0; j<nodes.isize(); j++) {
    //cout << "Iteration " << j << " of " << nodes.isize() << endl;
    if (!bCont)
      break;

    bCont = false;
    for (i=j; i<nodes.isize(); i++) {
    
      //cout << "NODE: " << i << endl;
      for (k=i+1; k<nodes.isize(); k++) {
	Node & f = nodes[i];
	Node & t = nodes[k];
	if (i != k) {
	  //cout << "FEED: " << i << endl;
	  if (t.Minimize(f, i)) {
	    bCont = true;
	  }
	}
	//cout << "k=" << k << " score: " << t.Score() << endl;
      }
    }
  }

  nodes[hi].Print();
  k = nodes[hi].Back();

  do {
    //cout << "Score " << nodes[k].Score() << " back: " << nodes[k].Back() << endl;
    nodes[k].Print();
    k = nodes[k].Back();
  } while (k != -1);
  

  
  
  return 0;
}
