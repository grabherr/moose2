#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"


int Run(const string & exec, const string & cmmd)
{
  string out = exec;
  //if (exec != "")
  //exec += "/";
  out += cmmd;
  cout << "Running " << out << endl;
  int ret = system(out.c_str());
  if (ret < 0) {
    cout << "ERROR: " << out << " died with exit code " << ret << endl;
    exit(-1);
  }
  return ret;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> repCmmd("-r","replicates file");
  commandArg<string> readsCmmd("-a","average read counts per sample", "0");
  commandArg<string> lenCmmd("-col","column that specifies the transcript length", "0");
  commandArg<string> wayCmmd("-w","waypoint gene file","");
  commandArg<string> firstCmmd("-f","first column with data (0-based)","1");
  commandArg<string> lastCmmd("-l","last column with data (0-based)","0");



  commandLineParser P(argc,argv);
  P.SetDescription("Runs a pipeline to re-estimate FPKM/RPKM values.");
  P.registerArg(fileCmmd);
  P.registerArg(repCmmd);
  P.registerArg(readsCmmd);
  P.registerArg(lenCmmd);
  P.registerArg(wayCmmd);
  P.registerArg(firstCmmd);
  P.registerArg(lastCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string repName = P.GetStringValueFor(repCmmd);
  string wayName = P.GetStringValueFor(wayCmmd);
  string col = P.GetStringValueFor(lenCmmd);
  string readCount = P.GetStringValueFor(readsCmmd);
  string first = P.GetStringValueFor(firstCmmd);
  string last = P.GetStringValueFor(lastCmmd);

  int i;

  const char * pExec = argv[0];
  cout << "Executing " << pExec << endl;

  char exec_dir[8192];
  strcpy(exec_dir, pExec);
  bool b = false;
  for (i = strlen(exec_dir)-1; i>=0; i--) {
    if (exec_dir[i] == '/') {
      exec_dir[i+1] = 0;
      b = true;
      break;
    }
  }
  if (!b)
    strcpy(exec_dir, "");

  
  //=============================================
  string cmmd = "RefGeneHMM -i ";
  cmmd += fileName;
  if (wayName != "") {
    cmmd += " -w ";
    cmmd += wayName;
  }
  cmmd += " -f ";
  cmmd += first;
  cmmd += " -l ";
  cmmd += last;

  cmmd += " > hmm_out";
  Run(exec_dir, cmmd);

  //=============================================
  cmmd = "RefGenePolyReg -w hmm_out -i ";
  cmmd += fileName;

  cmmd += " -f ";
  cmmd += first;
  cmmd += " -l ";
  cmmd += last;

  cmmd += " > normalized.out";
  Run(exec_dir, cmmd);

  
  //=============================================
  cmmd = "FitLaplaceNormal -i normalized.out ";
  cmmd += " -r ";
  cmmd += repName;
  cmmd += " -f ";
  cmmd += first;
  cmmd += " -l ";
  cmmd += last;

  cmmd += " > distribution.txt";
  Run(exec_dir, cmmd);


  //=============================================
  cmmd = "RefGeneSignificance -d distribution.txt -i ";
  cmmd += fileName;
  cmmd += " -r ";
  cmmd += repName;
  cmmd += " -a ";
  cmmd += readCount;
  cmmd += " -l ";
  cmmd += col;
 
  cmmd += " > significance.txt";
  Run(exec_dir, cmmd);

  cout << "SUCCESS: All done!!" << endl;


  return 0;
}
