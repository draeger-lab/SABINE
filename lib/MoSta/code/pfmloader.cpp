
#include "pfmloader.h"
#include <cstdlib>

using namespace std;

//constructor to load one or more PFMs from one or more files
// parameter:
//   sfn: filename (if prefix list: we assume filenames in sfn)
//   gc: gc content
//   bregularize: regularize with Rahmann 2005 method?
CPfmLoader::CPfmLoader(string sfn,double gc,bool bregularize,bool bsequences)
{
  if (sfn.substr(0,5)=="list:")
    {
      string sflist(sfn.substr(5));
      ifstream f;
      f.open(sflist.c_str());
      if (!f)
	{
	  cerr << "Error in pfmloader.cpp::CPfmLoader(): Cannot open file " 
	       << sflist 
	       << " for input.\n";
	  exit(1);
	}
      string s;
      while (getline(f,s)) 
	{
	  if (s!="")
	    push_back(CPfm(s,gc,bregularize,bsequences));
	}
      f.close();
    }
  else
    {
      ifstream f;
      f.open(sfn.c_str());
      while (!f.eof())
	{
	  push_back(CPfm(f,gc,bregularize,bsequences));
	  //remove empty lines
	  string s("");
	  int p;
	  p=f.tellg();
	  while (!(f.eof() || (s.length()>0)))
	    {
	      p=f.tellg();
	      getline(f,s);
	      if (s.length()>0)
		chomp(s);
	    }
	  if (!f.eof())
	    f.seekg(p);
	}
      f.close();
    }
  
}
