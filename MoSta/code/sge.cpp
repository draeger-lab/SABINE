
#include "sge.h"

using namespace std;

// parameter:
//   aoclient: obejct with sge_merge function implemented
//   csid: string with ID for jobs
//   n: number of jobs (-1 if not known)
CSGEMaster::CSGEMaster(ISGEClient* aoclient,const char* csid,int n) : sid(csid)
{
  oclient = aoclient;
  if (n>=0)
    vids.reserve(n);
  //get random number
  srand(time(0));
  inr=rand();
  //create directory
  stringstream s;
  s << "sgetemp/" << sid << "/" << inr;
  sdir = s.str();
  //remove directory if it exists
  stringstream srmdir;
  srmdir << "rm -R -f " << sdir;
  system(srmdir.str().c_str());
  //create directory
  string smkdir("mkdir -p ");
  smkdir += sdir;
  system(smkdir.c_str());
}

//submits a job
// parameter:
//   scmd: string with command, results are expected to be printed
//         on the command line
void CSGEMaster::submit(string scmd)
{
  int i = vids.size();

  //create shell script
  stringstream sfsh;
  sfsh << sdir << "/i" << i << ".sh";
  ofstream fsh(sfsh.str().c_str());
  fsh << "#!/bin/bash\n";
  fsh << "#$ -e " << sdir << "/i" << i << ".err\n";
  fsh << "#$ -o " << sdir << "/i" << i << ".stat\n";
  fsh << "ulimit -c 0\nhostname\necho -n 'SGE_START='\ndate +%s\n";
  fsh << scmd << " > " << sdir << "/i" << i;
  fsh << ".out\necho -n 'SGE_STOP='\ndate +%s\n";
  fsh.close();

  //submit to SGE engine
  stringstream sfjob;
  sfjob << sdir << "/i" << i << ".job";
  stringstream ssub;
  ssub << "submit2sge -p -q hpc64 -f " << sdir << "/i" << i << ".sh > " << sfjob.str();
  system(ssub.str().c_str());

  //read output
  ifstream fjob(sfjob.str().c_str());
  string sojob;
  getline(fjob,sojob);
  getline(fjob,sojob);
  vector<string> vr;
  if (SplitString(sojob,string(" "),vr,false)<3)
    {
      cerr << "Job \"" << scmd << "\" could not be submitted.\n";
      exit(0);
    }
  cerr << "Job " << i << " started with ID " << vr[2] << endl;
  vids.push_back(str2long(vr[2]));
}

//waits until all jobs are finished
void CSGEMaster::finish() throw(EWrongFormat,ESGEError)
{
  string sfjob(sdir+"/jobs.stat");
  //check whether jobs are still running
  string scmd("rm -f "+sfjob+";qstat | grep -e $USER > "+sfjob);
  long ssum=0;
  while (ssum!=-vids.size())
    {
      wait(5);
      system(scmd.c_str());
      vector<long> vrun;
      vector<bool> vstat;
      vrun.reserve(vids.size());
      vstat.reserve(vids.size());
      //open file and get IDs and states
      ifstream fjob(sfjob.c_str());
      string s;
      // don't need to skip first two lines because of grep
      //      getline(fjob,s);
      //      getline(fjob,s);
      while (getline(fjob,s))
	{
	  if (s!="")
	    {
	      vector<string> vr;
	      SplitString(s,string(" "),vr,false);
	      if (vr.size()<5)
		throw EWrongFormat(string("qstat returns line with wrong format: "+s));
	      
	      vrun.push_back(str2long(vr[0]));
	      vstat.push_back((vr[4].substr(0,1)=="E"));
	    }
	}
      fjob.close();
      //check which jobs are finished
      ssum=0;
      for (int i=0;i<vids.size();i++)
	{
	  if (vids[i]>=0)
	    {
	      bool brun=false;
	      for (int j=0;j<vrun.size();j++)
		{
		  brun = brun | (vids[i]==vrun[j]);
		  if (brun)
		    {
		      //check for error flag!
		      if (vstat[j])
			{
			  stringstream serr;
			  serr << "Error in job " << i << " with ID " << vids[i];
			  throw ESGEError(serr.str());
			}
		      break;
		    }
		}
	      //if job is finished: read results
	      if (!brun)
		{
		  stringstream sfile;
		  sfile << sdir << "/i" << i << ".out";
		  oclient->sge_merge(sfile.str(),i);
		  cerr << "Job " << i << " with ID " << vids[i] << " finished\n";
		  vids[i] = -1;
		}
	    }
	  ssum += vids[i];
	}
    }
  //all jobs are finished
}
