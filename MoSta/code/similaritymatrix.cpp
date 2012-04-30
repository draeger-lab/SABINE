#include "similaritymatrix.h"

using namespace std;

CSimilarityMatrix::CSimilarityMatrix(string &asflist,double &gc,string &astmethod,double &at,bool bdiagonal, bool bregularize) : sflist(asflist),stmethod(astmethod),igc(gc),t(at),bfilldiagonal(bdiagonal), vopfm(asflist,gc,bregularize,(astmethod=="nrwords"))
{
  if (vopfm.size()<2)
    {
      cerr << "Error in similaritymatrix.cpp::CSimilarityMatrix(): Matrix file contains less than 2 PFMs.\n";
      exit(1);
    }
}

//call-back function to merge similarity results
void CSimilarityMatrix::sge_merge(const string &sfn,int i) throw (EFileNotFound,EWrongFormat)
{
  //open output file
  vector<string> vfnmat;
  ifstream f;
  f.open(sfn.c_str());
  if ((!f) || f.eof())
    throw EFileNotFound();
  //iterate and put into matrix n
  string s;
  int k=0;
  while (getline(f,s))
    {
      if (s!="")
	{
	  //split tabs
	  vector<string> v;
	  if (SplitString(s,string("\t"),v,false)<6)
	    {
	      int i = str2int(v[0]);
	      int j = str2int(v[1]);

	      double s = str2double(v[2]);
	      double ssum = str2double(v[3]);
	      int imax = str2int(v[4]);
	      int bimaxp = str2int(v[5]);
	      m[i*nm+j] = s;
	      m[j*nm+i] = ssum;
	      mpos[i*nm+j] = imax;
	      mpos[j*nm+i] = bimaxp;
	      k++;
	    } else 
	    {
	      throw EWrongFormat();
	    }
	}
    }
  int nt = n;
  if (i==(int)(n/2)-1)
    nt = n/2;
  if (k<nt)
    {
      stringstream iserr;
      iserr << "Not sufficient number of lines in file ";
      iserr << i << " with name " << sfn;
      iserr << ": " << k << " instead of " << nt << ".\n";
      throw EWrongFormat(iserr.str());
    }
}

//computes the intial similarity matrix
// parameter:
//    kpartial: <-1 if complete matrix, -1: use sge, otherwise row kpartial and n-kpartial
//    bfilltwice: make matrix double size for clustering
void CSimilarityMatrix::fill_matrix(int akpartial,bool bfilltwice)
{
  kpartial = akpartial;
  int pfmoffset = max(kpartial,0);
  //create PFM objects
  valpha.reserve(vopfm.size());
  for (int i=pfmoffset;i<vopfm.size();i++)
    {
      vopfm[i].adjust_t(stmethod,t);
      valpha.push_back(vopfm[i].get_alpha());
    }

  //get space for matrix and initialize vectors
  n = valpha.size();
  //how much memory do we need?
  nm = n;
  if (bfilltwice)
    nm *= 2;

  //allocate memory
  m.reserve(nm*nm); 
  mpos.reserve(nm*nm);//save ipos in i,j and bpalindrom in j,i for i>j
  bused.reserve(nm);
  
  //initialize
  for (int i=0;i<nm;i++) 
    {
      if (i<n)
	bused.push_back(true);
      else
	bused.push_back(false);
      for (int j=0;j<nm;j++)
	{
	  m.push_back(0);
	  mpos.push_back(0);
	}
    }

  //fill diagonal matrix
  int jplus = bfilldiagonal ? 0 : 1;
  //compute similarity matrix
  if (kpartial<0)
    {
      if (kpartial==-1)
	{
	  //use SGE engine to compute similarity
	  int nmax = ((int)n/2);
	  CSGEMaster osge(this,"smat",nmax);
	  for (int i=0;i<nmax;i++)
	    {
	      stringstream istr;
	      istr << "sstat ";
	      istr << igc << " " << sflist << " " << stmethod << " " << t;
	      istr << " " << i << " " << bfilldiagonal;
	      osge.submit(istr.str());
	    }
	  //wait until all jobs have been processed
	  osge.finish();
	} else 
	{
	  int nc = n;
	  int ncc = nm;
          #pragma omp parallel default(shared)
	  {
	    //compute similarities
            #pragma omp for
	    for (int i=0;i<nc;i++)
	      {
                #pragma omp parallel shared(i,nc)
		{
                  #pragma omp for
		  for (int j=i+jplus;j<nc;j++)
		    {
		      //create object for similarity measure
		      CSimStat osstat(vopfm[i+pfmoffset],vopfm[j+pfmoffset],igc);
		      m[i*ncc+j] = osstat.smax;
		      m[j*ncc+i] = osstat.ssum;
		      mpos[i*ncc+j] = osstat.imax;
		      mpos[j*ncc+i] = (int)osstat.bimaxp;
		    }
		}
	      }
	  }
	}
    } else //partial execution
    {
      //consider first row
      int i=0;
      for (int j=jplus;j<n;j++)
	{
	  //create object for similarity measure
	  CSimStat osstat(vopfm[i+pfmoffset],vopfm[j+pfmoffset],igc);
	  m[i*nm+j] = osstat.smax;
	  m[j*nm+i] = osstat.ssum;
	  mpos[i*nm+j] = osstat.imax;
	  mpos[j*nm+i] = (int)osstat.bimaxp;
	}
      //get complementary row from back
      i = n-kpartial-2;
      if (i >0)
	{
	  for (int j=i+jplus;j<n;j++)
	    {
	      //create object for similarity measure
	      CSimStat osstat(vopfm[i+pfmoffset],vopfm[j+pfmoffset],igc);
	      m[i*nm+j] = osstat.smax;
	      m[j*nm+i] = osstat.ssum;
	      mpos[i*nm+j] = osstat.imax;
	      mpos[j*nm+i] = (int)osstat.bimaxp;
	    }
	}
    }
}

//function to printout Similarity Matrix object
ostream& operator<<(ostream& strm,const CSimilarityMatrix& obj)
{
  int n = obj.n;
  int nm = obj.nm;
  //fill diagonal matrix
  int jplus = obj.bfilldiagonal ? 0 : 1;
  int pfmoffset = max(obj.kpartial,0);
  //SGE client output?
  if (obj.kpartial<0)
    {
      for (int i=0;i<n;i++)
	{
	  for (int j=i+jplus;j<n;j++)
	    {
	      strm << obj.vopfm[i+pfmoffset].sid << "\t" << obj.vopfm[j+pfmoffset].sid;
	      strm << "\t" << obj.m[i*nm+j] << "\t" << obj.m[j*nm+i];
	      strm << "\t" << obj.mpos[i*nm+j] << "\t" << obj.mpos[j*nm+i];
	      strm << "\t" << obj.valpha[i];
	      strm << "\t" << obj.valpha[j] << "\n";
	    }
	}
    } else 
    {
      int i=0;
      for (int j=jplus;j<n;j++)
	{
	  strm << obj.kpartial+i << "\t" << obj.kpartial+j;
	  strm << "\t" << obj.m[i*nm+j] << "\t" << obj.m[j*nm+i];
	  strm << "\t" << obj.mpos[i*nm+j] << "\t" << obj.mpos[j*nm+i] << "\n";
	}
      //get complementary row from back
      i = n-obj.kpartial-2;
      if (i >0)
	{
	  for (int j=i+jplus;j<n;j++)
	    {
	      strm << obj.kpartial+i << "\t" << obj.kpartial+j;
	      strm << "\t" << obj.m[i*nm+j] << "\t" << obj.m[j*nm+i];
	      strm << "\t" << obj.mpos[i*nm+j] << "\t" << obj.mpos[j*nm+i] << "\n";
	    }
	}
    }
  return strm;
}

