#include "clustermatrix.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

CClusterMatrix::CClusterMatrix(string &asflist,double &gc,string &astmethod,double &at, bool absge, bool abloocv, bool abregularize) : CSimilarityMatrix(asflist,gc,astmethod,at,false,abregularize),bsge(absge),bloocv(abloocv)
{
}

//executes clustering: get similarity matrix and merge successivley.
// parameter:
//   p: percent quantile for smax_limit
void CClusterMatrix::do_cluster(double p)
{
  //call inherited method to fill matrix
  fill_matrix((int)bsge-2,true);

  //get highest similarity and quantiles
  double smax =0;
  int imax = -1;
  int jmax = -1;
  vector<double> vdistr;
  vdistr.reserve(n*(n-1)/2);
  for (int i=0;i<n;i++)
    {
      for (int j=i+1;j<n;j++)
	{
	  vdistr.push_back(m[i*nm+j]);
	  //save if maximum
	  if (m.at(i*nm+j)>smax)
	    {
	      smax = m[i*nm+j];
	      imax = i;
	      jmax = j;
	    }
	}
    }
  //sort distribution
  sort(vdistr.begin(),vdistr.end());
  //get p% quantile
  double smax_limit = max(vdistr[min((int)(p*vdistr.size()),(int)vdistr.size()-1)]-.00000001,0.0);
  cerr << "smax_limit=" << smax_limit << endl;

  //create file for matrix output
  string sfout(sflist);
  if (sflist.substr(0,5)=="list:")
    sfout = sflist.substr(5);
  stringstream sfmatrices;
  sfmatrices << sfout << ".matrices";
  ofstream fmatrices(sfmatrices.str().c_str());  

  //prepare member vector
  vector< vector<int> > vmembers;
  vmembers.reserve(n);
  int koffset = n;

  //FBPs for LOOCV
  vector<CPfm> vfbp_loocv;
  vfbp_loocv.reserve(n);
  for (int i=0;i<n;i++)
    vfbp_loocv.push_back(vopfm[i]);

  //print header
  cout << "matrixA\tmatrixB\tQA\tQB\ticA\ticB\tSmax\timax\tbimaxp\tQ\tic\n";
  //now merge maximum, get highest similarity and iterate k-1 times
  while (smax>smax_limit)
    {
      //merge maximum pair
      CPfm opfm_new(vopfm[imax],vopfm[jmax],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax]);

      cerr << "trying: " << vopfm[imax].sid << " with " << vopfm[jmax].sid << endl;

      //get members of new cluster
      vector<int> vm;
      if (imax<koffset)
	vm.push_back(imax);
      else
	{
	  for (int km=0;km<vmembers[imax-koffset].size();km++)
	    vm.push_back(vmembers[imax-koffset].at(km));
	}
      if (jmax<koffset)
	vm.push_back(jmax);
      else
	{
	  for (int km=0;km<vmembers[jmax-koffset].size();km++)
	    vm.push_back(vmembers[jmax-koffset].at(km));
	}
      
      //adjust matrix
      opfm_new.adjust_t(stmethod,t);
      //check similarity of members to new representative
      cerr << "\tmember\tSmax\n";
      int nmember = vm.size();
      vector<bool> vbs(nmember,false);
      #pragma omp parallel default(shared)
      {
        #pragma omp for
	for (int km=0;km<nmember;km++)
	  {
	    int i = vm[km];
	    CSimStat osstat(vopfm[i],opfm_new,igc);
	    cerr << "\t" << vopfm[i].sid << "\t" << osstat.smax << endl;
	    vbs[km] = osstat.smax>smax_limit;
	  }
      }
      //check whether criterium is fulfilled
      bool bnolimit=true;
      for (int km=0;km<nmember;km++)
	bnolimit = bnolimit && vbs[km];

      //check increase of information content
      if (bnolimit)
	{
	  cerr << "\t--> MERGED!\n";
	  //output
	  cout << vopfm[imax].sid << "\t" << vopfm[jmax].sid << "\t";
	  cout << vopfm[imax].get_beta() << "\t" << vopfm[jmax].get_beta() << "\t";
	  cout << vopfm[imax].get_ic() << "\t" << vopfm[jmax].get_ic() << "\t";
	  cout << smax << "\t" << mpos[imax*nm+jmax] << "\t" << mpos[jmax*nm+imax] << "\t";
	  
	  //add new matrix
	  int knew=n;
	  n++;
	  vopfm.push_back(opfm_new);
	  bused[knew] = true;

	  //add members
	  vmembers.push_back(vm);


	  //remove old matrices
	  bused[imax] = false;
	  bused[jmax] = false;

	  //output matrix into .matrices file
	  fmatrices << vopfm[knew];

	  //output beta, ic, ...
	  cout << vopfm[knew].get_beta() << "\t";
	  cout << vopfm[knew].get_ic() << endl;

	  //compute similarity to new pfm
          #pragma omp parallel default(shared)
	  {
            #pragma omp for
	    for (int i=0;i<knew;i++)
	      {
		if (bused[i])
		  {
		    CSimStat osstat(vopfm[i],vopfm[knew],igc);
		    m[i*nm+knew] = osstat.smax;
		    mpos[i*nm+knew] = osstat.imax;
		    mpos[knew*nm+i] = (int)osstat.bimaxp;
		  }
	      }
	  }

	  if (bloocv)
	    {
	      //update LOOCV FBPs
	      if (imax<koffset)
		{
		  CPfm opfm_fbp(vfbp_loocv[imax],vopfm[jmax],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax],true,false);
		  opfm_fbp.adjust_t(stmethod,t);
		  vfbp_loocv[imax] = opfm_fbp;
		}
	      else
		{
		  for (int km=0;km<vmembers[imax-koffset].size();km++)
		    {
		      int ii = vmembers[imax-koffset][km];
		      CPfm opfm_fbp(vfbp_loocv[ii],vopfm[jmax],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax]);
		      opfm_fbp.adjust_t(stmethod,t);
		      vfbp_loocv[ii] = opfm_fbp;
		    }
		}
	      if (jmax<koffset)
		{
		  CPfm opfm_fbp(vopfm[imax],vfbp_loocv[jmax],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax],false,true);
		  opfm_fbp.adjust_t(stmethod,t);
		  vfbp_loocv[jmax] = opfm_fbp;
		}
	      else
		{
		  for (int km=0;km<vmembers[jmax-koffset].size();km++)
		    {
		      int ii = vmembers[jmax-koffset][km];
		      CPfm opfm_fbp(vopfm[imax],vfbp_loocv[ii],mpos[imax*nm+jmax],(bool)mpos[jmax*nm+imax]);
		      opfm_fbp.adjust_t(stmethod,t);
		      vfbp_loocv[ii] = opfm_fbp;
		    }
		}
	    }

	} else
	{
	  //set smax to negative value since ic decreases
	  m[imax*nm+jmax] = smax_limit-1;
	}

      //NOW, find maximum
      smax = smax_limit;
      for (int i=0;i<n;i++)
	{
	  if (bused[i])
	    {
	      for (int j=i+1;j<n;j++)
		{
		  if (bused[j])
		    {
		      //save if maximum
		      if (m[i*nm+j]>smax)
			{
			  smax = m[i*nm+j];
			  imax = i;
			  jmax = j;
			}
		    }
		}
	    }
	}
    }
  fmatrices.close();

  //perform LOOCV
  if (bloocv)
    {
      int nsingle=0;
      int nright=0;
      int nwrong=0;
      int nright_nolimit=0;
      int nwrong_nolimit=0;
      //save wrong assignments
      string swrong("");

      //get for each PFM its cluster
      vector<int> vcluster(koffset,-1);
      for (int i=koffset;i<n;i++)
	{
	  for (int j=0;j<vmembers[i-koffset].size();j++)
	    vcluster[vmembers[i-koffset][j]] = i;
	}
      //      for (int i=0;i<koffset;i++)
      //	cout << vopfm[i].sid << "\t" << vcluster[i] << endl;
      //iterate over PFMs
      for (int i=0;i<koffset;i++)
	{
	  //ignore singletons
	  if (vcluster[i]>-1)
	    {
	      //compute similarity to own cluster
	      CSimStat osstat(vopfm[i],vfbp_loocv[i],igc);
	      double ownmax=osstat.smax;
	      cerr << vopfm[i].sid << "\t" << ownmax<< endl;
	      //compute all other similarities
	      int j=0;
	      double othermax=-1000;
	      while ((ownmax>othermax) && (j<n))
		{
		  if ((bused[j]) && (vcluster[i]!=j))
		    {
		      CSimStat oosstat(vopfm[i],vopfm[j],igc);
		      othermax=oosstat.smax;
		      cerr << "\t\t" << othermax << "\t" << vopfm[j].sid << endl;
		      if (othermax>=ownmax)
			{
			  cerr << vopfm[i].sid << " gets wrong classification: " << ownmax << " smaller than " << othermax << " with " << vopfm[j].sid << endl;
			  swrong += vopfm[i].sid + ":";
			}
		    }
		  j++;
		}
	      //still highest?
	      if (ownmax>othermax)
		{
		  nright++;
		  if (ownmax<smax_limit)
		    nright_nolimit++;
		} 
	      else
		{
		  nwrong++;
		  if (othermax<smax_limit)
		    nwrong_nolimit++;
		}
	    } 
	  else
	    {
	      nsingle++;
	    }
	}
      cout << "LOOCV\t" << swrong << "\t" << nright << "\t" << (nright+nwrong) << "\t"  << nsingle << "\t"  << nright_nolimit << "\t"  << nwrong_nolimit << "\t\t\t\t" << endl;
      cerr << "LOOCV:\t" << nright << " of " << (nright+nwrong) << " (" << ((double)nright/(nright+nwrong)*100) << "%) correctly classified, discarding " << nsingle << " singletons.\n\t" << nright_nolimit << "/" << nwrong_nolimit << " correctly/non-correctly classified do not reach the limit.\n";
    }
}

//prints transfac files for non-merged matrices
// parameter:
//   osout: output stream
void CClusterMatrix::print_pfms(ostream& osout)
{
  //iterate over matrices
  for (int i=0;i<bused.size();i++)
    {
      //print them if not merged
      if (bused[i])
	osout << vopfm[i];
    }
}
