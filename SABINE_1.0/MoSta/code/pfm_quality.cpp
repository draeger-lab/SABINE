#include "pfm_quality.h"

using namespace std;

// parameters:
//   opfm: PFM object to compute quality for
//   gc: gc content for background model
//   stmethod: threshold method, can also be 'optimize'
//   tp: threshold parameter
CPfmQuality::CPfmQuality(CPfm &opfm,double gc,string stmethod, double tp,bool bnr) : sid(opfm.sid)
{
  int n;
  //number of binding sites
  if (bnr)
    n = opfm.vbs.size();
  else
    n = opfm.vmbs.size();

  //get PFMs for binding site sequences
  vector<CPfm> vobs;
  //and weights
  vector<int> vw;
  vw.reserve(n);
  for (int i=0;i<n;i++)
    {
      if (bnr)
	{
	  CPfm obs(opfm.vbs[i],gc);
	  vobs.push_back(obs);
	  vw.push_back(opfm.vbs.vw.at(i));
	}
      else
	{
	  CPfm obs(opfm.vmbs[i],gc);
	  vobs.push_back(obs);
	  vw.push_back(opfm.vmbs.vw.at(i));
	}
    }

  //variance for binding sites
  double bsvar =0;
  double bsvar_nr =0;
  for (int i=0;i<n;i++)
    {
      for (int j=i;j<n;j++)
	{
	  //get variance
	  CSimStat obsvar(vobs[i],vobs[j],gc);
	  bsvar += (obsvar.ssum * (1+(int)(i!=j)))*vw[i]*vw[j];
	  bsvar_nr += obsvar.ssum * (1+(int)(i!=j));
	}
    }

  int tmax;
  if (stmethod=="optimize")
    {
      tmin = opfm.get_mint()+1;
      tmax = opfm.get_maxt();
    }
  else
    {
      opfm.adjust_t(stmethod,tp);
      tmin=opfm.get_t();
      tmax=tmin;
    }

  string st("threshold");
  bestt=tmin;

  //iterate over thresholds
  for (int t=tmin;t<=tmax;t++)
    {
      //      cerr << "t=" << t << endl;
      //adjust threshold for PFM
      opfm.adjust_t(st,(double)t);
      valpha.push_back(opfm.get_alpha());
      vbeta.push_back(opfm.get_beta());
      
      //check change, since, if no change quality won't change neither
      if ((t>tmin) && (valpha[t-tmin-1]==valpha[t-tmin]) && (vbeta[t-tmin-1]==vbeta[t-tmin]))
 	{
	  vquality.push_back(vquality[t-tmin-1]);
	  vquality_nr.push_back(vquality_nr[t-tmin-1]);
	  vsens.push_back(vsens[t-tmin-1]);
	  vspec.push_back(vspec[t-tmin-1]);
	  vprec.push_back(vprec[t-tmin-1]);
	}
      else
	{
	  //get variance
	  //for PFM
	  CSimStat ovar(opfm,opfm,gc);
	  
	  double q=0;
	  double q_nr=0;
	  double sw=0; //single word probability
	  double sj=0; //joint probability
	  //compute covariances
	  //iterate over given binding sites
	  for (int i=0;i<n;i++)
	    {
	      CSimStat osim(opfm,vobs[i],gc);
	      q += osim.ssum*vw[i];
	      q_nr += osim.ssum;
	      sw += vobs[i].get_alpha();
	      sj += osim.gABnull;
	    }
	  q /= pow((double)ovar.ssum*bsvar,1.0/2);
	  q_nr /= pow((double)ovar.ssum*bsvar_nr,1.0/2);
	  vquality.push_back(q);
	  vquality_nr.push_back(q_nr);
	  if (q>vquality[bestt-tmin])
	    bestt = t;
	  //sensitivity
	  vsens.push_back(sj/sw);
	  vspec.push_back((1.0-opfm.get_alpha()-sw+sj)/(1-sw));
	  vprec.push_back(sj/opfm.get_alpha());
	}
    }
}

//returns best threshold
int CPfmQuality::get_best_t()
{
  return(bestt);
}

//returns best quality
double CPfmQuality::get_best_q()
{
  return(vquality[bestt-tmin]);
}

//function to print results
ostream& operator<<(ostream& strm,const CPfmQuality& obj)
{
  for (int i=0;i<obj.vquality.size();i++)
    {
      strm << obj.sid << "\t"
	   << obj.vquality[i] << "\t"
	   << obj.vquality_nr[i] << "\t"
	   << obj.vsens[i] << "\t"
	   << obj.vspec[i] << "\t"
	   << obj.vprec[i] << "\t"
	   << i+obj.tmin << "\t"
	   << obj.valpha[i] << "\t"
	   << obj.vbeta[i] << endl;
    }
  return strm;
}

