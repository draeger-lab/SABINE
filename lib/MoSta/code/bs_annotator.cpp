

#include "bs_annotator.h"

using namespace std;

CBSAnnotator::CBSAnnotator(CPfm& opfm)
{
  //get parameters for annotation
  mpssm = opfm.mpssm;
  t = opfm.get_t();
  nlen = opfm.nlen;
  //and some for the output
  spfmid = opfm.sid;

  //make map to get indices for letters
  dcompl['A'] = 0;
  dcompl['C'] = 1;
  dcompl['G'] = 2;
  dcompl['T'] = 3;
  dcompl['N'] = 4;
}

int CBSAnnotator::annotate(const CSequence& oseq)
{
  //annotate 5'-3' strand
  _annotate(oseq,false);

  //annotate 3'-5' strand
  _annotate(oseq,true);
  return(vsid.size());
}

// annotates sequence
// parameter:
//   oseq: CSequence object to annotate
//   breverse: annotate reverse strand?
void CBSAnnotator::_annotate(const CSequence& oseq,bool breverse)
{
  //get sequence
  string sseq=oseq.get_sequence(breverse);
  string sseqid = oseq.get_sid();
  long int nseq = sseq.length();

  //sequence at least as long as motif?
  if (nseq>=nlen)
    {
      for (int i=0;i<nseq-nlen+1;i++)
	{
	  int iscore=0;
	  for (int j=0;j<nlen;j++)
	    {
	      int ibase = dcompl[sseq[i+j]];
	      if (ibase==4)
		{
		  i +=j;
		  iscore = t-1;
		  break;
		}
	      iscore += mpssm[j*4+ibase];
	    }
	  if (iscore>=t)
	    addhit(sseqid,i,nseq,breverse);
	}
    }
}

void CBSAnnotator::addhit(string& sid,long int ipos,long int n,bool breverse)
{
  //consider palindromic hits
  if (breverse)
    ipos = n-ipos;
  int istrand = 1;
  if (breverse)
      istrand = -1;

  //calibrate coordinates
  long int ipos_start = n-ipos;
  long int ipos_end = n-(ipos+istrand*(nlen-1));
  long int iseq_start = ipos;
  long int iseq_end = ipos+istrand*(nlen-1);

  //change coordinates for palindromic hits
  if (breverse)
    {
      //swap position
      long int k = ipos_start;
      ipos_start = ipos_end+1;
      ipos_end = k+1;
      //swap sequence position
      k = iseq_start;
      iseq_start = iseq_end-1;
      iseq_end = k-1;
    }

  //save hit
  vsid.push_back(sid);
  vstrand.push_back(istrand);
  vpos_start.push_back(ipos_start);
  vpos_end.push_back(ipos_end);
  vseq_start.push_back(iseq_start);
  vseq_end.push_back(iseq_end);
}

string CBSAnnotator::get_header()
{
  return(string("matrix\tgene\tstrand\tpos.start\tpos.end\tseq.start\tseq.end\n"));
}

//funtion to printout PFM in transfac format
ostream& operator<<(ostream& strm,const CBSAnnotator& obj)
{
  for (int i=0;i<obj.vsid.size();i++)
    {
      strm << obj.spfmid << "\t" << obj.vsid[i] << "\t" << obj.vstrand[i] << "\t";
      strm << obj.vpos_start[i]-1 << "\t" << obj.vpos_end[i]-1 << "\t";
      strm << obj.vseq_start[i]+1 << "\t" << obj.vseq_end[i]+1 << endl;
    }
  return strm;
}
