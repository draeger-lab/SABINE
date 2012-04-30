
#include "sequences.h"

using namespace std;

//initialize static variables

//create map for accepted letters
pair<char, char> a[] =
{
pair<char, char>('A', 'A'),
pair<char, char>('C', 'C'),
pair<char, char>('G', 'G'),
pair<char, char>('T', 'T'),
pair<char, char>('N', 'N'),
pair<char, char>('X', 'N'),
pair<char, char>(' ', 'N'),
pair<char, char>('-', 'N'),
};
map<char, char> CSequence::dtrans(a, a + sizeof(a) / sizeof(a[0]));

//create map for complement
pair<char, char> b[] =
{
pair<char, char>('A', 'T'),
pair<char, char>('C', 'G'),
pair<char, char>('G', 'C'),
pair<char, char>('T', 'A'),
pair<char, char>('N', 'N'),
};
map<char, char> CSequence::dcompl(b, b + sizeof(b) / sizeof(b[0]));

//definition of base transformation
// parameter:
//   cc: character to transform
//   brev: use the complement?
char CSequence::basetrans(char cc,bool brev) throw(EWrongFormat)
{
  char c = toupper(cc);
  //make map to get complement letters
 
  char cr='Z';

  if (c!='\\')
    {
      map<char, char>::iterator cur  = dtrans.find(c);

      //is letter allowed?
      if (cur == dtrans.end())
	{
	  stringstream serr("");
	  serr << "Sequence contains a prohibited letter: ";
	  serr << "'" << c << "'.\n";
	  throw(EWrongFormat(serr.str()));
	}
      cr = (*cur).second;
      
      if (brev)
	cr = dcompl[cr];
    }
  return(cr);
}

CSequence::CSequence(string& asid,const string& asseq) throw(EWrongFormat) : sseq(""),srseq(""),sid(asid),gc(0)
{
  
  //get gc content
  long int ngc=0;
  long int nat=0;

  //transform strings
  for (int i=0;i<asseq.length();i++)
    {
      
      char c = basetrans(asseq[i],false);
      if (c!='Z')
	sseq += c;

      char cr = basetrans(asseq[asseq.length()-i-1],true);
      if (c!='Z')
	srseq += cr;
      
      if ((c=='G') || (c=='C'))
	ngc++;
      if ((c=='A') || (c=='T'))
	nat++;
    }

  gc = (double)ngc/(ngc+nat);
}

string CSequence::get_sequence(bool breverse) const
{
  if (breverse)
    {
      return(srseq);
    } else 
    {
      return(sseq);
    }
}

double CSequence::get_gc()
{
  return(gc);
}

string CSequence::get_sid() const
{
  return(sid);
}

long int CSequence::get_length() const
{
  return(sseq.length());
}

//class CSequences

CSequences::CSequences(char* sf) throw(EFileNotFound,EWrongFormat) : fin(sf),beof(true)
{
  //could we open the file?
  if (!fin)
    throw EFileNotFound(string("CSequences::CSequences: File "+string(sf)+" could not be opened.\n"));
  
  //read first id
  string s;
  if (getline(fin,s))
    {
      if (s.substr(0,1)!=">")
	throw EWrongFormat(string("CSequences::CSsequences: File "+string(sf)+" is not in correct FASTA format - not starting with '>'.\n"));
      //name of next gene
      sid = s.substr(1);
      chomp(sid);
      beof=false;
    } else throw EWrongFormat(string("CSequences::CSequences: File "+string(sf)+" seems to be empty.\n"));  
}

CSequences::~CSequences()
{
  fin.close();
}

bool CSequences::eof()
{
  return(beof);
}

//make CSequence object containing next sequence
CSequence CSequences::next()
{
  //read next sequence from file
  string s;
  stringstream sseq (stringstream::in | stringstream::out);
  while (getline(fin,s) && (s.substr(0,1)!=">"))
    {
      chomp(s);
      sseq << s;
    }

  //create sequence object
  CSequence o(sid,sseq.str());

  //check whether we are at the end of the file
  if (s.substr(0,1)==">")
    {
      sid = s.substr(1);
      chomp(sid);
    } else 
    {
      beof = true;
    }
  return(o);
}

//adds a sequence to sequence set if not yet in the set
// parameter:
//   oseq: object containing sequence and possibly also Ns
void CNRSequences::add(CSequence oseq)
{
  int nmax = 100;

  string sseq(oseq.get_sequence());
  //create sequences without ambiguous letters
  vector<string> vs;
  vs.push_back("");
  for (int i=0;i<sseq.length();i++)
    {
      
      if (sseq[i]=='N')
	{
	  int nn = vs.size();
	  for (int j=0;j<nn;j++)
	    {
	      vs.push_back(vs[j]+"C");
	      vs.push_back(vs[j]+"G");
	      vs.push_back(vs[j]+"T");
	      vs[j] += "A";
	    }
	}
      else
	{
	  for (int j=0;j<vs.size();j++)
	    {
	      vs[j] += sseq[i];
	    }
	}
      if (vs.size()>nmax)
	break;
    }
  if (vs.size()<=nmax)
    {
      //only add non-redundant sequences
      for (int i=0;i<vs.size();i++)
	{
	  bool bfound=false;
	  for (int j=0;j<size();j++)
	    {
	      if (vs[i]==at(j).get_sequence())
		{
		  bfound = true;
		  vw[j]++;
		  break;
		}
	    }
	  if (!bfound)
	    {
	      string ssid(oseq.get_sid()+".1");
	      CSequence ooseq(ssid,vs[i]);
	      push_back(ooseq);
	      vw.push_back(1);
	    }
	  if (size()>nmax)
	    break;
	}
    }
}


//adds a sequence to sequence set if not yet in the set (where N=N)
// parameter:
//   oseq: object containing sequence and possibly also Ns
void CMSequences::add(CSequence oseq)
{
  string sseq(oseq.get_sequence());
  bool bfound=false;
  for (int j=0;j<size();j++)
    {
      if (sseq==at(j).get_sequence())
	{
	  bfound = true;
	  vw[j]++;
	  break;
	}
    }
  if (!bfound)
    {
      string ssid(oseq.get_sid()+".1");
      CSequence ooseq(ssid,sseq);
      push_back(ooseq);
      vw.push_back(1);
    }
}


//funtion to printout sequence in fasta format
ostream& operator<<(ostream& strm,const CSequence& obj)
{
  strm << ">" << obj.sid << endl;
  strm << obj.get_sequence() << endl;
  return strm;
}

