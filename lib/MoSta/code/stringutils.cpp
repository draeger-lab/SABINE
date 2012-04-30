#include "stringutils.h"

using namespace std;

//special trims a string (take first non-space part till next space)
void strim(string& s)
{
  vector<string> vr;
  if (SplitString(s,string(" "),vr,false)>0)
    s = vr[0];
}

//chomps a string
void chomp(string& s)
{
  while (s.substr(s.length()-1,1)=="\r")
    s = s.substr(0,s.length()-1);
}

int str2int(std::string &s)
{
  istringstream iss(s);
  int k;
  iss >> k;
  return(k);
}

double str2double(std::string &s)
{
  istringstream iss(s);
  double k;
  iss >> k;
  return(k);
}

long str2long(std::string &s)
{
  istringstream iss(s);
  long k;
  iss >> k;
  return(k);
}

int SplitString(const string& input, 
		const string& delimiter, vector<string>& results, 
		bool includeEmpties)
{
    int iPos = 0;
    int newPos = -1;
    int sizeS2 = (int)delimiter.size();
    int isize = (int)input.size();

    if( 
        ( isize == 0 )
        ||
        ( sizeS2 == 0 )
    )
    {
        return 0;
    }

    vector<int> positions;

    newPos = input.find (delimiter, 0);

    if( newPos < 0 )
    { 
        return 0; 
    }

    int numFound = 0;

    while( newPos >= iPos )
    {
        numFound++;
        positions.push_back(newPos);
        iPos = newPos;
        newPos = input.find (delimiter, iPos+sizeS2);
    }

    if( numFound == 0 )
    {
        return 0;
    }

    for( int i=0; i <= (int)positions.size(); ++i )
    {
        string s("");
        if( i == 0 ) 
        { 
	  s = input.substr( i, positions.at(i) ); 
        } else
	  {
	    int offset = positions.at(i-1) + sizeS2;
	    if( offset < isize )
	      {
		if( i == positions.size() )
		  {
		    s = input.substr(offset);
		  }
		else if( i > 0 )
		  {
		    s = input.substr( positions.at(i-1) + sizeS2, 
				      positions.at(i) - positions.at(i-1) - sizeS2 );
		  }
	      }
	  }
        if( includeEmpties || ( s.size() > 0 ) )
        {
            results.push_back(s);
        }
    }
    return numFound;
}
