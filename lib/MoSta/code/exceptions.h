#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <iostream>

//exception declaration
class EBase
{
 public:
  EBase(const std::string &ss) : s(ss) 
    {
      std::cerr << "Following Error occured: " << ss << std::endl;
    }
  EBase()
    {
      std::cerr << "Unkown Error occured.\n";
    }
  std::string s;
};

class EFileNotFound : public EBase
{
 public:
  EFileNotFound(const std::string &ss) : EBase(ss) {};
  EFileNotFound() {};
};

class EWrongFormat : public EBase
{
 public:
  EWrongFormat(const std::string &ss) : EBase(ss) {};
  EWrongFormat() {};
};

#endif
