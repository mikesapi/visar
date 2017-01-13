#include "ProgramInputParser.h"

#include <algorithm>

ProgramInputParser::ProgramInputParser(int argc, char **argv)
{
  for(int i = 1; i < argc; ++i)
  {
    tokens.push_back(argv[i]);
  }
}

bool ProgramInputParser::option_id_exists(const std::string& option) const
{
  std::vector<std::string>::const_iterator it = std::find(tokens.begin(), tokens.end(), option);

  if(it != tokens.end()) return true;
  else return false;
}

std::string ProgramInputParser::get_option(const std::string& option) const
{
  std::vector<std::string>::const_iterator it = std::find(tokens.begin(), tokens.end(), option);
  if(it != tokens.end())
  {
    if((++it != tokens.end()) && (it->compare(0, 1, "-", 1) != 0))
    {
      return *it;
    }
    else return "";
  }
  else return "";
}
