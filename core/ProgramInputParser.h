#ifndef VISAR_PROGRAMINPUTPARSER
#define VISAR_PROGRAMINPUTPARSER

// http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c

#include <vector>
#include <string>

class ProgramInputParser
{
private:
  std::vector<std::string> tokens;

public:
  ProgramInputParser(int argc, char **argv);

public:
  bool option_id_exists(const std::string& option) const;
  std::string get_option(const std::string& option) const;
};

#endif
