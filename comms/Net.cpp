#include "Net.h"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const net::Address& s)
{
  os << int(s.GetA()) << '.' << int(s.GetB()) << '.' << int(s.GetC()) << '.' << int(s.GetD()) << ':' << s.GetPort();
  return os;
}
