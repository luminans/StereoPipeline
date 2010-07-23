// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include "multiview.h"

using std::cout;
using std::endl;
using std::string;

namespace fs = boost::filesystem;
using namespace vw;
using namespace vw::cartography;

int main( int argc, char *argv[] ) {
  Options opts = parse_opts(argc, argv);
  cout << opts << endl;
  return 0;
}
