#ifndef RAAINPUTSETTINGS
#define RAAINPUTSETTINGS

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Settings{
  public:
  int nTrkMin = 185;
  int nTrkMax = 220;
  float ptMin = 0.4;
  float ptMax = 3; 

  float trkEtaCut = 2.4;
  float HFetaMin = 4.4;
  float HFetaMax = 5;

  float trkEtaGapMin = 0;
  float trkEtaGapMax = 1;

  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting setting.." << std::endl;
  return;
}

#endif
