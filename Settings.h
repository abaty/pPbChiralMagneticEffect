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

  static const int trkEtaGaps = 10;
  float etaGaps[trkEtaGaps+1] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};

  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting setting.." << std::endl;
  return;
}

#endif
