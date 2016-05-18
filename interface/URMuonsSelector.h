#ifndef URMuonsSelector_h
#define URMuonsSelector_h



#include "URAnalysis/AnalysisFW/interface/URGenericMuonsSelector.h"
// #include "URStreamer.h"



typedef URGenericMuonsSelector<Muon, &URStreamer::muons> URMuonsSelector;


#endif // URMuonsSelector_h
