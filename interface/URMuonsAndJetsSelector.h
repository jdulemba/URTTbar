#ifndef URMuonsAndJetsSelector_h
#define URMuonsAndJetsSelector_h



#include "URAnalysis/AnalysisFW/interface/URGenericMuonsAndJetsSelector.h"



typedef URGenericMuonsAndJetsSelector<Muon, &URStreamer::muons, Jet, &URStreamer::jets> URMuonsAndJetsSelector;



#endif // URMuonsAndJetsSelector_h
