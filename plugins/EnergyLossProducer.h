#ifndef EnergyLossProducer_H
#define EnergyLossProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace std;

namespace edm { class Run; class Event; class EventSetup; }
class TrackerGeometry;

class TFile;
class TH2F;

class EnergyLossProducer : public edm::EDProducer
{
public:
  explicit EnergyLossProducer(const edm::ParameterSet& ps);
  ~EnergyLossProducer();
  virtual void produce(edm::Event& ev, const edm::EventSetup& es);

private:
  void beginRun(edm::Run & run, const edm::EventSetup& es);
  void endJob();

  string trackProducer;
  double pixelToStripMultiplier, pixelToStripExponent;
  const TrackerGeometry * theTracker;

  TFile * resultFile;
  TH2F * hnor;
};
#endif

