#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VZero/interface/VZero.h"
#include "DataFormats/VZero/interface/VZeroFwd.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//#include "DataFormats/TrackReco/interface/TrackDeDxEstimate.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "QCDAnalysis/ChargedHadronSpectra/interface/Histograms.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"

#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"

// HF
//#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCounts.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include <fstream>
using namespace std;
using namespace reco;

/*****************************************************************************/
class QCDTrackAnalyzer : public edm::EDAnalyzer
{
 public:
   explicit QCDTrackAnalyzer(const edm::ParameterSet& pset);
   ~QCDTrackAnalyzer();
   virtual void beginJob(const edm::EventSetup& es);
   virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
   virtual void endJob();

 private:
   int getDetLayerId(const PSimHit& simHit);
   bool isAccepted(const TrackingParticle& simTrack);

   bool isPrimary(const edm::RefToBase<reco::Track> & recTrack);
   edm::RefToBase<reco::Track> getAssociatedRecTrack
     (const TrackingParticleRef & simTrack, int & nRec, bool hasToBePrimary);
   TrackingParticleRef getAssociatedSimTrack
     (const edm::RefToBase<reco::Track> & recTrack, int & nSim);

   float refitWithVertex(const reco::Track & recTrack);

   int processSimTracks();

   float getEnergyLoss(const reco::TrackRef & track);

   double truncate(double m);
   double getSigmaOfLogdEdx(double logde);
   double getLogdEdx(double bg);
   bool isCompatibleWithdEdx(double m, const reco::TrackRef & track);

   float getInvariantMass(const VZero & vZero, float m1, float m2);

   float getInvariantMass(const reco::Track * r1,
                          const reco::Track * r2,
                          float m1, float m2);

   void processVZeros();

   int processRecTracks(
     edm::Handle<reco::DeDxDataValueMap> elossCollection);

   const TrackerGeometry * theTracker;
   const TrackAssociatorByHits * theAssociatorByHits;
   const TransientTrackBuilder * theTTBuilder;
   const reco::BeamSpot* theBeamSpot;

//   TFile * resultFile; 
   TNtuple * trackSim;
   TNtuple * trackRec;
   TNtuple * vzeroRec;
   TNtuple * eventInfo;

   edm::Handle<edm::View<reco::Track> >    recCollection;
   edm::Handle<TrackingParticleCollection> simCollection;
   const reco::VZeroCollection * vZeros;
   const reco::DeDxDataValueMap * energyLoss;

   const reco::VertexCollection * vertices;

   reco::SimToRecoCollection simToReco;
   reco::RecoToSimCollection recoToSim;

   Histograms * histograms;

   string trackProducer;
//   string resultFileLabel;
   bool hasSimInfo;
   bool allRecTracksArePrimary;

   int proc, ntrk,nvtx;

   int ch;

   edm::View<reco::Track> oTrackCollection;
};

/*****************************************************************************/
QCDTrackAnalyzer::QCDTrackAnalyzer(const edm::ParameterSet& pset)
{
  trackProducer   = pset.getParameter<string>("trackProducer");
  hasSimInfo      = pset.getParameter<bool>("hasSimInfo");
  allRecTracksArePrimary = pset.getParameter<bool>("allRecTracksArePrimary");

  histograms = new Histograms(pset);
}

/*****************************************************************************/
QCDTrackAnalyzer::~QCDTrackAnalyzer()
{
}

/*****************************************************************************/
void QCDTrackAnalyzer::beginJob(const edm::EventSetup& es)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();
 
  // Get associator
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  es.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",
                                    theHitsAssociator);
  theAssociatorByHits =
   (const TrackAssociatorByHits*)theHitsAssociator.product();

  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  theTTBuilder = builder.product();
 
  // Root
  histograms->declareHistograms();
}

/*****************************************************************************/
void QCDTrackAnalyzer::endJob()
{
  histograms->writeHistograms();
}

/*****************************************************************************/
int QCDTrackAnalyzer::getDetLayerId(const PSimHit& simHit)
{
  int layerId;

  unsigned int id = simHit.detUnitId();

  if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel)
  {
    PXBDetId pid(id);
    layerId = pid.layer() - 1;
  }
  else
  {
    PXFDetId pid(id);
    layerId = 2 + pid.disk();
  }

  return layerId;
}

/*****************************************************************************/
bool QCDTrackAnalyzer::isAccepted
  (const TrackingParticle& simTrack_)
{
  TrackingParticle * simTrack = const_cast<TrackingParticle *>(&simTrack_);

  // How many pixel hits?
  const int nLayers = 5;
  vector<bool> filled(nLayers,false);

  std::vector<PSimHit> trackerPSimHit( simTrack->trackPSimHit(DetId::Tracker));
  
  for(std::vector<PSimHit>::const_iterator
        simHit = trackerPSimHit.begin();
        simHit!= trackerPSimHit.end(); simHit++)
  {

    if(simHit == trackerPSimHit.begin())
    if(simHit->particleType() != simTrack->pdgId())
      return false;

    unsigned int id = simHit->detUnitId();

    if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel ||
       theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelEndcap)
      filled[getDetLayerId(*simHit)] = true;
  }
  
  // Count the number of filled pixel layers
  int fLayers = 0;
  for(int i=0; i<nLayers; i++)
    if(filled[i] == true) fLayers++;
  
  // FIXME, should be in different layers, at least 2/3
  if(fLayers >= 3) return true;
              else return false;
}

/*****************************************************************************/
bool QCDTrackAnalyzer::isPrimary(const edm::RefToBase<reco::Track> & recTrack)
{
  if(allRecTracksArePrimary)
  {
    if(vertices->size() > 0)
    {
    // Look for the closest vertex in z
    float dzmin = -1.;
    for(reco::VertexCollection::const_iterator
        vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
    {
      float dz = fabs(recTrack->vertex().z() - vertex->position().z());
      if(vertex == vertices->begin() || dz < dzmin)
       dzmin = dz ;
    }

    // FIXME
    if(dzmin > 0.6) return false; // !!!
               else return true;
    }

    return true;
  }

  // Transverse impact paramter (0.2 cm or 5*sigma)
  double dt = fabs(recTrack->dxy(theBeamSpot->position()));
  double st = sqrt(recTrack->dxyError() * recTrack->dxyError() +
                  theBeamSpot->BeamWidth() * theBeamSpot->BeamWidth());

  if(dt > min(0.2, 5 * st)) return false;

  // Longitudinal impact paramter (0.2 cm or 5*sigma)
  // but only if there are vertices
  if(vertices->size() > 0)
  {
    // Look for the closest vertex in z
    float dzmin = -1.;
    for(reco::VertexCollection::const_iterator
        vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
    {
      float dz = fabs(recTrack->vertex().z() - vertex->position().z());
      if(vertex == vertices->begin() || dz < dzmin)
       dzmin = dz ;
    }

    if(dzmin > 0.3) return false;
// !!!   if(dzmin > min(0.2, 5 * recTrack->dzError())) return false;
  }

  return true;
}

/*****************************************************************************/
edm::RefToBase<reco::Track> QCDTrackAnalyzer::getAssociatedRecTrack
  (const TrackingParticleRef & simTrack, int & nRec, bool hasToBePrimary)
{
  edm::RefToBase<reco::Track> associatedRecTrack;
#ifndef NDEBUG
  vector<int> associatedRecId;
#endif

  float dmin = 1e+9;

  try
  {
    vector<pair<edm::RefToBase<reco::Track>, double> > recTracks = simToReco[simTrack];

    for(vector<pair<edm::RefToBase<reco::Track>,double> >::const_iterator
          it = recTracks.begin(); it != recTracks.end(); ++it)
    {
      edm::RefToBase<reco::Track> recTrack = it->first;
      float fraction                       = it->second;

      if(fraction > 0.5 && (!hasToBePrimary || isPrimary(recTrack)))
      {
#ifndef NDEBUG
        for(edm::View<reco::Track>::size_type j=0;
             j < recCollection.product()->size(); ++j)
        {
          edm::RefToBase<reco::Track> rt(recCollection, j);
          if(rt == recTrack) associatedRecId.push_back(j);
        }
#endif
        float d0 = fabs(recTrack->dxy(theBeamSpot->position()));

        if(d0 < dmin)
        { 
          dmin = d0;
          associatedRecTrack = recTrack;
        }
 
        nRec++;
      }
    }
  }
  catch (cms::Exception& event)
  { }

#ifndef NDEBUG
  if(nRec > 1)
  {
    sort(associatedRecId.begin(), associatedRecId.end());
    ostringstream o;

    for(vector<int>::const_iterator id = associatedRecId.begin();
                                    id!= associatedRecId.end(); id++)
      o << " #" << *id;

    LogTrace("MinBiasTracking")
      << " \033[22;32m" << "[TrackAnalyzer] multiple reco:"
      << o.str() << "\033[22;0m";
  }
#endif

  return associatedRecTrack;
}

/*****************************************************************************/
TrackingParticleRef QCDTrackAnalyzer::getAssociatedSimTrack
  (const edm::RefToBase<reco::Track> & recTrack, int & nSim)
{
  TrackingParticleRef associatedSimTrack;

  try
  {
    vector<pair<TrackingParticleRef, double> > simTracks = recoToSim[recTrack];

    for(vector<pair<TrackingParticleRef, double> >::const_iterator
          it = simTracks.begin(); it != simTracks.end(); ++it)
    {
      TrackingParticleRef simTrack = it->first;
      float fraction               = it->second;

      // If more than half is shared
      if(fraction > 0.5)
      {
        associatedSimTrack = simTrack; nSim++;
      }
    }
  }
  catch (cms::Exception& event)
  { }

  return associatedSimTrack;
}


/*****************************************************************************/
float QCDTrackAnalyzer::refitWithVertex
  (const reco::Track & recTrack)
{ 
  TransientTrack theTransientTrack = theTTBuilder->build(recTrack);
  
  // If there are vertices found
  if(vertices->size() > 0)
  { 
    float dzmin = -1.; 
    const reco::Vertex * closestVertex = 0;

    // Look for the closest vertex in z
    for(reco::VertexCollection::const_iterator
        vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
    {
      float dz = fabs(recTrack.vertex().z() - vertex->position().z());
      if(vertex == vertices->begin() || dz < dzmin)
      { dzmin = dz ; closestVertex = &(*vertex); }
    }

    // Get vertex position and error matrix
    GlobalPoint vertexPosition(closestVertex->position().x(),
                               closestVertex->position().y(),
                               closestVertex->position().z());

    float beamSize = theBeamSpot->BeamWidth();
    GlobalError vertexError(beamSize*beamSize, 0,
                            beamSize*beamSize, 0,
                            0,closestVertex->covariance(2,2));

    // Refit track with vertex constraint
    SingleTrackVertexConstraint stvc;
    pair<TransientTrack, float> result =
      stvc.constrain(theTransientTrack, vertexPosition, vertexError);

    return result.first.impactPointTSCP().pt();
  }
  else
    return recTrack.pt();
}

/*****************************************************************************/
int QCDTrackAnalyzer::processSimTracks()
{
  int ntrk = 0;

  for(TrackingParticleCollection::size_type i=0;
               i < simCollection.product()->size(); ++i)
  {
    const TrackingParticleRef simTrack(simCollection, i);

    SimTrack_t s;

    //
    s.ntrkr = ntrk; // ntrk

    // sim
    s.ids = simTrack->pdgId();               // ids
    s.prim = (simTrack->parentVertex()->position().T() == 0.); // prim

    s.etas = simTrack->eta();                 // etas
    s.pts  = simTrack->pt();                  // pts

    bool acc;
    
    if(simTrack->charge() != 0)
    {
      acc = isAccepted(*simTrack);
      
      // primary charged particles with |eta|<2.4
      if(s.prim && fabs(s.etas) < 2.4)
        ntrk++; 
    }
    else
    {
      acc = false;

      if(simTrack->decayVertices().size() == 1)
      if((simTrack->decayVertices()).at(0)->nDaughterTracks() == 2)
      { 
        const TrackingParticleRefVector& daughters =
          (simTrack->decayVertices()).at(0)->daughterTracks();

        acc = isAccepted(*(daughters.at(0))) &&
              isAccepted(*(daughters.at(1)));
      }
    }
    s.acc = acc; // acc

    // rec
    int nRec = 0;

    float charge = -9999; float eta    = -9999;
    float pt     = -9999; float dxy    = -9999;

    if(simTrack->charge() != 0)
    {
      edm::RefToBase<reco::Track> aRecTrack = 
        getAssociatedRecTrack(simTrack, nRec, true);
//        getAssociatedRecTrack(simTrack, nRec, true);

      if(nRec > 0)
      {
        charge = aRecTrack->charge();
        eta    = aRecTrack->eta();
        pt     = aRecTrack->pt();
        dxy    = aRecTrack->dxy(theBeamSpot->position());
      }

#ifndef DEBUG
     if(s.prim && nRec == 0)
       LogTrace("MinBiasTracking")
         << " \033[22;32m" << "[TrackAnalyzer] not reconstructed: #" << i
         << " (eta=" << s.etas << ", pt="  << s.pts << ")"
         << "\033[22;0m" << endl;
#endif

    } 
    else
    {
      if(simTrack->decayVertices().size() == 1)
      if((simTrack->decayVertices()).at(0)->nDaughterTracks() == 2)
      {
        const TrackingParticleRefVector& daughters = 
          (simTrack->decayVertices()).at(0)->daughterTracks();

        int nRec_1 = 0;
        edm::RefToBase<reco::Track> aRecTrack_0 =
          getAssociatedRecTrack(daughters.at(0), nRec_1, false);

        int nRec_2 = 0;
        edm::RefToBase<reco::Track> aRecTrack_1 =
          getAssociatedRecTrack(daughters.at(1), nRec_2, false);

        if(nRec_1 >= 1 && nRec_2 >= 1) // FIXME
        {
          edm::RefToBase<reco::Track> posRecTrack;
          edm::RefToBase<reco::Track> negRecTrack;

          if(daughters.at(0)->charge() > 0)
          { posRecTrack = aRecTrack_0; negRecTrack = aRecTrack_1; }
          else
          { posRecTrack = aRecTrack_1; negRecTrack = aRecTrack_0; }

          for(reco::VZeroCollection::const_iterator vZero = vZeros->begin();
                                                    vZero!= vZeros->end();
                                                    vZero++)
          if(&(*(edm::RefToBase<reco::Track>(vZero->positiveDaughter()))) ==
                                                &(*posRecTrack) &&
             &(*(edm::RefToBase<reco::Track>(vZero->negativeDaughter()))) == 
                                                &(*negRecTrack)) 
          {
            charge = 0.;

            GlobalVector momentum =  vZero->momenta().first +
                                     vZero->momenta().second;
          
            eta    = momentum.eta();
            pt     = momentum.perp();
            dxy    = vZero->impactMother();
            nRec++;
          }
        }  
      }
    }

    s.nrec = nRec;   // nrec

    // fill
    histograms->fillSimHistograms(s);
  }

  return ntrk;
}

/*****************************************************************************/
float QCDTrackAnalyzer::getEnergyLoss(const reco::TrackRef & track)
{
  const DeDxDataValueMap & eloss = *energyLoss;
  return eloss[track].dEdx();
}

/*****************************************************************************/
float QCDTrackAnalyzer::getInvariantMass(const VZero & vZero, float m1, float m2)
{
  GlobalVector p1 = vZero.momenta().first;
  GlobalVector p2 = vZero.momenta().second;

  double E1 = sqrt(p1*p1 + m1*m1);
  double E2 = sqrt(p2*p2 + m2*m2);

  return sqrt((E1+E2)*(E1+E2) - (p1+p2)*(p1+p2));
}

/*****************************************************************************/
float QCDTrackAnalyzer::getInvariantMass
  (const reco::Track * r1,
   const reco::Track * r2,
   float m1, float m2)
{
  GlobalVector p1(r1->momentum().x(),
                  r1->momentum().y(),
                  r1->momentum().z());
  GlobalVector p2(r2->momentum().x(),
                  r2->momentum().y(),
                  r2->momentum().z());

  double E1 = sqrt(p1*p1 + m1*m1);
  double E2 = sqrt(p2*p2 + m2*m2);

  return sqrt((E1+E2)*(E1+E2) - (p1+p2)*(p1+p2));
}


/*****************************************************************************/
double QCDTrackAnalyzer::getSigmaOfLogdEdx(double logde)
{
  return 0.3;
}

/****************************************************************************/
double QCDTrackAnalyzer::truncate(double m)
{
  const double s = 0.346;
  const double t = (ch == 0 ? 3.05 : 2.86);

  return m - s * sqrt(2/M_PI)*exp(-(m-t)*(m-t)/2/s/s)/(1+erf((t-m)/s/sqrt(2)));
}

/****************************************************************************/
double QCDTrackAnalyzer::getLogdEdx(double bg)
{
  const double a =  3.25 ;
  const double b =  0.288;
  const double c = -0.852;

  double beta = bg/sqrt(bg*bg + 1);
  double dedx = log( a/(beta*beta) + b * log(bg) + c );

   return dedx;
//  return truncate(dedx);
}

/****************************************************************************/
bool QCDTrackAnalyzer::isCompatibleWithdEdx(double m, const reco::TrackRef & track)
{  
  ch = (track->charge() > 0 ? 0 : 1);
   
  // no usable dE/dx if p > 2
  if(track->p() > 2) return true;

  double bg = track->p() / m;

  double theo = getLogdEdx(bg);

  // !!!!!!
  int nhitr = track->numberOfValidHits();

  double meas = log(getEnergyLoss(track));
  double sigm = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
   
  if(theo - 5 * sigm < meas && meas < theo + 5 * sigm)
    return true;
  else
    return false;
}

/*****************************************************************************/
void QCDTrackAnalyzer::processVZeros()
{
  const float mel = 0.511e-3;
  const float mpi = 0.13957018;
  const float mpr = 0.93827203;

  RecVzero_t v;

  for(reco::VZeroCollection::const_iterator vZero = vZeros->begin();
                                            vZero!= vZeros->end();
                                            vZero++)
  {
    GlobalVector momentum = vZero->momenta().first +
                            vZero->momenta().second;

    v.etar = momentum.eta();               // eta
    v.ptr  = momentum.perp();              // pt

    v.rhor = vZero->crossingPoint().rho(); // rho

//    result.push_back(vZero->armenterosPt());        // qt
//    result.push_back(vZero->armenterosAlpha());     // alpha

    if(isCompatibleWithdEdx(mel, ((vZero->positiveDaughter()))) &&
       isCompatibleWithdEdx(mel, ((vZero->negativeDaughter()))))
    {
      v.ima = getInvariantMass(*vZero, mel,mel); // gam
      histograms->fillVzeroHistograms(v, 0);
    }

    if(isCompatibleWithdEdx(mpi, ((vZero->positiveDaughter()))) &&
       isCompatibleWithdEdx(mpi, ((vZero->negativeDaughter()))))
    {
      v.ima = getInvariantMass(*vZero, mpi,mpi); // k0s
      histograms->fillVzeroHistograms(v, 1);
    }

    if(isCompatibleWithdEdx(mpr, ((vZero->positiveDaughter()))) &&
       isCompatibleWithdEdx(mpi, ((vZero->negativeDaughter()))))
    {
      v.ima = getInvariantMass(*vZero, mpr,mpi); // lam
      histograms->fillVzeroHistograms(v, 2);
    }

    if(isCompatibleWithdEdx(mpi, ((vZero->positiveDaughter()))) &&
       isCompatibleWithdEdx(mpr, ((vZero->negativeDaughter()))))
    {
      v.ima = getInvariantMass(*vZero, mpi,mpr); // ala
      histograms->fillVzeroHistograms(v, 3);
    }
  }
}

/*****************************************************************************/
int QCDTrackAnalyzer::processRecTracks
  (edm::Handle<reco::DeDxDataValueMap> elossCollection)
{
  int ntrk = 0;

  for(edm::View<reco::Track>::size_type i=0;
          i < recCollection.product()->size(); ++i)
  {
    RecTrack_t r;

    edm::RefToBase<reco::Track> recTrack(recCollection, i);

    r.ntrkr = ntrk;

    // rec
    r.charge = recTrack->charge();         // charge
    r.etar   = recTrack->eta();            // etar
    r.ptr    = refitWithVertex(*recTrack); // ptr  
    r.phir   = recTrack->phi();            // phir

    r.logpr = log(recTrack->p());
    r.nhitr = (*elossCollection.product())[recTrack.castTo<reco::TrackRef>()].numberOfMeasurements();


//    r.nhitr = recTrack->numberOfValidHits();           // nhitr

    r.prim  = isPrimary(recTrack);

    r.zr    = recTrack->dz(theBeamSpot->position());  // dzr

    r.logde = log((*elossCollection.product())[recTrack.castTo<reco::TrackRef>()].dEdx()); // log(dedx)

    if(r.prim && fabs(r.etas) < 2.4)
      ntrk++;

    // sim 
    if(hasSimInfo)
    {
    int nSim = 0;
    TrackingParticleRef aSimTrack =
      getAssociatedSimTrack(recTrack, nSim);

#ifndef DEBUG
  if(nSim == 0)
    LogTrace("MinBiasTracking")
      << " \033[22;32m" << "[TrackAnalyzer] fake track: #" << i << ","
      << " d0=" << recTrack->d0() << " cm" << "\033[22;0m";
#endif

    r.nsim = nSim; // nsim

    if(nSim > 0)
    {
      int parentId;

      if(aSimTrack->parentVertex()->nSourceTracks() == 0)
        parentId = 0;
      else
        parentId = (*(aSimTrack->parentVertex()->sourceTracks_begin()))->pdgId();

      r.ids    = aSimTrack->pdgId();      // ids
      r.parids = parentId;                // parids
      r.etas   = aSimTrack->eta();        // etas
      r.pts    = aSimTrack->pt();         // pts
    }
    }

    // fill
    histograms->fillRecHistograms(r);
  }

  return ntrk;
}

/*****************************************************************************/
void QCDTrackAnalyzer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  LogTrace("MinBiasTracking") << "[TrackAnalyzer]";

  // Get process id
  if(hasSimInfo)
  {
  edm::Handle<edm::HepMCProduct> hepEv;
  ev.getByType(hepEv);
  proc = hepEv->GetEvent()->signal_process_id();
  LogTrace("MinBiasTracking") << " [TrackAnalyzer] process = " << proc;

  ev.getByLabel("mergedtruth", simCollection);
  }
  else proc = 0;

  // Get HF info
/*
  edm::Handle<vector<L1GctJetCounts> >  jetCountDigi; 
  ev.getByLabel("hltGctDigis",          jetCountDigi);

  LogTrace("MinBiasTracking") << " [TrackAnalyzer] trigger(" 
    << jetCountDigi.product()->size()
    << ") :" 
    << " HF+("
    << jetCountDigi.product()->at(0).hfTowerCountPositiveEta() << ")"
    << " HF-("
    << jetCountDigi.product()->at(0).hfTowerCountNegativeEta() << ")";
*/

  // Get reconstructed tracks
  ev.getByLabel(trackProducer, recCollection);
  LogTrace("MinBiasTracking") << " [TrackAnalyzer] recTracks    = "
    << recCollection.product()->size();

  edm::View<reco::Track> rTrackCollection = *(recCollection.product());

  // Get reconstructed dE/dx
  edm::Handle<reco::DeDxDataValueMap>   elossCollection;
  ev.getByLabel("energyLoss", "energyLossStrHits", elossCollection);
  energyLoss = elossCollection.product();

  // Get reconstructed V0s
  edm::Handle<reco::VZeroCollection> vZeroCollection;
  ev.getByLabel("pixelVZeros",vZeroCollection);
  vZeros = vZeroCollection.product();

  // Get beamSpot
  edm::Handle<reco::BeamSpot>      beamSpotHandle;
  ev.getByLabel("offlineBeamSpot", beamSpotHandle);
  theBeamSpot = beamSpotHandle.product();

  LogTrace("MinBiasTracking")
    << " [TrackAnalyzer] beamSpot at " << theBeamSpot->position() 
    << " sigmaZ = "    << theBeamSpot->sigmaZ()
    << " BeamWidth = " << theBeamSpot->BeamWidth() << endl;

  // Get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  ev.getByLabel("pixelVertices",vertexCollection);
  vertices = vertexCollection.product();

  // Proocess eventInfo
  ntrk = recCollection.product()->size();
  nvtx = vertexCollection.product()->size();

  // Associators
  if(hasSimInfo)
  {
    LogTrace("MinBiasTracking") << " [TrackAnalyzer] associateSimToReco";
    simToReco =
    theAssociatorByHits->associateSimToReco(recCollection, simCollection,&ev);

    LogTrace("MinBiasTracking") << " [TrackAnalyzer] associateRecoToSim";
    recoToSim =
    theAssociatorByHits->associateRecoToSim(recCollection, simCollection,&ev);
  }

  // Analyze
  int prim_s_tracks = 0;
  if(hasSimInfo)
  {
    LogTrace("MinBiasTracking") << " [TrackAnalyzer] processSimTracks";
    prim_s_tracks = processSimTracks();
  }

  LogTrace("MinBiasTracking") << " [TrackAnalyzer] processRecTracks";
  int prim_r_tracks = processRecTracks(elossCollection);

  histograms->fillEventInfo(proc, prim_s_tracks,
                                  prim_r_tracks);

  LogTrace("MinBiasTracking") << " [TrackAnalyzer] processVZeros";
  processVZeros();
}

DEFINE_FWK_MODULE(QCDTrackAnalyzer);
