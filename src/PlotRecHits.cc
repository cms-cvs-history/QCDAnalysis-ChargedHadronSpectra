#include "QCDAnalysis/ChargedHadronSpectra/interface/PlotRecHits.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

#include <vector>

using namespace std;

/*****************************************************************************/
PlotRecHits::PlotRecHits
  (const edm::EventSetup& es, ofstream& file_) : file(file_)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> trackerHandle;
  es.get<TrackerDigiGeometryRecord>().get(trackerHandle);
  theTracker = trackerHandle.product();
}

/*****************************************************************************/
PlotRecHits::~PlotRecHits()
{
}

/*****************************************************************************/
void PlotRecHits::printPixelRecHit(const SiPixelRecHit * recHit)
{
  DetId id = recHit->geographicalId();

  // DetUnit
  double x = theTracker->idToDet(id)->surface().bounds().width() /2;
  double y = theTracker->idToDet(id)->surface().bounds().length()/2;
  double z = 0.;

  GlobalPoint p00 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x,-y,z));
  GlobalPoint p01 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x, y,z));
  GlobalPoint p10 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x,-y,z));
  GlobalPoint p11 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x, y,z));

  file << ", If[sd, {RGBColor[0.4,0.4,0.4], "
                 <<"Line[{{"<< p00.x()<<","<<p00.y()<<","<<p00.z()<<"*z}, "
                       <<"{"<< p01.x()<<","<<p01.y()<<","<<p01.z()<<"*z}, "
                       <<"{"<< p11.x()<<","<<p11.y()<<","<<p11.z()<<"*z}, "
                       <<"{"<< p10.x()<<","<<p10.y()<<","<<p10.z()<<"*z}, "
                       <<"{"<< p00.x()<<","<<p00.y()<<","<<p00.z()<<"*z}}]}]"
       << endl;
  
  // RecHit
  LocalPoint lpos; GlobalPoint p;
  
  lpos = LocalPoint(recHit->localPosition().x(),
                    recHit->localPosition().y(),
                    recHit->localPosition().z());

  p = theTracker->idToDet(id)->toGlobal(lpos);
  file << ", Point[{"<<p.x()<<","<<p.y()<<","<<p.z()<<"*z}]" << endl;

  // Cluster details
  SiPixelRecHit::ClusterRef const& cluster = recHit->cluster();
  vector<SiPixelCluster::Pixel> pixels = cluster->pixels();

  file << ", Text[StyleForm[\"r\", FontFamily->\"Helvetica\", URL -> \"";
  for(vector<SiPixelCluster::Pixel>::const_iterator
    pixel = pixels.begin(); pixel!= pixels.end(); pixel++)
  {
    file << "[" << int(pixel->x)
         << " " << int(pixel->y)
         << " " << int(pixel->adc/135) << "]";
  }
  file << "\"]";

  file << ", {"<< p.x()<<","<<p.y()<<","<<p.z()<<"*z}"
       << ", {0,1}]" << endl;
}

/*****************************************************************************/
void PlotRecHits::printStripRecHit(const SiStripRecHit2D * recHit)
{
  DetId id = recHit->geographicalId();

  // DetUnit
  double x = theTracker->idToDet(id)->surface().bounds().width() /2;
  double y = theTracker->idToDet(id)->surface().bounds().length()/2;
  double z = 0.;

  GlobalPoint p00 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x,-y,z));
  GlobalPoint p01 =  theTracker->idToDet(id)->toGlobal(LocalPoint(-x, y,z));
  GlobalPoint p10 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x,-y,z));
  GlobalPoint p11 =  theTracker->idToDet(id)->toGlobal(LocalPoint( x, y,z));

  file << ", If[sd, {RGBColor[0.6,0.6,0.6], "
                 <<"Line[{{"<<p00.x()<<","<<p00.y()<<","<<p00.z()<<"*z}, "
                       <<"{"<<p01.x()<<","<<p01.y()<<","<<p01.z()<<"*z}, "
                       <<"{"<<p11.x()<<","<<p11.y()<<","<<p11.z()<<"*z}, "
                       <<"{"<<p10.x()<<","<<p10.y()<<","<<p10.z()<<"*z}, "
                       <<"{"<<p00.x()<<","<<p00.y()<<","<<p00.z()<<"*z}}]}]"
       << endl;

  // RecHit
  LocalPoint lpos; GlobalPoint p;

  lpos = LocalPoint(recHit->localPosition().x(),
                 y, recHit->localPosition().z());
  p = theTracker->idToDet(id)->toGlobal(lpos);
  file << ", Line[{{"<<p.x()<<","<<p.y()<<","<<p.z()<<"*z}, {";

  lpos = LocalPoint(recHit->localPosition().x(),
                -y, recHit->localPosition().z());
  p = theTracker->idToDet(id)->toGlobal(lpos);
  file << ""<<p.x()<<","<<p.y()<<","<<p.z()<<"*z}}]" << endl;
}

/*****************************************************************************/
void PlotRecHits::printPixelRecHits(const edm::Event& ev)
{
  // Get pixel hit collections
  vector<edm::Handle<SiPixelRecHitCollection> > pixelColls;
  ev.getManyByType(pixelColls);

  for(vector<edm::Handle<SiPixelRecHitCollection> >::const_iterator
      pixelColl = pixelColls.begin();
      pixelColl!= pixelColls.end(); pixelColl++)
  {
    const SiPixelRecHitCollection* thePixelHits = (*pixelColl).product();

    for(SiPixelRecHitCollection::id_iterator
          id = thePixelHits->id_begin(); id!= thePixelHits->id_end(); id++)
    {
      SiPixelRecHitCollection::range range = thePixelHits->get(*id);

      // Take all hits
      for(SiPixelRecHitCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      if(recHit->isValid())
        printPixelRecHit(&(*recHit));
    }
  }
}

/*****************************************************************************/
void PlotRecHits::printStripRecHits(const edm::Event& ev)
{
  {
  // Get strip hit collections
  vector<edm::Handle<SiStripRecHit2DCollection> > stripColls;
  ev.getManyByType(stripColls);
  
  for(vector<edm::Handle<SiStripRecHit2DCollection> >::const_iterator
      stripColl = stripColls.begin();
      stripColl!= stripColls.end(); stripColl++)
  {
    const SiStripRecHit2DCollection* theStripHits = (*stripColl).product();
    
    for(SiStripRecHit2DCollection::id_iterator
          id = theStripHits->id_begin(); id!= theStripHits->id_end(); id++)
    {
      SiStripRecHit2DCollection::range range = theStripHits->get(*id);

      // Take all hits
      for(SiStripRecHit2DCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      if(recHit->isValid())
        printStripRecHit(&(*recHit));
    }
  } 
  }

  // Get matched strip hit collections
  {
  vector<edm::Handle<SiStripMatchedRecHit2DCollection> > stripColls;
  ev.getManyByType(stripColls);

  for(vector<edm::Handle<SiStripMatchedRecHit2DCollection> >::const_iterator
      stripColl = stripColls.begin();
      stripColl!= stripColls.end(); stripColl++)
  {
    const SiStripMatchedRecHit2DCollection* theStripHits = (*stripColl).product();

    for(SiStripMatchedRecHit2DCollection::id_iterator
          id = theStripHits->id_begin(); id!= theStripHits->id_end(); id++)
    {
      SiStripMatchedRecHit2DCollection::range range = theStripHits->get(*id);

      // Take all hits
      for(SiStripMatchedRecHit2DCollection::const_iterator
            recHit = range.first; recHit!= range.second; recHit++)
      {
        if(recHit->monoHit()->isValid())
          printStripRecHit((recHit->monoHit()));
        if(recHit->stereoHit()->isValid())
          printStripRecHit((recHit->stereoHit()));

        DetId id = recHit->geographicalId();
        LocalPoint lpos = recHit->localPosition();
        GlobalPoint p = theTracker->idToDet(id)->toGlobal(lpos);

        file << ", Point[{"<< p.x()<<","<<p.y()<<","<<p.z()<<"*z}]" << endl;
      }
    }
  }
  }
}

/*****************************************************************************/
void PlotRecHits::printRecHits(const edm::Event& ev)
{
  file << "AbsolutePointSize[5]";
  file << ", If[pr, {RGBColor[0.4,0.4,1.0]";
// FIXME
//  printPixelRecHits(ev);
  file << "}]";

  file << ", If[sr, {RGBColor[0.6,0.6,1.0]";
// FIXME
//  printStripRecHits(ev);
  file << "}]";
}

