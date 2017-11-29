//
// Created by damian on 6/22/17.
//

#include <FieldManager.h>
#include "KalmanFit.h"


namespace pixy_roimux {
    KalmanFit::KalmanFit(
            const RunParams &t_runParams,
            const std::string t_geoFileName,
            const bool t_initDisplay) :
            m_runParams(t_runParams),
            m_kalmanFitter(t_runParams.getKalmanUseRef(),
                           t_runParams.getKalmanDeltaPval(),
                           t_runParams.getKalmanDeltaWeight()),
            m_hits("genfit::mySpacepointDetectorHit"),
            m_measurementProducer(&m_hits),
            m_cov(6),
            m_posCov(3)
    {
        m_kalmanFitter.setMaxIterations(m_runParams.getKalmanMaxIterations());
        for (unsigned i = 0; i < 3; ++i) {
            double var = m_runParams.getKalmanPosErr().at(i) * m_runParams.getKalmanPosErr().at(i);
            m_posCov(i, i) = var;
            m_cov(i, i) = var;
            m_cov((i + 3), (i + 3)) = m_runParams.getKalmanMomErr().at(i) * m_runParams.getKalmanMomErr().at(i);
        }

        m_gaussRng.SetSeed(m_runParams.getKalmanRngSeed());
        new TGeoManager("genfitGeometry", "GENFIT geometry");
        TGeoManager::Import(t_geoFileName.c_str());
        genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 0.));
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        m_measurementFactory.addProducer(m_detId, &m_measurementProducer);

        if (t_initDisplay) {
            m_display = genfit::EventDisplay::getInstance();
        }
    }


    void KalmanFit::fitEvent(const Event &t_event)
    {
        m_hits.Clear();
        genfit::TrackCand trackCand;
        TVector3 trackPos(t_event.principalComponents.avePosition.at(0),
                          t_event.principalComponents.avePosition.at(1),
                          t_event.principalComponents.avePosition.at(2));
        TVector3 trackMom(t_event.principalComponents.eigenVectors.at(0).at(0),
                          t_event.principalComponents.eigenVectors.at(0).at(1),
                          t_event.principalComponents.eigenVectors.at(0).at(2));
        trackMom.SetMag(m_runParams.getKalmanMomMag());
        std::multimap<double, unsigned> hitOrderZ;
        unsigned hitId = 0;
        auto pcaId = t_event.pcaIds.cbegin();
        for (const auto &hitCandidates : t_event.hitCandidates) {
            unsigned hitCandidateId = 0;
            for (const auto &hit : hitCandidates) {
                TVector3 hitPos(hit.x, hit.y, hit.z);
                new(m_hits[hitId]) genfit::mySpacepointDetectorHit(hitPos, m_posCov);
                if (hitCandidateId == *pcaId){
                    hitOrderZ.insert({hit.z, hitId});
                }
                ++hitCandidateId;
                ++hitId;
            }
            ++pcaId;
        }
        if (hitOrderZ.empty()) {
            std::cerr << "WARNING: Event contains no usable hits, skipping..." << std::endl;
            return;
        }
        for (const auto& orderedHit : hitOrderZ) {
            trackCand.addHit(m_detId, orderedHit.second);
        }
        trackCand.setPosMomSeedAndPdgCode(trackPos, trackMom, m_runParams.getKalmanPdgCode());
        trackCand.setCovSeed(m_cov);

        genfit::Track fitTrack(trackCand, m_measurementFactory, new genfit::RKTrackRep(m_runParams.getKalmanPdgCode()));
        fitTrack.checkConsistency();
        try {
            m_kalmanFitter.processTrack(&fitTrack);
        }
        catch(genfit::Exception &exception) {
            std::cerr << exception.what();
            std::cerr << "Exception, next track." << std::endl;
            return;
        }
        fitTrack.checkConsistency();
        fitTrack.determineCardinalRep();
        try {
            fitTrack.getFittedState().Print();
        }
        catch(genfit::Exception &exception) {
            std::cerr << exception.what();
            std::cerr << "Exception, next track." << std::endl;
            return;
        }

        m_trackPtr = &fitTrack;
        m_eventId = t_event.eventId;
        m_tree->Fill();
        m_trackPtr = nullptr;

        if (m_display) {
            m_display->addEvent(&fitTrack);
        }
    }


    void KalmanFit::fit(
            const ChargeHits &t_chargeHits,
            const std::string t_treeFileName) {
        TFile treeFile(t_treeFileName.c_str(), "RECREATE");
        m_tree = std::unique_ptr<TTree>(new TTree("genfitTree", "genfitTree"));
        m_tree->Branch("Track", &m_trackPtr, 64000, 0);
        m_tree->Branch("eventId", &m_eventId, "eventId/i");
        for (const auto& event : t_chargeHits.getEvents()) {
            std::cout << "Fitting event number " << event.eventId << "...\n";
            if (event.pcaIds.empty()) {
                std::cerr << "WARNING: Event contains no PCA information, skipping..." << std::endl;
                continue;
            }
            else if (std::find(event.pcaIds.cbegin(), event.pcaIds.cend(), -1) != event.pcaIds.cend()) {
                std::cerr << "WARNING: Event contains PCA errors, skipping..." << std::endl;
                continue;
            }
            fitEvent(event);
        }
        m_tree->Write();
        m_tree.reset(nullptr);
        treeFile.Close();
    }
}