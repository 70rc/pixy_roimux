//
// Created by damian on 6/8/17.
//

#include "PrincipalComponentsCluster.h"


namespace pixy_roimux {
    void PrincipalComponentsCluster::analyseEvents(
            ChargeHits &t_chargeHits,
            const bool t_rejectOutliers,
            const bool t_rejectAmbiguities) {
        for (auto &&event : t_chargeHits.getEvents()) {
            std::cout << "Performing PCA for event number " << event.eventId << "...\n";
            event.pcaIds.resize(event.hitCandidates.size(), - 1);
            int err = analysis3D(event);
            if (err) {
                continue;
            }
            if (t_rejectAmbiguities) {
                rejectAmbiguities(event);
                analysis3D(event);
                if (t_rejectOutliers) {
                    unsigned iter = 0;
                    int totRejHits = 0;
                    int numRejHits;
                    int maxRejects = 0;
                    for (const auto &hitCandidates : event.hitCandidates) {
                        if (!hitCandidates.empty()) {
                            ++maxRejects;
                        }
                    }
                    maxRejects = static_cast<int>(0.4 * maxRejects);
                    do {
                        double maxRange = m_runParams.getPcaScaleFactor() * 0.5 *
                                (3. * sqrt(event.principalComponents.eigenValues.at(1)) +
                                 event.principalComponents.aveHitDoca);
                        numRejHits = rejectOutliers(event, maxRange);
                        totRejHits += numRejHits;
                        analysis3D(event);
                        ++iter;
                    } while ((iter <= m_runParams.getPcaMaxIterations()) && (numRejHits > 0) && (totRejHits < maxRejects));
                    std::cout << "Finished after " << iter << " iterations with " << totRejHits << " rejected hits.\n";
                }
            }
        }
    }


    int PrincipalComponentsCluster::analysis3D(Event &t_event) {
        std::vector<double> meanPos(3, 0.);
        unsigned numPairsInt = 0;

        auto pcaId = t_event.pcaIds.cbegin();
        for (const auto &hitCandidates : t_event.hitCandidates) {
            unsigned hitId = 0;
            for (const auto &hit : hitCandidates) {
                if ((*pcaId == hitId) || (*pcaId == -1)) {
                    meanPos.at(0) += hit.x;
                    meanPos.at(1) += hit.y;
                    meanPos.at(2) += hit.z;
                    ++numPairsInt;
                }
                ++hitId;
            }
            ++pcaId;
        }
        if (!numPairsInt) {
            std::cerr << "WARNING: No hits found for event " << t_event.eventId << '.' << std::endl;
            return 1;
        }

        double numPairs = static_cast<double>(numPairsInt);

        meanPos.at(0) /= numPairs;
        meanPos.at(1) /= numPairs;
        meanPos.at(2) /= numPairs;

        double xi2 = 0.;
        double xiyi = 0.;
        double xizi = 0.;
        double yi2 = 0.;
        double yizi = 0.;
        double zi2 = 0.;
        double weightSum = 0.;

        pcaId = t_event.pcaIds.cbegin();
        for (const auto &hitCandidates : t_event.hitCandidates) {
            unsigned hitId = 0;
            for (const auto &hit : hitCandidates) {
                if ((*pcaId == hitId) || (*pcaId == -1)) {
                    //double weight = hit.charge;
                    double weight = 1.;

                    double x = (hit.x - meanPos.at(0)) * weight;
                    double y = (hit.y - meanPos.at(1)) * weight;
                    double z = (hit.z - meanPos.at(2)) * weight;

                    weightSum += weight * weight;

                    xi2 += x * x;
                    xiyi += x * y;
                    xizi += x * z;
                    yi2 += y * y;
                    yizi += y * z;
                    zi2 += z * z;
                }
                ++hitId;
            }
            ++pcaId;
        }

        Eigen::Matrix3f sig;

        sig <<  xi2,    xiyi,   xizi,
                xiyi,   yi2,    yizi,
                xizi,   yizi,   zi2;

        sig *= 1. / weightSum;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

        if (eigenMat.info() == Eigen::ComputationInfo::Success) {
            using eigenValColPair = std::pair<float, size_t>;
            std::vector<eigenValColPair> eigenValColVec = {
                    eigenValColPair(eigenMat.eigenvalues()(0), 0),
                    eigenValColPair(eigenMat.eigenvalues()(1), 1),
                    eigenValColPair(eigenMat.eigenvalues()(2), 2)
            };

            std::sort(eigenValColVec.begin(), eigenValColVec.end(),
                      [](const eigenValColPair& left, const eigenValColPair& right){return left.first > right.first;});

            t_event.principalComponents.eigenValues = {
                    eigenValColVec.at(0).first,
                    eigenValColVec.at(1).first,
                    eigenValColVec.at(2).first
            };

            t_event.principalComponents.eigenVectors.clear();
            Eigen::Matrix3f eigenVecs(eigenMat.eigenvectors());
            for (const auto &pair : eigenValColVec) {
                std::vector<double> tempVec = {
                        eigenVecs(0, pair.second),
                        eigenVecs(1, pair.second),
                        eigenVecs(2, pair.second)
                };
                t_event.principalComponents.eigenVectors.emplace_back(tempVec);
            }

            t_event.principalComponents.numHitsUsed = numPairsInt;
            t_event.principalComponents.avePosition = meanPos;
        }
        else {
            std::cerr << "WARNING: PCA decompose failure for event " << t_event.eventId
                      << ", numPairs = " << numPairs << std::endl;
            return 1;
        }
        return 0;
    }


    double PrincipalComponentsCluster::computeDoca(
            const Hit3d &t_hit,
            const TVector3 &t_avePosition,
            const TVector3 &t_axisDirVec) {
        TVector3 clusPos(t_hit.x, t_hit.y, t_hit.z);
        TVector3 clusToHitVec = clusPos - t_avePosition;
        double arclenToPoca = clusToHitVec.Dot(t_axisDirVec);
        TVector3 docaPos = t_avePosition + arclenToPoca * t_axisDirVec;
        TVector3 docaPosToClusPos = clusPos - docaPos;

        return docaPosToClusPos.Mag();
    }


    void PrincipalComponentsCluster::rejectAmbiguities(Event &t_event) {
        TVector3 avePosition(t_event.principalComponents.avePosition.at(0),
                             t_event.principalComponents.avePosition.at(1),
                             t_event.principalComponents.avePosition.at(2));
        TVector3 axisDirVec(t_event.principalComponents.eigenVectors.at(0).at(0),
                             t_event.principalComponents.eigenVectors.at(0).at(1),
                             t_event.principalComponents.eigenVectors.at(0).at(2));
        std::vector<double> eventDocas;
        auto pcaId = t_event.pcaIds.begin();
        for (const auto &hitCandidates : t_event.hitCandidates) {
            if (!hitCandidates.empty()) {
                std::vector<double> hitDocas(hitCandidates.size());
                auto doca = hitDocas.begin();
                for (const auto &hit : hitCandidates) {
                    *doca = computeDoca(hit, avePosition, axisDirVec);
                    eventDocas.emplace_back(*doca);
                    ++doca;
                }
                *pcaId = static_cast<int>(std::distance(hitDocas.cbegin(),
                                                        std::min_element(hitDocas.cbegin(), hitDocas.cend())));
            }
            else {
                *pcaId = -3;
            }
            ++pcaId;
        }
        t_event.principalComponents.aveHitDoca =
                std::accumulate(eventDocas.cbegin(), eventDocas.cend(), static_cast<double>(0.))
                / static_cast<double>(eventDocas.size());
    }


    int PrincipalComponentsCluster::rejectOutliers(
            Event &t_event,
            const double maxDocaAllowed) {
        int numRejHits = 0;
        TVector3 avePosition(t_event.principalComponents.avePosition.at(0),
                             t_event.principalComponents.avePosition.at(1),
                             t_event.principalComponents.avePosition.at(2));
        TVector3 axisDirVec(t_event.principalComponents.eigenVectors.at(0).at(0),
                            t_event.principalComponents.eigenVectors.at(0).at(1),
                            t_event.principalComponents.eigenVectors.at(0).at(2));

        std::vector<double> docas;
        auto pcaId = t_event.pcaIds.begin();
        for (const auto &hitCandidates : t_event.hitCandidates) {
            if (*pcaId >= 0) {
                double doca = computeDoca(hitCandidates.at(static_cast<unsigned>(*pcaId)), avePosition, axisDirVec);
                docas.emplace_back(doca);
                if (doca > maxDocaAllowed) {
                    *pcaId = - 2;
                    ++numRejHits;
                }
            }
            ++pcaId;
        }
        t_event.principalComponents.aveHitDoca =
                std::accumulate(docas.cbegin(), docas.cend(), static_cast<double>(0.))
                / static_cast<double>(docas.size());

        return numRejHits;
    }
}
