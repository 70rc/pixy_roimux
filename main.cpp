#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "ChargeData.h"
#include "ChargeHits.h"
#include "NoiseFilter.h"
#include "PrincipalComponentsCluster.h"
#include "RunParams.h"
#include "KalmanFit.h"


int main(int argc, char** argv) {
    // Start time point for timer.
    auto clkStart = std::chrono::high_resolution_clock::now();

    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " runParamsFileName dataFileName rankingFileName geoFileName genfitTreeFileName csvBaseFileName [minRanking] [maxRanking]" << std::endl;
        exit(1);
    }
    const std::string runParamsFileName = std::string(argv[1]);
    const std::string dataFileName = std::string(argv[2]);
    const std::string rankingFileName = std::string(argv[3]);
    const std::string geoFileName = std::string(argv[4]);
    const std::string genfitTreeFileName = std::string(argv[5]);
    const std::string csvBaseFileName = std::string(argv[6]);
    int minRanking = 4;
    if (argc > 7) {
        minRanking = std::stoi(argv[7]);
    }
    int maxRanking = 4;
    if (argc > 8) {
        maxRanking = std::stoi(argv[8]);
    }
    const unsigned subrunId = 0;


    TFile rankingFile(rankingFileName.c_str(), "READ");
    if (!rankingFile.IsOpen()) {
        std::cerr << "ERROR: Failed to open ranking file " << rankingFileName << '!' << std::endl;
        exit(1);
    }
    TTree* rankingTree = nullptr;
    rankingFile.GetObject("RankingTime", rankingTree);
    if (!rankingTree) {
        std::cerr << "ERROR: Failed to find \"RankingTime\" tree in ranking file " << rankingFileName << '!' << std::endl;
        exit(1);
    }
    int ranking;
    rankingTree->SetBranchAddress("Ranking", &ranking);
    std::vector<unsigned> eventIds;
    for (unsigned event = 0; event < rankingTree->GetEntries(); ++event) {
        rankingTree->GetEvent(event);
        if ((ranking >= minRanking) && (ranking <= maxRanking)) {
            eventIds.push_back(event);
        }
    }
    rankingFile.Close();

    // Create the VIPER runParams containing all the needed run parameters.
    const pixy_roimux::RunParams runParams(runParamsFileName);

    // Load events directly from ROOT file.
    std::cout << "Extracting chargeData...\n";
    pixy_roimux::ChargeData chargeData(dataFileName, eventIds, subrunId, runParams);

    // Noise filter
    std::cout << "Filtering chargeData...\n";
    pixy_roimux::NoiseFilter noiseFilter(runParams);
    noiseFilter.filterData(chargeData);

    // Find the chargeHits.
    std::cout << "Initialising hit finder...\n";
    pixy_roimux::ChargeHits chargeHits(chargeData, runParams);
    // Set up TSpectrum.
    std::cout << "Running hit finder...\n";
    chargeHits.findHits();

    std::cout << "Initialising principle components analysis...\n";
    pixy_roimux::PrincipalComponentsCluster principalComponentsCluster(runParams);
    std::cout << "Running principle components analysis...\n";
    principalComponentsCluster.analyseEvents(chargeHits);

    std::cout << "Initialising Kalman Fitter...\n";
    pixy_roimux::KalmanFit kalmanFit(runParams, geoFileName, true);
    std::cout << "Running Kalman Fitter...\n";
    kalmanFit.fit(chargeHits, genfitTreeFileName);

    // Write chargeHits of events in eventIds vector to CSV files so we can plot them with viper3Dplot.py afterwards.
    unsigned nHitCandidates = 0;
    unsigned nAmbiguities = 0;
    unsigned nUnmatchedPixelHits = 0;
    // Loop through events.
    for (const auto& event : chargeHits.getEvents()) {
        std::cout << "Writing event number " << event.eventId << " to file...\n";
        // Compose CSV filename and open file stream.
        const std::string csvEventBaseFileName = csvBaseFileName + "_event" + std::to_string(event.eventId);
        const std::string csvHitsFileName = csvEventBaseFileName + "_hits.csv";
        std::ofstream csvHitsFile(csvHitsFileName, std::ofstream::out);
        csvHitsFile << "X,Y,Z,Q,A" << std::endl;
        // viperEvents.hitCandidates is a vector of vectors so we need two loops.
        // Outer loop. Actually loops through pixel chargeHits of current event.
        auto pcaId = event.pcaIds.cbegin();
        for (const auto& hitCandidates : event.hitCandidates) {
            // Inner loop. Loops through all hit candidates of current pixel hit.
            int hitId = 0;
            for (const auto& hit : hitCandidates) {
                int reject = 0;
                if (*pcaId == - 2) {
                    reject = 1;
                }
                else if (hitId != *pcaId) {
                    reject = 2;
                }
                // Append coordinates and charge to file.
                csvHitsFile << hit.x << ',' << hit.y << ',' << hit.z << ',' << hit.chargeInt << ',' << reject << std::endl;
                ++hitId;
            }
            ++pcaId;
        }
        // Close CSV file.
        csvHitsFile.close();
        const std::string csvPcaFileName = csvEventBaseFileName + "_pca.csv";
        std::ofstream csvPcaFile(csvPcaFileName, std::ofstream::out);
        if (!event.principalComponents.avePosition.empty()) {
            csvPcaFile << event.principalComponents.avePosition.at(0) << ','
                       << event.principalComponents.avePosition.at(1) << ','
                       << event.principalComponents.avePosition.at(2) << std::endl;
        }
        else {
            csvPcaFile << "0,0,0" << std::endl;
        }
        for (const auto &eigenVector : event.principalComponents.eigenVectors) {
            csvPcaFile << eigenVector.at(0) << ','
                       << eigenVector.at(1) << ','
                       << eigenVector.at(2) << std::endl;
        }
        csvPcaFile.close();
        // Calculate some stats.
        // Loop over all pixel chargeHits using the pixel to ROI hit runParams.
        for (const auto &candidateRoiHitIds : event.pixel2roi) {
            // Add number of ROI hit candidates for current pixel hit.
            nHitCandidates += candidateRoiHitIds.size();
            // If there's no ROI hit candidates, increment the unmatched counter.
            if (candidateRoiHitIds.empty()) {
                ++nUnmatchedPixelHits;
            }
                // If there's more than one ROI hit candidate, increment the ambiguity counter.
            else if (candidateRoiHitIds.size() > 1) {
                ++nAmbiguities;
            }
        }
    }
    unsigned long nEvents = chargeHits.getEvents().size();
    float averageHitCandidates = static_cast<float>(nHitCandidates) / static_cast<float>(nEvents);
    float averageAmbiguities = static_cast<float>(nAmbiguities) / static_cast<float>(nEvents);
    float averageUnmatchedPixelHits = static_cast<float>(nUnmatchedPixelHits) / static_cast<float>(nEvents);
    const std::string statsFileName = csvBaseFileName + "_stats.txt";
    std::ofstream statsFile(statsFileName, std::ofstream::out);
    statsFile << "Number of events processed: " << nEvents << std::endl;
    statsFile << "Average number of hit candidates per event: " << averageHitCandidates << std::endl;
    statsFile << "Average number of ambiguities per event: " << averageAmbiguities << std::endl;
    statsFile << "Average number of unmatched pixel chargeHits per event: " << averageUnmatchedPixelHits << std::endl;
    statsFile.close();

    std::cout << "Done.\n";

    // Stop time point for timer.
    auto clkStop = std::chrono::high_resolution_clock::now();
    // Calculate difference between timer start and stop.
    auto clkDuration = std::chrono::duration_cast<std::chrono::milliseconds>(clkStop - clkStart);
    std::cout << "Elapsed time for " << chargeHits.getEvents().size() << " processed events is: "
              << clkDuration.count() << "ms\n";

    kalmanFit.openEventDisplay();

    return 0;
}
