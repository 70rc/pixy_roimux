#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cerrno>
#include <getopt.h>
#include "TFile.h"
#include "TTree.h"
#include "ChargeData.h"
#include "ChargeHits.h"
#include "NoiseFilter.h"
#include "PrincipalComponentsCluster.h"
#include "RunParams.h"
#include "KalmanFit.h"


void printUsage(const struct option *options) {
    std::cerr << "Usage: "
              << program_invocation_short_name
              << " --"  << options[0].name  << "=FILE"
              << " --"  << options[1].name  << "=FILE"
              << " --"  << options[2].name  << "=FILE"
              << " --"  << options[3].name  << "=DIR"
              << " [--" << options[4].name  << "=NUM]"
              << " [--" << options[5].name  << "=NUM]"
              << " [--" << options[6].name  << "=FILE]"
              << " [--" << options[7].name  << "=NUM]"
              << " [--" << options[8].name  << "=NUM]"
              << " [--" << options[9].name  << "=NUM]"
              << " [--" << options[10].name << "]"
              << " [--" << options[11].name << "]"
              << " [--" << options[12].name << "]"
              << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[0].val) << ", --" << options[0].name
              << "=FILE\tRead run parameters from FILE." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[1].val) << ", --" << options[1].name
              << "=FILE\t\tRead data from FILE." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[2].val) << ", --" << options[2].name
              << "=FILE\t\tRead geometry from FILE." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[3].val) << ", --" << options[3].name
              << "=DIR\tSave results to DIR." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[4].val) << ", --" << options[4].name
              << "=NUM\t\tProcess only NUM events." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[5].val) << ", --" << options[5].name
              << "=NUM\t\tStart with event number NUM." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[6].val) << ", --" << options[6].name
              << "=FILE\tRead ranking tree from FILE." <<std::endl;
    std::cerr << "\t-" << static_cast<char>(options[7].val) << ", --" << options[7].name
              << "=NUM\tOnly process events with ranking >= NUM." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[8].val) << ", --" << options[8].name
              << "=NUM\tOnly process events with ranking <= NUM." << std::endl;
    std::cerr << "\t-" << static_cast<char>(options[9].val) << ", --" << options[9].name
              << "=NUM\tSet subrun ID to NUM." << std::endl;
    std::cerr << "\t--" << options[10].name
              << "\t\tEnable GENFIT event display." << std::endl;
    std::cerr << "\t--" << options[11].name
              << "\t\tDisable principal components analysis." << std::endl;
    std::cerr << "\t--" << options[12].name
              << "\t\tDisable Kalman Filter." << std::endl;
}


int main(int argc, char** argv) {
    std::string runParamsFileName;
    std::string dataFileName;
    std::string geoFileName;
    std::string outputPath;
    unsigned nEvents = 0;
    unsigned firstEvent = 0;
    std::string rankingFileName;
    int minRanking = 4;
    int maxRanking = 4;
    unsigned subrunId = 0;
    int display = 0;
    int pca = 1;
    int kalman = 1;

    static struct option long_options[] =
            {
                    {"params",      required_argument,  nullptr,    'p'},
                    {"data",        required_argument,  nullptr,    'i'},
                    {"geo",         required_argument,  nullptr,    'g'},
                    {"output",      required_argument,  nullptr,    'o'},
                    {"nevt",        required_argument,  nullptr,    'n'},
                    {"first",       required_argument,  nullptr,    'f'},
                    {"ranktree",    required_argument,  nullptr,    'r'},
                    {"minrank",     required_argument,  nullptr,    'm'},
                    {"maxrank",     required_argument,  nullptr,    'M'},
                    {"subrun",      required_argument,  nullptr,    's'},
                    {"display",     no_argument,        &display,   1},
                    {"nopca",       no_argument,        &pca,       0},
                    {"nokalman",    no_argument,        &kalman,    0},
                    {nullptr, 0, nullptr, 0}
            };
    while(true) {
        int c;
        // getopt_long stores the option index here.
        int option_index = 0;
        c = getopt_long(argc, argv, "p:i:g:o:n:f:r:m:M:s:", long_options, &option_index);
        // Detect the end of the options.
        if (c == - 1) {
            break;
        }
        switch(c) {
            case 0:
                // Handled by long_options.
                break;
            case 'p':
                runParamsFileName = optarg;
                break;
            case 'i':
                dataFileName = optarg;
                break;
            case 'g':
                geoFileName = optarg;
                break;
            case 'o':
                outputPath = std::string(optarg) + '/';
                break;
            case 'n':
                nEvents = static_cast<unsigned>(std::stoul(optarg));
                break;
            case 'f':
                firstEvent = static_cast<unsigned>(std::stoul(optarg));
                break;
            case 'r':
                rankingFileName = optarg;
                break;
            case 'm':
                minRanking = std::stoi(optarg);
                break;
            case 'M':
                maxRanking = std::stoi(optarg);
                break;
            case 's':
                subrunId = static_cast<unsigned>(std::stoul(optarg));
                break;
            case '?':
                printUsage(long_options);
                exit(1);
            default:
                std::cerr << "Parsing error!" << std::endl;
                exit(1);
        }
    }
    if (runParamsFileName.empty() || dataFileName.empty() || geoFileName.empty() || outputPath.empty()) {
        printUsage(long_options);
        exit(1);
    }

    // Start time point for timer.
    auto clkStart = std::chrono::high_resolution_clock::now();

    std::vector<unsigned> eventIds;
    if (!rankingFileName.empty()) {
        TFile rankingFile(rankingFileName.c_str(), "READ");
        if (!rankingFile.IsOpen()) {
            std::cerr << "ERROR: Failed to open ranking file " << rankingFileName << '!' << std::endl;
            exit(1);
        }
        TTree *rankingTree = nullptr;
        rankingFile.GetObject("RankingTime", rankingTree);
        if (!rankingTree) {
            std::cerr << "ERROR: Failed to find \"RankingTime\" tree in ranking file " << rankingFileName << '!'
                      << std::endl;
            exit(1);
        }
        int ranking;
        rankingTree->SetBranchAddress("Ranking", &ranking);
        for (unsigned event = firstEvent; event < rankingTree->GetEntries(); ++event) {
            rankingTree->GetEvent(event);
            if ((ranking >= minRanking) && (ranking <= maxRanking)) {
                eventIds.emplace_back(event);
            }
            if ((eventIds.size() >= nEvents) && (nEvents > 0)) {
                break;
            }
        }
        rankingFile.Close();
        if ((eventIds.size() < nEvents) || eventIds.empty()) {
            std::cerr << "ERROR: Failed to find " << ((nEvents == 0) ? "any" : std::to_string(nEvents))
                      << " events with " << minRanking << " <= ranking <= " << maxRanking << '!' << std::endl;
            exit(1);
        }
    }
    else {
        for (unsigned event = firstEvent; event < (firstEvent + nEvents); ++event) {
            eventIds.emplace_back(event);
        }
    }

    // Create the VIPER runParams containing all the needed run parameters.
    const pixy_roimux::RunParams runParams(runParamsFileName);

    // Load events directly from ROOT file.
    std::cout << "Extracting chargeData...\n";
    pixy_roimux::ChargeData chargeData = (eventIds.empty()
                                          ? pixy_roimux::ChargeData(dataFileName, subrunId, runParams)
                                          : pixy_roimux::ChargeData(dataFileName, eventIds, subrunId, runParams));

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
    if (pca) {
        std::cout << "Running principle components analysis...\n";
        principalComponentsCluster.analyseEvents(chargeHits);
    }

    std::cout << "Initialising Kalman Fitter...\n";
    pixy_roimux::KalmanFit kalmanFit(runParams, geoFileName, static_cast<bool>(display));
    if (kalman) {
        std::cout << "Running Kalman Fitter...\n";
        std::ostringstream genfitTreeFileName;
        genfitTreeFileName << outputPath << "genfit.root";
        kalmanFit.fit(chargeHits, genfitTreeFileName.str());
    }

    // Write chargeHits of events in eventIds vector to CSV files so we can plot them with viper3Dplot.py afterwards.
    unsigned nHitCandidates = 0;
    unsigned nAmbiguities = 0;
    unsigned nUnmatchedPixelHits = 0;
    // Loop through events.
    for (const auto& event : chargeHits.getEvents()) {
        std::cout << "Writing event number " << event.eventId << " to file...\n";
        // Compose CSV filename and open file stream.
        std::ostringstream csvEventBaseFileName;
        csvEventBaseFileName << outputPath << "event" << event.eventId;
        std::ostringstream csvHitsFileName;
        csvHitsFileName << csvEventBaseFileName.str() << "_hits.csv";
        std::ofstream csvHitsFile(csvHitsFileName.str(), std::ofstream::out);
        csvHitsFile << "X,Y,Z,Q,A" << std::endl;
        // viperEvents.hitCandidates is a vector of vectors so we need two loops.
        // Outer loop. Actually loops through pixel chargeHits of current event.
        auto pcaId = event.pcaIds.cbegin();
        for (const auto& hitCandidates : event.hitCandidates) {
            // Inner loop. Loops through all hit candidates of current pixel hit.
            int hitId = 0;
            for (const auto &hit : hitCandidates) {
                int reject = -1;
                if (pcaId < event.pcaIds.cend()) {
                    if (*pcaId == hitId) {
                        reject = 0;
                    }
                    else if (*pcaId == -2) {
                        reject = 1;
                    }
                    else if (*pcaId >= 0) {
                        reject = 2;
                    }
                }
                // Append coordinates and charge to file.
                csvHitsFile << hit.x << ',' << hit.y << ',' << hit.z << ',' << hit.chargeInt << ',' << reject << std::endl;
                ++hitId;
            }
            ++pcaId;
        }
        // Close CSV file.
        csvHitsFile.close();
        std::ostringstream csvPcaFileName;
        csvPcaFileName << csvEventBaseFileName.str() << "_pca.csv";
        std::ofstream csvPcaFile(csvPcaFileName.str(), std::ofstream::out);
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
    const unsigned long nProcessedEvents = chargeHits.getEvents().size();
    const double averageHitCandidates = static_cast<double>(nHitCandidates) / static_cast<double>(nProcessedEvents);
    const double averageAmbiguities = static_cast<double>(nAmbiguities) / static_cast<double>(nProcessedEvents);
    const double averageUnmatchedPixelHits = static_cast<double>(nUnmatchedPixelHits)
                                             / static_cast<double>(nProcessedEvents);
    std::ostringstream statsFileName;
    statsFileName << outputPath << "stats.txt";
    std::ofstream statsFile(statsFileName.str(), std::ofstream::out);
    statsFile << "Number of events processed: " << nProcessedEvents << std::endl;
    statsFile << "Average number of hit candidates per event: " << averageHitCandidates << std::endl;
    statsFile << "Average number of ambiguities per event: " << averageAmbiguities << std::endl;
    statsFile << "Average number of unmatched pixel chargeHits per event: " << averageUnmatchedPixelHits << std::endl;
    statsFile.close();

    std::cout << "Done.\n";

    // Stop time point for timer.
    auto clkStop = std::chrono::high_resolution_clock::now();
    // Calculate difference between timer start and stop.
    auto clkDuration = std::chrono::duration_cast<std::chrono::milliseconds>(clkStop - clkStart);
    std::cout << "Elapsed time for " << nProcessedEvents << " processed events is: "
              << clkDuration.count() << "ms\n";

    if (display) {
        kalmanFit.openEventDisplay();
    }

    return 0;
}
