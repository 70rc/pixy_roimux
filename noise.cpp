#include <array>
#include <chrono>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "ChargeData.h"
#include "ChargeHits.h"
#include "NoiseFilter.h"
#include "PrincipalComponentsCluster.h"
#include "RunParams.h"
#include "KalmanFit.h"


void analyseNoise(const pixy_roimux::ChargeData &data,
                  std::ofstream &noiseParamsFile,
                  const pixy_roimux::RunParams &runParams,
                  const std::string &histoNamePrefix) {
    const unsigned nAdcBits = 14;
    const std::array<unsigned, 2> nChannels{runParams.getNPixels(), runParams.getNRois()};
    const std::array<std::set<unsigned>, 2> blindedChannels = {std::set<unsigned>{9}, std::set<unsigned>{}};
    std::array<TH2D, 2> channelNoiseHistos;
    std::array<TH1D, 2> noiseHistos;
    for (unsigned channelType = 0; channelType < 2; ++channelType) {
        std::ostringstream histoName;
        histoName << histoNamePrefix << pixy_roimux::ChannelTypeString.at(channelType) << "ChannelNoiseHisto";
        channelNoiseHistos.at(channelType) = TH2D(histoName.str().c_str(), histoName.str().c_str(), (1 << nAdcBits),
                                                  - (1 << (nAdcBits - 1)), (1 << (nAdcBits - 1)),
                                                  nChannels.at(channelType), 0, nChannels.at(channelType));
        histoName.str("");
        histoName << histoNamePrefix << pixy_roimux::ChannelTypeString.at(channelType) << "NoiseHisto";
        noiseHistos.at(channelType) = TH1D(histoName.str().c_str(), histoName.str().c_str(), (1 << nAdcBits),
                                           - (1 << (nAdcBits - 1)), (1 << (nAdcBits - 1)));
    }

    auto eventId = data.getEventIds().cbegin();
    for (const auto &event : data.getReadoutHistos()) {
        std::cout << "Analysing event number " << *eventId << "...\n";
        for (unsigned channelType = 0; channelType < 2; ++channelType) {
            for (unsigned channel = 0; channel < nChannels.at(channelType); ++channel) {
                const bool blinded = (blindedChannels.at(channelType).find(channel)
                                      != blindedChannels.at(channelType).cend());
                for (unsigned sample = 0; sample < runParams.getNSamples(); ++sample) {
                    double amplitude = event.at(channelType).GetBinContent((sample + 1), (channel + 1));
                    channelNoiseHistos.at(channelType).Fill(amplitude, channel);
                    if (!blinded) {
                        noiseHistos.at(channelType).Fill(amplitude);
                    }
                }
            }
        }
        ++eventId;
    }

    for (unsigned channelType = 0; channelType < 2; ++channelType) {
        channelNoiseHistos.at(channelType).Write();
        noiseHistos.at(channelType).Write();
    }

    for (unsigned channelType = 0; channelType < 2; ++channelType){
        for (unsigned channel = 0; channel <= nChannels.at(channelType); ++channel) {
            std::ostringstream channelString;
            channelString << pixy_roimux::ChannelTypeString.at(channelType) << ' ';
            std::unique_ptr<TH1D> histo;
            bool blinded = false;
            if (channel == nChannels.at(channelType)) {
                channelString << 'S';
                histo = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(noiseHistos.at(channelType).Clone()));
            }
            else {
                channelString << channel;
                histo = std::unique_ptr<TH1D>(channelNoiseHistos.at(channelType).ProjectionX(
                        "channelHisto", (channel + 1), (channel + 1)));
                blinded = (blindedChannels.at(channelType).find(channel) != blindedChannels.at(channelType).cend());
            }
            const double mean = histo->GetMean();
            const double stdDev = histo->GetStdDev();
            TF1 gauss("gauss", "gaus");
            //histo->GetXaxis()->SetRangeUser((mean - stdDev), (mean + stdDev));
            histo->Fit(&gauss, "QN");
            noiseParamsFile << channelString.str() << ":\t" << mean << ", " << stdDev << ", "
                            << gauss.GetParameter(1) << ", " << gauss.GetParameter(2);
            if (blinded) {
                noiseParamsFile << "\tBLINDED";
            }
            noiseParamsFile << std::endl;
        }
    }
}


int main(int argc, char** argv) {
    // Start time point for timer.
    auto clkStart = std::chrono::high_resolution_clock::now();

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " runParamsFileName dataFileName outputPath [nEvents]" << std::endl;
        exit(1);
    }
    const std::string runParamsFileName(argv[1]);
    const std::string dataFileName(argv[2]);
    const std::string outputPath = std::string(argv[3]) + '/';
    const unsigned nEvents = ((argc > 4) ? static_cast<unsigned>(std::stoul(argv[4])) : 0);
    const unsigned subrunId = 0;
    std::vector<unsigned> eventIds;
    for (unsigned event = 0; event < nEvents; ++event) {
        eventIds.emplace_back(event);
    }

    // Create the VIPER runParams containing all the needed run parameters.
    const pixy_roimux::RunParams runParams(runParamsFileName);

    // Load events directly from ROOT file.
    std::cout << "Extracting chargeData...\n";
    pixy_roimux::ChargeData chargeData = (eventIds.empty()
                                          ? pixy_roimux::ChargeData(dataFileName, subrunId, runParams)
                                          : pixy_roimux::ChargeData(dataFileName, eventIds, subrunId, runParams));

    std::ostringstream rootFileName;
    rootFileName << outputPath << "noise.root";
    TFile rootFile(rootFileName.str().c_str(), "RECREATE");
    if (!rootFile.IsOpen()) {
        std::cerr << "ERROR: Failed to open ROOT file " << rootFileName << '!' << std::endl;
        exit(1);
    }
    std::ostringstream noiseParamFileName;
    noiseParamFileName << outputPath << "noiseParams.csv";
    std::ofstream noiseParamsFile(noiseParamFileName.str(), std::ofstream::out);
    noiseParamsFile << "Noise Parameters" << std::endl;
    noiseParamsFile << "================" << std::endl;
    noiseParamsFile << "Histogram Mean, Histogram Standard Deviation, Gaussian Mu, Gaussian Sigma" << std::endl;

    std::cout << "Analysing noise of unfiltered data...\n";
    noiseParamsFile << std::endl;
    noiseParamsFile << "Unfiltered data:" << std::endl;
    noiseParamsFile << "----------------" << std::endl;
    analyseNoise(chargeData, noiseParamsFile, runParams, "unfiltered");

    // Noise filter
    std::cout << "Filtering chargeData...\n";
    pixy_roimux::NoiseFilter noiseFilter(runParams);
    noiseFilter.filterData(chargeData);

    std::cout << "Analysing noise of filtered data...\n";
    noiseParamsFile << std::endl;
    noiseParamsFile << "Filtered data:" << std::endl;
    noiseParamsFile << "--------------" << std::endl;
    analyseNoise(chargeData, noiseParamsFile, runParams, "filtered");

    rootFile.Close();
    noiseParamsFile.close();

    std::cout << "Done.\n";

    // Stop time point for timer.
    auto clkStop = std::chrono::high_resolution_clock::now();
    // Calculate difference between timer start and stop.
    auto clkDuration = std::chrono::duration_cast<std::chrono::milliseconds>(clkStop - clkStart);
    std::cout << "Elapsed time for " << chargeData.getReadoutHistos().size() << " processed events is: "
              << clkDuration.count() << "ms\n";

    return 0;
}
