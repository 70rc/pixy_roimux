#include <chrono>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"
#include "ChargeData.h"
#include "ChargeHits.h"
#include "NoiseFilter.h"
#include "PrincipalComponentsCluster.h"
#include "RunParams.h"
#include "KalmanFit.h"


void plotRawHistos(const pixy_roimux::ChargeData &data,
                   const pixy_roimux::RunParams &runParams,
                   TCanvas &canvas,
                   const std::string &outputPath,
                   const std::string &nameSuffix) {
    canvas.SetRightMargin(.15);

    auto eventId = data.getEventIds().cbegin();
    for (const auto &histos : data.getReadoutHistos()) {
        for (unsigned channelType = 0; channelType < 2; ++channelType) {
            TH2D histoCopy("histoCopy", "histoCopy",
                           histos.at(channelType).GetNbinsX(), 0, histos.at(channelType).GetNbinsX(),
                           histos.at(channelType).GetNbinsY(), 0, histos.at(channelType).GetNbinsY());
            for (unsigned channel = 0; channel < histos.at(channelType).GetNbinsY(); ++channel) {
                for (unsigned sample = 0; sample < histos.at(channelType).GetNbinsX(); ++sample) {
                    histoCopy.SetBinContent((sample + 1), (channel + 1),
                                            histos.at(channelType).GetBinContent((sample + 1), (channel + 1)));
                }
            }
            histoCopy.Scale(runParams.getAdcLsb());
            histoCopy.GetXaxis()->SetLimits(0., (histoCopy.GetNbinsX() * runParams.getSampleTime()));
            histoCopy.GetXaxis()->SetTitle("Time [us]");
            histoCopy.GetYaxis()->SetTitle("Channel #");
            histoCopy.GetZaxis()->SetTitle("Amplitude [mV]");
            std::ostringstream histoString;
            histoString << "Event " << *eventId << ' ' << nameSuffix << ' '
                        << pixy_roimux::ChannelTypeString.at(channelType) << " Raw Data";
            histoCopy.SetTitle(histoString.str().c_str());
            canvas.cd();
            histoCopy.Draw("colz");
            histoString.str("");
            histoString << outputPath << "event" << *eventId << "_raw" << nameSuffix
                        << pixy_roimux::ChannelTypeString.at(channelType) << ".pdf";
            canvas.Print(histoString.str().c_str());
        }
        ++eventId;
    }

}


void plotHitHisto(const std::array<TH2S, 2> &rawHistos,
                  const std::array<std::vector<std::array<double, 2>>, 2> &noiseParams,
                  const pixy_roimux::Event &event,
                  const pixy_roimux::RunParams &runParams,
                  TCanvas &canvas,
                  const unsigned hitId,
                  const unsigned channelType,
                  const std::string &rejectionString,
                  const std::string &plotFileName) {
    canvas.SetRightMargin(.1);
    TLegend legend(.7, .7, .9, .9);
    TLine line;
    const unsigned lineWidth = 1;
    line.SetLineWidth(lineWidth);
    TBox box;

    const pixy_roimux::Hit2d &hit = ((channelType == pixy_roimux::kPixel)
                                     ? event.pixelHits.at(hitId) : event.roiHits.at(hitId));
    auto histo = std::shared_ptr<TH1D>(rawHistos.at(channelType).ProjectionX(
            "channelHisto", (hit.channel + 1), (hit.channel + 1)));
    unsigned xMin = static_cast<unsigned>(
            std::max((static_cast<int>(hit.firstSample) - static_cast<int>(hit.posPulseWidth)), 0));
    unsigned xMax = std::min((hit.lastSample + hit.posPulseWidth), static_cast<unsigned>(histo->GetNbinsX()));
    histo->GetXaxis()->SetRange(xMin, xMax);
    histo->Scale(runParams.getAdcLsb());
    histo->GetXaxis()->SetLimits(0, (histo->GetNbinsX() * runParams.getSampleTime()));
    histo->GetXaxis()->SetTitle("Time [us]");
    histo->GetYaxis()->SetTitle("Amplitude [mV]");
    histo->SetLineColor(kBlue);
    histo->SetLineWidth(lineWidth);
    std::ostringstream histoString;
    histoString << "Event " << event.eventId << ' ' << pixy_roimux::ChannelTypeString.at(channelType) << ' '
                << hit.channel << " Hit " << hitId << rejectionString;
    histo->SetTitle(histoString.str().c_str());
    canvas.cd();
    histo->Draw("hist");
    canvas.Update();
    double y1 = gPad->GetUymin();
    double y2 = gPad->GetUymax();
    double x1 = hit.firstSample * runParams.getSampleTime();
    double x2 = x1;
    line.SetLineColor(kGreen);
    line.SetLineStyle(2);
    line.DrawLine(x1, y1, x2, y2);
    x1 = hit.posPeakSample * runParams.getSampleTime();
    x2 = x1;
    line.SetLineColor(kGreen);
    line.SetLineStyle(1);
    line.DrawLine(x1, y1, x2, y2);
    if (channelType == pixy_roimux::kRoi) {
        x1 = hit.zeroCrossSample * runParams.getSampleTime();
        x2 = x1;
        line.SetLineColor(kOrange);
        line.SetLineStyle(1);
        line.DrawLine(x1, y1, x2, y2);
        x1 = hit.negPeakSample * runParams.getSampleTime();
        x2 = x1;
        line.SetLineColor(kRed);
        line.SetLineStyle(1);
        line.DrawLine(x1, y1, x2, y2);
        line.SetLineColor(kRed);
        line.SetLineStyle(2);
    }
    else {
        line.SetLineColor(kGreen);
        line.SetLineStyle(2);
    }
    x1 = hit.lastSample * runParams.getSampleTime();
    x2 = x1;
    line.DrawLine(x1, y1, x2, y2);
    x1 = gPad->GetUxmin();
    x2 = gPad->GetUxmax();
    double mean = noiseParams.at(channelType).at(hit.channel).at(0);
    double sigma = noiseParams.at(channelType).at(hit.channel).at(1);
    if (channelType == pixy_roimux::kPixel) {
        y1 = (mean + runParams.getDiscSigmaPixelLead() * sigma) * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kGreen);
        line.SetLineStyle(2);
        line.DrawLine(x1, y1, x2, y2);
        y1 = (mean + std::max((runParams.getDiscSigmaPixelPeak() * sigma), runParams.getDiscAbsPixelPeak()))
             * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kGreen);
        line.SetLineStyle(1);
        line.DrawLine(x1, y1, x2, y2);
        y1 = (mean + runParams.getDiscSigmaPixelTrail() * sigma) * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kGreen);
        line.SetLineStyle(2);
        line.DrawLine(x1, y1, x2, y2);
    }
    else {
        y1 = (mean + runParams.getDiscSigmaRoiPosLead() * sigma) * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kGreen);
        line.SetLineStyle(2);
        line.DrawLine(x1, y1, x2, y2);
        y1 = (mean + std::max((runParams.getDiscSigmaRoiPosPeak() * sigma), runParams.getDiscAbsRoiPosPeak()))
            * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kGreen);
        line.SetLineStyle(1);
        line.DrawLine(x1, y1, x2, y2);
        y1 = (mean - std::max((runParams.getDiscSigmaRoiNegPeak() * sigma), runParams.getDiscAbsRoiNegPeak()))
            * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kRed);
        line.SetLineStyle(1);
        line.DrawLine(x1, y1, x2, y2);
        y1 = (mean - runParams.getDiscSigmaRoiNegTrail() * sigma) * runParams.getAdcLsb();
        y2 = y1;
        line.SetLineColor(kRed);
        line.SetLineStyle(2);
        line.DrawLine(x1, y1, x2, y2);
    }
    legend.Clear();
    line.SetLineColor(kBlack);
    line.SetLineStyle(1);
    legend.AddEntry(line.Clone(), "Peak sample/threshold", "L");
    line.SetLineColor(kBlack);
    line.SetLineStyle(2);
    legend.AddEntry(line.Clone(), "Edge sample/threshold", "L");
    box.SetFillColor(kGreen);
    legend.AddEntry(box.Clone(), "Positive pulse", "F");
    if (channelType == pixy_roimux::kRoi) {
        box.SetFillColor(kOrange);
        legend.AddEntry(box.Clone(), "Zero crossing", "F");
        box.SetFillColor(kRed);
        legend.AddEntry(box.Clone(), "Negative pulse", "F");
    }
    legend.Draw();
    canvas.Print(plotFileName.c_str());
}


int main(int argc, char** argv) {
    // Start time point for timer.
    auto clkStart = std::chrono::high_resolution_clock::now();

    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " runParamsFileName dataFileName eventID geoFileName outputPath" << std::endl;
        exit(1);
    }
    const std::string runParamsFileName(argv[1]);
    const std::string dataFileName(argv[2]);
    const std::vector<unsigned> eventIds{static_cast<unsigned>(std::stoul(argv[3]))};
    const std::string geoFileName(argv[4]);
    const std::string outputPath = std::string(argv[5]) + "/";
    const unsigned subrunId = 0;

    gStyle->SetOptStat(0);

    // Create the VIPER runParams containing all the needed run parameters.
    const pixy_roimux::RunParams runParams(runParamsFileName);

    TCanvas canvas("canvas", "canvas", 1920, 1080);

    // Load events directly from ROOT file.
    std::cout << "Extracting chargeData...\n";
    pixy_roimux::ChargeData chargeData(dataFileName, eventIds, subrunId, runParams);

    std::cout << "Plotting unfiltered raw histograms...\n";
    plotRawHistos(chargeData, runParams, canvas, outputPath, "Unfiltered");

    // Noise filter
    std::cout << "Filtering chargeData...\n";
    pixy_roimux::NoiseFilter noiseFilter(runParams);
    noiseFilter.filterData(chargeData);

    std::cout << "Plotting filtered raw histograms...\n";
    plotRawHistos(chargeData, runParams, canvas, outputPath, "Filtered");

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


    std::cout << "Plotting hit histograms...\n";
    auto rawHistos = chargeData.getReadoutHistos().cbegin();
    auto noiseParams = chargeData.getNoiseParams().cbegin();
    for (const auto &event : chargeHits.getEvents()) {
        auto pcaId = event.pcaIds.cbegin();
        for (const auto &hitCandidates : event.hitCandidates) {
            if (hitCandidates.empty()) {
                ++pcaId;
                continue;
            }
            int acceptedRoiHitId;
            if (*pcaId < 0) {
                acceptedRoiHitId = *pcaId;
            }
            else {
                acceptedRoiHitId = hitCandidates.at(*pcaId).roiHitId;
            }
            std::set<unsigned> rejectedRoiHitIds;
            unsigned hitId = 0;
            for (const auto &hitCandidate : hitCandidates) {
                if (hitId != *pcaId) {
                    rejectedRoiHitIds.insert(hitCandidate.roiHitId);
                }
                ++hitId;
            }
            unsigned pixelHitId = hitCandidates.at(0).pixelHitId;
            std::ostringstream plotFileName;
            plotFileName << outputPath << "event" << event.eventId << "_pixel" << event.pixelHits.at(pixelHitId).channel
                         << "_hit" << pixelHitId << ".pdf";
            canvas.Print(std::string(plotFileName.str() + '[').c_str());
            plotHitHisto(*rawHistos, *noiseParams, event, runParams, canvas, pixelHitId, pixy_roimux::kPixel, "",
                         plotFileName.str());
            for (const auto &roiHitId : event.pixel2roi.at(pixelHitId)) {
                std::string rejectionString;
                if (roiHitId == acceptedRoiHitId) {
                    rejectionString = " PCA Accepted";
                }
                else if (rejectedRoiHitIds.find(roiHitId) != rejectedRoiHitIds.cend()){
                    if (acceptedRoiHitId >= 0) {
                        rejectionString = " PCA Rejected Ambiguity";
                    }
                    else if (acceptedRoiHitId == -1) {
                        rejectionString = " PCA Unassigned";
                    }
                    else if (acceptedRoiHitId == - 2) {
                        rejectionString = " PCA Rejected Hit";
                    }
                    else {
                        rejectionString = " PCA Failure";
                    }
                }
                else {
                    rejectionString = " Duplicate 3D Hit";
                }
                plotHitHisto(*rawHistos, *noiseParams, event, runParams, canvas, roiHitId, pixy_roimux::kRoi,
                             rejectionString, plotFileName.str());
            }
            ++pcaId;
            canvas.Print(std::string(plotFileName.str() + ']').c_str());
        }
        ++rawHistos;
        ++noiseParams;
    }

    std::cout << "Initialising Kalman Fitter...\n";
    pixy_roimux::KalmanFit kalmanFit(runParams, geoFileName, true);
    std::cout << "Running Kalman Fitter...\n";
    std::ostringstream genfitTreeFileName;
    genfitTreeFileName << outputPath << "genfit.root";
    kalmanFit.fit(chargeHits, genfitTreeFileName.str());

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
    unsigned long nEvents = chargeHits.getEvents().size();
    float averageHitCandidates = static_cast<float>(nHitCandidates) / static_cast<float>(nEvents);
    float averageAmbiguities = static_cast<float>(nAmbiguities) / static_cast<float>(nEvents);
    float averageUnmatchedPixelHits = static_cast<float>(nUnmatchedPixelHits) / static_cast<float>(nEvents);
    std::ostringstream statsFileName;
    statsFileName << outputPath << "stats.txt";
    std::ofstream statsFile(statsFileName.str(), std::ofstream::out);
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
