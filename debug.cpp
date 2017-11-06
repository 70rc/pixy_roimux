#include <chrono>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>
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


void copyRawHisto(const TH2S &histo,
                  TH2D &histoCopy) {
    for (unsigned channel = 0; channel < histo.GetNbinsY(); ++channel) {
        for (unsigned sample = 0; sample < histo.GetNbinsX(); ++sample) {
            histoCopy.SetBinContent((sample + 1), (channel + 1), histo.GetBinContent((sample + 1), (channel + 1)));
        }
    }
}


void formatRawHisto(TH2D &histo,
                    const pixy_roimux::RunParams &runParams){
    histo.Scale(runParams.getAdcLsb());
    histo.GetXaxis()->SetLimits((0 * runParams.getSampleTime()), (histo.GetNbinsX() * runParams.getSampleTime()));
    histo.GetXaxis()->SetTitle(std::string("Sample [us]").c_str());
    histo.GetYaxis()->SetTitle(std::string("Channel #").c_str());
    histo.GetZaxis()->SetTitle(std::string("Amplitude [mV]").c_str());
}


void formatHitHisto(const std::shared_ptr<TH1D> &histo,
                    const pixy_roimux::Hit2d &hit,
                    const pixy_roimux::RunParams &runParams){
    unsigned xMin = static_cast<unsigned>(
            std::max((static_cast<int>(hit.firstSample) - static_cast<int>(hit.posPulseWidth)), 0));
    unsigned xMax = std::min((hit.lastSample + hit.posPulseWidth), static_cast<unsigned>(histo->GetNbinsX()));
    histo->GetXaxis()->SetRange(xMin, xMax);
    histo->Scale(runParams.getAdcLsb());
    histo->GetXaxis()->SetLimits(0, (histo->GetNbinsX() * runParams.getSampleTime()));
    histo->GetXaxis()->SetTitle(std::string("Sample [us]").c_str());
    histo->GetYaxis()->SetTitle(std::string("Amplitude [mV]").c_str());
}


void drawThresholds(TLine &line,
                    const pixy_roimux::Hit2d &hit,
                    const pixy_roimux::RunParams &runParams,
                    const bool bipolar = false) {
    double y1 = gPad->GetUymin();
    double y2 = gPad->GetUymax();
    double x = hit.firstSample * runParams.getSampleTime();
    line.SetLineColor(kGreen);
    line.SetLineStyle(2);
    line.DrawLine(x, y1, x, y2);
    x = hit.posPeakSample * runParams.getSampleTime();
    line.SetLineColor(kGreen);
    line.SetLineStyle(1);
    line.DrawLine(x, y1, x, y2);
    if (bipolar) {
        x = hit.zeroCrossSample * runParams.getSampleTime();
        line.SetLineColor(kOrange);
        line.SetLineStyle(1);
        line.DrawLine(x, y1, x, y2);
        x = hit.negPeakSample * runParams.getSampleTime();
        line.SetLineColor(kRed);
        line.SetLineStyle(1);
        line.DrawLine(x, y1, x, y2);
        line.SetLineColor(kRed);
        line.SetLineStyle(2);
    }
    else {
        line.SetLineColor(kGreen);
        line.SetLineStyle(2);
    }
    x = hit.lastSample * runParams.getSampleTime();
    line.DrawLine(x, y1, x, y2);
}


int main(int argc, char** argv) {
    // Start time point for timer.
    auto clkStart = std::chrono::high_resolution_clock::now();

    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " runParamsFileName dataFileName eventID geoFileName outputPath" << std::endl;
        exit(1);
    }
    const std::string runParamsFileName = std::string(argv[1]);
    const std::string dataFileName = std::string(argv[2]);
    const std::vector<unsigned> eventIds = {static_cast<unsigned>(std::stoul(argv[3]))};
    const std::string geoFileName = std::string(argv[4]);
    const std::string outputPath = std::string(argv[5]) + "/";
    const unsigned subrunId = 0;

    gStyle->SetOptStat(0);

    // Create the VIPER runParams containing all the needed run parameters.
    const pixy_roimux::RunParams runParams(runParamsFileName);

    const std::string th2DrawOpts("colz");
    const std::string th1DrawOpts("hist");
    TCanvas canvas("canvas", "canvas", 1920, 1080);
    TLegend legend(.7, .7, .9, .9);
    TLine line;
    line.SetLineWidth(1);
    TBox box;

    // Load events directly from ROOT file.
    std::cout << "Extracting chargeData...\n";
    pixy_roimux::ChargeData chargeData(dataFileName, eventIds, subrunId, runParams);

    std::cout << "Plotting unfiltered raw histograms...\n";
    canvas.SetRightMargin(.15);
    auto eventId = eventIds.cbegin();
    for (const auto &event : chargeData.getReadoutHistos()) {
        TH2D histoCopy("histoCopy", "histoCopy",
                       event.first.GetNbinsX(), 0, event.first.GetNbinsX(),
                       event.first.GetNbinsY(), 0, event.first.GetNbinsY());
        copyRawHisto(event.first, histoCopy);
        formatRawHisto(histoCopy, runParams);
        histoCopy.SetTitle(std::string("Event " + std::to_string(*eventId) + " Unfiltered Pixel Raw Data").c_str());
        canvas.cd();
        histoCopy.Draw(th2DrawOpts.c_str());
        canvas.Print(std::string(outputPath + "event" + std::to_string(*eventId) + "_rawUnfilteredPixel.pdf").c_str());
        histoCopy = TH2D("histoCopy", "histoCopy",
                         event.second.GetNbinsX(), 0, event.second.GetNbinsX(),
                         event.second.GetNbinsY(), 0, event.second.GetNbinsY());
        copyRawHisto(event.second, histoCopy);
        formatRawHisto(histoCopy, runParams);
        histoCopy.SetTitle(std::string("Event " + std::to_string(*eventId) + " Unfiltered ROI Raw Data").c_str());
        canvas.cd();
        histoCopy.Draw(th2DrawOpts.c_str());
        canvas.Print(std::string(outputPath + "event" + std::to_string(*eventId) + "_rawUnfilteredROI.pdf").c_str());
        ++eventId;
    }

    // Noise filter
    std::cout << "Filtering chargeData...\n";
    pixy_roimux::NoiseFilter noiseFilter(runParams);
    noiseFilter.filterData(chargeData);

    std::cout << "Plotting filtered raw histograms...\n";
    canvas.SetRightMargin(.15);
    eventId = eventIds.cbegin();
    for (const auto &event : chargeData.getReadoutHistos()) {
        TH2D histoCopy("histoCopy", "histoCopy",
                       event.first.GetNbinsX(), 0, event.first.GetNbinsX(),
                       event.first.GetNbinsY(), 0, event.first.GetNbinsY());
        copyRawHisto(event.first, histoCopy);
        formatRawHisto(histoCopy, runParams);
        histoCopy.SetTitle(std::string("Event " + std::to_string(*eventId) + " Filtered Pixel Raw Data").c_str());
        canvas.cd();
        histoCopy.Draw(th2DrawOpts.c_str());
        canvas.Print(std::string(outputPath + "event" + std::to_string(*eventId) + "_rawFilteredPixel.pdf").c_str());
        histoCopy = TH2D("histoCopy", "histoCopy",
                         event.second.GetNbinsX(), 0, event.second.GetNbinsX(),
                         event.second.GetNbinsY(), 0, event.second.GetNbinsY());
        copyRawHisto(event.second, histoCopy);
        formatRawHisto(histoCopy, runParams);
        histoCopy.SetTitle(std::string("Event " + std::to_string(*eventId) + " Filtered ROI Raw Data").c_str());
        canvas.cd();
        histoCopy.Draw(th2DrawOpts.c_str());
        canvas.Print(std::string(outputPath + "event" + std::to_string(*eventId) + "_rawFilteredROI.pdf").c_str());
        ++eventId;
    }

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
    canvas.SetRightMargin(.1);
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
            const pixy_roimux::Hit2d &pixelHit = event.pixelHits.at(pixelHitId);
            std::string plotFileName = outputPath;
            plotFileName += "event";
            plotFileName += std::to_string(event.eventId);
            plotFileName += "_pixel";
            plotFileName += std::to_string(pixelHit.channel);
            plotFileName += "_hit";
            plotFileName += std::to_string(pixelHitId);
            plotFileName += ".pdf";
            canvas.Print(std::string(plotFileName + '[').c_str());
            auto pixelHisto = std::shared_ptr<TH1D>(rawHistos->first.ProjectionX(
                    "pixelChannelHisto", (pixelHit.channel + 1), (pixelHit.channel + 1)));
            formatHitHisto(pixelHisto, pixelHit, runParams);
            pixelHisto->SetLineColor(kBlue);
            pixelHisto->SetTitle(std::string("Event " + std::to_string(event.eventId)
                                             + " Pixel " + std::to_string(pixelHit.channel)
                                             + " Hit " + std::to_string(pixelHitId)).c_str());
            canvas.cd();
            pixelHisto->Draw(th1DrawOpts.c_str());
            canvas.Update();
            drawThresholds(line, pixelHit, runParams);
            double x1 = gPad->GetUxmin();
            double x2 = gPad->GetUxmax();
            double mean = noiseParams->first.at(pixelHit.channel).first;
            double sigma = noiseParams->first.at(pixelHit.channel).second;
            double y = (mean + runParams.getDiscSigmaPixelLead() * sigma) * runParams.getAdcLsb();
            line.SetLineColor(kGreen);
            line.SetLineStyle(2);
            line.DrawLine(x1, y, x2, y);
            y = (mean + std::max((runParams.getDiscSigmaPixelPeak() * sigma), runParams.getDiscAbsPixelPeak()))
                * runParams.getAdcLsb();
            line.SetLineColor(kGreen);
            line.SetLineStyle(1);
            line.DrawLine(x1, y, x2, y);
            y = (mean + runParams.getDiscSigmaPixelTrail() * sigma) * runParams.getAdcLsb();
            line.SetLineColor(kGreen);
            line.SetLineStyle(2);
            line.DrawLine(x1, y, x2, y);
            legend.Clear();
            line.SetLineColor(kBlack);
            line.SetLineStyle(1);
            legend.AddEntry(line.Clone(), "Peak sample/threshold", "L");
            line.SetLineColor(kBlack);
            line.SetLineStyle(2);
            legend.AddEntry(line.Clone(), "Edge sample/threshold", "L");
            box.SetFillColor(kGreen);
            legend.AddEntry(box.Clone(), "Positive pulse", "F");
            legend.Draw();
            canvas.Print(plotFileName.c_str());
            for (const auto &roiHitId : event.pixel2roi.at(pixelHitId)) {
                std::string rejectionString;
                if (roiHitId == acceptedRoiHitId) {
                    rejectionString = "PCA Accepted";
                }
                else if (rejectedRoiHitIds.find(roiHitId) != rejectedRoiHitIds.cend()){
                    if (acceptedRoiHitId >= 0) {
                        rejectionString = "PCA Rejected Ambiguity";
                    }
                    else if (acceptedRoiHitId == -1) {
                        rejectionString = "PCA Unassigned";
                    }
                    else if (acceptedRoiHitId == - 2) {
                        rejectionString = "PCA Rejected Hit";
                    }
                    else {
                        rejectionString = "PCA Failure";
                    }
                }
                else {
                    rejectionString = "Duplicate 3D Hit";
                }
                const pixy_roimux::Hit2d &roiHit = event.roiHits.at(roiHitId);
                auto roiHisto = std::shared_ptr<TH1D>(rawHistos->second.ProjectionX(
                        "roiChannelHisto", (roiHit.channel + 1), (roiHit.channel + 1)));
                formatHitHisto(roiHisto, roiHit, runParams);
                roiHisto->SetLineColor(kBlue);
                roiHisto->SetTitle(std::string("Event " + std::to_string(event.eventId)
                                               + " ROI " + std::to_string(roiHit.channel)
                                               + " Hit " + std::to_string(roiHitId)
                                               + ' ' + rejectionString).c_str());
                canvas.cd();
                roiHisto->Draw(th1DrawOpts.c_str());
                canvas.Update();
                drawThresholds(line, roiHit, runParams, true);
                x1 = gPad->GetUxmin();
                x2 = gPad->GetUxmax();
                mean = noiseParams->second.at(roiHit.channel).first;
                sigma = noiseParams->second.at(roiHit.channel).second;
                y = (mean + runParams.getDiscSigmaRoiPosLead() * sigma) * runParams.getAdcLsb();
                line.SetLineColor(kGreen);
                line.SetLineStyle(2);
                line.DrawLine(x1, y, x2, y);
                y = (mean + std::max((runParams.getDiscSigmaRoiPosPeak() * sigma), runParams.getDiscAbsRoiPosPeak()))
                    * runParams.getAdcLsb();
                line.SetLineColor(kGreen);
                line.SetLineStyle(1);
                line.DrawLine(x1, y, x2, y);
                y = (mean - std::max((runParams.getDiscSigmaRoiNegPeak() * sigma), runParams.getDiscAbsRoiNegPeak()))
                    * runParams.getAdcLsb();
                line.SetLineColor(kRed);
                line.SetLineStyle(1);
                line.DrawLine(x1, y, x2, y);
                y = (mean - runParams.getDiscSigmaRoiNegTrail() * sigma) * runParams.getAdcLsb();
                line.SetLineColor(kRed);
                line.SetLineStyle(2);
                line.DrawLine(x1, y, x2, y);
                legend.Clear();
                line.SetLineColor(kBlack);
                line.SetLineStyle(1);
                legend.AddEntry(line.Clone(), "Peak sample/threshold", "L");
                line.SetLineColor(kBlack);
                line.SetLineStyle(2);
                legend.AddEntry(line.Clone(), "Edge sample/threshold", "L");
                box.SetFillColor(kGreen);
                legend.AddEntry(box.Clone(), "Positive pulse", "F");
                box.SetFillColor(kOrange);
                legend.AddEntry(box.Clone(), "Zero crossing", "F");
                box.SetFillColor(kRed);
                legend.AddEntry(box.Clone(), "Negative pulse", "F");
                legend.Draw();
                canvas.Print(plotFileName.c_str());
            }
            ++pcaId;
            canvas.Print(std::string(plotFileName + ']').c_str());
        }
        ++rawHistos;
        ++noiseParams;
    }

    std::cout << "Initialising Kalman Fitter...\n";
    pixy_roimux::KalmanFit kalmanFit(runParams, geoFileName, true);
    std::cout << "Running Kalman Fitter...\n";
    const std::string genfitTreeFileName = outputPath + "genfit.root";
    kalmanFit.fit(chargeHits, genfitTreeFileName);

    // Write chargeHits of events in eventIds vector to CSV files so we can plot them with viper3Dplot.py afterwards.
    unsigned nHitCandidates = 0;
    unsigned nAmbiguities = 0;
    unsigned nUnmatchedPixelHits = 0;
    // Loop through events.
    for (const auto& event : chargeHits.getEvents()) {
        std::cout << "Writing event number " << event.eventId << " to file...\n";
        // Compose CSV filename and open file stream.
        const std::string csvEventBaseFileName = outputPath + "event" + std::to_string(event.eventId);
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
    const std::string statsFileName = outputPath + "stats.txt";
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
