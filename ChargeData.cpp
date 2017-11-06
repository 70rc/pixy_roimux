//
// Created by damian on 6/3/17.
//

#include "ChargeData.h"


namespace pixy_roimux {
    ChargeData::ChargeData(
            const std::string t_rootFileName,
            const unsigned t_subrunId,
            const pixy_roimux::RunParams &t_runParams) :
            m_subrunId(t_subrunId),
            m_runParams(t_runParams) {
        int eventId = -1;
        const TH2S *indHisto = nullptr;
        const TH2S *colHisto = nullptr;

        // Open the ROOT file read-only.
        TFile rootFile(t_rootFileName.c_str(), "READ");
        // Loop through file.
        do {
            // Don't add the histo to the DAQ histo vector in the first pass as we don't have one yet.
            if (eventId >= 0) {
                m_daqHistos.emplace_back(std::pair<const TH2S *, const TH2S *>(indHisto, colHisto));
                m_eventIds.emplace_back(static_cast<unsigned>(eventId));
            }
            ++eventId;
            // Get the histogram pointers from the ROOT file using the GetObject method which performs basic sanity checks.
            std::cout << "Reading event number " << eventId << "...\n";
            std::string histoName;
            histoName = "Ind_" + std::to_string(eventId);
            rootFile.GetObject(histoName.c_str(), indHisto);
            histoName = "Col_" + std::to_string(eventId);
            rootFile.GetObject(histoName.c_str(), colHisto);
        }
            // Check whether we got something from the file and leave the loop if we didn't.
        while (indHisto && colHisto);

        // Convert the DAQ histos to readout histos.
        convertHistos();

        // Clear the DAQ histo vector as we don't need them anymore.
        m_daqHistos.clear();
        // After pointers are cleared, we can close the ROOT file.
        rootFile.Close();
    }


    ChargeData::ChargeData(
            const std::string t_rootFileName,
            const std::vector<unsigned> &t_eventIds,
            const unsigned t_subrunId,
            const pixy_roimux::RunParams &t_runParams) :
            m_eventIds(t_eventIds),
            m_subrunId(t_subrunId),
            m_runParams(t_runParams) {
        // Preallocate the DAQ histo vector for speed.
        m_daqHistos.resize(m_eventIds.size());
        // Open ROOT file read-only.
        TFile rootFile(t_rootFileName.c_str(), "READ");
        // DAQ histo vector iterator.
        auto daqHisto = m_daqHistos.begin();
        // Loop over event IDs.
        for (const auto &eventId : m_eventIds) {
            // Get the histogram pointers from the ROOT file using the GetObject method which performs basic sanity checks.
            std::cout << "Reading event number " << eventId << "...\n";
            std::string histoName;
            histoName = "Ind_" + std::to_string(eventId);
            const TH2S *indHisto = nullptr;
            rootFile.GetObject(histoName.c_str(), indHisto);
            histoName = "Col_" + std::to_string(eventId);
            const TH2S *colHisto = nullptr;
            rootFile.GetObject(histoName.c_str(), colHisto);
            // If we got something from the file, add its pointer to the DAQ histo vector.
            if (indHisto && colHisto) {
                *daqHisto = std::pair<const TH2S *, const TH2S *>(indHisto, colHisto);
            }
                // Else die.
            else {
                std::cerr << "ERROR: Failed to load event ID " << eventId
                          << " from file " << t_rootFileName << '!' << std::endl;
                exit(1);
            }
            // Increment the DAQ histo vector iterator.
            ++daqHisto;
        }

        // Convert the DAQ histos to readout histos.
        convertHistos();

        // Clear the DAQ histo vector as we don't need them anymore.
        m_daqHistos.clear();
        // After pointers are cleared, we can close the ROOT file.
        rootFile.Close();
    }


    ChargeData::ChargeData(
            const std::vector<std::pair<const TH2S *, const TH2S *>> &t_daqHistos,
            const std::vector<unsigned> &t_eventIds,
            const unsigned t_subrunId,
            const pixy_roimux::RunParams &t_runParams) :
            m_daqHistos(t_daqHistos),
            m_eventIds(t_eventIds),
            m_subrunId(t_subrunId),
            m_runParams(t_runParams) {
        // As we get all the necessary data passed to the constructor, we just need to convert the DAQ histos to readout
        // histos.
        convertHistos();
    }


    ChargeData::ChargeData(
            const TH2S *const t_indHisto,
            const TH2S *const t_colHisto,
            const unsigned t_eventId,
            const unsigned t_subrunId,
            const pixy_roimux::RunParams &t_runParams) :
            m_subrunId(t_subrunId),
            m_runParams(t_runParams) {
        // Fill the DAQ histos passed to the constructor into the DAQ histo vector, as the convertHistos() method reads from
        // this.
        m_daqHistos.emplace_back(std::pair<const TH2S *, const TH2S *>(t_indHisto, t_colHisto));
        m_eventIds.emplace_back(t_eventId);

        // Convert the DAQ histos to readout histos.
        convertHistos();

        // Clear the DAQ histo vector as we don't need them anymore.
        m_daqHistos.clear();
    }


    void ChargeData::convertHistos() {
        // Preallocate the readout histo vector for speed.
        m_readoutHistos.resize(m_daqHistos.size());

        // DAQ histo vector const iterator.
        auto daqHisto = m_daqHistos.cbegin();
        // Readout histo vector iterator.
        auto readoutHisto = m_readoutHistos.begin();
        // Loop over events.
        for (const auto &eventId : m_eventIds) {
            std::cout << "Converting event number " << eventId << "...\n";
            const std::string pixelHistoName = "pixelHisto_" + std::to_string(eventId);
            const std::string roiHistoName = "roiHisto_" + std::to_string(eventId);
            // Create the readout histos.
            *readoutHisto = std::pair<TH2S, TH2S>(
                    TH2S(pixelHistoName.c_str(), pixelHistoName.c_str(),
                         m_runParams.getNSamples(), 0, m_runParams.getNSamples(), m_runParams.getNPixels(), 0,
                         m_runParams.getNPixels()),
                    TH2S(roiHistoName.c_str(), roiHistoName.c_str(),
                         m_runParams.getNSamples(), 0, m_runParams.getNSamples(), m_runParams.getNRois(), 0,
                         m_runParams.getNRois()));
            // Assign the histo pointers.
            // DAQ histo pointers for simplicity.
            const TH2S &indHisto = *(daqHisto->first);
            const TH2S &colHisto = *(daqHisto->second);
            // Readout histo pointers for simplicity.
            TH2S &pixelHisto = readoutHisto->first;
            TH2S &roiHisto = readoutHisto->second;
            // Pay attention to the histo bin numbering!!! Loops (as everything else) start at 0, histos start at 1!!!
            // We use the pixel2daq and roi2daq methods of the viperMap to match the desired readout channels to DAQ
            // channels. Samples are copied 1 by 1 up to m_nSamples.
            // Loop through pixels.
            for (unsigned pixel = 0; pixel < m_runParams.getNPixels(); ++pixel) {
                // Loop through samples.
                for (unsigned sample = 0; sample < m_runParams.getNSamples(); ++sample) {
                    // Current pixel is in first DAQ histo (indHisto).
                    if (m_runParams.pixel2daq(pixel) < (m_runParams.getNChans() / 2)) {
                        // Fill pixel histo.
                        pixelHisto.SetBinContent((sample + 1), (pixel + 1), indHisto.GetBinContent(
                                (sample + 1), (m_runParams.pixel2daq(pixel) + 1)));
                    }
                        // Current pixel is in second DAQ histo (colHisto).
                    else {
                        // Fill ROI histo.
                        pixelHisto.SetBinContent((sample + 1), (pixel + 1), colHisto.GetBinContent(
                                (sample + 1), (m_runParams.pixel2daq(pixel) - m_runParams.getNChans() / 2 + 1)));
                    }
                }
            }
            // Loop through ROIs.
            for (unsigned roi = 0; roi < m_runParams.getNRois(); ++roi) {
                // Loop through samples.
                for (unsigned sample = 0; sample < m_runParams.getNSamples(); ++sample) {
                    // Fill ROI histo.
                    roiHisto.SetBinContent((sample + 1), (roi + 1), colHisto.GetBinContent(
                            (sample + 1), (m_runParams.roi2daq(roi) - m_runParams.getNChans() / 2 + 1)));
                }
            }
            // Increment DAQ and readout histo vector iterators.
            ++daqHisto;
            ++readoutHisto;
        }
    }
}