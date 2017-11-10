//
// Created by damian on 6/3/17.
//

#include "ChargeHits.h"


namespace pixy_roimux {
    void ChargeHits::find2dHits(
            const TH2S &t_histo,
            const std::vector<std::array<double, 2>> & t_noiseParams,
            std::vector<Hit2d> &t_hits,
            std::multimap<unsigned, unsigned> &t_hitOrderLead,
            std::multimap<unsigned, unsigned> &t_hitOrderTrail,
            unsigned &t_nMissed,
            const bool t_bipolar,
            const double t_discSigmaPosLead,
            const double t_discSigmaPosPeak,
            const double t_discAbsPosPeak,
            const double t_discSigmaPosTrail,
            const double t_discSigmaNegLead,
            const double t_discSigmaNegPeak,
            const double t_discAbsNegPeak,
            const double t_discSigmaNegTrail) {
        // Clear the hit vector and maps from potential old data.
        t_hits.clear();
        t_hitOrderLead.clear();
        t_hitOrderTrail.clear();
        // Index of the hits vector needed for the ordered maps.
        unsigned hitId = 0;
        // Loop over all channels of the input histo.
        // Pay attention to bin numbers!!! Loops (and everything else) start at 0, histos start at 1!!!
        auto noiseParams = t_noiseParams.cbegin();
        for (unsigned channel = 0; channel < t_histo.GetNbinsY(); ++channel) {
            // Get the histo of a single channel using the ProjectionX method of TH2.
            auto channelHisto = std::shared_ptr<TH1D>(
                    t_histo.ProjectionX("channelHisto", (channel + 1), (channel + 1)));
            //if (noiseParams.first < 1.) {
            //    noiseParams.first = 1.;
            //}
            double noiseBaseline    = noiseParams->at(0);
            double thrPosLead       = noiseParams->at(0) + t_discSigmaPosLead * noiseParams->at(1);
            double thrPosPeak       = noiseParams->at(0) + std::max(t_discSigmaPosPeak * noiseParams->at(1), t_discAbsPosPeak);
            double thrPosTrail      = noiseParams->at(0) + t_discSigmaPosTrail * noiseParams->at(1);
            //double thrNegLead       = noiseParams->first - t_discSigmaNegLead * noiseParams->second;
            double thrNegPeak       = noiseParams->at(0) - std::max(t_discSigmaNegPeak * noiseParams->at(1), t_discAbsNegPeak);
            double thrNegTrail      = noiseParams->at(0) - t_discSigmaNegTrail * noiseParams->at(1);
            //std::cout << "noiseBaseline " << noiseBaseline << std::endl;
            //std::cout << "thrPosLead " << thrPosLead << std::endl;
            //std::cout << "thrPosPeak " << thrPosPeak << std::endl;
            //std::cout << "thrPosTrail " << thrPosTrail << std::endl;
            //std::cout << "thrNegPeak " << thrNegPeak << std::endl;
            //std::cout << "thrNegTrail " << thrNegTrail << std::endl;
            unsigned posPeakSample = static_cast<unsigned>(channelHisto->GetMaximumBin()) - 1;
            int posPeakValue = static_cast<int>(channelHisto->GetBinContent(posPeakSample + 1));
            while (posPeakValue >= thrPosPeak) {
                // Start and end of the pulse.
                int firstSample;
                int zeroCrossSample;
                int negPeakSample;
                int negPeakValue = 0;
                int lastSample;
                bool foundFirstSample = false;
                bool foundZeroCrossSample = false;
                bool crossedThrNegPeak = false;
                bool foundLastSample = false;
                // Loop to find the first/last sample of the pulse using constant fraction discrimination.
                // We only search in the specified range.
                for (int sampleOffset = 1; sampleOffset <= m_runParams.getDiscRange(); ++sampleOffset) {
                    // First sample. Only look if we haven't found it yet.
                    if (!foundFirstSample) {
                        firstSample = posPeakSample - sampleOffset;
                        // Check whether we've crossed the constant fraction threshold.
                        if ((firstSample >= 0) &&
                            (channelHisto->GetBinContent(firstSample + 1) < thrPosLead)) {
                            foundFirstSample = true;
                        }

                    }
                    // Last sample. Only look if we haven't found it yet.
                    if (!foundLastSample) {
                        lastSample = posPeakSample + sampleOffset;
                        // Check whether we've crossed the constant fraction threshold.
                        if ((lastSample < channelHisto->GetNbinsX()) &&
                            (channelHisto->GetBinContent(lastSample + 1) < thrPosTrail)) {
                            foundLastSample = true;
                        }
                    }
                }
                if (t_bipolar && foundLastSample) {
                    foundLastSample = false;
                    for (int sample = lastSample;
                         (sample <= (posPeakSample + 3 * m_runParams.getDiscRange())) && (sample < channelHisto->GetNbinsX());
                         ++sample) {
                        int binContent = static_cast<int>(channelHisto->GetBinContent(sample + 1));
                        if (!foundZeroCrossSample) {
                            if (binContent < noiseBaseline) {
                                zeroCrossSample = sample;
                                foundZeroCrossSample = true;
                            }
                        } else if (!crossedThrNegPeak) {
                            if (binContent <= thrNegPeak) {
                                crossedThrNegPeak = true;
                            }
                        }
                        // Last sample. Only look if we haven't found it yet.
                        if (!foundLastSample) {
                            lastSample = sample;
                            // Check whether we've crossed the constant fraction threshold.
                            if (crossedThrNegPeak &&
                                (binContent > thrNegTrail)) {
                                foundLastSample = true;
                            } else if (binContent < negPeakValue) {
                                negPeakSample = sample;
                                negPeakValue = binContent;
                            }
                        }
                    }
                }
                // If we detected both the rising and the falling edge, build a 2D hit.
                if (foundFirstSample && foundLastSample) {
                    Hit2d hit;
                    hit.channel = channel;
                    hit.firstSample = static_cast<unsigned>(firstSample);
                    hit.lastSample = static_cast<unsigned>(lastSample);
                    hit.posPeakSample = posPeakSample;
                    hit.posPulseHeight = posPeakValue;
                    if (t_bipolar) {
                        hit.zeroCrossSample = static_cast<unsigned>(zeroCrossSample);
                        hit.negPeakSample = static_cast<unsigned>(negPeakSample);
                        hit.negPulseHeight = negPeakValue;
                        hit.posPulseWidth = hit.zeroCrossSample - hit.firstSample;
                        hit.negPulseWidth = hit.lastSample - hit.zeroCrossSample + 1;
                    } else {
                        hit.zeroCrossSample = 0;
                        hit.negPeakSample = 0;
                        hit.negPulseHeight = 0;
                        hit.posPulseWidth = hit.lastSample - hit.firstSample + 1;
                        hit.negPulseWidth = 0;
                        //std::cout << "posPulseWidth " << hit.posPulseWidth << std::endl;
                    }
                    hit.pulseIntegral = 0;
                    hit.pulseRaw.clear();
                    hit.pulseRaw.resize(hit.posPulseWidth + hit.negPulseWidth);
                    // Raw pulse data vector iterator.
                    auto pulseRaw = hit.pulseRaw.begin();
                    // Loop over all pulse samples to store and integrate them.
                    for (unsigned sample = hit.firstSample; sample <= hit.lastSample; ++sample) {
                        const int pulseData = static_cast<int>(channelHisto->GetBinContent(sample + 1));
                        hit.pulseIntegral += pulseData;
                        *pulseRaw = pulseData;
                        // Increment the raw pulse data vector iterator.
                        ++pulseRaw;
                    }
                    // Push the hit to the hits vector.
                    t_hits.push_back(hit);
                    // Insert the hit ID into the hit order maps. The key is the sample where the signal rises/falls
                    // above/below the constant fraction.
                    t_hitOrderLead.insert({firstSample, hitId});
                    t_hitOrderTrail.insert({lastSample, hitId});
                    // Increment the hit ID.
                    ++hitId;
                } else {
                    ++t_nMissed;
                }
                if (firstSample < 0) {
                    firstSample = 0;
                }
                if (lastSample >= channelHisto->GetNbinsX()) {
                    lastSample = channelHisto->GetNbinsX() - 1;
                }
                for (unsigned sample = firstSample; sample <= lastSample; ++sample) {
                    channelHisto->SetBinContent((sample + 1), noiseBaseline);
                }
                posPeakSample = static_cast<unsigned>(channelHisto->GetMaximumBin()) - 1;
                posPeakValue = static_cast<int>(channelHisto->GetBinContent(posPeakSample + 1));
            }
            ++noiseParams;
        }
    }


    void ChargeHits::find3dHits(Event &t_event) {
        // Clear the match vectors from potential old data.
        t_event.pixel2roi.clear();
        t_event.roi2pixel.clear();
        // Preallocate the match vectors.
        t_event.pixel2roi.resize(t_event.pixelHits.size());
        t_event.roi2pixel.resize(t_event.roiHits.size());
        // Loop through ROI hits, sorted by rising pulse edge.
        for (const auto &roiHitOrderEntry : t_event.roiHitOrderLead) {
            // Get the ROI hit ID from the map entry.
            const unsigned roiHitId = roiHitOrderEntry.second;
            // Loop over pixel hits, sorted by falling pulse edge.
            // First loop over all pixel hits falling in between the risng and the falling edge of the current ROI hit.
            for (auto pixelHitOrderEntry = t_event.pixelHitOrderTrail.lower_bound(
                    t_event.roiHits.at(roiHitId).firstSample);
                 pixelHitOrderEntry != t_event.pixelHitOrderTrail.upper_bound(t_event.roiHits.at(roiHitId).lastSample);
                 ++pixelHitOrderEntry) {
                // Get the pixel hit ID from the map entry.
                const unsigned pixelHitId = pixelHitOrderEntry->second;
                // Append matches to the match vectors.
                // Because we're currently inside the ROI pulse, we're sure this is an actual match.
                t_event.pixel2roi.at(pixelHitId).push_back(roiHitId);
                t_event.roi2pixel.at(roiHitId).push_back(pixelHitId);
            }
            // Now loop from the end of the ROI pulse until twice the peak finding range m_discRange after the end of the
            // pulse. This is the maximum length a pixel pulse can have. Thus, outside this range, a match to this ROI hit
            // is not possible.
            for (auto pixelHitOrderEntry = t_event.pixelHitOrderTrail.lower_bound(
                    t_event.roiHits.at(roiHitId).lastSample);
                 pixelHitOrderEntry !=
                 t_event.pixelHitOrderTrail.upper_bound(t_event.roiHits.at(roiHitId).lastSample + 2 * m_runParams.getDiscRange());
                 ++pixelHitOrderEntry) {
                // Get the pixel hit ID from the map entry.
                const unsigned pixelHitId = pixelHitOrderEntry->second;
                // Because we're no longer inside the ROI pulse, we need to check whether there's an actual overlap between
                // the pixel pulse and the ROI pulse.
                if (t_event.pixelHits.at(pixelHitId).firstSample <= t_event.roiHits.at(roiHitId).lastSample) {
                    // If they actually overlap, append the match to the match vectors.
                    t_event.pixel2roi.at(pixelHitId).push_back(roiHitId);
                    t_event.roi2pixel.at(roiHitId).push_back(pixelHitId);
                }
            }
        }
    }


    void ChargeHits::buildHitCandidates(Event &t_event) {
        // Clear the hit candidate vector from potential old data.
        t_event.hitCandidates.clear();
        // Preallocate the hit candidate vector.
        t_event.hitCandidates.resize(t_event.pixelHits.size());
        // Pixel hits vector const iterator.
        auto pixelHit = t_event.pixelHits.cbegin();
        unsigned pixelHitId = 0;
        // Hit candidates vector iterator.
        auto hitCandidate = t_event.hitCandidates.begin();
        // Loop over all pixel hits in the pixel to ROI map.
        for (const auto &candidateRoiHitIds : t_event.pixel2roi) {
            // Get the pixel ID from the pixel hits vector.
            const unsigned pixelId = pixelHit->channel;
            // Check for duplicate 3dHits
            std::set<unsigned> roiIds;
            // Loop over all ROI hit IDs matched to the current pixel hit.
            for (const auto &candidateRoiHitId : candidateRoiHitIds) {
                // Get the ROI ID from the roi hits vector of the event using the ROI hit ID from the pixel to ROI map.
                const unsigned roiId = t_event.roiHits.at(candidateRoiHitId).channel;
                if (roiIds.find(roiId) != roiIds.cend()) {
                    continue;
                }
                roiIds.insert(roiId);
                // Build the 3D hit using the coordinates and calibration constants from the RunParams.
                Hit3d hit;
                // Calculate x,y,z in the units of pixel pitch and drift speed times drift time.
                // Origin in the center.
                hit.x = static_cast<float>((m_runParams.getRoiCoor(roiId, 0) + m_runParams.getPixelCoor(pixelId, 0))
                                           * m_runParams.getPixelPitch() + m_runParams.getTpcOrigin().at(0));
                hit.y = static_cast<float>((m_runParams.getRoiCoor(roiId, 1) + m_runParams.getPixelCoor(pixelId, 1))
                                           * m_runParams.getPixelPitch() + m_runParams.getTpcOrigin().at(1));
                hit.z = static_cast<float>((pixelHit->posPeakSample - m_runParams.getAnodeSample())
                                           * m_runParams.getSampleTime() * m_runParams.getDriftSpeed()
                                           + m_runParams.getTpcOrigin().at(2));
                // Calculate charge in C.
                hit.chargeInt = static_cast<float>(pixelHit->pulseIntegral * m_runParams.getSampleTime()
                                                   * (m_runParams.getAdcLsb() / m_runParams.getPreampTransimpedance()));
                hit.chargePeak = static_cast<float>(pixelHit->posPulseHeight
                                                    * (m_runParams.getAdcLsb() / m_runParams.getPreampGain()));
                hit.pixelHitId = pixelHitId;
                hit.roiHitId = candidateRoiHitId;
                hitCandidate->push_back(hit);
            }
            // Increment the pixel hits and 3D hit candidates vector iterators.
            ++pixelHit;
            ++pixelHitId;
            ++hitCandidate;
        }
    }


    void ChargeHits::findHits(const bool t_bipolarRoiHits) {
        // Clear events vector in case there's old data in it.
        m_events.clear();
        // Preallocate fHits for speed.
        m_events.resize(m_chargeData.getReadoutHistos().size());
        // Raw data vector const iterator. We'll read the raw data from there event by event.
        // Gives us access to an array of histos containing the pixel (0) and the ROI (1) data.
        auto eventData = m_chargeData.getReadoutHistos().cbegin();
        auto eventNoise = m_chargeData.getNoiseParams().cbegin();
        // Events vector iterator. We'll store the events we built from the raw data in there.
        auto event = m_events.begin();
        // Loop over all events using the event IDs vector.
        for (const auto &eventId : m_chargeData.getEventIds()) {
            std::cout << "Processing event number " << eventId << ":\n";
            std::cout << "Running 2D hit finder...\n";
            unsigned nMissedPixelHits = 0;
            unsigned nMissedRoiHits = 0;
            // Find pixel hits.
            find2dHits(eventData->at(kPixel),
                       eventNoise->at(kPixel),
                       event->pixelHits,
                       event->pixelHitOrderLead,
                       event->pixelHitOrderTrail,
                       nMissedPixelHits,
                       false,
                       m_runParams.getDiscSigmaPixelLead(),
                       m_runParams.getDiscSigmaPixelPeak(),
                       m_runParams.getDiscAbsPixelPeak(),
                       m_runParams.getDiscSigmaPixelTrail(),
                       0,
                       0,
                       0,
                       0);
            std::cout << "Found " << event->pixelHits.size() << " pixel hits.\n";
            std::cout << "Missed " << nMissedPixelHits << " pixel hits.\n";
            // Find ROI hits.
            find2dHits(eventData->at(kRoi),
                       eventNoise->at(kRoi),
                       event->roiHits,
                       event->roiHitOrderLead,
                       event->roiHitOrderTrail,
                       nMissedRoiHits,
                       t_bipolarRoiHits,
                       m_runParams.getDiscSigmaRoiPosLead(),
                       m_runParams.getDiscSigmaRoiPosPeak(),
                       m_runParams.getDiscAbsRoiPosPeak(),
                       m_runParams.getDiscSigmaRoiPosTrail(),
                       m_runParams.getDiscSigmaRoiNegLead(),
                       m_runParams.getDiscSigmaRoiNegPeak(),
                       m_runParams.getDiscAbsRoiNegPeak(),
                       m_runParams.getDiscSigmaRoiNegTrail());
            std::cout << "Found " << event->roiHits.size() << " ROI hits.\n";
            std::cout << "Missed " << nMissedRoiHits << " ROI hits.\n";

            std::cout << "Running 3D hit finder...\n";
            // Search for matches between pixel and ROI 2D hits.
            find3dHits(*event);
            // Build 3D hit candidates from the matches.
            buildHitCandidates(*event);
            // Some statistics.
            // Number of hit candidates for this event.
            unsigned nHitCandidates = 0;
            // Number of ambiguous hit candidates for this event.
            unsigned nAmbiguities = 0;
            // Number of pixel hits for this event that couldn't be matched to any ROI hits.
            unsigned nUnmatchedPixelHits = 0;
            // Loop over all pixel hits using the pixel to ROI hit map.
            for (const auto &hitCandidates : event->hitCandidates) {
                // Add number of ROI hit candidates for current pixel hit.
                nHitCandidates += hitCandidates.size();
                // If there's no ROI hit candidates, increment the unmatched counter.
                if (hitCandidates.empty()) {
                    ++nUnmatchedPixelHits;
                }
                    // If there's more than one ROI hit candidate, increment the ambiguity counter.
                else if (hitCandidates.size() > 1) {
                    ++nAmbiguities;
                }
            }
            std::cout << "Found " << nHitCandidates << " 3D hit candidates.\n";
            std::cout << "Found " << nAmbiguities << " ambiguities.\n";
            std::cout << "Failed to match " << nUnmatchedPixelHits << " pixel hits.\n";

            // Store run, subrun and event ID to the event struct.
            event->runId = m_runParams.getRunId();
            event->subrunId = m_chargeData.getSubrunId();
            event->eventId = eventId;

            // Increment raw data and events vector iterators.
            ++eventData;
            ++eventNoise;
            ++event;
        }
    }
}
