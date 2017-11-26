//
// Created by damian on 6/3/17.
//

#ifndef PIXY_ROIMUX_CHARGEHITS_H
#define PIXY_ROIMUX_CHARGEHITS_H


#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "TH1D.h"
#include "TH2S.h"
#include "TSpectrum.h"
#include "ChargeData.h"
#include "Event.h"
#include "NoiseFilter.h"
#include "RunParams.h"


namespace pixy_roimux {
    /// \brief Class implementing a hit finder generating 3D spacepoints.
    ///
    /// As input it needs a ChargeData object containing the raw data and a RunParams object for the pixel and ROI
    /// coordinates and the peak finder thresholds. For the pulse detection, the maximum of a channel histogram is
    /// first compared against a threshold in terms of RMS of the noise. In case it is above threshold, the pulse is
    /// scanned until the leading and trailing edge have fallen below a respective threshold. In case of a bipolar pulse,
    /// the scan goes on to find the zero crossing, the negative peak and then the trailing edge crossing another
    /// threshold. If all thresholds were found, the pulse is stored and deleted from the original histogram. If not,
    /// the scanned range is deleted from the original histogram. Then, the whole process starts again by finding the
    /// next new maximum in the histogram. The process stops when the new maximum of the histogram is below the peak
    /// threshold. After this, the found pixel pulses are searched for overlaps with ROI pulses and all matches are
    /// written to both a vector that contains the indices of all matched ROI pulses for a pixel pulse and vice versa.
    /// Besides, the raw pulse, its integral and 3D hit candidates are calculated. All the data is stored in a vector
    /// containing one Event struct per event. They can be accessed by const reference or by reference.
    class ChargeHits {
    public:

        /// \brief Constructor.
        /// \param t_chargeData raw data.
        /// \param t_runParams run parameters.
        ChargeHits(
                const ChargeData &t_chargeData,
                const RunParams &t_runParams) :
                m_chargeData(t_chargeData),
                m_runParams(t_runParams) {
        }

        /// \brief Run the hit finder.
        void findHits();

        /// \brief Get the vector containing all the hits of all events.
        /// \return events.
        std::vector<Event> &getEvents() {
            return m_events;
        }

        /// \brief Get the vector containing all the hist of all events as const reference.
        /// \return events.
        const std::vector<Event> &getEvents() const {
            return m_events;
        }


    private:

        /// \brief Private method used internally to find all the 2D hits in a histrogram.
        /// \param t_histo 2D histogram containing the raw data.
        /// \param t_noiseParams noise parameters used for thresholds.
        /// \param t_hits vector to store the hits.
        /// \param t_hitOrderLead hit indices ordered by leading pulse edge.
        /// \param t_hitOrderTrail hit indices ordered by trailing pulse edge.
        /// \param t_nMissed number of missed pulses for which a above threshold was detected but not all thresholds.
        /// \param t_bipolar search for bipolar pulses (used for ROI pulses).
        /// \param t_discSigmaPosLead threshold in noise RMS for the discrimination of the leading edge of the positive pulse.
        /// \param t_discSigmaPosPeak threshold in noise RMS for the discrimination of the peak of the positive pulse.
        /// \param t_discAbsPosPeak absolute threshold for the discrimination of the peak of a positive pulse
        /// \param t_discSigmaPosTrail threshold in noise RMS for the discrimination of the trailing edge of the positive pulse.
        /// \param t_discSigmaNegLead threshold in noise RMS for the discrimination of the leading edge of the negative pulse.
        /// \param t_discSigmaNegPeak threshold in noise RMS for the discrimination of the peak of the negative pulse.
        /// \param t_discAbsNegPeak absolute threshold for the discrimination of the peak of a negative pulse
        /// \param t_discSigmaNegTrail threshold in noise RMS for the discrimination of the trailing edge of the negative pulse.
        void find2dHits(
                const TH2S &t_histo,
                const std::vector<std::array<double, 2>> &t_noiseParams,
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
                const double t_discSigmaNegTrail
        );

        /// \brief Private method used internally to find the 3D hits using the 2D hits found in both readout histos by
        /// the 2D hit finder.
        ///
        /// The method reads the 2D hits from the Event struct passed by reference and writes back two vectors containing
        /// all the matched ROI(pixel) hits for each pixel(ROI) hit. A match occurs if a pixel and a ROI pulse overlap.
        /// This method only performs the matching. BuildHitCandidates is called afterwards to build the 3D hit
        /// candidates containing 3D coordinates and reconstructed charge.
        /// \param t_event event struct.
        void find3dHits(Event &t_event);

        /// \brief Find hit candidates.
        ///
        /// This method builds Hit3ds using the matches found by Find3dHits. It reads the vectors mapping pixels to ROIs
        /// from the Event struct passed by reference, builds 3D hits and writes them back to the Event struct.
        /// \param t_event event struct.
        void buildHitCandidates(Event &t_event);

        /// \brief Vector containing all events.
        std::vector<Event> m_events;

        /// \brief ChargeData object containing the raw data.
        const ChargeData &m_chargeData;

        /// \brief Run parameters.
        const RunParams &m_runParams;
    };
}


#endif //PIXY_ROIMUX_CHARGEHITS_H
