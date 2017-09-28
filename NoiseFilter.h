//
// Created by damian on 6/6/17.
//

#ifndef PIXY_ROIMUX_NOISEFILTER_H
#define PIXY_ROIMUX_NOISEFILTER_H


#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "TF1.h"
#include "TH1D.h"
#include "TH2S.h"
#include "ChargeData.h"
#include "RunParams.h"


namespace pixy_roimux {
    /// \brief Class implementing a filter to reduce common mode noise in the raw data.
    ///
    /// Filters the histograms contained in a ChargeData object. First, an amplitude spectrum is created for each
    /// channel by filling the data into a new ROOT TH1S histogram. Then, it is fitted with a Gaussian. For a reliable
    /// fit, the range is constrained to the mean +- 1 RMS of the spectrum. After this, for each sample, an average is
    /// formed over all channels excluding values above a threshold in standard deviations of the Gaussian specified in the
    /// RunParams. Finally, the average is subtracted from all channels. The method to generate the mean and standard
    /// deviation of the noise from a 1D channel histogram is static so it can also be used without prior instantiation.
    class NoiseFilter {
    public:

        /// \brief Constructor
        /// \param t_runParams run parameters.
        explicit NoiseFilter(const RunParams &t_runParams) : m_runParams(t_runParams) {};

        /// \brief Compute mean and standard deviation of the noise by fitting a Gaussian to the amplitude spectrum of
        /// the 1D histogram of a single channel.
        ///
        /// Static so it can be used without prior instantiation.
        /// \param t_channelHisto TH1D of a single channel.
        /// \param t_fit fit Gaussian to the amplitude spectrum. If false, mean and RMS of the TH1D are used instead.
        /// \return pair containing mean and standard deviation of the amplitude spectrum.
        static std::pair<double, double> computeNoiseParams(std::shared_ptr<TH1D> t_channelHisto,
                                                            const bool t_fit);

        /// \brief Filter data.
        /// \param t_data ChargeData object containing the histograms to be filtered.
        void filterData(ChargeData &t_data);


    private:

        /// \brief Filter a single TH2S from a dataset.
        /// \param t_histo 2D histogram to be filtered.
        void filterHisto(TH2S &t_histo);

        /// \brief Run parameters.
        const RunParams &m_runParams;
    };
}


#endif //PIXY_ROIMUX_NOISEFILTER_H
