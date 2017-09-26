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
    ///
    /// Loud noises.
    ///
    class NoiseFilter {
    public:

        ///
        /// Constructor
        ///
        explicit NoiseFilter(const RunParams &t_runParams) : m_runParams(t_runParams) {};

        ///
        /// Compute mean and standard deviation of the noise by fitting a Gaussian to the amplitude distribution of the
        /// 2D histogram of a single channel.
        ///
        static std::pair<double, double> computeNoiseParams(std::shared_ptr<TH1D> t_channelHisto,
                                                            const bool t_fit);

        ///
        /// Filter data.
        ///
        void filterData(ChargeData &t_data);


    private:

        ///
        /// Filter a single TH2S from a dataset.
        ///
        void filterHisto(TH2S &t_histo);

        ///
        /// Run parameters.
        ///
        const RunParams &m_runParams;
    };
}


#endif //PIXY_ROIMUX_NOISEFILTER_H
