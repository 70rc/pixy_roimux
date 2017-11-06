//
// Created by damian on 6/6/17.
//

#include "NoiseFilter.h"


namespace pixy_roimux{
    std::pair<double, double> NoiseFilter::computeNoiseParams(std::shared_ptr<TH1D> t_channelHisto,
                                                              const bool t_fit) {
        int histoMin = static_cast<int>(t_channelHisto->GetBinContent(t_channelHisto->GetMinimumBin())) - 1;
        int histoMax = static_cast<int>(t_channelHisto->GetBinContent(t_channelHisto->GetMaximumBin())) + 1;
        unsigned nBins = static_cast<unsigned>(histoMax - histoMin);
        TH1S amplitudeSpectrum("amplitudeSpectrum", "amplitudeSpectrum", nBins, histoMin, histoMax);
        for (unsigned sample = 0; sample < t_channelHisto->GetNbinsX(); ++sample) {
            amplitudeSpectrum.Fill(t_channelHisto->GetBinContent(sample + 1));
        }
        if (t_fit) {
            TF1 gauss("gauss", "gaus");
            amplitudeSpectrum.GetXaxis()->SetRangeUser(amplitudeSpectrum.GetMean() - 1 * amplitudeSpectrum.GetStdDev(),
                                                       amplitudeSpectrum.GetMean() + 1 * amplitudeSpectrum.GetStdDev());
            amplitudeSpectrum.Fit(&gauss, "QN");
            //std::cout << "mean stats: " << amplitudeSpectrum.GetMean() << "\tMean gauss: " << gauss.GetParameter(1)
            //          << std::endl;
            //std::cout << "stddev stats: " << amplitudeSpectrum.GetStdDev() << "\tstddev gauss: "
            //          << gauss.GetParameter(2) << std::endl;
            return {gauss.GetParameter(1), gauss.GetParameter(2)};
        }
        else {
            return {amplitudeSpectrum.GetMean(), amplitudeSpectrum.GetStdDev()};
        }
    }


    void NoiseFilter::filterHisto(
            TH2S &t_histo,
            std::vector<std::pair<double, double>> &t_noiseParams) {
        unsigned nSamples = static_cast<unsigned>(t_histo.GetNbinsX());
        unsigned nChannels = static_cast<unsigned>(t_histo.GetNbinsY());
        std::vector<std::pair<double, double>> thresholds(nChannels);
        for (unsigned channel = 0; channel < nChannels; ++channel) {
            auto channelHisto = std::shared_ptr<TH1D>(t_histo.ProjectionX("channelHisto", (channel + 1), (channel + 1)));
            //std::cout << "channel: " << channel << std::endl;
            std::pair<double, double> noiseParams = computeNoiseParams(channelHisto, true);
            thresholds.at(channel).first = noiseParams.first - m_runParams.getNoiseFilterSigma() * noiseParams.second;
            thresholds.at(channel).second = noiseParams.first + m_runParams.getNoiseFilterSigma() * noiseParams.second;
        }
        for (unsigned sample = 0; sample < nSamples; ++sample) {
            double commonModeNoise = 0.;
            unsigned nCommonModeChannels = 0;
            for (unsigned channel = 0; channel < nChannels; ++channel) {
                int binContent = static_cast<int>(t_histo.GetBinContent((sample + 1), (channel + 1)));
                if ((binContent >= thresholds.at(channel).first) && binContent <= thresholds.at(channel).second) {
                    commonModeNoise += binContent;
                    ++nCommonModeChannels;
                }
            }
            commonModeNoise /= static_cast<double>(nCommonModeChannels);
            commonModeNoise = round(commonModeNoise);
            for (unsigned channel = 0; channel < nChannels; ++channel) {
                t_histo.SetBinContent((sample + 1), (channel + 1),
                                    (t_histo.GetBinContent((sample + 1), (channel + 1)) - commonModeNoise));
            }
        }
        t_noiseParams.resize(nChannels);
        unsigned channel = 0;
        for (auto &&noiseParams : t_noiseParams) {
            auto channelHisto = std::shared_ptr<TH1D>(t_histo.ProjectionX("channelHisto", (channel + 1), (channel + 1)));
            noiseParams = computeNoiseParams(channelHisto, true);
            ++channel;
        }
    }


    void NoiseFilter::filterData(ChargeData &t_data) {
        t_data.getNoiseParams().clear();
        t_data.getNoiseParams().resize(t_data.getReadoutHistos().size());
        auto eventId = t_data.getEventIds().cbegin();
        auto noiseParams = t_data.getNoiseParams().begin();
        for (auto &&histos : t_data.getReadoutHistos()) {
            std::cout << "Filtering event number " << *eventId << "...\n";
            filterHisto(histos.first, noiseParams->first);
            filterHisto(histos.second, noiseParams->second);
            ++eventId;
            ++noiseParams;
        }
    }
}
