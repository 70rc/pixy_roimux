//
// Created by damian on 6/3/17.
//

#ifndef PIXY_ROIMUX_RUNPARAMS_H
#define PIXY_ROIMUX_RUNPARAMS_H


#include <array>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/rapidjson.h"


namespace pixy_roimux {
    /// \brief Class that reads all the run parameters from a JSON file and makes them available to all the other reco
    /// classes.
    ///
    /// Namely, the map from DAQ channels to readout channels and vice versa, and the mechanical coordinates of the pixels
    /// and the regions of interest (ROI). All the getter methods are const and thus an instance of this class can be passed
    /// by const reference to other objects requiring it. DAQ channels are numbered from 0 to 63 with channels 0 through 31
    /// corresponding to the "Ind_x" histogram in the raw data and channels 32 through 63 corresponding to the "Col_x"
    /// histogram in the raw data. Caution has to be applied not to confuse those with the actual pixel collection channels
    /// and ROI induction channels as the naming of the histograms is hardcoded in the DAQ driver. Additionally, this class
    /// also stores all the run parameters necessary for reconstruction such as calibration constants.
    class RunParams {
    public:

        /// \brief Constructor.
        /// \param t_runParamsFileName JSON file name.
        explicit RunParams(const std::string t_runParamsFileName);

        /// \brief Convert DAQ channel to pixel channel.
        /// \param t_daqChan DAQ channel.
        /// \return index of the pixel.
        unsigned daq2pixel(const unsigned t_daqChan) const {
            return m_daq2readout.at(t_daqChan);
        }

        /// \brief Convert pixel channel to DAQ channel.
        /// \param t_pixelInd index of the pixel.
        /// \return DAQ channel.
        unsigned pixel2daq(const unsigned t_pixelInd) const {
            return m_readout2daq.at(t_pixelInd);
        }

        /// \brief Convert DAQ channel to ROI channel.
        /// \param t_daqChan DAQ channel.
        /// \return index of the ROI.
        unsigned daq2roi(const unsigned t_daqChan) const {
            return m_daq2readout.at(t_daqChan) - m_nPixels;
        }

        /// \brief Convert ROI channel to DAQ channel.
        /// \param t_roiInd index of the ROI.
        /// \return DAQ channel.
        unsigned roi2daq(const unsigned t_roiInd) const {
            return m_readout2daq.at(t_roiInd + m_nPixels);
        }

        /// \brief Get pixel coordinates in units of pixel pitch.
        ///
        /// These are just relative offsets within one ROI. To get absolute coordinates, this number will be added to
        /// the according ROI coordinates obtained from GetPixelCoor.
        /// \param t_pixelInd index of the pixel
        /// \param t_dim 0 for X, 1 for Y.
        /// \return pixel coordinate.
        int getPixelCoor(
                const unsigned t_pixelInd,
                const unsigned t_dim)
        const {
            return m_pixelCoor.at(t_pixelInd).at(t_dim);
        }

        /// \brief Get ROI coordinates in units of pixel pitch.
        ///
        /// To get absolut pixel coordinates, this number will be added to the pixel coordinate obtained from
        /// GetPixelCoor.
        /// \param t_roiInd index of the ROI
        /// \param t_dim 0 for X, 1 for Y.
        /// \return ROI coordinate.
        int getRoiCoor(
                const unsigned t_roiInd,
                const unsigned t_dim)
        const {
            return m_roiCoor.at(t_roiInd).at(t_dim);
        }

        /// \brief Get the run ID.
        /// \return run ID.
        unsigned getRunId() const {
            return m_runId;
        }

        /// \brief Get the number of pixels.
        /// \return number of pixels.
        unsigned getNPixels() const {
            return m_nPixels;
        }

        /// \brief Get the number of ROIs.
        /// \return number of ROIs.
        unsigned getNRois() const {
            return m_nRois;
        }

        /// \brief Get the total number of readout channels.
        ///
        /// Divided by 2 this gives the number of channels of the two DAQ histograms ("Ind_x" and "Col_x").
        /// \return number of channels.
        unsigned getNChans() const {
            return m_nChans;
        }

        /// \brief Get absolute coordinates of the TPC origin.
        ///
        /// These are used to shift the origin of the TPC. They are added to the pixel and drift coordinates.
        /// \return TPC origin.
        std::vector<double> getTpcOrigin() const {
            return m_tpcOrigin;
        };

        /// \brief Get the pixel pitch in cm.
        /// \return pixel pitch.
        double getPixelPitch() const {
            return m_pixelPitch;
        }

        /// \brief Get the drift length in cm.
        /// \return drift length.
        double getDriftLength() const {
            return m_driftLength;
        }

        /// \brief Get the sample time in us.
        /// \return sample time.
        double getSampleTime() const {
            return m_sampleTime;
        }

        /// \brief Get the drift speed in cm/us.
        /// \return drift speed.
        double getDriftSpeed() const {
            return m_driftSpeed;
        }

        /// \brief Get the location of the anode in histogram samples.
        /// \return anode sample.
        unsigned getAnodeSample() const {
            return m_anodeSample;
        }

        /// \brief Get the ADC LSB in mV.
        /// \return ADC LSB.
        double getAdcLsb() const {
            return m_adcLsb;
        }

        /// \brief Get the preamplifier gain in mV/fC.
        /// \return preamplifier gain.
        double getPreampGain() const {
            return m_preampGain;
        }

        /// \brief Get the number of samples to process.
        /// \return number of samples.
        unsigned getNSamples() const {
            return m_nSamples;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the noise filter.
        /// \return threshold.
        double getNoiseFilterSigma() const {
            return m_noiseFilterSigma;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the leading edge of a pixel pulse.
        /// \return threshold.
        double getDiscSigmaPixelLead() const {
            return m_discSigmaPixelLead;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the peak of a pixel pulse.
        /// \return threshold.
        double getDiscSigmaPixelPeak() const {
            return m_discSigmaPixelPeak;
        }

        /// \brief Get the absolute threshold for the discrimination of the peak of a pixel pulse.
        /// \return threshold.
        double getDiscAbsPixelPeak() const {
            return m_discAbsPixelPeak;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a pixel pulse.
        /// \return threshold.
        double getDiscSigmaPixelTrail() const {
            return m_discSigmaPixelTrail;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the leading edge of a positive ROI
        /// pulse.
        /// \return threshold.
        double getDiscSigmaRoiPosLead() const {
            return m_discSigmaRoiPosLead;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the peak of a positive ROI pulse.
        /// \return threshold.
        double getDiscSigmaRoiPosPeak() const {
            return m_discSigmaRoiPosPeak;
        }

        /// \brief Get the absolute threshold for the discrimination of the peak of a positive ROI pulse.
        /// \return threshold.
        double getDiscAbsRoiPosPeak() const {
            return m_discAbsRoiPosPeak;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a positive ROI
        /// pulse.
        /// \return threshold.
        double getDiscSigmaRoiPosTrail() const {
            return m_discSigmaRoiPosTrail;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the leading edge of a negative ROI
        /// pulse.
        /// \return threshold.
        double getDiscSigmaRoiNegLead() const {
            return m_discSigmaRoiNegLead;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the peak of a negative ROI pulse.
        /// \return threshold.
        double getDiscSigmaRoiNegPeak() const {
            return m_discSigmaRoiNegPeak;
        }

        /// \brief Get the absolute threshold for the discrimination of the peak of a negative ROI pulse.
        /// \return threshold.
        double getDiscAbsRoiNegPeak() const {
            return m_discAbsRoiNegPeak;
        }

        /// \brief Get the threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a negative ROI
        /// pulse.
        /// \return threshold.
        double getDiscSigmaRoiNegTrail() const {
            return m_discSigmaRoiNegTrail;
        }

        /// \brief Get the range in samples within which a discriminated egde/peak needs to be found.
        /// \return discriminator range.
        unsigned getDiscRange() const {
            return m_discRange;
        }

        /// \brief Get the scale factor for the principal components analysis.
        /// \return scale factor.
        double getPcaScaleFactor() const {
            return m_pcaScaleFactor;
        }

        /// \brief Get the maximum number of iterations for the principal components analysis.
        /// \return maximum number of iterations.
        unsigned getPcaMaxIterations() const {
            return m_pcaMaxIterations;
        }

        /// \brief Get the position error for the Kalman fitter.
        /// \return position error.
        std::vector<double> getKalmanPosErr() const {
            return m_kalmanPosErr;
        }

        /// \brief Get the start momentum magnitude for the Kalman fitter.
        /// \return start momentum magnitude.
        double getKalmanMomMag() const {
            return m_kalmanMomMag;
        }

        /// \brief Get the momentum error for the Kalman fitter.
        /// \return momentum error.
        std::vector<double> getKalmanMomErr() const {
            return m_kalmanMomErr;
        }

        /// \brief Get the seed for the RNG used by the Kalman fitter.
        /// \return RNG seed.
        unsigned getKalmanRngSeed() const {
            return m_kalmanRngSeed;
        }

        /// \brief Get the maximum number of iterations for the Kalman fitter.
        /// \return maximum number of iterations.
        unsigned getKalmanMaxIterations() const {
            return m_kalmanMaxIterations;
        }

        /// \brief Get useRefKalman for the Kalman fitter.
        /// \return useRefKalman.
        bool getKalmanUseRef() const {
            return m_kalmanUseRef;
        }

        /// \brief Get deltaPval for the Kalman fitter.
        /// \return deltaPval.
        double getKalmanDeltaPval() const {
            return m_kalmanDeltaPval;
        }

        /// \brief Get deltaWeight for the Kalman fitter.
        /// \return deltaWeight.
        double getKalmanDeltaWeight() const {
            return m_kalmanDeltaWeight;
        }

        /// \brief Get the PDG code of the Kalman fitter particle hypothesis.
        /// \return PDG code.
        int getKalmanPdgCode() const {
            return m_kalmanPdgCode;
        }


    private:

        /// \brief Get a member from the JSON file and perform some basic checks.
        /// \param t_memberName name of the member.
        /// \param t_memberType type of the member. Will be checked against.
        /// \param t_arraySize size in case t_memberType is kArrayType.
        /// \param t_arrayType type of the array elements in case t_memberType is kArrayType.
        /// \return the value read from the JSON file.
        const rapidjson::Value & getJsonMember(
                const std::string t_memberName,
                const rapidjson::Type t_memberType,
                const unsigned t_arraySize = 0,
                const rapidjson::Type t_arrayType = rapidjson::kNullType);

        /// \brief Strings associated with the rapisjson::Type enumerator.
        const std::array<std::string, 7> m_jsonTypes = {{"Null", "False", "True", "Object", "Array", "String", "Number"}};

        /// \brief rapidjson Document storing all the data read from the JSON file.
        rapidjson::Document m_jsonDoc;

        /// \brief Run ID.
        unsigned m_runId;

        /// \brief Number of pixels.
        unsigned m_nPixels;

        /// \brief Number of ROIs.
        unsigned m_nRois;

        /// \brief Total number of readout channels.
        unsigned m_nChans;

        /// \brief Absolute coordinates of the TPC origin.
        std::vector<double> m_tpcOrigin;

        /// \brief Pixel pitch in cm.
        double m_pixelPitch;

        /// \brief Drift length in cm.
        double m_driftLength;

        /// \brief Sample time in us.
        double m_sampleTime;

        /// \brief Drift speed in cm/us.
        double m_driftSpeed;

        /// \brief Location of the anode in histogram samples.
        unsigned m_anodeSample;

        /// \brief Analog-to-digital converter least significant bit in mV.
        ///
        /// Used to convert the value recorded by the ADC to a voltage.
        double m_adcLsb;

        /// \brief Preamp gain in mV/fC.
        ///
        /// Used to convert the voltage recorded by the ADC to charge.
        double m_preampGain;

        /// \brief Array mapping DAQ channels to readout channels.
        std::vector<unsigned> m_daq2readout;

        /// \brief Array mapping readout channels to DAQ channels.
        std::vector<unsigned> m_readout2daq;

        /// \brief Array containing the 2D pixel coordinates.
        std::vector<std::vector<int>> m_pixelCoor;

        /// \brief Array containing the 2D ROI coordinates.
        std::vector<std::vector<int>> m_roiCoor;

        /// \brief Number of samples to process.
        unsigned m_nSamples;

        /// \brief Threshold in sigma of noise Gaussian for the noise filter.
        double m_noiseFilterSigma;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the leading edge of a pixel pulse.
        double m_discSigmaPixelLead;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the peak of a pixel pulse.
        double m_discSigmaPixelPeak;

        /// \brief Absolute threshold for the discrimination of the peak of a pixel pulse.
        double m_discAbsPixelPeak;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a pixel pulse.
        double m_discSigmaPixelTrail;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the leading edge of a positive ROI pulse.
        double m_discSigmaRoiPosLead;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the peak of a positive ROI pulse.
        double m_discSigmaRoiPosPeak;

        /// \brief Absolute threshold for the discrimination of the peak of a positive ROI pulse.
        double m_discAbsRoiPosPeak;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a positive ROI pulse.
        double m_discSigmaRoiPosTrail;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the leading edge of a negative ROI pulse.
        double m_discSigmaRoiNegLead;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the peak of a negative ROI pulse.
        double m_discSigmaRoiNegPeak;

        /// \brief Absolute threshold for the discrimination of the peak of a negative ROI pulse.
        double m_discAbsRoiNegPeak;

        /// \brief Threshold in sigma of noise Gaussian for the discrimination of the trailing edge of a negative ROI pulse.
        double m_discSigmaRoiNegTrail;

        /// \brief Range in samples within which a discriminated egde/peak needs to be found.
        unsigned m_discRange;

        /// \brief Scale factor for the principal components analysis.
        double m_pcaScaleFactor;

        /// \brief Maximum number of iterations for the principal components analysis.
        unsigned m_pcaMaxIterations;

        /// \brief Position error for the Kalman fitter.
        std::vector<double> m_kalmanPosErr;

        /// \brief Start momentum magnitude for the Kalman fitter.
        double m_kalmanMomMag;

        /// \brief Momentum error for the Kalman fitter.
        std::vector<double> m_kalmanMomErr;

        /// \brief Seed for the RNG used by the Kalman fitter.
        unsigned m_kalmanRngSeed;

        /// \brief Maximum number of iterations for the Kalman fitter.
        unsigned m_kalmanMaxIterations;

        /// \brief Use a reference track for the Kalman fitter.
        bool m_kalmanUseRef;

        /// \brief deltaPval for the Kalman fitter.
        double m_kalmanDeltaPval;

        /// \brief deltaWeight for the Kalman fitter.
        double m_kalmanDeltaWeight;

        /// \brief PDG code of the Kalman fitter particle hypothesis.
        int m_kalmanPdgCode;
    };
}


#endif //PIXY_ROIMUX_RUNPARAMS_H
