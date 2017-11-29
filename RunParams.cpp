//
// Created by damian on 6/3/17.
//

#include "RunParams.h"


namespace pixy_roimux {
    RunParams::RunParams(const std::string t_runParamsFileName) {
        FILE *runParamsFile = fopen(t_runParamsFileName.c_str(), "r");
        char readBuffer[65536];
        rapidjson::FileReadStream jsonStream(runParamsFile, readBuffer, sizeof(readBuffer));
        m_jsonDoc.ParseStream(jsonStream);
        fclose(runParamsFile);
        if (m_jsonDoc.HasParseError() || !m_jsonDoc.IsObject()) {
            std::cerr << "ERROR: Failed to parse run parameter file " << t_runParamsFileName << '!' << std::endl;
            exit(1);
        }

        m_runId                     = getJsonMember("runId", rapidjson::kNumberType).GetUint();
        m_nPixels                   = getJsonMember("nPixels", rapidjson::kNumberType).GetUint();
        m_nRois                     = getJsonMember("nRois", rapidjson::kNumberType).GetUint();
        m_nChans                    = m_nPixels + m_nRois;
        m_pixelPitch                = getJsonMember("pixelPitch", rapidjson::kNumberType).GetDouble();
        m_driftLength               = getJsonMember("driftLength", rapidjson::kNumberType).GetDouble();
        m_sampleTime                = getJsonMember("sampleTime", rapidjson::kNumberType).GetDouble();
        m_driftSpeed                = getJsonMember("driftSpeed", rapidjson::kNumberType).GetDouble();
        m_anodeSample               = getJsonMember("anodeSample", rapidjson::kNumberType).GetUint();
        m_adcLsb                    = getJsonMember("adcLsb", rapidjson::kNumberType).GetDouble();
        m_preampGain                = getJsonMember("preampGain", rapidjson::kNumberType).GetDouble();
        m_preampTransimpedance      = getJsonMember("preampTransimpedance", rapidjson::kNumberType).GetDouble();
        m_nSamples                  = getJsonMember("nSamples", rapidjson::kNumberType).GetUint();
        m_noiseFilterSigma          = getJsonMember("noiseFilterSigma", rapidjson::kNumberType).GetDouble();
        m_discSigmaPixelLead        = getJsonMember("discSigmaPixelLead", rapidjson::kNumberType).GetDouble();
        m_discSigmaPixelPeak        = getJsonMember("discSigmaPixelPeak", rapidjson::kNumberType).GetDouble();
        m_discAbsPixelPeak          = getJsonMember("discAbsPixelPeak", rapidjson::kNumberType).GetDouble();
        m_discSigmaPixelTrail       = getJsonMember("discSigmaPixelTrail", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiPosLead       = getJsonMember("discSigmaRoiPosLead", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiPosPeak       = getJsonMember("discSigmaRoiPosPeak", rapidjson::kNumberType).GetDouble();
        m_discAbsRoiPosPeak         = getJsonMember("discAbsRoiPosPeak", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiPosTrail      = getJsonMember("discSigmaRoiPosTrail", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiNegLead       = getJsonMember("discSigmaRoiNegLead", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiNegPeak       = getJsonMember("discSigmaRoiNegPeak", rapidjson::kNumberType).GetDouble();
        m_discAbsRoiNegPeak         = getJsonMember("discAbsRoiNegPeak", rapidjson::kNumberType).GetDouble();
        m_discSigmaRoiNegTrail      = getJsonMember("discSigmaRoiNegTrail", rapidjson::kNumberType).GetDouble();
        m_discRange                 = getJsonMember("discRange", rapidjson::kNumberType).GetUint();
        m_detectBipolarRoiPulses    = getJsonMember("detectBipolarRoiPulses", rapidjson::kTrueType).GetBool();
        m_breakAtUsedSamples        = getJsonMember("breakAtUsedSamples", rapidjson::kTrueType).GetBool();
        m_pcaScaleFactor            = getJsonMember("pcaScaleFactor", rapidjson::kNumberType).GetDouble();
        m_pcaMaxIterations          = getJsonMember("pcaMaxIterations", rapidjson::kNumberType).GetUint();
        m_pcaUseChargeWeights       = getJsonMember("pcaUseChargeWeights", rapidjson::kNumberType).GetUint();
        m_kalmanRngSeed             = getJsonMember("kalmanRngSeed", rapidjson::kNumberType).GetUint();
        m_kalmanMaxIterations       = getJsonMember("kalmanMaxIterations", rapidjson::kNumberType).GetUint();
        m_kalmanUseRef              = getJsonMember("kalmanUseRef", rapidjson::kTrueType).GetBool();
        m_kalmanDeltaPval           = getJsonMember("kalmanDeltaPval", rapidjson::kNumberType).GetDouble();
        m_kalmanDeltaWeight         = getJsonMember("kalmanDeltaWeight", rapidjson::kNumberType).GetDouble();
        m_kalmanPdgCode             = getJsonMember("kalmanPdgCode", rapidjson::kNumberType).GetInt();
        m_kalmanMomMag              = getJsonMember("kalmanMomMag", rapidjson::kNumberType).GetDouble();

        m_tpcOrigin = std::vector<double>(3);
        auto jsonArrayItr = getJsonMember("tpcOrigin", rapidjson::kArrayType, 3, rapidjson::kNumberType).Begin();
        for (auto &&component : m_tpcOrigin) {
            component = jsonArrayItr->GetDouble();
            ++jsonArrayItr;
        }
        m_daq2readout = std::vector<unsigned>(m_nChans);
        jsonArrayItr = getJsonMember("daq2readout", rapidjson::kArrayType, m_nChans, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_daq2readout) {
            channel = jsonArrayItr->GetUint();
            ++jsonArrayItr;
        }
        m_readout2daq = std::vector<unsigned>(m_nChans);
        jsonArrayItr = getJsonMember("readout2daq", rapidjson::kArrayType, m_nChans, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_readout2daq) {
            channel = jsonArrayItr->GetUint();
            ++jsonArrayItr;
        }
        m_pixelCoor = std::vector<std::vector<int>>(m_nPixels, std::vector<int>(2));
        jsonArrayItr = getJsonMember("pixelCoorX", rapidjson::kArrayType, m_nPixels, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_pixelCoor) {
            channel.at(0) = jsonArrayItr->GetInt();
            ++jsonArrayItr;
        }
        jsonArrayItr = getJsonMember("pixelCoorY", rapidjson::kArrayType, m_nPixels, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_pixelCoor) {
            channel.at(1) = jsonArrayItr->GetInt();
            ++jsonArrayItr;
        }
        m_roiCoor = std::vector<std::vector<int>>(m_nRois, std::vector<int>(2));
        jsonArrayItr = getJsonMember("roiCoorX", rapidjson::kArrayType, m_nRois, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_roiCoor) {
            channel.at(0) = jsonArrayItr->GetInt();
            ++jsonArrayItr;
        }
        jsonArrayItr = getJsonMember("roiCoorY", rapidjson::kArrayType, m_nRois, rapidjson::kNumberType).Begin();
        for (auto &&channel : m_roiCoor) {
            channel.at(1) = jsonArrayItr->GetInt();
            ++jsonArrayItr;
        }
        m_kalmanPosErr = std::vector<double>(3);
        jsonArrayItr = getJsonMember("kalmanPosErr", rapidjson::kArrayType, 3, rapidjson::kNumberType).Begin();
        for (auto &&component : m_kalmanPosErr) {
            component = jsonArrayItr->GetDouble();
            ++jsonArrayItr;
        }
        m_kalmanMomErr = std::vector<double>(3);
        jsonArrayItr = getJsonMember("kalmanMomErr", rapidjson::kArrayType, 3, rapidjson::kNumberType).Begin();
        for (auto &&component : m_kalmanMomErr) {
            component = jsonArrayItr->GetDouble();
            ++jsonArrayItr;
        }
    }


    const rapidjson::Value & RunParams::getJsonMember(
            const std::string t_memberName,
            const rapidjson::Type t_memberType,
            const unsigned t_arraySize,
            const rapidjson::Type t_arrayType)
    {
        if (!m_jsonDoc.HasMember(t_memberName.c_str())) {
            std::cerr << "ERROR: Entry \"" << t_memberName << "\" in run parameter file not found!" << std::endl;
            exit(1);
        }
        rapidjson::Value &member = m_jsonDoc[t_memberName.c_str()];
        if ((t_memberType == rapidjson::kTrueType) || (t_memberType == rapidjson::kFalseType)) {
            if ((member.GetType() != rapidjson::kTrueType) && (member.GetType() != rapidjson::kFalseType)) {
                std::cerr << "ERROR: Entry \"" << t_memberName << "\" in run parameter file has wrong type!"
                          << std::endl;
                std::cerr << "Expected " << m_jsonTypes.at(rapidjson::kTrueType)
                          << " or " << m_jsonTypes.at(rapidjson::kFalseType)
                          << ", got " << m_jsonTypes.at(member.GetType()) << '.' << std::endl;
                exit(1);
            }
        }
        else if (member.GetType() != t_memberType) {
            std::cerr << "ERROR: Entry \"" << t_memberName << "\" in run parameter file has wrong type!"
                      << std::endl;
            std::cerr << "Expected " << m_jsonTypes.at(t_memberType) << ", got " << m_jsonTypes.at(member.GetType())
                      << '.' << std::endl;
            exit(1);
        }
        if (member.GetType() == rapidjson::kArrayType) {
            if (member.Size() != t_arraySize) {
                std::cerr << "ERROR: Size mismatch for array \"" << t_memberName << "\" in run parameter file!"
                          << std::endl;
                std::cerr << "Expected " << t_arraySize << ", got " << member.Size() << '.' << std::endl;
                exit(1);
            }
            for (const auto &value : member.GetArray()) {
                if (value.GetType() != t_arrayType) {
                    std::cerr << "ERROR: Type mismatch in array \"" << t_memberName << "\" in run parameter file!"
                              << std::endl;
                    std::cerr << "Expected " << m_jsonTypes.at(t_arrayType) << ", got "
                              << m_jsonTypes.at(value.GetType()) << '.' << std::endl;
                    exit(1);
                }
            }
        }
        return member;
    }
}