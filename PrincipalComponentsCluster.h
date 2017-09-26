//
// Created by damian on 6/8/17.
//

#ifndef PIXY_ROIMUX_PRINCIPALCOMPONENTSCLUSTER_H
#define PIXY_ROIMUX_PRINCIPALCOMPONENTSCLUSTER_H


#include <string>
#include <functional>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include "TVector3.h"
#include "ChargeHits.h"
#include "Event.h"
#include "RunParams.h"


namespace pixy_roimux {
    class PrincipalComponentsCluster {
    public:
        explicit PrincipalComponentsCluster(const RunParams &t_runParams) : m_runParams(t_runParams) {};

        void analyseEvents(
                ChargeHits &t_chargeHits,
                const bool t_rejectOutliers = true,
                const bool t_rejectAmbiguities = true);


    private:
        int analysis3D(Event &t_event);

        double computeDoca(
                const Hit3d &t_hit,
                const TVector3 &t_avePosition,
                const TVector3 &t_axisDirVec);

        void rejectAmbiguities(Event &t_event);

        int rejectOutliers(
                Event &t_event,
                const double maxDocaAllowed);


        const RunParams &m_runParams;
    };
}


#endif //PIXY_ROIMUX_PRINCIPALCOMPONENTSCLUSTER_H
