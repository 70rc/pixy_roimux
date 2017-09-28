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
    /// \brief Class implementing a principal component analysis (PCA).
    ///
    /// A 3D PCA calculates the three orthogonal eigen vectors of a 3D spacepoint cloud and their eigen values.
    /// The longest eigen vector is used as an estimate for the direction of a track or shower. After a first iteration,
    /// ambiguities in the 3D spacepoints are resolved by selecting the ones closest to this estimate. In further
    /// iterations, outliers are rejected that are outside a cylinder around the track estimate. The radius of the
    /// cylinder is defined by the second largest eigen value, the average distance of closest approach (doca) to the
    /// track estimate and a scale factor defined in the RunParams. With the reduced hit sample, the outlier
    /// rejection procedure is repeated. The process stops when no more outliers are found or the maximum number of
    /// iterations specified in the RunParams is reached.
    class PrincipalComponentsCluster {
    public:
        /// \brief Constructor
        /// \param t_runParams run parameters.
        explicit PrincipalComponentsCluster(const RunParams &t_runParams) : m_runParams(t_runParams) {};

        /// \brief Run iterative PCA on a collection of events.
        /// \param t_chargeHits ChargeHits object containing the 3D spacepoints.
        /// \param t_rejectOutliers perform iterative outlier rejection.
        /// \param t_rejectAmbiguities perform ambiguity rejection.
        void analyseEvents(
                ChargeHits &t_chargeHits,
                const bool t_rejectOutliers = true,
                const bool t_rejectAmbiguities = true);


    private:
        /// \brief Perform a single PCA on a single event.
        /// \param t_event Event struct.
        /// \return error code indicating success or failure of PCA.
        int analysis3D(Event &t_event);

        /// \brief Compute the distance of closest approach (doca) of a 3D hit to a vector.
        /// \param t_hit 3D hit.
        /// \param t_avePosition position of the vector.
        /// \param t_axisDirVec direction of the vector.
        /// \return doca.
        double computeDoca(
                const Hit3d &t_hit,
                const TVector3 &t_avePosition,
                const TVector3 &t_axisDirVec);

        /// \brief Reject ambiguities of a single event.
        /// \param t_event Event struct.
        void rejectAmbiguities(Event &t_event);

        /// \brief Iteratively reject outliers of a single event.
        /// \param t_event Event struct.
        /// \param maxDocaAllowed hits further away from the track estimate are rejected.
        /// \return number of rejected hits.
        int rejectOutliers(
                Event &t_event,
                const double maxDocaAllowed);


        /// \brief run parameters.
        const RunParams &m_runParams;
    };
}


#endif //PIXY_ROIMUX_PRINCIPALCOMPONENTSCLUSTER_H
