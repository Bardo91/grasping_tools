#ifndef SPHERICALPRIOR_H
#define SPHERICALPRIOR_H

//#include "GPISPrior.h"
//
//namespace gpis {
//    class SphericalPrior: public GPISPrior
//    {
//    public:
//        // constructor
//        SphericalPrior(const int numParts,
//                       const int numDims,
//                       const GPISParam gpisParam,
//                       const arma::mat &covMat);
//        // public methods
//        const arma::vec& location() const { return location_; }
//        void computeAuxFPlusVec(const std::vector<Part> &parts);
//        void samplePriorParameters(const std::vector<Part> &parts,
//                                           const arma::uvec &assocParts,
//                                           const int numAssocParts,
//                                           const arma::mat &fullCovMat);
//        void printInfo();
//    private:
//        arma::vec location_;
//        double radius_;
//        double twoTimesRadius_;
//        double radiusSquared_;
//    };
//}   // namespace gpis

#endif // SPHERICALPRIOR_H
