#ifndef FINITEPLANEPRIOR_H
#define FINITEPLANEPRIOR_H
//
//#include "GPISPrior.h"
//
//namespace gpis {
//    class FinitePlanePrior: public GPISPrior
//    {
//    public:
//        // constructor
//        FinitePlanePrior(const int numParts,
//                         const int numDims,
//                         const GPISParam gpisParam,
//                         const arma::mat &covMat);
//        // public methods
//        const arma::vec& planeParams() const {return planeParams_;}
//        const arma::vec& location() const {return location_;}
//        void computeAuxFPlusVec(const std::vector<Part> &parts);
//        void samplePriorParameters(const std::vector<Part> &parts,
//                                           const arma::uvec &assocParts,
//                                           const int numAssocParts,
//                                           const arma::mat &fullCovMat);
//        double computeGeometricAssocProb(const Part &part);
//        void printInfo();
//    private:
//        arma::vec planeParams_;
//        arma::vec locInPlane_;
//        double stDevInPlane_;
//        arma::mat planeTrafoMat_;
//        arma::vec location_;
//        void computeLocation();
//
//    };
//}   // namespace gpis

#endif // INFPLANEPRIOR_H
