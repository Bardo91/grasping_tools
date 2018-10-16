//#include "InfPlanePrior.h"
//#include "../../Utils/utilsProb.h"
//
//using namespace std;
//using namespace arma;
//
//namespace gpis {
//    // constructor
//    InfPlanePrior::InfPlanePrior(const int numParts,
//                                 const int numDims,
//                                 const GPISParam gpisParam,
//                                 const arma::mat &covMat)
//    : GPISPrior(numParts, numDims, gpisParam, covMat),
//      planeParams_(numDimsPlusOne_)
//    {
//        planeParams_.zeros();
//    }
//
//    // public methods
//    void InfPlanePrior::computeAuxFPlusVec(const std::vector<Part> &parts)
//    {
//        // FPlus needs to be computed for each part in the whole scene (for this object).
//        // FPlus is a vector that packs together all parts. For each part, there is
//        // one value and numDims gradient components. So the length of FPlus is numParts * (1 + numDims).
//        // These values represent the "GP error" for each part with respect to this object.
//
//        // loop through all parts
//        int partIndex = -1;
//        for (std::vector<Part>::const_iterator iter = parts.begin();
//                iter != parts.end(); iter++)
//        {
//            partIndex++;
//
//            ///////////////////////////////////////////////////////////////////////////////////////
//            auxFPlusVec_(partIndex * numDimsPlusOne_) =
//                    - planeParams_(0) * iter->location(0)
//                    - planeParams_(1) * iter->location(1)
//                    - planeParams_(2) * iter->location(2)
//                    - planeParams_(3);
//
//            for (int d = 0; d < numDims_; ++d)
//            {
//                auxFPlusVec_(partIndex * numDimsPlusOne_ + d + 1) =
//                        iter->surfaceNormals(d) - planeParams_(d);
//            }
//        }
//    }
//
//    void InfPlanePrior::samplePriorParameters(const vector<Part> &parts,
//                                              const uvec &assocParts,
//                                              const int numAssocParts,
//                                              const mat &fullCovMat)
//    {
//        // update the objects. at this point the parts are associated to some objects.
//        // we look at this one object and find the parts associated with it.
//        // based on the parts we sample the parameters of that object.
//
//
//        mat rHSMat = zeros<mat>(numDimsPlusOne_ * numAssocParts, numDimsPlusOne_);
//        vec muTilde = zeros<vec>(numDimsPlusOne_ * numAssocParts);
//        for (int n = 0; n < numAssocParts; ++n)
//        {
//
//            rHSMat.row((n * numDimsPlusOne_)) = join_vert(parts[assocParts(n)].location(), vec{1}).t();
//            for (int d = 0; d < numDims_; ++d) //TODO: Could be vectorised
//            {
//                rHSMat.at((n * numDimsPlusOne_) + d + 1, d) = 1;
//                muTilde((n * numDimsPlusOne_) + d + 1) = parts[assocParts(n)].surfaceNormals(d);
//            }
//        }
//        mat P = solve(trimatl(matrixData_.lCholDec()), rHSMat);
//        mat Q = solve(trimatl(matrixData_.lCholDec()), muTilde);// TODO: recycle solution
//        vec B = trans(P) * Q;
//
//        mat planeParamPrecMat = trans(P) * P;
//        mat planeParamCovMat = inv(planeParamPrecMat);
//        vec planeParamMean = planeParamCovMat * B;
//
//        // resample the parameters of the object
//        vec nVar(numDimsPlusOne_);
////        multiVarGaussSample(nVar);
//        nVar.zeros();
//        planeParams_  = planeParamMean + chol(planeParamCovMat, "lower") * nVar;
//    }
//
//
//    void InfPlanePrior::printInfo()
//    {
//        cout << "This is an infinite plane GPISObject" << endl;
//        cout << "Plane parameters:" << planeParams_.t() << endl;
//    }
//
//}   // namespace gpis
