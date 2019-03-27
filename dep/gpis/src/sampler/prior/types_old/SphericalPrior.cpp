//#include "SphericalPrior.h"
//#include "../../Utils/utilsProb.h"
//
//using namespace std;
//using namespace arma;
//
//namespace gpis {
//    // constructor
//    SphericalPrior::SphericalPrior(const int numParts,
//                                   const int numDims,
//                                   const GPISParam gpisParam,
//                                   const arma::mat &covMat)
//    : GPISPrior(numParts, numDims, gpisParam, covMat),
//    location_(numDims_),
//    radius_(gpisParam.radius()),
//    twoTimesRadius_(2.0 * radius_),
//    radiusSquared_(pow(radius_, 2))
//    {
//        location_.zeros();
//    }
//
//    // public methods
//    void SphericalPrior::computeAuxFPlusVec(const std::vector<Part> &parts)
//    {
//        // FPlus needs to be computed for each part in the whole scene (for this object).
//        // FPlus is a vector that packs together all parts. For each part, there is
//        // one value and numDims gradient components. So the length of FPlus is numParts * (1 + numDims).
//        // These values represent the "GP error" for each part with respect to this object.
//        vec diffVec(numDims_);
//
//        // loop through all parts
//        int partIndex = -1;
//        for (std::vector<Part>::const_iterator iter = parts.begin();
//        iter != parts.end(); iter++)
//        {
//            partIndex++;
//            ///////////////////////////////////////////////////////////////////////////////////////
//            // for each part, we need to compute four numbers and write them into the Fplus vector
//
//            // this is a distance for each dimension
//            diffVec = iter->location() - location_;
//
//            // dot product
//            auxFPlusVec_(partIndex * numDimsPlusOne_) = 1.0 / twoTimesRadius_ *
//                (as_scalar(trans(diffVec)*diffVec) - radiusSquared_);
//
//            // now the gradient components
//            for (int d = 0; d < numDims_; ++d)
//            {
//                auxFPlusVec_(partIndex * numDimsPlusOne_ + d + 1) =
//                    iter->surfaceNormals(d) - 1 / radius_ * diffVec[d];
//            }
//        }
//    }
//
//
//    void SphericalPrior::samplePriorParameters(const vector<Part> &parts,
//                                               const uvec &assocParts,
//                                               const int numAssocParts,
//                                               const mat &fullCovMat)
//    {
//        // update the objects. at this point the parts are associated to some objects.
//        // we look at this one object and find the parts associated with it.
//        // based on the parts we sample the location of that object.
//
//        // compute the mean and covariance matrix for the location distribution for this object
//        mat eyeRepMat = zeros<mat>(numDims_ * numAssocParts, numDims_);
//        uvec derIndices(numDims_ * numAssocParts);
//        vec muTilde(numDims_ * numAssocParts);
//
//        for (int n = 0; n < numAssocParts; ++n)
//        {
//            for (int d = 0; d < numDims_; ++d)
//            {
//                derIndices.at((n * numDims_) + d) = assocParts(n) * (numDimsPlusOne_)+ d + 1;
//
//                eyeRepMat.at((n * numDims_) + d, d) = 1;
//                muTilde((n * numDims_) + d) = parts[assocParts(n)].location()(d) -
//                    radius_ * parts[assocParts(n)].surfaceNormals(d);
//
//            }
//        }
//        mat lChol = chol(fullCovMat.submat(derIndices, derIndices), "lower");
//
//        mat P = solve(trimatl(lChol), eyeRepMat);
//        mat Q = solve(trimatl(lChol), muTilde);
//        vec B = trans(P) * Q / radiusSquared_;
//
//        mat locPrecMat = trans(P) * P / radiusSquared_;
//        mat locCovMat = inv(locPrecMat);
//        vec locMean = locCovMat * B;
//
//        // resample the location of the object
//        vec nVar(numDims_);
//        multiVarGaussSample(nVar);
//        location_ = locMean + chol(locCovMat, "lower") * nVar;
//
//    }
//
//
//    void SphericalPrior::printInfo()
//    {
//        cout << "This is a spherical GPISObject" << endl;
//        cout << "Location:" << location_.t() << endl;
//    }
//
//
//
//
//}   // namespace gpis
