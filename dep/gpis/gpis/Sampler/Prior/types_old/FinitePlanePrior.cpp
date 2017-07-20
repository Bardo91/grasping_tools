//#include "FinitePlanePrior.h"
//#include "../../Utils/utilsProb.h"
//
//using namespace std;
//using namespace arma;
//
//namespace gpis {
//    // constructor
//    FinitePlanePrior::FinitePlanePrior(const int numParts,
//                                       const int numDims,
//                                       const GPISParam gpisParam,
//                                       const arma::mat &covMat)
//        : GPISPrior(numParts, numDims, gpisParam, covMat),
//        planeParams_(numDimsPlusOne_),
//        locInPlane_(2),
//        stDevInPlane_(gpisParam.radius()),
//        planeTrafoMat_(2,3),
//        location_(3)
//    {
//        planeParams_.zeros();
//        locInPlane_.zeros();
//        planeTrafoMat_.zeros();
//        location_.zeros();
//    }
//
//    // public methods
//    void FinitePlanePrior::computeAuxFPlusVec(const std::vector<Part> &parts)
//    {
//        // FPlus needs to be computed for each part in the whole scene (for this object).
//        // FPlus is a vector that packs together all parts. For each part, there is
//        // one value and numDims gradient components. So the length of FPlus is numParts * (1 + numDims).
//        // These values represent the "GP error" for each part with respect to this object.
//
//        // loop through all parts
//        int partIndex = -1;
//        for (std::vector<Part>::const_iterator iter = parts.begin();
//             iter != parts.end(); iter++)
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
//    void FinitePlanePrior::samplePriorParameters(const vector<Part> &parts,
//                                                 const uvec &assocParts,
//                                                 const int numAssocParts,
//                                                 const mat &fullCovMat)
//    {
//        // based on the parts associated to this object we sample the parameters of that prior.
//
//        // Sample the inifinite plane parameters
//        mat rHSMat = zeros<mat>(numDimsPlusOne_ * numAssocParts, numDimsPlusOne_);
//        vec muTilde = zeros<vec>(numDimsPlusOne_ * numAssocParts);
//        mat locationMat = zeros<mat>(numDims_, numAssocParts);
//        for (int n = 0; n < numAssocParts; ++n)
//        {
//            vec loc = parts[assocParts(n)].location();
//            locationMat.col(n) = loc;
//            rHSMat.row((n * numDimsPlusOne_)) = join_vert(loc, vec{1}).t();
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
//        multiVarGaussSample(nVar);
//
//        planeParams_  = planeParamMean + chol(planeParamCovMat, "lower") * nVar;
//
//        // sample the in-plane centroid
//        // Compute in-plane x-,y-coordinates
//        double nXYMag = norm(planeParams_(span(0,1)));
//        double nXYZMag = norm(planeParams_(span(0,2)));
//        planeTrafoMat_ = {{-planeParams_(1)/nXYMag, planeParams_(0)/nXYMag, 0},
//                       {-planeParams_(0)*planeParams_(2)/(nXYMag * nXYZMag),
//                        -planeParams_(1)*planeParams_(2)/(nXYMag * nXYZMag),
//                        nXYMag/nXYZMag}};
//        mat inPlaneXY = planeTrafoMat_ * locationMat;
//        // Compute in-plane mean and variance
//        double prec = 1/(stDevInPlane_ * stDevInPlane_);
//        vec inPlaneXMuAndPrec = productOf1DGaussians(inPlaneXY.row(0), prec);
//        vec inPlaneYMuAndPrec = productOf1DGaussians(inPlaneXY.row(1), prec);
//        // Sample in-plane centroid
//        vec nVarXY(2);
//        multiVarGaussSample(nVarXY);
//        locInPlane_(0) = inPlaneXMuAndPrec(0) + 1/sqrt(inPlaneXMuAndPrec(1)) * nVarXY(0);
//        locInPlane_(1) = inPlaneYMuAndPrec(0) + 1/sqrt(inPlaneYMuAndPrec(1)) * nVarXY(1);
//        computeLocation();
//    }
//
//    double FinitePlanePrior::computeGeometricAssocProb(const Part &part)
//    {
//        vec inPlaneXY = planeTrafoMat_ * part.location();
//        mat covMat(2,2,fill::eye);
//        covMat *= stDevInPlane_ * stDevInPlane_;
//        return computeGaussianPdf(inPlaneXY, locInPlane_, covMat);
//    }
//
//    void FinitePlanePrior::computeLocation()
//    {
//        // transform location in plane back into global coordinate system
//        location_ = planeTrafoMat_.t() * locInPlane_ - planeParams_(3) * planeParams_.head(3);
//    }
//
//    void FinitePlanePrior::printInfo()
//    {
//        cout << "This is a finite plane GPISObject" << endl;
//        cout << "Plane centroid:" << location() << endl;
//    }
//
//}   //namespace gpis
