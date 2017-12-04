///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "ContactPoint.h"

#include <armadillo>
#include <cassert>

using namespace arma;

namespace grasping_tools {
	//-----------------------------------------------------------------------------------------------------------------
	ContactPoint::ContactPoint(const arma::colvec3 &_position, const arma::mat33 &_covPos, const arma::colvec3 &_n, const arma::mat33 &_covNormal, const eContactTypes _contactType, double _maxForce, double _friction, double _spinFriction) {
		mPosition = _position;
		mPositionCovariance = _covPos;
		
		assert((norm(_n) - 1.0) < 1e-10);
		mNormal = _n/norm(_n);

		mNormalCovariance = _covNormal;

		mContactType = _contactType;
		
		mMu = _friction;
		mMus = _spinFriction;
		
		mMaxForce = _maxForce;
		computeContactFrame();
	}

	//-----------------------------------------------------------------------------------------------------------------
	colvec3 ContactPoint::position() const {
		return mPosition;
	}

	arma::mat33 ContactPoint::posCov() const{
		return mPositionCovariance;
	}

	//-----------------------------------------------------------------------------------------------------------------
	colvec3 ContactPoint::normal() const {
		return mNormal;
	}

	arma::mat33 ContactPoint::norCov() const {
		return mNormalCovariance;
	}

	//-----------------------------------------------------------------------------------------------------------------
	mat33 ContactPoint::frame() {
		if (!mInitContactFrame) {	// Lazy initialization
			computeContactFrame();
		}

		return mContactFrame;
	}

	//-----------------------------------------------------------------------------------------------------------------
	eContactTypes ContactPoint::contactType() const {
		return mContactType;
	}

	//-----------------------------------------------------------------------------------------------------------------
	mat ContactPoint::selectionMatrix() const {
		switch (mContactType) {
		case eContactTypes::PFL:
			return B_PFL;
		case eContactTypes::PWF:
			return B_PWF;
		case eContactTypes::SFC:
			return B_SFC;
		default:
			assert(false);
			return mat();
		}
	}

	//-----------------------------------------------------------------------------------------------------------------
	double ContactPoint::maxForce() const {
		return mMaxForce;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double ContactPoint::friction() const {
		return mMu;
	}

	//-----------------------------------------------------------------------------------------------------------------
	double ContactPoint::spinFriction() const {
		return mMus;
	}

	//-----------------------------------------------------------------------------------------------------------------
	arma::mat ContactPoint::aprxWrenchCone(unsigned _pointsPerContact) {
		if (!mComputedAppWrenchCone) {
			if (!mComputedAppForceCone) {
				forceCone(_pointsPerContact);
			}
			mWrenchCone = arma::mat(6, 0);

			for (unsigned i = 0; i < mForceCone.n_cols; i++) {
				arma::colvec wrench(6);
				wrench.rows(0, 2) = mForceCone.col(i);							// Force
				wrench.rows(3, 5) = arma::cross(mPosition, mForceCone.col(i));	// Moments 
				wrench.rows(3, 5) += (i%2 == 0? -1:1)*mMus*mMaxForce*normal();  // 666 Temporarly solution to create a fully dimensional CH
				mWrenchCone.insert_cols(mWrenchCone.n_cols, wrench);
			}


			mComputedAppWrenchCone = true;
		}
		return mWrenchCone;
	}

	//-----------------------------------------------------------------------------------------------------------------
	arma::mat ContactPoint::forceCone(unsigned _pointsPerContact) {
		if (!mComputedAppForceCone) {
			const double C_PI = 3.141598;
			mForceCone = arma::mat(3, 0);
			for (unsigned i = 0; i < _pointsPerContact; i++) {

				arma::colvec force = -mMaxForce*mContactFrame*arma::colvec(
				{ 1,
					mMu*cos(2 * C_PI * i / _pointsPerContact),
					mMu*sin(2 * C_PI * i / _pointsPerContact),
				});

				mForceCone.insert_cols(mForceCone.n_cols, force);
			}
			mComputedAppForceCone = true;
		}

		return mForceCone;
	}

	//-----------------------------------------------------------------------------------------------------------------
	ContactPoint ContactPoint::sample() {
		auto newPos = mPosition + mNormal*arma::randn()*mPositionCovariance(0,0);
		auto newNormal = mNormal + arma::colvec({	arma::randn()*mNormalCovariance(0,0),
													arma::randn()*mNormalCovariance(1,1),
													arma::randn()*mNormalCovariance(2,2)});
		return ContactPoint(newPos, mPositionCovariance,
							newNormal, mNormalCovariance, 
							mContactType, mMaxForce, mMu, mMus);
	}

	//-----------------------------------------------------------------------------------------------------------------
	void ContactPoint::computeContactFrame() {
		if (!mInitContactFrame) {
			if (abs(dot(mNormal, C_AXIS_X)) < 0.8) {
				mContactFrame.col(0) = mNormal;
				mContactFrame.col(2) = cross(mNormal, C_AXIS_X);
				mContactFrame.col(1) = cross(mContactFrame.col(2), mNormal);
			}
			else {
				mContactFrame.col(0) = mNormal;
				mContactFrame.col(1) = cross(mNormal, C_AXIS_Y);
				mContactFrame.col(2) = cross(mContactFrame.col(1), mNormal);
			}
		}
		mInitContactFrame = true;
	}
}