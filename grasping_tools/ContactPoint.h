///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef GPISGRASPING_CONTACTPOINT_H_
#define GPISGRASPING_CONTACTPOINT_H_

#include <armadillo>

namespace gpisGrasping {
	// Contact types.
	enum eContactTypes {PFL, PWF, SFC};

	/// Selection matrix of contact type. Punctual friction less (r = 1).
	const arma::colvec6 B_PFL = { { 0, 0, 1, 0, 0, 0 } };

	/// Selection matrix of contact type. Punctual with friction (r = 3).
	const arma::mat B_PWF = {	{ 1, 0, 0 },
								{ 0, 1, 0 },
								{ 0, 0, 1 },
								{ 0, 0, 0 },
								{ 0, 0, 0 },
								{ 0, 0, 0 } };

	/// Selection matrix of contact type. Soft contact (r = 4).
	const arma::mat B_SFC = {	{ 1, 0, 0, 0 },
								{ 0, 1, 0, 0 },
								{ 0, 0, 1, 0 },
								{ 0, 0, 0, 0 },
								{ 0, 0, 0, 0 },
								{ 0, 0, 0, 1 } };

	const arma::colvec3 C_AXIS_X = { 1, 0, 0 };
	const arma::colvec3 C_AXIS_Y = { 0, 1, 0 };
	const arma::colvec3 C_AXIS_Z = { 0, 0, 1 };

	class   ContactPoint {
	public:
		ContactPoint(const arma::colvec3 &_position, const arma::mat33 &_covPos, const arma::colvec3 &_n, const arma::mat33 &_covNormal, const eContactTypes _contactType, double _maxForce, double _friction = 0, double _spinFriction = 0);

		arma::colvec3	position()	const;
		arma::mat33		posCov()	const;
		arma::colvec3	normal()	const;
		arma::mat33		norCov()	const;
		arma::mat33		frame();

		eContactTypes	contactType()		const;
		arma::mat		selectionMatrix()	const;

		double			maxForce()		const;
		double			friction()		const;
		double			spinFriction()	const;

		/// Compute the 6-dimensional cone resulting of the union of the friction cones of the contact points.
		/// \param _pointsPerContact: number of points used to approximate the each of the contact friction crones
		arma::mat aprxWrenchCone(unsigned _pointsPerContact);

		arma::mat forceCone(unsigned _pointsPerContact);

		/// Sample new cp using covariances
		ContactPoint sample();
	private:
		void computeContactFrame();

	private:
		arma::colvec3	mPosition;	// Position of the contact point.
		arma::mat33		mPositionCovariance;
		arma::colvec3	mNormal;			// Normal of the contact point (pointing inward).
		arma::mat33		mNormalCovariance;

		arma::mat33		mContactFrame;
		bool			mInitContactFrame = false;

		eContactTypes mContactType;	// The contact type determine the selection matrix (H) that select how the contact affect to the object.
									// In 3d: 
									//		-> Punctual friction less: H is 6x1.
									//		-> Punctual with friction: H is 6x3.
									//		-> Soft contact: H is 6x4.

		double mMaxForce;	// Maximum force that can be applied by this finger.

		double mMu;			// Tangential friction coefficient.
		double mMus;		// Spin friction coefficient.

		bool mComputedAppWrenchCone = false;
		arma::mat mWrenchCone;
		bool mComputedAppForceCone = false;
		arma::mat mForceCone;

	};
}	//	gpisGrasping

#endif	//	GPISGRASPING_CONTACTPOINT_H_