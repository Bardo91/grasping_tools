///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <grasping_tools/hands/GripperHandTwoFingers.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/io/ply_io.h>

#include <grasping_tools/mathTools.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <pcl/PCLPointCloud2.h>

namespace grasping_tools {
	GripperHandTwoFingers::GripperHandTwoFingers(double _aperture) {
		mAperture = _aperture;
        //mModel = "/home/bardo91/programing/gpisGrasping/projects/grasping/hands/gripper/";
        //
        //if (pcl::io::loadPolygonFileSTL(mModel + "base.stl", mMeshBase) == 0) {
        //    PCL_ERROR("Failed to load STL file\n");
        //}
        //if (pcl::io::loadPolygonFileSTL(mModel + "leftGripper.stl", mLeftGripper) == 0) {
        //    PCL_ERROR("Failed to load STL file\n");
        //}
        //if (pcl::io::loadPolygonFileSTL(mModel + "rightGripper.stl", mRightGripper) == 0) {
        //    PCL_ERROR("Failed to load STL file\n");
        //}
	}

	pcl::PolygonMesh GripperHandTwoFingers::generateHandModel(const Grasp & _grasp, double _rotation) {
        // auto cps = _grasp.contactPoints();
		// pcl::PointCloud<pcl::PointXYZ> finalMeshCloud;
		pcl::PolygonMesh finalMesh;

        // if(!mLoadedModel){
        //     return finalMesh;
        // }

		// // Compute aperture
		// arma::vec vAperture = cps[1].position() - cps[0].position();
		// double separation = arma::norm(vAperture);
		// pcl::PointCloud<pcl::PointXYZ> cloudGripper1;
		// pcl::fromPCLPointCloud2(mLeftGripper.cloud, cloudGripper1);

		// pcl::PointCloud<pcl::PointXYZ> cloudGripper2;
		// pcl::fromPCLPointCloud2(mRightGripper.cloud, cloudGripper2);

		// pcl::PointCloud<pcl::PointXYZ> cloudGripper3;
		// pcl::fromPCLPointCloud2(mMeshBase.cloud, cloudGripper3);

		// Eigen::Matrix4f t = Eigen::Matrix4f::Identity();
		// t(1, 3) = -0.36;
		// pcl::transformPointCloud(cloudGripper3, cloudGripper3, t);

		// t = Eigen::Matrix4f::Identity();
		// t(0, 3) = -separation / 2;
		// pcl::transformPointCloud(cloudGripper1, cloudGripper1, t);
		// t(0, 3) = separation / 2;
		// pcl::transformPointCloud(cloudGripper2, cloudGripper2, t);

		// finalMeshCloud.points.insert(finalMeshCloud.end(), cloudGripper1.begin(), cloudGripper1.end());
		// finalMesh.polygons.insert(finalMesh.polygons.end(), mLeftGripper.polygons.begin(), mLeftGripper.polygons.end());
		// finalMeshCloud.points.insert(finalMeshCloud.end(), cloudGripper2.begin(), cloudGripper2.end());
		// finalMesh.polygons.insert(finalMesh.polygons.end(), mRightGripper.polygons.begin(), mRightGripper.polygons.end());
		// for (unsigned i = 0; i < mRightGripper.polygons.size(); i++) {
		// 	finalMesh.polygons[i + mLeftGripper.polygons.size()].vertices[0] += mLeftGripper.cloud.width*mLeftGripper.cloud.height;
		// 	finalMesh.polygons[i + mLeftGripper.polygons.size()].vertices[1] += mLeftGripper.cloud.width*mLeftGripper.cloud.height;
		// 	finalMesh.polygons[i + mLeftGripper.polygons.size()].vertices[2] += mLeftGripper.cloud.width*mLeftGripper.cloud.height;
		// }

		// int prevPoints = finalMeshCloud.size();
		// int prevPolygons = finalMesh.polygons.size();
		// finalMeshCloud.points.insert(finalMeshCloud.end(), cloudGripper3.begin(), cloudGripper3.end());
		// finalMesh.polygons.insert(finalMesh.polygons.end(), mMeshBase.polygons.begin(), mMeshBase.polygons.end());

		// for (unsigned i = 0; i < mMeshBase.polygons.size(); i++) {
		// 	finalMesh.polygons[i + prevPolygons].vertices[0] += prevPoints;
		// 	finalMesh.polygons[i + prevPolygons].vertices[1] += prevPoints;
		// 	finalMesh.polygons[i + prevPolygons].vertices[2] += prevPoints;
		// }


		// // Rotate pointcloud of mesh
		// t = Eigen::Matrix4f::Identity();
		// auto rot = rotUnitaryVectors({ 1,0,0 }, arma::vec(cps[1].position() - cps[0].position()));
		// t.block<3, 3>(0, 0) = Eigen::Matrix3d(&rot(0, 0)).cast<float>();
		// pcl::transformPointCloud(finalMeshCloud, finalMeshCloud, t);

		// // Set extra rotation
		// t = Eigen::Matrix4f::Identity();
		// arma::vec cp1 = cps[0].position();
		// arma::vec cp2 = cps[1].position();
		// Eigen::Vector3f axis = { float(cp2[0] - cp1[0]) ,float(cp2[1] - cp1[1]),float(cp2[2] - cp1[2]) };
		// axis /= axis.norm();
		// auto rotExtra = Eigen::AngleAxisf(_rotation, axis)*
		// 			Eigen::AngleAxisf(0, Eigen::Vector3f::UnitX());
		// t.block<3, 3>(0, 0) = rotExtra.matrix();
		// pcl::transformPointCloud(finalMeshCloud, finalMeshCloud, t);

		// // Translate pointcloud of mesh
		// t = Eigen::Matrix4f::Identity();
		// t(0, 3) = cps[0].position()[0] + vAperture[0]/2;
		// t(1, 3) = cps[0].position()[1] + vAperture[1]/2;
		// t(2, 3) = cps[0].position()[2] + vAperture[2]/2;
		// pcl::transformPointCloud(finalMeshCloud, finalMeshCloud, t);


		// pcl::toPCLPointCloud2(finalMeshCloud, finalMesh.cloud);

		return finalMesh;
	}
}
