/*v1*************************************************************************************/
/* This is developed at Carnegie Mellon University in collaboration with Autel Robotics */
/*                                                                                      */
/* PI:                                                                                  */
/* George Kantor                                                                        */
/*                                                                                      */
/* Authors:                                                                             */ 
/* Weizhao Shao                                                                         */
/* Cong Li                                                                              */
/* Srinivasan Vijayarangan                                                              */
/*                                                                                      */
/* Please refer to the contract document for details on license/copyright information.  */
/****************************************************************************************/
#ifndef LIDARSCANREGISTRATION_H
#define LIDARSCANREGISTRATION_H

#include <cmath>
#include <vector>
#include <fstream>

#include <time_based_retriever.h>
#include <opencv/cv.h>
#include <nav_msgs/Odometry.h>
#include <opencv/cv.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl_ros/transforms.h>
#include <pcl_ros/point_cloud.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <time.h>
#include <string>
#include <yaml-cpp/yaml.h>
#include <queue>
#include <unordered_set>

typedef pcl::PointCloud< pcl::PointXYZI > PointCloudI;
typedef pcl::PointXYZI PointType;

class LidarScanRegistration
{
public:
	LidarScanRegistration(int frameCount, int initFrameCount, Eigen::Matrix4f input_matrix, bool auto_delete,
  	int N_SCANS);

	void inline laserCloudHandler(const sensor_msgs::PointCloud2ConstPtr& laserCloudMsg);

	void inline OdometryCallback(const nav_msgs::Odometry::ConstPtr& odom);

	void inline transform_to_end(PointType pi, PointType& po, Eigen::Matrix4f dtrans);

	void inline myhomogeneousToVect(float* x, Eigen::Matrix4f H);

	void inline getScanInd(const sensor_msgs::PointCloud2ConstPtr& laserCloudMsg,
	std::vector<pcl::PointCloud<PointType> >& laserCloudScans, int &pointsNum);

	Eigen::Matrix3f  inline myangleAxisToEigenRot(float *x);

	Eigen::Matrix4f inline myvectToHomogeneous(float *x);

	void inline get_feature_points(pcl::PointCloud<PointType>::Ptr laserCloud,
		                           std::vector<int> scanStartInd, std::vector<int> scanEndInd, std::vector<int> pointInd,
		                           std::vector<float> pointCurvature, pcl::PointCloud<PointType>& cornerPointsSharp,
  	                               pcl::PointCloud<PointType>& surfPointsFlat,
  	                               std::unordered_set<int> selectedPoint);

    void inline get_dewarped_points(pcl::PointCloud<PointType> point_cloud,
		                            double cur_time, Eigen::Matrix4f cur_close_odometry,
                                pcl::PointCloud< pcl::PointXYZI >::Ptr dewarped_edge_points);

	void inline get_dewarped_all_points(pcl::PointCloud<PointType>::Ptr point_cloud,
		                            double cur_time, Eigen::Matrix4f cur_close_odometry,
                                pcl::PointCloud< pcl::PointXYZI >::Ptr dewarped_edge_points);
    void inline get_dewarped_points_const_velo(pcl::PointCloud<PointType> point_cloud,
                                                                  double cur_time, Eigen::Matrix4f prev_prev_odom, Eigen::Matrix4f prev_odom, 
                                                                  pcl::PointCloud< pcl::PointXYZI >::Ptr dewarped_points);
	void refinedPoseCallback(const nav_msgs::Odometry::ConstPtr &odom);
protected:
	int frameCount;
	int initFrameCount;
	double odom_timestamp_pre;
	double scanPeriod;
	int N_SCANS;
	// segment one scan line to small regions, each region extracts pre-defined number of edge points and surface
	int lineSegNum;
	int maxEdgeNum;
	int maxSurfNum;
	double edgeThreshold;
	bool use_vio_odom;



	TimeBasedRetriever<Eigen::Matrix4f> odoMsg;

	Eigen::Matrix4f trans_odom_prev;
	Eigen::Matrix4f dtrans;
	Eigen::Matrix4f prev_close_odometry;
	Eigen::Matrix4f prev_odom;
	Eigen::Matrix4f prev_prev_odom;

	pcl::PointCloud<PointType> cornerPointsSharp_prev;
	pcl::PointCloud<PointType> surfPointsFlat_prev;

	pcl::PointCloud<PointType> cornerPointsSharp_prev_prev;
	pcl::PointCloud<PointType> surfPointsFlat_prev_prev;

	pcl::PointCloud<PointType>::Ptr laserCloud_prev;
	pcl::PointCloud<PointType>::Ptr laserCloud_prev_prev;


	std::vector<float> pointCurvature;
	std::vector<int> pointInd;
	std::unordered_set<int> selectedPoint;

	std::vector<std::vector<int> > lineIndArr;
	std::vector<std::vector<float> > timePortionOfScan;
	//std::vector<std::priority_queue<std::pair<float, int> > > surfaceFeatureInd;
	//std::vector<std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int> >,
	 //std::greater<std::pair<float, int> > > > lineFeatureInd;

	 std::vector<std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int> >,
	 std::greater<std::pair<float, int> > > > surfaceFeatureInd;
	 std::vector<std::priority_queue<std::pair<float, int> > > lineFeatureInd;
    float transform[6] = {0};
	double laser_timestamp_pre = 0;
	double time_count = 0;

	clock_t start, end;
};

#endif // LIDARSCANREGISTRATION_H
