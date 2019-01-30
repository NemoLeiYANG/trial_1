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
#ifndef SCAN_TO_MAP_MATCHING
#define SCAN_TO_MAP_MATCHING
#include <ros/ros.h>
#include <ros/callback_queue.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <geometry_msgs/PointStamped.h>
#include <nav_msgs/Path.h>

#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/io/pcd_io.h>

#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>

#include <Eigen/Geometry>
#include <Eigen/Core>

#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <fstream>
#include <yaml-cpp/yaml.h>


typedef pcl::PointCloud< pcl::PointXYZI > PointCloudI;
typedef pcl::PointCloud< pcl::PointNormal > PointCloudN;

typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef pcl::PointXYZI PointType;

struct intpoint {
    int x, y, z;
    intpoint(int x, int y, int z): x(x), y(y), z(z) {}
    friend bool operator<(const intpoint &l, const intpoint &r) {
        return l.x<r.x || (l.x==r.x && l.y<r.y) || (l.x==r.x && l.y==r.y && l.z<r.z);
    }
};



class scan_to_map_matching {
public:
	scan_to_map_matching(bool flag_dense_edge_received, bool flag_dense_plane_received, 
		bool flag_odom_received, bool flag_dewarped_cloud_received);

	float stdDeviation(Eigen::VectorXf &x);

	void matrixToVec(float *x, Eigen::Matrix4f H);

    void DenseEdgePointsCallback(const pcl::PointCloud< pcl::PointXYZI >::ConstPtr& cloud);

	void DensePlanePointsCallback(const pcl::PointCloud< pcl::PointXYZI >::ConstPtr& cloud);

	void OdometryCallback(const nav_msgs::Odometry::ConstPtr& odom);

	void VelodyenDewarpedCloudCallback(const pcl::PointCloud< pcl::PointXYZI >::ConstPtr& cloud);

    void RelocalizeCallback(const nav_msgs::Odometry::ConstPtr &odom);

	Eigen::Matrix4f vecToMatrix(float *x);

	intpoint get_map_cell(const pcl::PointXYZI p);

	void process(scan_to_map_matching &lidar);

	void initialize();

	void update_feature();

	void update_cells();

	void update_neighbor();

	void compute_edge_jacobian(edgePlaneJacobian lidar_jacobian, const float _correspondence_thresh);

	void compute_plane_jacobian(edgePlaneJacobian lidar_jacobian, const float _correspondence_thresh);

    void correctMap(const nav_msgs::Path::ConstPtr &path);


protected:
	bool flag_dense_edge_received;
	bool flag_dense_plane_received;
	bool flag_odom_received;
	bool flag_dewarped_cloud_received;
    bool flag_relocalize_received = false;
	double map_cell_size;
	nav_msgs::Odometry::Ptr received_odom;
	Eigen::Matrix4f trans_odom_get;
    Eigen::Matrix4f relocalize_odom;
	float transform[6] = {0};
	double received_edge_timestamp;
	double received_plane_timestamp;
	double received_odom_timestamp;
	double received_dewarped_timestamp;
    double received_relocalize_timestamp;
	double map_edge_voxel_size;
	double map_plane_voxel_size;
	double full_map_voxel_size;
    int mapMatchingScan;
    double lambda;
    double min_error;
    double _correspondence_thresh_;
    double bisWeight;
    double edge_timestamp;
    double plane_timestamp;
    double odom_timestamp;
    double dewarped_timestamp;
    double init_distance;
    int received_counter = 0;
    int matched_edge_num, matched_plane_num;
    std::string PCD_outdir;
    int save_interval;


	std::map<intpoint, PointCloudI::Ptr> _global_edges, _global_planes;
    // pose map to store the pose of the lidar odometry corresponding to each lidar time stamp
    std::map<double, Eigen::Matrix4f> poseMap;
    std::set<intpoint> map_cells;
    std::set<intpoint> edges_map_cells;
    std::set<intpoint> planes_map_cells;
    std::set<int> invalidEdgeFeatureInd;
    std::set<int> invalidSurfFeatureInd;

    pcl::PointXYZI input_point;
    pcl::KdTreeFLANN<pcl::PointXYZI> edges_map_kdtree;
    pcl::KdTreeFLANN<pcl::PointXYZI> planes_map_kdtree;

     // transformation matrix from lidar coordinate system to world coordinate system
    Eigen::Matrix4f trans_map;
    // the previous transformation from lidar coordinate system to world coordinate system(VIN's measurement)
    Eigen::Matrix4f trans_odom_pre;
    // the current transformation from lidar coordinate system to world frame (VIN's measurement)
    Eigen::Matrix4f trans_odom_cur;
    // relative transformation of current scan (VIN's measurement)
    Eigen::Matrix4f dtrans_odom;

    PointCloudI::Ptr dense_edges;
    pcl::PointCloud<PointType>::Ptr received_dense_edge_points;
    pcl::PointCloud<PointType>::Ptr received_dense_planar_points;
    pcl::PointCloud<PointType>::Ptr received_dewarped_cloud;
    PointCloudI::Ptr dense_planes;
    PointCloudI::Ptr dense_edges_trans;
    PointCloudI::Ptr dense_planes_trans;
    PointCloudI::Ptr point_residuals;
    PointCloudN::Ptr jacobians;
    PointCloudI::Ptr edges_map_near;
    PointCloudI::Ptr planes_map_near;
    PointCloudI::Ptr edges_map_full;
    PointCloudI::Ptr planes_map_full;
    PointCloudI::Ptr dewarped_cloud;
    PointCloudI::Ptr edes_map_all;
    PointCloudI::Ptr planes_map_all;
    PointCloudI::Ptr dewarped_points_all;

    pcl::ApproximateVoxelGrid<pcl::PointXYZI> edge_voxel_filter;
    pcl::ApproximateVoxelGrid<pcl::PointXYZI> plane_voxel_filter;
    pcl::ApproximateVoxelGrid<pcl::PointXYZI> full_map_voxel_filter;

};
#endif