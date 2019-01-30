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
#include "edge_plane_fns.h"
#include "scan_to_map_matching.h"
#include <pcl/io/ply_io.h>
#include "std_msgs/String.h"
#include <sstream>
#include <pcl/visualization/cloud_viewer.h>
#include <Eigen/Eigenvalues>
std::ofstream myfile;

#define PI 3.141592
nav_msgs::Path path;
bool flag = false;
bool startingFlag = true;
nav_msgs::Odometry::ConstPtr curOdom;
std::map<long, pcl::PointCloud<PointType>::ConstPtr> edgeMap;
std::map<long, pcl::PointCloud<PointType>::ConstPtr> planeMap;
ros::Publisher correct_finish;
bool refineMap = false;
std::pair<PointCloudI::Ptr, PointCloudI::Ptr> edge_match;
std::pair<PointCloudI::Ptr, PointCloudI::Ptr> plane_match;
bool correct = false;
ros::Publisher pub_registered_cloud; 
std::string pointCloudSavePath = "/home/cmuautel/Data/lidar_scan/";
std::string all_map = "/home/cmuautel/Data/all_map/";
int seq = 0;

scan_to_map_matching::scan_to_map_matching(bool flag_dense_edge_received, bool flag_dense_plane_received,
                                           bool flag_odom_received, bool flag_dewarped_cloud_received) : flag_dense_edge_received(flag_dense_edge_received),
                                                                                                         flag_dense_plane_received(flag_dense_plane_received),
                                                                                                         flag_odom_received(flag_odom_received),
                                                                                                         flag_dewarped_cloud_received(flag_dewarped_cloud_received)
{
    trans_odom_get.setIdentity();
    trans_map.setIdentity();
    trans_odom_pre.setIdentity();
    trans_odom_cur.setIdentity();
    dtrans_odom.setIdentity();
    received_dense_edge_points = pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>());
    received_dense_planar_points = pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>());
    received_dewarped_cloud = pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>());
    // load yaml file
    std::string yaml_file = "/home/cong/catkin_ws/src/lidar_odometry/config/indoor.yaml";
    YAML::Node config = YAML::LoadFile(yaml_file);
    // Do Map Matching Once in the Number of scans
    mapMatchingScan = config["map_matching_scan"].as<int>();
    std::cout << "mapMatchingScan: " << mapMatchingScan << std::endl;
    lambda = config["Lambda"].as<double>();
    std::cout << "lambda: " << lambda << std::endl;
    min_error = config["min_error"].as<double>();
    std::cout << "min_error: " << min_error << std::endl;
    _correspondence_thresh_ = config["_correspondence_thresh_"].as<double>();
    std::cout << "correspondence_thresh: " << _correspondence_thresh_ << std::endl;
    bisWeight = config["bis_weight"].as<double>();
    std::cout << "bisWeight: " << bisWeight << std::endl;
    map_cell_size = config["map_cell_size"].as<double>();
    std::cout << "map cell size: " << map_cell_size << std::endl;
    map_edge_voxel_size = config["map_edge_voxel_size"].as<double>();
    std::cout << "map edge voxel size: " << map_edge_voxel_size << std::endl;
    map_plane_voxel_size = config["map_plane_voxel_size"].as<double>();
    std::cout << "map plane voxel size: " << map_plane_voxel_size << std::endl;
    full_map_voxel_size = config["full_map_voxel_size"].as<double>();
    std::cout << "full_map_voxel_size: " << full_map_voxel_size << std::endl;
    PCD_outdir = config["PCD_outdir"].as<std::string>();
    std::cout << "PCD_outdir: " << PCD_outdir << std::endl;
    save_interval = config["save_interval"].as<int>();
    std::cout << "save_interval: " << save_interval << std::endl;
    init_distance = config["init_distance"].as<double>();
    std::cout << "init_distance: " << init_distance << std::endl;
}

/*
  transform the transformation matrix to vector 
*/
void scan_to_map_matching::matrixToVec(float *x, Eigen::Matrix4f H)
{
    //ceres method is column major
    // 1 4 7
    // 2 5 8
    // 3 6 9

    float r[9];
    r[0] = H(0, 0);
    r[1] = H(1, 0);
    r[2] = H(2, 0);
    r[3] = H(0, 1);
    r[4] = H(1, 1);
    r[5] = H(2, 1);
    r[6] = H(0, 2);
    r[7] = H(1, 2);
    r[8] = H(2, 2);

    x[3] = H(0, 3);
    x[4] = H(1, 3);
    x[5] = H(2, 3);

    ceres::RotationMatrixToAngleAxis(r, x);
}

/*
   transform the 1 by 6 vector to transformation matrix
*/
Eigen::Matrix4f scan_to_map_matching::vecToMatrix(float *x)
{
    Eigen::Matrix4f H(Eigen::Matrix4f::Identity());

    Eigen::Matrix3f rot;

    float r[9];
    ceres::AngleAxisToRotationMatrix(x, r);
    //ceres method is column major
    // 1 4 7
    // 2 5 8
    // 3 6 9
    H(0, 0) = r[0];
    H(1, 0) = r[1];
    H(2, 0) = r[2];
    H(0, 1) = r[3];
    H(1, 1) = r[4];
    H(2, 1) = r[5];
    H(0, 2) = r[6];
    H(1, 2) = r[7];
    H(2, 2) = r[8];
    H(0, 3) = x[3];
    H(1, 3) = x[4];
    H(2, 3) = x[5];

    return H;
}
// subscribe edge feature points
void scan_to_map_matching::DenseEdgePointsCallback(const pcl::PointCloud<pcl::PointXYZI>::ConstPtr &cloud)
{
    received_edge_timestamp = pcl_conversions::fromPCL(cloud->header.stamp).toSec();
    *received_dense_edge_points = *cloud;
    flag_dense_edge_received = true;
    //add the edge point to the map
    double time = cloud->header.stamp / 1000000.0;
    long time1 = (long)(time * 10000);
    edgeMap[time1] = cloud;
}

// subscribe surface feature points
void scan_to_map_matching::DensePlanePointsCallback(const pcl::PointCloud<pcl::PointXYZI>::ConstPtr &cloud)
{
    received_plane_timestamp = pcl_conversions::fromPCL(cloud->header.stamp).toSec();

    *received_dense_planar_points = *cloud;
    flag_dense_plane_received = true;
    //add plane points to the map
    double time = cloud->header.stamp / 1000000.0;
    long time1 = (long)(time * 10000);
    planeMap[time1] = cloud;
}

// subscribe odometry information
void scan_to_map_matching::OdometryCallback(const nav_msgs::Odometry::ConstPtr &odom)
{
    received_odom_timestamp = odom->header.stamp.toSec();
    curOdom = odom;

    Eigen::Quaternionf q_b;
    Eigen::Matrix3f rot_b2o; // body to origin
    q_b.x() = odom->pose.pose.orientation.x;
    q_b.y() = odom->pose.pose.orientation.y;
    q_b.z() = odom->pose.pose.orientation.z;
    q_b.w() = odom->pose.pose.orientation.w;
    rot_b2o = q_b.toRotationMatrix(); //Eigen::Matrix<float, 3, 3>::Identity();
    trans_odom_get.matrix().topLeftCorner<3, 3>() = rot_b2o;
    trans_odom_get(0, 3) = odom->pose.pose.position.x;
    trans_odom_get(1, 3) = odom->pose.pose.position.y;
    trans_odom_get(2, 3) = odom->pose.pose.position.z; // odom: k to 0 transform

    flag_odom_received = true;
}

// subscribe dewarped mapping points
void scan_to_map_matching::VelodyenDewarpedCloudCallback(const pcl::PointCloud<pcl::PointXYZI>::ConstPtr &cloud)
{
    received_dewarped_timestamp = pcl_conversions::fromPCL(cloud->header.stamp).toSec();

    *received_dewarped_cloud = *cloud;
    flag_dewarped_cloud_received = true;
}

// subscribe relocalization odometry information
void scan_to_map_matching::RelocalizeCallback(const nav_msgs::Odometry::ConstPtr &odom)
{
    received_relocalize_timestamp = odom->header.stamp.toSec();

    Eigen::Quaternionf q_b;
    Eigen::Matrix3f rot_b2o; // body to origin
    q_b.x() = odom->pose.pose.orientation.x;
    q_b.y() = odom->pose.pose.orientation.y;
    q_b.z() = odom->pose.pose.orientation.z;
    q_b.w() = odom->pose.pose.orientation.w;
    rot_b2o = q_b.toRotationMatrix(); //Eigen::Matrix<float, 3, 3>::Identity();
    relocalize_odom.matrix().topLeftCorner<3, 3>() = rot_b2o;
    relocalize_odom(0, 3) = odom->pose.pose.position.x;
    relocalize_odom(1, 3) = odom->pose.pose.position.y;
    relocalize_odom(2, 3) = odom->pose.pose.position.z; // odom: k to 0 transform

    flag_relocalize_received = true;
}

// subscribe refined pose graph information to correct map
void scan_to_map_matching::correctMap(const nav_msgs::Path::ConstPtr &path)
{
    clock_t start = clock();
    refineMap = true;
    for (int i = 0; i < path->poses.size(); i++)
    {
        double curTime = path->poses[i].header.stamp.toSec();
        long curTime1 = (long)(curTime * 10000);
        if (edgeMap.count(curTime1) > 0 && planeMap.count(curTime1) > 0)
        {

            //get the corrsponding point cloud
            geometry_msgs::Pose curPose = path->poses[i].pose;
            *received_dense_edge_points = *(edgeMap[curTime1]);
            *received_dense_planar_points = *(planeMap[curTime1]);

            Eigen::Matrix4f trans_map;
            Eigen::Quaternionf q_b;
            Eigen::Matrix3f rot_b2o; // body to origin
            q_b.x() = curPose.orientation.x;
            q_b.y() = curPose.orientation.y;
            q_b.z() = curPose.orientation.z;
            q_b.w() = curPose.orientation.w;

            rot_b2o = q_b.toRotationMatrix(); //Eigen::Matrix<float, 3, 3>::Identity();

            trans_map.matrix().topLeftCorner<3, 3>() = rot_b2o;
            trans_map(0, 3) = curPose.position.x;
            trans_map(1, 3) = curPose.position.y;
            trans_map(2, 3) = curPose.position.z; // odom: k to 0 transform

            Eigen::Matrix4f T_L_I;
            T_L_I << 1, 0, 0, -0.011356,
                0, -1, 0, -0.002352,
                0, 0, -1, -0.08105,
                0, 0, 0, 1;
            trans_map = trans_map * T_L_I.inverse();
            pcl::PointCloud<PointType>::Ptr cur_edge_trans = pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>());
            pcl::PointCloud<PointType>::Ptr cur_plane_trans = pcl::PointCloud<PointType>::Ptr(new pcl::PointCloud<PointType>());
            pcl::transformPointCloud(*received_dense_edge_points, *cur_edge_trans, trans_map);
            pcl::transformPointCloud(*received_dense_edge_points, *cur_plane_trans, trans_map);
            *dense_edges_trans += *cur_edge_trans;
            *dense_planes_trans += *cur_plane_trans;
        }
    }
    *dewarped_points_all = *dense_edges_trans + *dense_planes_trans;
    pub_registered_cloud.publish(dewarped_points_all);
    //pcl::io::savePCDFile(PCD_outdir, *dewarped_points_all, true); // Binary format*/
    _global_edges.clear();
    _global_planes.clear();
    update_cells();
    std_msgs::String msg;
    std::stringstream ss;
    ss << "pose correct finished";
    msg.data == ss.str();
    correct_finish.publish(msg);
    refineMap = false;
    clock_t timeElapsed = clock() - start;
    unsigned msElapsed = timeElapsed / (CLOCKS_PER_SEC / 1000);
    std::cout << "time: " << msElapsed << " ms\n";
    correct = true;
}

// calculate the cube index of the point
intpoint scan_to_map_matching::get_map_cell(const pcl::PointXYZI curPlaneFeature)
{

    int x = floor(curPlaneFeature.x / map_cell_size),
        y = floor(curPlaneFeature.y / map_cell_size),
        z = floor(curPlaneFeature.z / map_cell_size);
    return intpoint(x, y, z);
}

// intilize variables
void scan_to_map_matching::initialize()
{
    dense_edges = PointCloudI::Ptr(new PointCloudI);
    dense_planes = PointCloudI::Ptr(new PointCloudI);
    dense_edges_trans = PointCloudI::Ptr(new PointCloudI);
    dense_planes_trans = PointCloudI::Ptr(new PointCloudI);
    point_residuals = PointCloudI::Ptr(new PointCloudI);
    jacobians = PointCloudN::Ptr(new PointCloudN);
    edges_map_near = PointCloudI::Ptr(new PointCloudI);
    planes_map_near = PointCloudI::Ptr(new PointCloudI);
    edges_map_full = PointCloudI::Ptr(new PointCloudI);
    planes_map_full = PointCloudI::Ptr(new PointCloudI);
    dewarped_cloud = PointCloudI::Ptr(new PointCloudI);
    edes_map_all = PointCloudI::Ptr(new PointCloudI);
    planes_map_all = PointCloudI::Ptr(new PointCloudI);
    dewarped_points_all = PointCloudI::Ptr(new PointCloudI);
    edges_map_near->header.frame_id = "world";
    planes_map_near->header.frame_id = "world";
    edges_map_full->header.frame_id = "world";
    planes_map_full->header.frame_id = "world";
    dewarped_cloud->header.frame_id = "world";
    planes_map_all->header.frame_id = "world";
    edes_map_all->header.frame_id = "world";
    dewarped_points_all->header.frame_id = "world";
}

void scan_to_map_matching::update_feature()
{
    *dense_edges = *received_dense_edge_points;
    *dense_planes = *received_dense_planar_points;
    *dewarped_cloud = *received_dewarped_cloud;

    edge_timestamp = received_edge_timestamp;
    plane_timestamp = received_plane_timestamp;
    odom_timestamp = received_odom_timestamp;
    dewarped_timestamp = received_dewarped_timestamp;
}

void scan_to_map_matching::update_cells()
{

    map_cells.clear();

    // update map_cells, edges_map_cells and _global_edges
    for (unsigned int i = 0; i < dense_edges_trans->size(); i++)
    {

        pcl::PointXYZI curEdgeFeature = dense_edges_trans->at(i);
        intpoint map_cell_ind = get_map_cell(curEdgeFeature);
        map_cells.insert(map_cell_ind);
        if (!_global_edges.count(map_cell_ind))
        {
            _global_edges[map_cell_ind] = PointCloudI::Ptr(new PointCloudI);
            edges_map_cells.insert(map_cell_ind);
        }
        _global_edges[map_cell_ind]->push_back(curEdgeFeature);
    }

    // filter global edge points
    for (auto iter = map_cells.begin();
         iter != map_cells.end();
         iter++)
    {

        intpoint map_cell_ind = *iter;
        edge_voxel_filter.setInputCloud(_global_edges[map_cell_ind]);
        edge_voxel_filter.filter(*(_global_edges[map_cell_ind]));
    }

    map_cells.clear();

    // update map_cells, planes_map_cells and _global_planes
    for (unsigned int i = 0; i < dense_planes->size(); i++)
    {

        pcl::PointXYZI curPlaneFeature = dense_planes_trans->at(i);
        intpoint map_cell_ind = get_map_cell(curPlaneFeature);
        map_cells.insert(map_cell_ind);
        if (!_global_planes.count(map_cell_ind))
        {
            _global_planes[map_cell_ind] = PointCloudI::Ptr(new PointCloudI);
            planes_map_cells.insert(map_cell_ind);
        }
        _global_planes[map_cell_ind]->push_back(curPlaneFeature);
    }

    // filter global plane points
    for (auto iter = map_cells.begin();
         iter != map_cells.end();
         iter++)
    {
        intpoint map_cell_ind = *iter;

        plane_voxel_filter.setInputCloud(_global_planes[map_cell_ind]);
        plane_voxel_filter.filter(*(_global_planes[map_cell_ind]));
    }
}

void scan_to_map_matching::update_neighbor()
{
    edges_map_near->clear();
    planes_map_near->clear();

    map_cells.clear();

    // find cubes containing edge feature points
    for (unsigned int i = 0; i < dense_edges_trans->size(); i++)
    {
        pcl::PointXYZI curEdgeFeature = dense_edges_trans->at(i);
        intpoint map_cell_ind = get_map_cell(curEdgeFeature);
        map_cells.insert(map_cell_ind);
    }

    // update neighbor map of edge feature points
    for (auto iter = map_cells.begin();
         iter != map_cells.end();
         iter++)
    {
        intpoint map_cell_ind = *iter;
        if (_global_edges.count(map_cell_ind))
        {
            *edges_map_near += *_global_edges[map_cell_ind];
        }
    }

    // find cubes containing plane feature points
    map_cells.clear();
    for (unsigned int i = 0; i < dense_planes_trans->size(); i++)
    {
        pcl::PointXYZI curPlaneFeature = dense_planes_trans->at(i);
        intpoint map_cell_ind = get_map_cell(curPlaneFeature);
        map_cells.insert(map_cell_ind);
    }

    // ipdate neighbor map of plane feature points
    for (auto iter = map_cells.begin();
         iter != map_cells.end();
         iter++)
    {
        intpoint map_cell_ind = *iter;
        if (_global_planes.count(map_cell_ind))
        {
            *planes_map_near += *_global_planes[map_cell_ind];
        }
    }

    // update edges and planes kdtree
    //cout << "retrieved edge point size: " << edges_map_near->size() << endl;
    //cout << "retrived plane point size: " << planes_map_near->size() << endl;
    if (edges_map_near->size() == 0 || planes_map_near->size() == 0)
    {
        flag = true;
        return;
    }
    edges_map_kdtree.setInputCloud(edges_map_near);
    planes_map_kdtree.setInputCloud(planes_map_near);
}

// calculate jacobian function for edge feature
void scan_to_map_matching::compute_edge_jacobian(edgePlaneJacobian lidar_jacobian,
                                                 const float _correspondence_thresh_)
{
    int valid_edge =0, nonvalid_edge =0;
    for (unsigned int i = 0; i < dense_edges_trans->size(); i++)
    {
        // previously, we only retrieve 2 closest points for each edge point,
        // now we select 5 closest points for each edge point to form covariance matrix.
        pcl::PointXYZI search(dense_edges_trans->at(i));
        std::vector<int> ids_edge(5);
        std::vector<float> dist2(5);
        // use kd tree to search two closest points
        if (edges_map_kdtree.nearestKSearch(search, 5, ids_edge, dist2) <= 0 || // dist2 is squared distance
            dist2[0] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[1] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[2] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[3] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[4] > _correspondence_thresh_ * _correspondence_thresh_) 
        {
            // <= 0 failed
            continue;
        }


        Eigen::Vector3f vc(0.0, 0.0, 0.0);

        for (int j = 0; j < 5; j++) {
            vc[0] += edges_map_near->at(ids_edge[j]).x;
            vc[1] += edges_map_near->at(ids_edge[j]).y;
            vc[2] += edges_map_near->at(ids_edge[j]).z;
        }

        vc /= 5.0;

        Eigen::Matrix3f tmpMat;
        tmpMat.setZero();

        for (int j = 0; j < 5; j++) {
            Eigen::Vector3f tmp(0.0, 0.0, 0.0);
            tmp[0] = edges_map_near->at(ids_edge[j]).x - vc[0];
            tmp[1] = edges_map_near->at(ids_edge[j]).y - vc[1];
            tmp[2] = edges_map_near->at(ids_edge[j]).z - vc[2];

            tmpMat(0, 0) += tmp[0] * tmp[0];
            tmpMat(0, 1) += tmp[0] * tmp[1];
            tmpMat(0, 2) += tmp[0] * tmp[2];
            tmpMat(1, 1) += tmp[1] * tmp[1];
            tmpMat(1, 2) += tmp[1] * tmp[2];
            tmpMat(2, 2) += tmp[2] * tmp[2];
        }

        tmpMat /= 5.0;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> esolver(tmpMat);
        Eigen::Matrix<float, 1, 3> matD1;
        Eigen::Matrix3f matV1;

        matD1 = esolver.eigenvalues().real();
        matV1 = esolver.eigenvectors().real();

        float x0, y0, z0, x1, y1, z1, x2, y2, z2;


        if (matD1(0, 2) <= 3 * matD1(0, 1)) {
            //std::cout << matD1(0, 2) << matD1(0, 1) << std::endl;
            nonvalid_edge++;
            //std::cout << "***** this is not the edge point *****" << std::endl;
            continue;
        }
        else {
            //std::cout << matD1(0, 2) << matD1(0, 1) << std::endl;
            //std::cout << "!!!!! this is the edge points !!!!!" << std::endl;
            valid_edge++;
            x0 = dense_edges_trans->at(i).x;
            y0 = dense_edges_trans->at(i).y;
            z0 = dense_edges_trans->at(i).z;
            x1 = vc[0] + 0.1 * matV1(0, 2);
            y1 = vc[1] + 0.1 * matV1(1, 2);
            z1 = vc[2] + 0.1 * matV1(2, 2);
            x2 = vc[0] - 0.1 * matV1(0, 2);
            y2 = vc[1] - 0.1 * matV1(1, 2);
            z2 = vc[2] - 0.1 * matV1(2, 2);
        }

        // compute the distance of transformed point to the line formed by two closest points
        float de = sqrt(pow(fabs((x0 - x1) * (y0 - y2) - (x0 - x2) * (y0 - y1)), 2.0) +
                        pow(fabs((x0 - x1) * (z0 - z2) - (x0 - x2) * (z0 - z1)), 2.0) +
                        pow(fabs((y0 - y1) * (z0 - z2) - (y0 - y2) * (z0 - z1)), 2.0)) *
                   1.0 / sqrt(pow(fabs(x1 - x2), 2.0) + pow(fabs(y1 - y2), 2.0) + pow(fabs(z1 - z2), 2.0));

        if (!std::isnan(de))
        {

            input_point = dense_edges->at(i);
            input_point.intensity = de;

            // update residuals
            point_residuals->push_back(input_point);

            // compute Jacobian
            pcl::PointNormal temp_jacobian = lidar_jacobian.calc_edge_jacobians(input_point.x, input_point.y,
                                                                                 input_point.z, x1, y1, z1, x2, y2, z2, x0, y0, z0, transform);
            jacobians->push_back(temp_jacobian);
        }
    }
    std::cout << "valid_edge number: " << valid_edge << std::endl;
    std::cout << "nonvalid_edge number: " << nonvalid_edge << std::endl;

    matched_edge_num = point_residuals->points.size();
}

// compute jacobian function for surface feature points
void scan_to_map_matching::compute_plane_jacobian(edgePlaneJacobian lidar_jacobian,
                                                  const float _correspondence_thresh_)
{
    int valid_suf =0, nonvalid_surf =0;
    for (unsigned int i = 0; i < dense_planes_trans->size(); i++)
    {
        bool valid_flag = true;

        pcl::PointXYZI search(dense_planes_trans->at(i));
        std::vector<int> ids_plane(5);
        std::vector<float> dist2(5);
        // previously, use the kdtree to retrieve the 3 most closing points,
        // now retrieve 5 most closing points to form covariance matrix
        if (planes_map_kdtree.nearestKSearch(search, 3, ids_plane, dist2) <= 0 || // dist2 is squared distance
            dist2[0] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[1] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[2] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[3] > _correspondence_thresh_ * _correspondence_thresh_ ||
            dist2[4] > _correspondence_thresh_ * _correspondence_thresh_)
        {
            // <= 0 means failed
            continue;
        }

        Eigen::Matrix<float, 5, 3> mat1;
        Eigen::Vector3f mat2;
        mat2.setZero();

        for (int j = 0; j < 5; j++) {
            mat1(j, 0) = planes_map_near->at(ids_plane[j]).x;
            mat1(j, 1) = planes_map_near->at(ids_plane[j]).y;
            mat1(j, 2) = planes_map_near->at(ids_plane[j]).z;
        }

        Eigen::Matrix<float, 5, 1> mat3;
        mat3.setConstant(-1);

        mat2 = mat1.colPivHouseholderQr().solve(mat3);

        float p1 = mat2(0, 0);
        float p2 = mat2(1, 0);
        float p3 = mat2(2, 0);
        float p4 = 1;

        float ps = sqrt(p1 * p1 + p2 * p2 + p3 * p3);
        p1 /= ps;
        p2 /= ps;
        p3 /= ps;
        p4 /= ps;


        bool planeValid = true;
        for (int j = 0; j < 5; j++) {
            if (fabs(p1 * planes_map_near->at(ids_plane[j]).x +
                        p2 *planes_map_near->at(ids_plane[j]).y +
                        p3 * planes_map_near->at(ids_plane[j]).z + p4) > 0.2)
               {
                  //std::cout << "***** this is not the surface point *****" << std::endl;
                  valid_flag = false;
                  break;
            }
            else {
                //std::cout << "!!!! this is not the surface point !!!!" << std::endl;
            }
        }
        if (!valid_flag) {
            nonvalid_surf++;
            continue;
        }
        valid_suf++;

        float x0 = dense_planes_trans->at(i).x;
        float y0 = dense_planes_trans->at(i).y;
        float z0 = dense_planes_trans->at(i).z;
        float x1 = planes_map_near->at(ids_plane[0]).x;
        float y1 = planes_map_near->at(ids_plane[0]).y;
        float z1 = planes_map_near->at(ids_plane[0]).z;
        float x2 = planes_map_near->at(ids_plane[1]).x;
        float y2 = planes_map_near->at(ids_plane[1]).y;
        float z2 = planes_map_near->at(ids_plane[1]).z;
        float x3 = planes_map_near->at(ids_plane[2]).x;
        float y3 = planes_map_near->at(ids_plane[2]).y;
        float z3 = planes_map_near->at(ids_plane[2]).z;

        // compute the distance to the plane
        float dh = (z0 - z1) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) * 1.0 / sqrt(pow(fabs((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2.0) + pow(fabs((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2.0) + pow(fabs((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2.0)) - (y0 - y1) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) * 1.0 / sqrt(pow(fabs((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2.0) + pow(fabs((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2.0) + pow(fabs((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2.0)) + (x0 - x1) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)) * 1.0 / sqrt(pow(fabs((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2.0) + pow(fabs((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2.0) + pow(fabs((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2.0));

        if (!std::isnan(dh))
        {

            input_point = dense_planes->at(i);
            input_point.intensity = dh;

            point_residuals->push_back(input_point);

            // compute Jacobian
            pcl::PointNormal temp_jacobian = lidar_jacobian.calc_plane_jacobians(input_point.x, input_point.y, input_point.z,
                                                                                  x1, y1, z1, x2, y2, z2, x3, y3, z3, transform);

            jacobians->push_back(temp_jacobian);
        }
    }
    std::cout << "valid_surf number: " << valid_suf << std::endl;
    std::cout << "nonvalid_surf number: " << nonvalid_surf << std::endl;

    matched_plane_num = point_residuals->size() - matched_edge_num;
}

void scan_to_map_matching::process(scan_to_map_matching &lidarOdom)
{
    ros::NodeHandle nh;
    pub_registered_cloud = nh.advertise<PointCloudI>("/velodyne_registered_cloud", 1000);

    edgePlaneJacobian lidar_jacobian;


    ros::Subscriber sub_dense_edges = nh.subscribe("/edge_points", 1000,
                                                   &scan_to_map_matching::DenseEdgePointsCallback, &lidarOdom);

    ros::Subscriber sub_dense_planes = nh.subscribe("/planar_points", 1000,
                                                    &scan_to_map_matching::DensePlanePointsCallback, &lidarOdom);

    ros::Subscriber sub_odometry = nh.subscribe("/wheel_velocity_vehicle_odom",
                                                1000, &scan_to_map_matching::OdometryCallback, &lidarOdom);

    ros::Subscriber sub_dewarped_cloud = nh.subscribe("/velodyne_dewarped_cloud", 1000,
                                                      &scan_to_map_matching::VelodyenDewarpedCloudCallback, &lidarOdom);

    ros::Subscriber sub_relocalize_odom = nh.subscribe("/relocalization_odom", 1000, &scan_to_map_matching::RelocalizeCallback, &lidarOdom);

    ros::Subscriber sub_refined_pose_graph = nh.subscribe("/loop_closure/global_pose_graph", 1000, &scan_to_map_matching::correctMap, &lidarOdom);

    ros::Publisher pub_lidar_odom = nh.advertise<nav_msgs::Odometry>("/lidar_camera_odom", 1000);
    ros::Publisher pub_edges_map_near = nh.advertise<PointCloudI>("/edges_map_near", 1000);
    ros::Publisher pub_planes_map_near = nh.advertise<PointCloudI>("/planes_map_near", 1000);
    ros::Publisher lidar_path_pub = nh.advertise<nav_msgs::Path>("/lidarPath", 1000);
    ros::Publisher feature_points = nh.advertise<PointCloudI>("/feature_points", 1000);
    correct_finish = nh.advertise<std_msgs::String>("/pose_correct", 1000);

    nav_msgs::Odometry::Ptr odom;
    nav_msgs::Odometry lidar_odom;
    lidar_odom.header.frame_id = "world";

    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header.frame_id = "world";
    path.header.frame_id = "world";

    //PointCloudI::Ptr dense_edges(new PointCloudI);
    lidarOdom.initialize();

    // L-M algorithm parameters
    float prev_error = 0.0, curr_error = 0.0;

    //when the robot moves greater than the distance, the map updates
    float map_update_distance = 0.01;

    Vector6f curr_delta;
    Eigen::Vector3f last_update_position;
    last_update_position.x() = 0.0;
    last_update_position.y() = 0.0;
    last_update_position.z() = 0.0;

    ros::Rate rate(100);

    edge_voxel_filter.setLeafSize(map_edge_voxel_size, map_edge_voxel_size, map_edge_voxel_size);
    plane_voxel_filter.setLeafSize(map_plane_voxel_size, map_plane_voxel_size, map_plane_voxel_size);
    full_map_voxel_filter.setLeafSize(full_map_voxel_size, full_map_voxel_size, full_map_voxel_size);

    tf::TransformBroadcaster lidar_br;
    tf::Transform lidar_tf;

    while (ros::ok())
    {

        lidarOdom.update_feature();

        if (flag_dense_edge_received && flag_dense_plane_received && flag_odom_received && flag_dewarped_cloud_received &&
            fabs(odom_timestamp - edge_timestamp) < 0.001 && fabs(odom_timestamp - plane_timestamp) < 0.001 && fabs(odom_timestamp - dewarped_timestamp) < 0.001)
        {
            // if the relocalization message received
            if (false)
            {
                flag_relocalize_received = false;

                // correct current pose
                auto iter = poseMap.find(received_relocalize_timestamp);
                double timeStamp = iter->first;
                // retrieve the previous pose before optimization
                Eigen::Matrix4f loop_prev_pose = iter->second;
                Eigen::Matrix4f refinedPose;
                // update poseMap using pose after optimization
                poseMap[timeStamp] = relocalize_odom;
                while (iter != poseMap.end())
                {
                    timeStamp = iter->first;
                    Eigen::Matrix4f prevPose = iter->second;
                    // get the refinedPose using the corrected pose
                    refinedPose = relocalize_odom * loop_prev_pose.inverse() * prevPose;
                    // update pose map
                    poseMap[timeStamp] = refinedPose;
                    iter++;
                }
                // update transmap
                trans_map = refinedPose;
                continue;
            }
            if (refineMap)
            {
                continue;
            }
            flag = false;
            received_counter++;

            trans_odom_cur = trans_odom_get;

            flag_dense_edge_received = false;
            flag_dense_plane_received = false;
            flag_odom_received = false;
            bool is_degenerate = false;

            // Initialization, add first scan to map
            if ((received_counter == 1 || sqrt(pow(trans_odom_cur(0, 3), 2.0) + pow(trans_odom_cur(1, 3), 2.0) + pow(trans_odom_cur(2, 3), 2.0)) < init_distance) && startingFlag)
            {

                dtrans_odom = trans_odom_pre.inverse() * trans_odom_cur;
                // update the transformation from lidar to world (initial value to be optimized)
                trans_map = trans_map * dtrans_odom;
                // update previous odometry information
                trans_odom_pre = trans_odom_cur;
                // transform edge points from lidar frame to world frame
                pcl::transformPointCloud(*dense_edges, *dense_edges_trans, trans_map);

                // transform plane points from lidar frame to world frame
                pcl::transformPointCloud(*dense_planes, *dense_planes_trans, trans_map);

                lidarOdom.update_cells();

                last_update_position.x() = 0.0;
                last_update_position.y() = 0.0;
                last_update_position.z() = 0.0;
                //*dewarped_cloud = *dense_edges_trans + *dense_planes_trans;

                // publish odometry information using Vio's estimation
                pub_lidar_odom.publish(curOdom);

                // get path
                pose_stamped.header.stamp = curOdom->header.stamp;
                pose_stamped.pose = curOdom->pose.pose;
                path.header.stamp = curOdom->header.stamp;
                path.poses.push_back(pose_stamped);
                lidar_path_pub.publish(path);

                // save current odometry information in the pose map
                poseMap[odom_timestamp] = trans_map;
                std::cout << "odom_timestamp: " << odom_timestamp << std::endl;
            }

            else
            {
                clock_t filter = clock();
                startingFlag = false;
                // calculate relative transformation of current scan from VIN's measurement
                dtrans_odom = trans_odom_pre.inverse() * trans_odom_cur;
                // update the transformation from lidar to world (initial value to be optimized)
                trans_map = trans_map * dtrans_odom;
                // update previous odometry information
                trans_odom_pre = trans_odom_cur;

                // transform the trans_map to 1-6 vector
                lidarOdom.matrixToVec(transform, trans_map);



                // filter the receiving edge feature points
                edge_voxel_filter.setInputCloud(dense_edges);
                edge_voxel_filter.filter(*dense_edges);

                // filter the receiving plane feature points
                plane_voxel_filter.setInputCloud(dense_planes);
                plane_voxel_filter.filter(*dense_planes);


                // transform edge points from lidar frame to world frame
                pcl::transformPointCloud(*dense_edges, *dense_edges_trans, trans_map);

                // transform plane points from lidar frame to world frame
                pcl::transformPointCloud(*dense_planes, *dense_planes_trans, trans_map);

                // renew neighbor-map

                lidarOdom.update_neighbor();
                clock_t filterTime = clock() - filter;
                unsigned msElapsed2 = filterTime / (CLOCKS_PER_SEC / 1000);
                std::cout << "filter feature points time: " << msElapsed2 << " ms\n";

                if (flag)
                {
                    continue;
                }

                clock_t start = clock();

                // LM optimization
                for (int iter = 0; iter < 20; iter++)
                {
                    point_residuals->clear();
                    jacobians->clear();

                    clock_t start1 = clock();

                    // compute residuals and jacobians of error function related to edge feature
                    lidarOdom.compute_edge_jacobian(lidar_jacobian, _correspondence_thresh_);

                    // compute residuals and jacobians of error function related to plane feature
                    lidarOdom.compute_plane_jacobian(lidar_jacobian, _correspondence_thresh_);

                    clock_t jacobian_time = clock() - start1;
                    unsigned msElapsed1 = jacobian_time / (CLOCKS_PER_SEC / 1000);

                    std::cout << "jacobian time: " << msElapsed1 << std::endl;

                    // total matched points number
                    int opt_points_size = point_residuals->size();

                    // for vanilla optimizer
                    pcl::PointNormal temp_jacobian;

                    Eigen::VectorXf d(opt_points_size);
                    Eigen::VectorXf residuals(opt_points_size);
                    Eigen::MatrixXf J(opt_points_size, 6);
                    Eigen::MatrixXf JtJ;
                    Eigen::VectorXf diag;
                    Eigen::DiagonalMatrix<float, Eigen::Dynamic> diagonal;
                    Eigen::MatrixXf temp1;
                    Eigen::MatrixXf temp2;

                    // copy residual before adding weight
                    for (int i = 0; i < opt_points_size; i++)
                    {
                        input_point = point_residuals->points[i];
                        temp_jacobian = jacobians->points[i];
                        float d2 = input_point.intensity;
                        residuals(i) = d2;
                    }

                    // adaptive bisquare weight, not working...
                    int num1 = 0;
                    int num2 = 0;

                    for (int i = 0; i < opt_points_size; i++)
                    {
                        input_point = point_residuals->points[i];
                        temp_jacobian = jacobians->points[i];

                        float d2 = input_point.intensity;

                        // account for bisquare weight, adjust for different dataset

                        float weight;
                        if (fabs(d2) <= bisWeight)
                        {
                            //std::cout << fabs(d2) << std::endl;
                            weight = pow(1 - pow(d2 / bisWeight, 2), 2);
                            num1++;
                        }
                        else
                        {
                            //std::cout << fabs(d2) << std::endl;
                            //std::cout << "weight 0" << std::endl;
                            weight = 0;
                            num2++;
                        }
                        

                        J(i, 0) = weight * temp_jacobian.x;
                        J(i, 1) = weight * temp_jacobian.y;
                        J(i, 2) = weight * temp_jacobian.z;
                        J(i, 3) = weight * temp_jacobian.normal_x;
                        J(i, 4) = weight * temp_jacobian.normal_y;
                        J(i, 5) = weight * temp_jacobian.normal_z;
                        d(i) = weight * d2;
                    }

                    JtJ = J.transpose() * J;
                    diag = JtJ.diagonal();
                    diagonal.diagonal() = diag;
                    temp1 = lambda * diagonal;

                    // Identity Matrix works better
                    Eigen::Matrix<float, 6, 6> temp_I6(Eigen::Matrix<float, 6, 6>::Identity());

                    // temp1 = lambda * temp_I6;
                    temp2 = JtJ + temp1;

                    //LM optimization is of the form
                    // x <- x - (J'J + lambda * diag(J'J))^-1 J'residuals
                    curr_delta = temp2.inverse() * J.transpose() * d;

                    Eigen::Matrix<float, 6, 6> matP;

                    if (iter == 0) {
                        is_degenerate = false;
                        Eigen::Matrix<float, 6, 6> matV;
                        Eigen::Matrix<float, 6, 6> matV2;
                        Eigen::SelfAdjointEigenSolver< Eigen::Matrix<float, 6, 6> > esolver(JtJ);
                        Eigen::Matrix<float, 1, 6> matE = esolver.eigenvalues().real();
                        matV = esolver.eigenvectors().real();
                        matV2 = matV;
                        for (int i = 0; i < 6; i++) {              
                            std::cout << "eigen value: " << matE(0, i) << std::endl;
                            if (matE(0, i) < 10) {
                                is_degenerate = true;
                                for (int j = 0; j < 6; j++) {
                                    matV2(i, j) = 0;
                                }
                                is_degenerate = true;
                            }
                            else {
                                break;
                            }
                        }
                        matP = matV.inverse() * matV2;    
                    }

                    if (is_degenerate) {
                        Vector6f matX2;
                        for (int i = 0; i < 6; i++) {
                            matX2(0, i) = curr_delta(0, i);
                        }
                        curr_delta = matP * matX2;
                    }

                    curr_error = 0;
                    for (int i = 0; i < opt_points_size; i++)
                    {
                        curr_error += pow(residuals(i), 2);
                    }
                    curr_error = sqrt(curr_error);

                    //std::cout << "curr_error:" << curr_error << std::endl;

                    if (iter == 1)
                    {
                        // lambda = 1e-3 * diag.mean();
                        lambda = 0.001;

                        //for the first iteration we just accept the update
                        for (int i = 0; i < 6; i++)
                        {
                            transform[i] -= curr_delta(i);
                        }
                        trans_map = vecToMatrix(transform);
                    }
                    else
                    {
                        if (curr_error > prev_error)
                        {

                            //current error is larger, lets get aggressive and set lambda larger
                            lambda = lambda * 10;

                            //and ignore this update
                        }
                        else
                        {

                            //current error is smaller, lets scale back lambda
                            lambda = lambda / 10;

                            //and accept this udpate
                            for (int i = 0; i < 6; i++)
                            {
                                transform[i] -= curr_delta(i);
                            }

                            // trans_map = util::pose_mat2vec(transform);
                            trans_map = vecToMatrix(transform);
                        }

                        //exit criteria
                        if (fabs(curr_error - prev_error) < min_error)
                        {
                            break;
                        }
                    }

                    //let update the prev pointers
                    prev_error = curr_error;
                    std::cout << "curr_error: " << curr_error / opt_points_size << std::endl;

                    // update the input points with new transform
                    pcl::transformPointCloud(*dense_edges, *dense_edges_trans, trans_map);
                    pcl::transformPointCloud(*dense_planes, *dense_planes_trans, trans_map);
                }

                // compute the distance traveled from last update
                float traveled_distance = sqrt(pow(last_update_position.x() - transform[3], 2) + pow(last_update_position.y() - transform[4], 2) +
                                               pow(last_update_position.z() - transform[5], 2));

                // determine whether updating map
                /*if (abs(last_update_position.z() - transform[5]) > 0.1)
                {
                    std::cout << "drop this frame" << std::endl;
                    continue;
                }*/
                // store current pose in the poseMap
                //std::cout << "odom_timestamp: " << odom_timestamp << std::endl;
                poseMap[odom_timestamp] = trans_map;
                if (traveled_distance > map_update_distance)
                {
                    last_update_position.x() = transform[3];
                    last_update_position.y() = transform[4];
                    last_update_position.z() = transform[5];
                    lidarOdom.update_cells();
                }

                // publish odom
                Eigen::Matrix3f rot_b2o = trans_map.matrix().topLeftCorner<3, 3>(); // body to origin
                Eigen::Quaternionf q_b(rot_b2o);
                auto euler = q_b.toRotationMatrix().eulerAngles(2,1,0);
                std::cout << "euler :" << euler << std::endl;

                std::cout << "map matching pose x: " << transform[3] << ", y: " << transform[4] << ", z:" << transform[5] << std::endl;

                myfile << transform[3] << " " << transform[4] << " " << transform[5] << " " << std::endl;

                lidar_odom.header.stamp.fromSec(odom_timestamp);
                lidar_odom.pose.pose.orientation.x = q_b.x();
                lidar_odom.pose.pose.orientation.y = q_b.y();
                lidar_odom.pose.pose.orientation.z = q_b.z();
                lidar_odom.pose.pose.orientation.w = q_b.w();
                lidar_odom.pose.pose.position.x = transform[3];
                lidar_odom.pose.pose.position.y = transform[4];
                lidar_odom.pose.pose.position.z = transform[5];

                pub_lidar_odom.publish(lidar_odom);

                // get path
                pose_stamped.header.stamp = lidar_odom.header.stamp;
                pose_stamped.pose = lidar_odom.pose.pose;
                path.header.stamp = lidar_odom.header.stamp;
                path.poses.push_back(pose_stamped);
                lidar_path_pub.publish(path);

                // publish registered cloud
                pcl::transformPointCloud(*dewarped_cloud, *dewarped_cloud, trans_map);
                //*dewarped_points_all += *dewarped_cloud;
                *dewarped_points_all += *dewarped_cloud;
                *dewarped_cloud = *dense_edges_trans + *dense_planes_trans;
                //feature_points.publish(dewarped_cloud);

                // save the individual scan
               /*std::string outDir = pointCloudSavePath + std::to_string(seq) + ".pcd";
                pcl::io::savePCDFile(outDir, *dewarped_cloud, true); 
                std::string outDir1 = all_map + std::to_string(seq) + ".pcd";
                full_map_voxel_filter.setInputCloud(dewarped_points_all);
                full_map_voxel_filter.filter(*dewarped_points_all);
                pcl::io::savePCDFile(outDir1, *dewarped_points_all, true);*/
                seq++;
                if (received_counter % save_interval == 0 || correct == true)
                {
                    correct = false;
                    /*full_map_voxel_filter.setInputCloud(dewarped_points_all);
                    full_map_voxel_filter.filter(*dewarped_points_all);
                    pub_registered_cloud.publish(dewarped_points_all);
                    feature_points.publish(dewarped_cloud);*/
                    full_map_voxel_filter.setInputCloud(dewarped_points_all);
                    full_map_voxel_filter.filter(*dewarped_points_all);
                    pcl::io::savePCDFile(PCD_outdir, *dewarped_points_all, true); // Binary format*/
                }
                // broadcast tf
                tf::Quaternion q_tf(q_b.x(), q_b.y(), q_b.z(), q_b.w());
                ros::Time time_tf;
                time_tf.fromSec(odom_timestamp);
                lidar_tf.setOrigin(tf::Vector3(transform[3], transform[4], transform[5]));
                lidar_tf.setRotation(q_tf);
                lidar_br.sendTransform(tf::StampedTransform(lidar_tf, time_tf, "world", "lidar"));
                clock_t timeElapsed = clock() - start;
                unsigned msElapsed = timeElapsed / (CLOCKS_PER_SEC / 1000);
                std::cout << "total time for optimization: " << msElapsed << " ms\n";
            }
        }

        // rate.sleep();
        // ros::spinOnce();
        ros::getGlobalCallbackQueue()->callOne(ros::WallDuration(0));
    }
}

int main(int argc, char **argv)
{
    cout << std::fixed;
    cout.precision(17);
    myfile.open("/home/cong/odom.txt");
    ros::init(argc, argv, "scan_to_map_matching");
    scan_to_map_matching lidarOdom(false, false, false, false);
    lidarOdom.process(lidarOdom);
    return 0;
}
