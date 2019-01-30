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
#include "lidarScanRegistration.h"

// define ros::publisher
ros::Publisher edge_points_pub;
ros::Publisher planar_points_pub;
ros::Publisher dewarped_cloud_pub;
ros::Publisher odom_vehicle_pub;
ros::Publisher sparse_points_pub;
std::string odomTopic;
std::string lidarTopic;

int publishCount = 1;
int Skipcount = 0;
pcl::PointCloud<PointType>::Ptr cornerPointsSharp_total(new PointCloudI);
pcl::PointCloud<PointType>::Ptr surfPointsFlat_total(new PointCloudI);

// construction function
LidarScanRegistration::LidarScanRegistration(int _frameCount, int _initFrameCount,
                                             Eigen::Matrix4f _Identity_matrix, bool _auto_delete,
                                             int _N_SCANS) : frameCount(_frameCount),
                                                             initFrameCount(_initFrameCount),
                                                             odoMsg(_auto_delete),
                                                             N_SCANS(_N_SCANS),
                                                             dtrans(_Identity_matrix),
                                                             prev_close_odometry(_Identity_matrix)
{
  pcl::PointCloud<PointType>::Ptr laserCloud_prev_prev(new pcl::PointCloud<PointType>());
  // load yaml file
  std::string yaml_file = "/home/cong/catkin_ws/src/lidar_odometry/config/indoor.yaml";
  YAML::Node config = YAML::LoadFile(yaml_file);
  lidarTopic = config["lidar_topic"].as<std::string>();
  std::cout << "lidarTopic: " << lidarTopic << std::endl;
  odomTopic = config["odom_topic"].as<std::string>();
  std::cout << "odomTopic: " << odomTopic << std::endl;
  lineSegNum = config["line_seg_num"].as<int>();
  std::cout << "lineSegNum: " << lineSegNum << std::endl;
  maxEdgeNum = config["max_edge_num"].as<int>();
  std::cout << "maxEdgeNum: " << maxEdgeNum << std::endl;
  maxSurfNum = config["max_surf_num"].as<int>();
  std::cout << "maxSurfNum: " << maxSurfNum << std::endl;
  edgeThreshold = config["edge_threshold"].as<double>();
  std::cout << "edgeThreshold: " << edgeThreshold << std::endl;
  scanPeriod = config["scan_period"].as<double>();
  std::cout << "scanPeriod: " << scanPeriod << std::endl;
}

// transform the lidar point to the end of scan
void inline LidarScanRegistration::transform_to_end(PointType pi, PointType &po, Eigen::Matrix4f dtrans)
{
  po.x = dtrans(0, 0) * pi.x + dtrans(0, 1) * pi.y + dtrans(0, 2) * pi.z + dtrans(0, 3);
  po.y = dtrans(1, 0) * pi.x + dtrans(1, 1) * pi.y + dtrans(1, 2) * pi.z + dtrans(1, 3);
  po.z = dtrans(2, 0) * pi.x + dtrans(2, 1) * pi.y + dtrans(2, 2) * pi.z + dtrans(2, 3);
  po.intensity = pi.intensity;
}

// get the line ind of each lidar point for each scan
void inline LidarScanRegistration::getScanInd(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg,
                                              std::vector<pcl::PointCloud<PointType>> &laserCloudScans, int &pointsNum)
{
  pcl::PointCloud<pcl::PointXYZ> curLidarScan;
  // get input laser
  pcl::fromROSMsg(*laserCloudMsg, curLidarScan);
  std::vector<int> index;
  // remove invalid point
  pcl::removeNaNFromPointCloud(curLidarScan, curLidarScan, index);
  pointsNum = curLidarScan.points.size();
  // get the start orientation and end orientation of point cloud
  float startAngle = -atan2(curLidarScan.points[0].y, curLidarScan.points[0].x);
  float endAngle = -atan2(curLidarScan.points[pointsNum - 1].y,
                          curLidarScan.points[pointsNum - 1].x);
  // adjust the pitchAngle
  while (endAngle - startAngle > 3 * M_PI)
  {
    endAngle -= 2 * M_PI;
  }

  while (endAngle - startAngle < M_PI)
  {
    endAngle += 2 * M_PI;
  }

  bool hasScanedFirstPart = false;
  int validPointNum = pointsNum;
  PointType point;
  lineIndArr.clear();
  timePortionOfScan.clear();
  std::vector<int> tmp;
  tmp.resize(0);
  std::vector<float> tmp1;
  tmp1.resize(0);
  lineIndArr.resize(N_SCANS, tmp);
  timePortionOfScan.resize(N_SCANS, tmp1);
  for (unsigned int i = 0; i < pointsNum; i++)
  {
    point.x = curLidarScan.points[i].x;
    point.y = curLidarScan.points[i].y;
    point.z = curLidarScan.points[i].z;

    // calculate the pitch pitchAngle
    float pitchAngle = atan(curLidarScan.points[i].z / sqrt(pow(curLidarScan.points[i].x, 2) + pow(curLidarScan.points[i].y, 2))) * 180 / M_PI;
    int lineInd;
    // get the lineInd of this point
    if (pitchAngle <= 0)
    {
      lineInd = int(pitchAngle - 0.5) + (N_SCANS - 1);
    }
    else
    {
      lineInd = int(pitchAngle + 0.5);
    }
    // if lineInd is out of boundary, discard this point
    if (lineInd < 0 || lineInd >= N_SCANS)
    {
      validPointNum--;
      continue;
    }

    // calculate the yaw pitchAngle to determine the time of this point in the whole scan
    float yawAngle = -atan2(point.y, point.x);
    if (hasScanedFirstPart)
    {
      yawAngle += 2 * M_PI;
    }
    // for the first half scanning points

    while (yawAngle < startAngle)
    {
      yawAngle += 2 * M_PI;
      if (yawAngle - startAngle > M_PI)
      {
        hasScanedFirstPart = true;
      }
    }

    // get the realTime of this point
    float portionOfScan = (yawAngle - startAngle) / (endAngle - startAngle);
    // save the lineInd and realTime of this point in the intensity channel of point(Integer part: lineInd, decimal part: realTime)
    point.intensity = lineInd + scanPeriod * portionOfScan;
    // save this point in the laserCloudScans
    laserCloudScans[lineInd].push_back(point);
    lineIndArr[lineInd].emplace_back(lineInd);
    timePortionOfScan[lineInd].emplace_back(portionOfScan);
  }
  pointsNum = validPointNum; //*/
}

// extract edge points and surface points from point cloud
void inline LidarScanRegistration::get_feature_points(pcl::PointCloud<PointType>::Ptr laserCloud,
                                                      std::vector<int> startPointInd, std::vector<int> endPointInd, std::vector<int> pointInd,
                                                      std::vector<float> pointCurvature, pcl::PointCloud<PointType> &cornerPointsSharp,
                                                      pcl::PointCloud<PointType> &surfPointsFlat,
                                                      std::unordered_set<int> selectedPoint)
{
  // Iterate all scan lines (16 for this project)
  for (int i = 0; i < N_SCANS; i++)
  {
    // segment one scan line to 6 regions, and extract fixed number feature points
    for (int j = 0; j < lineSegNum; j++)
    {
      std::priority_queue<std::pair<float, int>> tmp;
      std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<std::pair<float, int>>> tmp1;
      // get the starting ind and the end ind
      int lowerBound = (startPointInd[i] * (lineSegNum - j) + endPointInd[i] * j) / lineSegNum;
      int upperBound = (startPointInd[i] * (lineSegNum - 1 - j) + endPointInd[i] * (j + 1)) / lineSegNum - 1;

      // add the points to the priority queue
      for (int k = lowerBound; k <= upperBound; k++)
      {
        tmp.push({pointCurvature[pointInd[k]], pointInd[k]});
        tmp1.push({pointCurvature[pointInd[k]], pointInd[k]});
      }

      // extract the edge feature for each line sement, pointcloud is sorted,
      //the last has the largest curvature
      int curEdgeNum = 0, curSurfNum = 0;

      while (!tmp.empty() && curEdgeNum < maxEdgeNum)
      {
        std::pair<float, int> curFeatureInd = tmp.top();
        tmp.pop();
        int ind = curFeatureInd.second;
        // This point shouldn't be picked or close to previous picked points
        // the curvature of edge point should be largen than 0.1
        if (!selectedPoint.count(ind) && pointCurvature[ind] > edgeThreshold)
        {
          cornerPointsSharp.push_back(laserCloud->points[ind]);
          curEdgeNum++;
        }
        else
        {
          continue;
        }
        // do not pick up too close neighbor points
        selectedPoint.insert(ind);
        // check the neighboring points
        float curX = laserCloud->points[ind].x;
        float curY = laserCloud->points[ind].y;
        float curZ = laserCloud->points[ind].z;
        for (int l = -5; l <= 5; l++)
        {
          float diffX = laserCloud->points[ind + l].x - curX;
          float diffY = laserCloud->points[ind + l].y - curY;
          float diffZ = laserCloud->points[ind + l].z - curZ;
          float dis = sqrt(pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2));
          if (dis > 0.25)
          {
            continue;
          }

          selectedPoint.insert(ind + l);
        }
      }

      // extract surface points
      while (!tmp1.empty() && curSurfNum < maxSurfNum)
      {
        std::pair<float, int> curSurfInd = tmp1.top();
        tmp1.pop();
        int ind = curSurfInd.second;
        if (!selectedPoint.count(ind) &&
            pointCurvature[ind] < edgeThreshold)
        {

          surfPointsFlat.push_back(laserCloud->points[ind]);
          curSurfNum++;

          selectedPoint.insert(ind);
          float curX = laserCloud->points[ind].x;
          float curY = laserCloud->points[ind].y;
          float curZ = laserCloud->points[ind].z;
          for (int l = -5; l <= 5; l++)
          {
            float diffX = laserCloud->points[ind + l].x - curX;
            float diffY = laserCloud->points[ind + l].y - curY;
            float diffZ = laserCloud->points[ind + l].z - curZ;
            float dis = sqrt(pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2));
            if (dis > 0.25)
            {
              continue;
            }

            selectedPoint.insert(ind + l);
          }
        }
      }
    }
  }
}

// dewarp featuer points
void inline LidarScanRegistration::get_dewarped_points(pcl::PointCloud<PointType> point_cloud,
                                                       double cur_time, Eigen::Matrix4f cur_close_odometry,
                                                       PointCloudI::Ptr dewarped_points)
{
  for (unsigned int i = 0; i < point_cloud.size(); i++)
  {
    PointType point;
    point = point_cloud.points[i];
    bool cur_odometry_success = false;
    double point_retrieve_time = 0;
    // the lidar points belong to the previous preivous sacn, so we should eliminate the offset(two scan time)
    double point_time = cur_time - 2 * scanPeriod + (point.intensity - int(point.intensity));
    ros::Time cur_point_time = (ros::Time)point_time;
    // retrieve the odometry information
    Eigen::Matrix4f cur_point_odometry =
        odoMsg.GetClosestEntry(cur_point_time, cur_odometry_success, point_retrieve_time);
    // calculate transformation from this point to the end of scan
    dtrans = cur_close_odometry.inverse() * cur_point_odometry;
    // calculate time difference
    double df = point_time - point_retrieve_time;
    // transform this point to the end of scan
    transform_to_end(point, point, dtrans);
    if (cur_odometry_success && fabs(df) < 0.005)
    {
      dewarped_points->push_back(point);
    }
  }
}

// dewarp all lidar points
void inline LidarScanRegistration::get_dewarped_all_points(pcl::PointCloud<PointType>::Ptr point_cloud,
                                                           double cur_time, Eigen::Matrix4f cur_close_odometry,
                                                           PointCloudI::Ptr dewarped_points)
{
  for (unsigned int i = 0; i < point_cloud->size(); i++)
  {
    PointType point;
    point = point_cloud->points[i];
    bool cur_odometry_success = false;
    double point_retrieve_time = 0;
    double point_time = cur_time - 2 * scanPeriod + (point.intensity - int(point.intensity));
    ros::Time cur_point_time = (ros::Time)point_time;
    Eigen::Matrix4f cur_point_odometry =
        odoMsg.GetClosestEntry(cur_point_time, cur_odometry_success, point_retrieve_time);
    dtrans = cur_close_odometry.inverse() * cur_point_odometry;
    double df = point_time - point_retrieve_time;
    transform_to_end(point, point, dtrans);
    if (cur_odometry_success && fabs(df) < 0.005)
    {
      dewarped_points->push_back(point);
    }
  }
}

// process lidar point clouds
// The Laser cloud is directly picked up from the topic as a costant pointer to the laserCloudMsg object 
// Constant pointer prevents the laser cloud handles from changing anything 
void inline LidarScanRegistration::laserCloudHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudMsg)
{
  // all dewarped points
  PointCloudI::Ptr dewarped_cloud(new PointCloudI);
  // dewarped edge points
  PointCloudI::Ptr dewarped_edge_points(new PointCloudI);
  // dewarped surface points
  PointCloudI::Ptr dewarped_flat_points(new PointCloudI);
  // all feature points
  PointCloudI::Ptr sparse_points(new PointCloudI);

  std::vector<int> startPointInd(N_SCANS, 0);
  std::vector<int> endPointInd(N_SCANS, 0);

  std::vector<pcl::PointCloud<PointType>> laserCloudScans(N_SCANS);
  pointCurvature.clear();
  pointCurvature.resize(50000, 0);
  selectedPoint.clear();
  pointInd.clear();
  pointInd.resize(50000, 0);

  int pointsNum = 0;

  // assigen each lidar point to each sacn line
  getScanInd(laserCloudMsg, laserCloudScans, pointsNum);

  double timeScanCur = laserCloudMsg->header.stamp.toSec();

  pcl::PointCloud<PointType>::Ptr laserCloud(new pcl::PointCloud<PointType>());
  for (int i = 0; i < N_SCANS; i++)
  {
    *laserCloud += laserCloudScans[i];
  }
  int scanCount = -1;
  // calculate the curvature of each lidar point
  for (int i = 5; i < pointsNum - 5; i++)
  {
    float diffX = -10 * laserCloud->points[i].x;
    float diffY = -10 * laserCloud->points[i].y;
    float diffZ = -10 * laserCloud->points[i].z;
    // compare with neighboring 10 lidar point of the same scan line (2d feature!)
    for (int j = -5; j <= 5; j++)
    {
      if (j == 0)
      {
        continue;
      }
      diffX += laserCloud->points[i + j].x;
      diffY += laserCloud->points[i + j].y;
      diffZ += laserCloud->points[i + j].z;
    }
    // assign initial value
    pointCurvature[i] = (pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2));
    // pointInd is used for sorting point cloud of the same scan line according to the curvature
    pointInd[i] = i;
    // intensity stores the line number and realtime of the point
    // Intensity = (integer.decimal), integer part stores the lineInd information
    // Decimal part stores the real time information(0 - 1 of the whole scan)
    if (int(laserCloud->points[i].intensity) != scanCount)
    {
      scanCount = int(laserCloud->points[i].intensity);
      // get the starInd and EndInd of each scan line
      if (scanCount > 0 && scanCount < N_SCANS)
      {
        startPointInd[scanCount] = i + 5;
        endPointInd[scanCount - 1] = i - 5;
      }
    }
  }
  startPointInd[0] = 5;
  endPointInd.back() = pointsNum - 5;
  /* remove points on local planar surfaces that are roughly parallel
   to the laser beams (case 1)
  avoid points that are on boundary of occlude regions (case2)
  If the point is for these cases, we set the cloudNeighborPicked of this point as 1*/
  for (int i = pointsNum - 5; i > 5; i--)
  {
    float diffX = laserCloud->points[i].x - laserCloud->points[i - 1].x;
    float diffY = laserCloud->points[i].y - laserCloud->points[i - 1].y;
    float diffZ = laserCloud->points[i].z - laserCloud->points[i - 1].z;
    float diff = pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2);
    float point1X = laserCloud->points[i - 1].x;
    float point1Y = laserCloud->points[i - 1].y;
    float point1Z = laserCloud->points[i - 1].z;
    float point2X = laserCloud->points[i].x;
    float point2Y = laserCloud->points[i].y;
    float point2Z = laserCloud->points[i].z;
    float depth1 = sqrt(pow(point1X, 2) + pow(point1Y, 2) + pow(point1Z, 2));
    float depth2 = sqrt(pow(point2X, 2) + pow(point2Y, 2) + pow(point2Z, 2));
    // deal with case1
    if (diff > 0.1)
    {

      if (depth1 > depth2)
      {
        diffX = point2X - point1X * depth2 / depth1;
        diffY = point2Y - point1Y * depth2 / depth1;
        diffZ = point2Z - point1Z * depth2 / depth1;

        if (sqrt(pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2)) / depth2 < 0.1)
        {
          for (int j = -5; j <= 0; j++)
          {
            selectedPoint.insert(i + j);
          }
        }
      }
      else
      {
        diffX = point2X * depth1 / depth2 - point1X;
        diffY = point2Y * depth1 / depth2 - point1Y;
        diffZ = point2Z * depth1 / depth2 - point1Z;

        if (sqrt(pow(diffX, 2) + pow(diffY, 2) + pow(diffZ, 2)) / depth1 < 0.1)
        {
          for (int j = 1; j <= 6; j++)
          {
            selectedPoint.insert(i + j);
          }
        }
      }
    }
    // deal with case2

    float diffX2 = point2X - point1X;
    float diffY2 = point2Y - point1Y;
    float diffZ2 = point2Z - point1Z;
    float diff2 = pow(diffX2, 2) + pow(diffY, 2) + pow(diffZ, 2);
    if (diff > 0.0002 * depth2 && diff2 > 0.0002 * depth2)
    {
      selectedPoint.insert(i);
    }
  }

  pcl::PointCloud<PointType> cornerPointsSharp;
  pcl::PointCloud<PointType> surfPointsFlat;

  get_feature_points(laserCloud, startPointInd, endPointInd, pointInd, pointCurvature,
                     cornerPointsSharp,
                     surfPointsFlat, selectedPoint);

  bool odometry_success = false;
  double cur_time = laserCloudMsg->header.stamp.toSec() - scanPeriod;
  ros::Time cur_laser_time = (ros::Time)cur_time;
  double retrived_time = 0;
  Eigen::Matrix4f cur_close_odometry =
      odoMsg.GetClosestEntry(cur_laser_time, odometry_success, retrived_time);

  double df1 = cur_time - retrived_time;

  //dewarp the laser cloud
  Eigen::Matrix4f dtrans(Eigen::Matrix4f::Identity());
  if (odometry_success && fabs(df1) < 0.003)
  {
    start = clock();

    // dewarping edge points
    get_dewarped_points(cornerPointsSharp_prev_prev, cur_time, cur_close_odometry, dewarped_edge_points);

    // dewarping surface points
    get_dewarped_points(surfPointsFlat_prev_prev, cur_time, cur_close_odometry, dewarped_flat_points);
    end = clock();
    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    //std::cout << "cpu_time_used: " << cpu_time_used << std::endl;

    // update the timestamp
    dewarped_edge_points->header.stamp = (cur_time)*1000000;
    dewarped_flat_points->header.stamp = (cur_time)*1000000;
    dewarped_cloud->header.stamp = (cur_time)*1000000;

    dewarped_edge_points->header.frame_id = "world";
    dewarped_flat_points->header.frame_id = "world";
    dewarped_cloud->header.frame_id = "world";

    *cornerPointsSharp_total += *dewarped_edge_points;
    *surfPointsFlat_total += *dewarped_flat_points;
    cornerPointsSharp_total->header.stamp = dewarped_edge_points->header.stamp;
    surfPointsFlat_total->header.stamp = dewarped_edge_points->header.stamp;
    cornerPointsSharp_total->header.frame_id = "world";
    surfPointsFlat_total->header.frame_id = "world";

    Eigen::Matrix3f mat = cur_close_odometry.matrix().topLeftCorner<3, 3>();
    Eigen::Quaternionf q_b(mat);
    nav_msgs::Odometry velocity_vehicle_odom;

    //velocity_vehicle_odom.header.stamp = (ros::Time) time_count;
    velocity_vehicle_odom.header.stamp = (ros::Time)(cur_time);
    velocity_vehicle_odom.header.frame_id = "/vehicle_init";
    velocity_vehicle_odom.child_frame_id = "/vehicle";
    velocity_vehicle_odom.pose.pose.orientation.x = q_b.x();
    velocity_vehicle_odom.pose.pose.orientation.y = q_b.y();
    velocity_vehicle_odom.pose.pose.orientation.z = q_b.z();
    velocity_vehicle_odom.pose.pose.orientation.w = q_b.w();
    velocity_vehicle_odom.pose.pose.position.x = cur_close_odometry(0, 3);
    velocity_vehicle_odom.pose.pose.position.y = cur_close_odometry(1, 3);
    velocity_vehicle_odom.pose.pose.position.z = cur_close_odometry(2, 3);

    prev_close_odometry = cur_close_odometry;

    time_count++;
    Skipcount++;
    //std::cout << "dewarped_cloud_size:" << dewarped_cloud->size() << std::endl;
    //std::cout << "dewarped_flat_points size:" << dewarped_flat_points->size() << std::endl;
    //std::cout << "dewarped_edge_points size:" << dewarped_edge_points->size() << std::endl;
    //
    if (time_count > initFrameCount && Skipcount % publishCount == 0)
    {
      // filter the map points
      //pcl::VoxelGrid<PointType> mapFilter;
      //mapFilter.setInputCloud(laserCloud_prev_prev);
      //mapFilter.setLeafSize(0.01, 0.01, 0.01);
      //mapFilter.filter(*laserCloud_prev_prev);
      // dewarp all points
      get_dewarped_points(*laserCloud_prev_prev, cur_time, cur_close_odometry, dewarped_cloud);

      std::cout << "edge size " << dewarped_edge_points->size() << std::endl;
      std::cout << "plan size: " << dewarped_flat_points->size() << std::endl;

      edge_points_pub.publish(cornerPointsSharp_total);
      std::cout << "corner size: " << cornerPointsSharp_total->size() << std::endl;
      planar_points_pub.publish(surfPointsFlat_total);
      std::cout << "surface size: " << surfPointsFlat_total->size() << std::endl;
      dewarped_cloud_pub.publish(dewarped_cloud);
      *sparse_points = *dewarped_edge_points + *dewarped_flat_points;
      sparse_points_pub.publish(sparse_points);
      odom_vehicle_pub.publish(velocity_vehicle_odom);
      cornerPointsSharp_total->clear();
      surfPointsFlat_total->clear();
    }
  }

  //save the feature point
  cornerPointsSharp_prev_prev = cornerPointsSharp_prev;
  surfPointsFlat_prev_prev = surfPointsFlat_prev;

  cornerPointsSharp_prev = cornerPointsSharp;
  surfPointsFlat_prev = surfPointsFlat;

  laserCloud_prev_prev = laserCloud_prev;
  laserCloud_prev = laserCloud;
  cur_time = laserCloudMsg->header.stamp.toSec() - 0.202;
  cur_laser_time = (ros::Time)cur_time;
  odoMsg.DeletePrevious(cur_laser_time);
}

void LidarScanRegistration::OdometryCallback(const nav_msgs::Odometry::ConstPtr &odom)
{
  // Because we use the odometry information to dewarp the previous previous point clouds, so we should wait at least 2 lidar scans
  frameCount++;
  if (frameCount < initFrameCount)
  {
    return;
  }

  Eigen::Matrix4f trans_odom_cur(Eigen::Matrix4f::Identity());
  Eigen::Quaternionf q_b;
  Eigen::Matrix3f rot_b2o;       // body to origin
  Eigen::Matrix3f rot_lidar2imu; // lidar to imu
  Eigen::Quaternionf q_lidar2imu;
  Eigen::Matrix4f lidar2imu(Eigen::Matrix4f::Identity());
  Eigen::Matrix4f imu2lidar(Eigen::Matrix4f::Identity());
  //extrinsic, first platform
  /*q_lidar2imu.x() = 0;
  q_lidar2imu.y() = 0;
  q_lidar2imu.z() = 0.70711;
  q_lidar2imu.w() = 0.70711;
  rot_lidar2imu = q_lidar2imu.toRotationMatrix();
  lidar2imu.matrix().topLeftCorner<3,3>() = rot_lidar2imu;
  lidar2imu(0,3) = -0.02803;
  lidar2imu(1,3) = 0.03496;
  lidar2imu(2,3) = -0.087869;*/
  //extrinsic
  q_lidar2imu.x() = 1; //0;
  q_lidar2imu.y() = 0; //0;
  q_lidar2imu.z() = 0; //0.70711;
  q_lidar2imu.w() = 0; //0.70711;
  rot_lidar2imu = q_lidar2imu.toRotationMatrix();
  lidar2imu.matrix().topLeftCorner<3, 3>() = rot_lidar2imu;
  lidar2imu(0, 3) = -0.011356;
  lidar2imu(1, 3) = -0.002352;
  lidar2imu(2, 3) = -0.08105;
  imu2lidar = lidar2imu.inverse();

  q_b.x() = odom->pose.pose.orientation.x;
  q_b.y() = odom->pose.pose.orientation.y;
  q_b.z() = odom->pose.pose.orientation.z;
  q_b.w() = odom->pose.pose.orientation.w;

  rot_b2o = q_b.toRotationMatrix();
  trans_odom_cur.matrix().topLeftCorner<3, 3>() = rot_b2o;
  trans_odom_cur(0, 3) = odom->pose.pose.position.x;
  trans_odom_cur(1, 3) = odom->pose.pose.position.y;
  trans_odom_cur(2, 3) = odom->pose.pose.position.z;
  // transform the pose in lidar coordinate system to lidar coordinate system
  trans_odom_cur = trans_odom_cur * imu2lidar;
  odoMsg.AddEntry(odom->header.stamp, trans_odom_cur);
}

int main(int argc, char **argv)
{

  ros::init(argc, argv, "scanRegistration");

  ros::NodeHandle nh;

  // lidarScanRegistration(frame_count, init_frame count, identity matrix, bool message, total frames)
  // Does this mean it takes only 16 frames at a time?
  LidarScanRegistration lidarScanRegistration(0, 2,
                                              Eigen::Matrix4f::Identity(), false, 16);
  // Takes into two objects as reference with callback functions 
  // Things like Lidar topic name are loaded from a YAML file 

  // For subscribers, its an option for you to mention the type name of messages or not 
  ros::Subscriber subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2>(lidarTopic, 100, &LidarScanRegistration::laserCloudHandler, &lidarScanRegistration);
  
  // Define custom msg name in order to prevent random message types                                                                       
  ros::Subscriber subOdometry = nh.subscribe(odomTopic, 1000, &LidarScanRegistration::OdometryCallback, &lidarScanRegistration);
  // Advertise all these messages 
  edge_points_pub = nh.advertise<PointCloudI>("/edge_points", 100);
  planar_points_pub = nh.advertise<PointCloudI>("/planar_points", 100);
  dewarped_cloud_pub = nh.advertise<PointCloudI>("/velodyne_dewarped_cloud", 100);
  odom_vehicle_pub = nh.advertise<nav_msgs::Odometry>("wheel_velocity_vehicle_odom", 100);
  sparse_points_pub = nh.advertise<PointCloudI>("/sparse_points", 100);
  ros::spin();
  return 0;
}
