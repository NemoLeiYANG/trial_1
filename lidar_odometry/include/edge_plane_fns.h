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
#ifndef EDGE_PLANE_FNS_H
#define EDGE_PLANE_FNS_H
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl_ros/transforms.h>
#include <pcl_ros/point_cloud.h>

class edgePlaneJacobian {
public:
	inline pcl::PointNormal calc_edge_jacobians(float x00, float y00, float z00, float x1, float y1,
	 float z1, float x2, float y2, float z2, float x0, float y0, float z0, float *transform) {
	// angle axis
    float rx = transform[0];
    float ry = transform[1];
    float rz = transform[2];

    // det of cross product
    float a012 = sqrt(((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1))
             * ((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
             + ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1))
             * ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) 
             + ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1))
             * ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1)));

    float l12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));

    //compute the derivative of distance error function to x0, y0, and z0
    float d_x0 = ((y1 - y2)*((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
           + (z1 - z2)*((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1))) / a012 / l12;

    float d_y0 = -((x1 - x2)*((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
           - (z1 - z2)*((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1))) / a012 / l12;

    float d_z0 = -((x1 - x2)*((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) 
         + (y1 - y2)*((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1))) / a012 / l12;

    // compute the derivative of x0, y0, z0 to rx, ry, rz, tx, ty, tz, respectively

    float d_x0_rx = -x00*((rx*(ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_x0_ry = -x00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((ry*-2.0)/(rx*rx+ry*ry+rz*rz)+(ry*ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+ry*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((ry*ry)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_x0_rz = -x00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rz*-2.0)/(rx*rx+ry*ry+rz*rz)+(rz*rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(ry*ry)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+z00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((rz*rz)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_rx = -y00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rx*-2.0)/(rx*rx+ry*ry+rz*rz)+(rx*rx*rx)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(rx*rx)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((rx*rx)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_ry = -y00*(((rx*rx)*ry*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+ry*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_rz = -y00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rz*-2.0)/(rx*rx+ry*ry+rz*rz)+(rz*rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(rx*rx)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+z00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+x00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((rz*rz)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_rx = -z00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rx*-2.0)/(rx*rx+ry*ry+rz*rz)+(rx*rx*rx)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(rx*rx)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((rx*rx)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_ry = -z00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((ry*-2.0)/(rx*rx+ry*ry+rz*rz)+(ry*ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(rx*rx)*ry*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+x00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((ry*ry)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_rz = -z00*(((rx*rx)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(ry*ry)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    // chain rules
    float df_drx = d_x0 * d_x0_rx + d_y0 * d_y0_rx + d_z0 * d_z0_rx;

    float df_dry = d_x0 * d_x0_ry + d_y0 * d_y0_ry + d_z0 * d_z0_ry;

    float df_drz = d_x0 * d_x0_rz + d_y0 * d_y0_rz + d_z0 * d_z0_rz;

    float df_dtx = d_x0 * 1.0;

    float df_dty = d_y0 * 1.0;

    float df_dtz = d_z0 * 1.0;
    
    pcl::PointNormal edge_jacobians;
    edge_jacobians.x = df_drx;
    edge_jacobians.y = df_dry;
    edge_jacobians.z = df_drz;
    edge_jacobians.normal_x = df_dtx;
    edge_jacobians.normal_y = df_dty;
    edge_jacobians.normal_z = df_dtz;

    return edge_jacobians;

    }
	inline pcl::PointNormal calc_plane_jacobians(float x00, float y00, float z00, float x1, float y1,
	 float z1, float x2, float y2, float z2, float x3, float y3, float z3, float *transform) {
	float rx = transform[0];
    float ry = transform[1];
    float rz = transform[2];

    float d_x0 = (y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1);
    float d_y0 = (z2 - z1) * (x3 - x1) - (z3 - z1) * (x2 - x1);
    float d_z0 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

    float ps = sqrt(d_x0 * d_x0 + d_y0 * d_y0 + d_z0 * d_z0);
    d_x0 /= ps;
    d_y0 /= ps;
    d_z0 /= ps;


    float d_x0_rx = -x00*((rx*(ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_x0_ry = -x00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((ry*-2.0)/(rx*rx+ry*ry+rz*rz)+(ry*ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+ry*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((ry*ry)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_x0_rz = -x00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rz*-2.0)/(rx*rx+ry*ry+rz*rz)+(rz*rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(ry*ry)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((ry*ry)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+z00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((rz*rz)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_rx = -y00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rx*-2.0)/(rx*rx+ry*ry+rz*rz)+(rx*rx*rx)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(rx*rx)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((rx*rx)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_ry = -y00*(((rx*rx)*ry*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+ry*(rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(ry*ry)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+z00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_y0_rz = -y00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rz*-2.0)/(rx*rx+ry*ry+rz*rz)+(rz*rz*rz)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(rx*rx)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(rz*rz)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+z00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+x00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((rz*rz)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_rx = -z00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((rx*-2.0)/(rx*rx+ry*ry+rz*rz)+(rx*rx*rx)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+rx*(ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+rx*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(rx*rx)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)-(rx*rx)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+((rx*rx)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_ry = -z00*((cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*((ry*-2.0)/(rx*rx+ry*ry+rz*rz)+(ry*ry*ry)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(rx*rx)*ry*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+y00*(-(rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*ry*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*ry*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+(ry*ry)*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+x00*(-sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz)+(ry*ry)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)-((ry*ry)*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+rx*ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*ry*rz*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float d_z0_rz = -z00*(((rx*rx)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0+(ry*ry)*rz*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)+rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*((rx*rx)/(rx*rx+ry*ry+rz*rz)+(ry*ry)/(rx*rx+ry*ry+rz*rz))*1.0/sqrt(rx*rx+ry*ry+rz*rz))+x00*(-(rx*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)-(ry*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)+ry*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+rx*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0)+y00*(-(ry*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0))/(rx*rx+ry*ry+rz*rz)+(rx*rz*cos(sqrt(rx*rx+ry*ry+rz*rz)))/(rx*rx+ry*ry+rz*rz)-rx*rz*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*sin(sqrt(rx*rx+ry*ry+rz*rz))*1.0/pow(rx*rx+ry*ry+rz*rz,3.0/2.0)+ry*(rz*rz)*(cos(sqrt(rx*rx+ry*ry+rz*rz))-1.0)*1.0/pow(rx*rx+ry*ry+rz*rz,2.0)*2.0);

    float df_drx = d_x0 * d_x0_rx + d_y0 * d_y0_rx + d_z0 * d_z0_rx;

    float df_dry = d_x0 * d_x0_ry + d_y0 * d_y0_ry + d_z0 * d_z0_ry;

    float df_drz = d_x0 * d_x0_rz + d_y0 * d_y0_rz + d_z0 * d_z0_rz;

    float df_dtx = d_x0 * 1.0;

    float df_dty = d_y0 * 1.0;

    float df_dtz = d_z0 * 1.0;
    pcl::PointNormal plane_jacobians;
    plane_jacobians.x = df_drx;
    plane_jacobians.y = df_dry;
    plane_jacobians.z = df_drz;
    plane_jacobians.normal_x = df_dtx;
    plane_jacobians.normal_y = df_dty;
    plane_jacobians.normal_z = df_dtz;

    return plane_jacobians;
    }
};
#endif