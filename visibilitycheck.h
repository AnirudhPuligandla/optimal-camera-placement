#ifndef VISIBILITYCHECK_H
#define VISIBILITYCHECK_H

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 200

#include <cl2.hpp>
#include <memory>
#include <algorithm>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <QVector>
#include <QVector3D>
#include <QDebug>
#include <QFile>
#include <QDataStream>

class VisibilityCheck
{
private:
    // Structure to define the view frustum
    struct Frustum
    {
        /* Far plane
         * Although 4 points can be used, we use only 3 because three points sufficient to define the plane
         */
        QVector<QVector3D> far;
        // Top, bottom, left and right planes
        QVector<QVector3D> top, bottom, left, right;
    };

    // Structure for cameras
    struct cameraSpecs
    {
        // horizontal and vertical FoV angles
        float hFOV, vFOV;
        // depth of field
        float z_f;
        // horizontal and vertical half FoV values in voxels
        int hFar, vFar;
    };

    // private data members to setup input data
    QVector<QVector3D> camPos, backgroundVox;
    QVector<QVector<QVector3D>> camDir, camUp, camRight;
    // OpenCL program
    cl::Program camProgram;
    // one plus 6 camera specs
    cameraSpecs onePlus6;

public:
    // General Constructor
    VisibilityCheck(QVector<QVector3D> &boundaryVoxels, QVector<QVector3D> &boundaryNormals, QVector<QVector3D> &backgroundVoxels);
    /* Function to initialize GPU and associated parameters
     * Also initialize jernel strings
     */
    void setupGPU();
    /* Function to calculate camera up and right vectors
     * this function will also calculate different orientations for each camera position
     */
    void setupCameras();
    /* Funtion to setup the g matrix and write into a file
     */
    void setup_g();
};

#endif // VISIBILITYCHECK_H
