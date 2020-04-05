#include "visibilitycheck.h"

VisibilityCheck::VisibilityCheck(QVector<QVector3D> &boundaryVoxels, QVector<QVector3D> &boundaryNormals, QVector<QVector3D> &backgroundVoxels)
{
    // Initialize input data
    camPos.clear();
    camDir.clear();
    backgroundVox.clear();
    camUp.clear();
    camRight.clear();

    int numCamPos = boundaryVoxels.size();
    int numRotations = 9;
    camDir.resize(numRotations);
    camUp.resize(numRotations);
    camRight.resize(numRotations);
    for(int i = 0; i < numRotations; i++)
    {
        camDir[i].resize(numCamPos);
        camUp[i].resize(numCamPos);
        camRight[i].resize(numCamPos);
    }

    camPos = boundaryVoxels;
    camDir[0] = boundaryNormals;
    backgroundVox = backgroundVoxels;

    /* Define camera specifications considering One plus 6 camera
     * the horizontal field of view angle is 63.1 deg
     * aspect ratio is 16:9
     * depth of field is decided according to the 12m radius (approx. 190 voxels)
     * based on calculations hFar*2 estimates at approx. 223 voxels rounded to 224 (multiple of 16)
     * from aspect ratio it can be said that vertical FoV angle is roughly 35.5 deg
     * and vFar*2 will be 126 voxels
     */
    onePlus6.hFOV = 63.1f;
    onePlus6.vFOV = 35.5f;
    onePlus6.z_f = 190;
    onePlus6.hFar = 112;
    onePlus6.vFar = 63;

    // Initialize the kernel strings and other GPU parameters
    setupGPU();

    // generate rotations for all camera positions
    setupCameras();

    // call the funtion to setup the g matrix and write it to file
    setup_g();

    qDebug() << "Visibility checks completed!";
}

void VisibilityCheck::setupGPU()
{
    /* GPU device selection
     *
     * Filter for a 2.0 platform and set it as default
     */
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    cl::Platform plat;
    //qDebug()<<"No. of platforms: "<<platforms.size()<<"\n";
    for (auto &p : platforms)
    {
        std::string platver = p.getInfo<CL_PLATFORM_VERSION>();
        //std::string platvend = p.getInfo<CL_PLATFORM_VENDOR>();
        std::cout << platver << "\n";
        //std::cout << platvend << "\n";
        if(platver.find("OpenCL 2.") != std::string::npos)
            plat = p;
    }
    if(plat() == 0)
    {
        qDebug()<< "No OpenCL 2.0 platform found.\n";
    }

    cl::Platform newP = cl::Platform::setDefault(plat);
    if(newP != plat)
    {
        qDebug()<<"Error setting default platform.\n";
    }

    /* Create a openCL kernel for calculating camera up and right vectors from camDir (view direction)
     *
     * involves the following steps:
     * 1) Assume x-z plane faces in +ve y-direction
     * 2) get direction of line intersecting between surface plane and x-z plane
     *      -- given by camDir x (0,1,0)
     * 3) get camera up vector
     *      -- given by camRight x camUp
     *
     * Inputs : camera view direction (surface normal) at voxel 'camDirection'
     * Outputs : camera up and right direction vectors 'camUp', 'camRight'
     */
    std::string cameraKernel{
        "kernel void defineCamera(global const float3 *camDirection, global float3 *camUp, global float3 *camRight){"
        "   int iterId = get_global_id(0);"
        "   float3 xzNormal = (float3)(0.0f, 1.0f, 0.0f);"
        "   camRight[iterId] = normalize(cross(camDirection[iterId], xzNormal));"
        "   camUp[iterId] = normalize(cross(camRight[iterId], camDirection[iterId]));"
    "}"};

    /* OpenCL kernel to calculate orientation vector after rotation about an 'axis' by an 'angle'
     *
     * steps involved:
     * 1) construct the rotation matrix -> get the matrix multiplication output with the 'axis' vector
     *
     * Inputs: primary orientation 'camDirection'
     *         axis of rotation 'axis'
     *         'angle' of rotation in radians
     * Outputs: direction vector after rotation 'rotatedDirection'
     */
    std::string camPoseKernel{
        "kernel void calcOrientations(global const float3 *camDirection, global const float3 *Axis, float angle, global float3 *rotatedDirection){"
        "   int iterId = get_global_id(0);"
        "   float3 axis = Axis[iterId];"
        "   float r11 = cos(angle) + pown(axis.s0, 2) * (1 - cos(angle));"
        "   float r12 = axis.s0 * axis.s1 * (1 - cos(angle)) - axis.s2 * sin(angle);"
        "   float r13 = axis.s0 * axis.s2 * (1 - cos(angle)) + axis.s1 * sin(angle);"
        "   float r21 = axis.s1 * axis.s0 * (1 - cos(angle)) + axis.s2 * sin(angle);"
        "   float r22 = cos(angle) + pown(axis.s1, 2) * (1 - cos(angle));"
        "   float r23 = axis.s1 * axis.s2 * (1 - cos(angle)) - axis.s0 * sin(angle);"
        "   float r31 = axis.s2 * axis.s0 * (1 - cos(angle)) - axis.s1 * sin(angle);"
        "   float r32 = axis.s2 * axis.s1 * (1 - cos(angle)) + axis.s0 * sin(angle);"
        "   float r33 = cos(angle) + pown(axis.s2, 2) * (1 - cos(angle));"
        "   rotatedDirection[iterId].s0 = r11 * camDirection[iterId].s0 + r12 * camDirection[iterId].s1 + r13 * camDirection[iterId].s2;"
        "   rotatedDirection[iterId].s1 = r21 * camDirection[iterId].s0 + r22 * camDirection[iterId].s1 + r23 * camDirection[iterId].s2;"
        "   rotatedDirection[iterId].s2 = r31 * camDirection[iterId].s0 + r32 * camDirection[iterId].s1 + r33 * camDirection[iterId].s2;"
    "}"};

    // Kernel to calculate frustum at each camera location
    std::string frustumKernel{
        "kernel void frustumCalc(global const float3 *camPoints, global const float3 *camDirection, global const float3 *camUp, global const float3 *camRight, global float16 *camFrustum, global float3 *farPlanePoint, float3 camParams){"
        "  int iterId = get_global_id(0);"
        // Field of view parameter calculation
        "  float3 farCenter = camPoints[iterId] + camDirection[iterId]*camParams.z;"
        "  float3 far_tl = farCenter + camUp[iterId]*camParams.y - camRight[iterId]*camParams.x;"
        "  float3 far_bl = far_tl - camUp[iterId]*2*camParams.y;"
        "  float3 far_br = far_bl + camRight[iterId]*2*camParams.x;"
        "  float3 far_tr = far_br + camUp[iterId]*2*camParams.y;"
        "  float3 farNormal = normalize(cross(far_bl - far_tl, far_tr - far_tl));"
        "  float3 leftNormal = normalize(cross(far_bl - camPoints[iterId], far_tl - camPoints[iterId]));"
        "  float3 bottomNormal = normalize(cross(far_br - camPoints[iterId], far_bl - camPoints[iterId]));"
        "  float3 rightNormal = normalize(cross(far_tr - camPoints[iterId], far_br - camPoints[iterId]));"
        "  float3 topNormal = normalize(cross(far_tl - camPoints[iterId], far_tr - camPoints[iterId]));"
        "  camFrustum[iterId] = (float16)(farNormal, leftNormal, bottomNormal, rightNormal, topNormal, 0.0f);"
        "  farPlanePoint[iterId] = far_tl;"
        //"  printf(\"Processed id %d \", iterId);"
        //"printf(\" -----%d------\", count);"
    "}"};

    // Kernel to check control points for each camera location
    std::string pointCheckKernel{
        "kernel void checkControlPoints(global const float3 *contPoints, const float16 camFrustum, float3 camPoint, float3 far_tl, global int *outMatG){"
        "  int iterId = get_global_id(0);"
        "  float3 farNormal = (float3)(camFrustum.s0, camFrustum.s1, camFrustum.s2);"
        "  float3 leftNormal = (float3)(camFrustum.s3, camFrustum.s4, camFrustum.s5);"
        "  float3 bottomNormal = (float3)(camFrustum.s6, camFrustum.s7, camFrustum.s8);"
        "  float3 rightNormal = (float3)(camFrustum.s9, camFrustum.sA, camFrustum.sB);"
        "  float3 topNormal = (float3)(camFrustum.sC, camFrustum.sD, camFrustum.sE);"
        "  int check = 0;"
        "  if(dot(contPoints[iterId] - far_tl, farNormal) >= 0)"
        "    if(dot(contPoints[iterId] - camPoint, leftNormal) >= 0)"
        "      if(dot(contPoints[iterId] - camPoint, bottomNormal) >= 0)"
        "        if(dot(contPoints[iterId] - camPoint, rightNormal) >= 0)"
        "          if(dot(contPoints[iterId] - camPoint, topNormal) >= 0)"
        "            check = 1;"
        "  outMatG[iterId] = check;"
    "}"};

    // Newer Opencl simpler string interface
    std::vector<std::string> programStrings{cameraKernel, camPoseKernel, frustumKernel, pointCheckKernel};

    camProgram = cl::Program(programStrings);

#if defined (CL_HPP_ENABLE_EXCEPTIONS)
    try{
        camProgram.build("-cl-std=CL2.0");
    }
    catch(...){
        // Print build info for all devices
        cl_int buildErr = CL_SUCCESS;
        auto buildInfo = camProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(&buildErr);
        for(auto &pair : buildInfo)
        {
            std::cerr << pair.second << std::endl << std::endl;
        }
    }
#else
    cl_int buildErr = camProgram.build("-cl-std=CL2.0");
    if(buildErr != CL_SUCCESS)
    {
        std::cerr << "Build error: "<< buildErr << "\n";
    }
#endif
}

void VisibilityCheck::setupCameras()
{
    // create SVM data structure type
    cl::SVMAllocator<int, cl::SVMTraitCoarse<>> svmAlloc;
    qDebug() << "Max alloc size: " << svmAlloc.max_size() << " bytes\n";

    // SVM allocations for passing into kernels
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camDirection(svmAlloc);         // input svm buffer containing camera view direction vectors for all orientations
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> rotationAxisUp(svmAlloc);         // input svm buffer containing the up rotation axes for corresponding position
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> rotationAxisRight(svmAlloc);         // input svm buffer containing the right rotation axes for corresponding position
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camDirectionRotated(svmAlloc);  // output svm buffer to hold camera view direction after rotation

    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camUpVectors(svmAlloc);         // output svm buffers to hold camera up vectors for all view directions
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camRightVectors(svmAlloc);      // output svm buffer to hold camera right vectors

    // initialize kernels from strings
    auto cameraKernel =
            cl::KernelFunctor<
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&
            >(camProgram, "defineCamera");

    auto camPoseKernel =
            cl::KernelFunctor<
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                cl_float,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&
            >(camProgram, "calcOrientations");

    // get the angle step for horizontal and vertical rotations
    float angleH = ((90 - onePlus6.hFOV/2)/2)*static_cast<float>((M_PI/180));
    float angleV = ((90 - onePlus6.vFOV/2)/2)*static_cast<float>((M_PI/180));

    /* Create a list to change the angle of rotation in steps of angleH and angleV
     * We use steps of 4 in each horizaontal and vertical directions to obtain 9 rotations in total for each camera position
     * ***NOTE***
     * When changing number of rotations here, change it also in the constructor
     * ***NOTE***
     */
    QVector<int> stepList = {-2,-1,1,2};

    //----------First compute the up and right vectors for the primary view direction--------------------------------

    // initialize the svm buffers
    int numCamPos = camPos.size();
    camUpVectors.clear();
    camRightVectors.clear();
    camDirection.resize(numCamPos);
    camUpVectors.resize(numCamPos);
    camRightVectors.resize(numCamPos);
    // initialize primary orientation data
    for(int i = 0; i < numCamPos; i++)
    {
        camDirection[i].x = camDir[0][i].x();
        camDirection[i].y = camDir[0][i].y();
        camDirection[i].z = camDir[0][i].z();
    }

    // Unmap the input vector
    cl::unmapSVM(camDirection);
    cl::unmapSVM(camUpVectors);
    cl::unmapSVM(camRightVectors);

    // Call the function camera definition kernel
    cl_int error;
    cameraKernel(cl::EnqueueArgs(cl::NDRange(numCamPos)), camDirection, camUpVectors, camRightVectors, error);

    // Grab the output SVM vectors
    cl::mapSVM(camUpVectors);
    cl::mapSVM(camRightVectors);

    // Save the outputs into global containers
    for(int i = 0; i < numCamPos; i++)
    {
        camUp[0][i].setX(camUpVectors[i].x);
        camUp[0][i].setY(camUpVectors[i].y);
        camUp[0][i].setZ(camUpVectors[i].z);

        camRight[0][i].setX(camRightVectors[i].x);
        camRight[0][i].setY(camRightVectors[i].y);
        camRight[0][i].setZ(camRightVectors[i].z);
    }

    //-------------------------Now calculate the 8 rotations and corresponding up and right vectors---------------------

    // Counter to keep track of global rotations
    int counter = 1;
    cl_float angle;

    // intitalize the rotation axes with the up and right vectors of the primary orientation
    rotationAxisUp.clear();
    rotationAxisRight.clear();
    rotationAxisUp.resize(numCamPos);
    rotationAxisRight.resize(numCamPos);
    for(int i = 0; i < numCamPos; i++)
    {
        rotationAxisUp[i].x = camUp[0][i].x();
        rotationAxisUp[i].y = camUp[0][i].y();
        rotationAxisUp[i].z = camUp[0][i].z();

        rotationAxisRight[i].x = camRight[0][i].x();
        rotationAxisRight[i].y = camRight[0][i].y();
        rotationAxisRight[i].z = camRight[0][i].z();
    }

    // Iterate over the number of steps
    foreach (int val, stepList)
    {
        //-------------Calculate rotaion in horizontal direction-----------------------------
        // Compute the angle for horizontal rotations
        angle = val * angleH;
        // initialize svm containers
        camDirectionRotated.clear();
        camDirectionRotated.resize(numCamPos);
        // unmap input vectors for direction calculating kernel
        cl::unmapSVM(camDirection);
        cl::unmapSVM(rotationAxisUp);
        cl::unmapSVM(camDirectionRotated);
        // Call the kernel function to calculate rotated directions
        cl_int error;
        camPoseKernel(cl::EnqueueArgs(cl::NDRange(numCamPos)), camDirection, rotationAxisUp, angle, camDirectionRotated, error);
        // grab the output rotated vectors
        cl::mapSVM(camDirectionRotated);
        // Save into global containers
        for(int i = 0; i < numCamPos; i++)
        {
            camDir[counter][i].setX(camDirectionRotated[i].x);
            camDir[counter][i].setY(camDirectionRotated[i].y);
            camDir[counter][i].setZ(camDirectionRotated[i].z);
        }

        // ----------------Calculate the up and right vectors for new rotation-------------------
        // re-initialize the output vectors
        camUpVectors.clear();
        camRightVectors.clear();
        camUpVectors.resize(numCamPos);
        camRightVectors.resize(numCamPos);
        // unmap the svm vectors
        cl::unmapSVM(camDirectionRotated);
        cl::unmapSVM(camUpVectors);
        cl::unmapSVM(camRightVectors);
        // Call the function camera definition kernel
        cameraKernel(cl::EnqueueArgs(cl::NDRange(numCamPos)), camDirectionRotated, camUpVectors, camRightVectors, error);
        // grab the output vectors
        cl::mapSVM(camUpVectors);
        cl::mapSVM(camRightVectors);
        // save into global containers
        for(int i = 0; i < numCamPos; i++)
        {
            camUp[counter][i].setX(camUpVectors[i].x);
            camUp[counter][i].setY(camUpVectors[i].y);
            camUp[counter][i].setZ(camUpVectors[i].z);

            camRight[counter][i].setX(camRightVectors[i].x);
            camRight[counter][i].setY(camRightVectors[i].y);
            camRight[counter][i].setZ(camRightVectors[i].z);
        }

        // increment the counter
        counter++;

        //---------------Calculate rotation in vertical direction-----------------------------
        // compute angle for vertical rotations
        angle = val * angleV;
        // Initialize svm containers
        camDirectionRotated.clear();
        camDirectionRotated.resize(numCamPos);
        // unmap input vectors for direction calculating kernel
        cl::unmapSVM(camDirection);
        cl::unmapSVM(rotationAxisRight);
        cl::unmapSVM(camDirectionRotated);
        // Call the kernel function to calculate rotated directions
        camPoseKernel(cl::EnqueueArgs(cl::NDRange(numCamPos)), camDirection, rotationAxisRight, angle, camDirectionRotated, error);
        // grab the output rotated vectors
        cl::mapSVM(camDirectionRotated);
        // Save into global containers
        for(int i = 0; i < numCamPos; i++)
        {
            camDir[counter][i].setX(camDirectionRotated[i].x);
            camDir[counter][i].setY(camDirectionRotated[i].y);
            camDir[counter][i].setZ(camDirectionRotated[i].z);
        }

        // ----------------Calculate the up and right vectors for new rotation-------------------
        // re-initialize the output vectors
        camUpVectors.clear();
        camRightVectors.clear();
        camUpVectors.resize(numCamPos);
        camRightVectors.resize(numCamPos);
        // unmap the svm vectors
        cl::unmapSVM(camDirectionRotated);
        cl::unmapSVM(camUpVectors);
        cl::unmapSVM(camRightVectors);
        // Call the function camera definition kernel
        cameraKernel(cl::EnqueueArgs(cl::NDRange(numCamPos)), camDirectionRotated, camUpVectors, camRightVectors, error);
        // grab the output vectors
        cl::mapSVM(camUpVectors);
        cl::mapSVM(camRightVectors);
        // save into global containers
        for(int i = 0; i < numCamPos; i++)
        {
            camUp[counter][i].setX(camUpVectors[i].x);
            camUp[counter][i].setY(camUpVectors[i].y);
            camUp[counter][i].setZ(camUpVectors[i].z);

            camRight[counter][i].setX(camRightVectors[i].x);
            camRight[counter][i].setY(camRightVectors[i].y);
            camRight[counter][i].setZ(camRightVectors[i].z);
        }

        // increment the counter
        counter++;
    }
}

void VisibilityCheck::setup_g()
{
    // Create the required input and output svm data structures
    cl::SVMAllocator<int, cl::SVMTraitCoarse<>> svmAlloc;

    // SVM allocations for frustum calculation kernel
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camPoints(svmAlloc);        // input vector containing camera locations
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camDirection(svmAlloc);     // input vector containing camera directions
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camUpVectors(svmAlloc);     // input svm buffer containing camera up vectors
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camRightVectors(svmAlloc);  // input svm buffer containing camera right vectors
    std::vector<cl_float16, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> camFrustum(svmAlloc);      // output vector to store the calculated frustum
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> farPlanePoint(svmAlloc);    // output vector to store far plane top left point coordinates

    // SVM allocations for control points checking kernel
    std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> contPoints(svmAlloc);
    std::vector<cl_int, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>> outMatG(svmAlloc);

    // Create the kernel functors
    auto frustumCalcKernel =
            cl::KernelFunctor<
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float16, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                cl_float3
            >(camProgram, "frustumCalc");

    auto pointCheckKernel =
            cl::KernelFunctor<
                std::vector<cl_float3, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&,
                cl_float16,
                cl_float3,
                cl_float3,
                std::vector<cl_int, cl::SVMAllocator<int, cl::SVMTraitCoarse<>>>&
            >(camProgram, "checkControlPoints");

    // Setup the camera parameters (hfar, vfar and z_f) to pass it to the kernel function
    cl_float3 camParams;
    camParams.x = onePlus6.hFar;
    camParams.y = onePlus6.vFar;
    camParams.z = onePlus6.z_f;

    // Open a file to write the calculated g matrix to binary files
    QString fileName("/media/anirudh/Data/Documents/PhD/Qt_projects/virtual_machine/MIP_optimizer/gmat32Down/g" + QString::number(1) + ".dat");
    QFile file(fileName);
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);

    int numOrientations = camDir.size();
    int numCamPos = camPos.size();
    int numBackPos = backgroundVox.size();
    //qDebug() << numOrientations;

    // initialize svm vectors for calls to background point checking kernel (done outside the loop because it is constant set of points)
    contPoints.clear();
    contPoints.resize(numBackPos);
    for(int i = 0; i < numBackPos; i++)
    {
        contPoints[i].x = backgroundVox[i].x();
        contPoints[i].y = backgroundVox[i].y();
        contPoints[i].z = backgroundVox[i].z();
    }

    // Loop over the number of orientations
    for(int k = 0; k < numOrientations; k++)
    {
        // initialize the input svm buffers for kernel call to calculate the camera frustums
        camPoints.clear();
        camPoints.resize(numCamPos);
        camDirection.clear();
        camDirection.resize(numCamPos);
        camUpVectors.clear();
        camUpVectors.resize(numCamPos);
        camRightVectors.clear();
        camRightVectors.resize(numCamPos);
        // initialize output svm buffer for fustum calculation kernel
        camFrustum.clear();
        camFrustum.resize(numCamPos);
        farPlanePoint.clear();
        farPlanePoint.resize(numCamPos);
        // Populate the input buffers with data
        for(int i = 0; i < numCamPos; i++)
        {
            // camera positions
            camPoints[i].x = camPos[i].x();
            camPoints[i].y = camPos[i].y();
            camPoints[i].z = camPos[i].z();
            // camera view directions
            camDirection[i].x = camDir[k][i].x();
            camDirection[i].y = camDir[k][i].y();
            camDirection[i].z = camDir[k][i].z();
            // camera up vectors
            camUpVectors[i].x = camUp[k][i].x();
            camUpVectors[i].y = camUp[k][i].y();
            camUpVectors[i].z = camUp[k][i].z();
            // camera right vectors
            camRightVectors[i].x = camRight[k][i].x();
            camRightVectors[i].y = camRight[k][i].y();
            camRightVectors[i].z = camRight[k][i].z();
        }
        // unmap all svm vectors
        cl::unmapSVM(camPoints);
        cl::unmapSVM(camDirection);
        cl::unmapSVM(camUpVectors);
        cl::unmapSVM(camRightVectors);
        cl::unmapSVM(camFrustum);
        cl::unmapSVM(farPlanePoint);

        // Call the kernel function to calculate camera view frustums
        cl_int error;
        frustumCalcKernel(
                    cl::EnqueueArgs(cl::NDRange(numCamPos)),
                    camPoints,
                    camDirection,
                    camUpVectors,
                    camRightVectors,
                    camFrustum,
                    farPlanePoint,
                    camParams,
                    error
                    );

        // unmap output svm vectors
        cl::unmapSVM(camFrustum);
        cl::unmapSVM(farPlanePoint);

        // loop over the number of camera positions to check each gainst all the background points
        for(int i = 0; i < numCamPos; i++)
        {
            // initialize the output svm buffer
            outMatG.clear();
            outMatG.resize(numBackPos);
            // unmap svm
            cl::unmapSVM(contPoints);
            cl::unmapSVM(outMatG);

            // call the kernel function to check all the background points against each camera position and it's frustum
            cl_int error;
            pointCheckKernel(
                        cl::EnqueueArgs(cl::NDRange(numBackPos)),
                        contPoints,
                        camFrustum[i],
                        camPoints[i],
                        farPlanePoint[i],
                        outMatG,
                        error
                        );

            // unmap output svm
            cl::unmapSVM(outMatG);

            // Write the g values variable into the file
            for (int j = 0; j < numBackPos; j++)
            {
                if(outMatG[j] == 1)
                {
                    out << (qint32)k << (qint32)i << (qint32)j << (bool)outMatG[j];
                }
            }
        }
    }
    // close the file after processing all the camera positions
    file.close();
}
