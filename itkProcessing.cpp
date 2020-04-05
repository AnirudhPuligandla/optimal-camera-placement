#include "itkProcessing.h"

void ITKProcessing::vtkToItk(vtkSmartPointer<vtkStructuredPoints> &vtkData, ImageType::ConstPointer &itkData)
{
    // VTK image to ITK image conversion filter
    using FilterType = itk::VTKImageToImageFilter<ImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(vtkData);

    try {
        filter->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "ERROR: " << error << std::endl;
    }

    itkData = filter->GetOutput();
}

void ITKProcessing::itkToVtk(ImageType::ConstPointer &itkData, vtkSmartPointer<vtkStructuredPoints> &vtkData)
{
    // ITK image to VTK image conversion filter
    using FilterType = itk::ImageToVTKImageFilter<ImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(itkData);
    try
    {
        filter->Update();
    }
    catch(itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }

    vtkData->DeepCopy(filter->GetOutput());
}

void ITKProcessing::floatVtkToItk(vtkSmartPointer<vtkStructuredPoints> &vtkData, GradientImageType::ConstPointer &itkData)
{
    // VTK image to ITK image conversion filter
    using FilterType = itk::VTKImageToImageFilter<GradientImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(vtkData);

    try {
        filter->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "ERROR: " << error << std::endl;
    }

    itkData = filter->GetOutput();
}

void ITKProcessing::floatItkToVtk(GradientImageType::ConstPointer &itkData, vtkSmartPointer<vtkStructuredPoints> &vtkData)
{
    // ITK image to VTK image conversion filter
    using FilterType = itk::ImageToVTKImageFilter<GradientImageType>;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(itkData);
    try
    {
        filter->Update();
    }
    catch(itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }

    vtkData->DeepCopy(filter->GetOutput());
}

void ITKProcessing::volumeDilation(vtkSmartPointer<vtkStructuredPoints> &voxData)
{
    // Convert vtk voxel image to itk image
    ImageType::ConstPointer itkData;
    vtkToItk(voxData, itkData);

    // We want to dilate the object (represented with the maximum intensity in the dataset)
    // Get maximum intensity value
    unsigned int maxVal = voxData->GetScalarRange()[1];
    qDebug() << "maximum intensity value: " << maxVal;

    // Create a 3D box structuring element
    using StructuringElementType = itk::FlatStructuringElement<3>;
    StructuringElementType::RadiusType elementRadius;
    elementRadius.Fill(1);              // Use radius=1 for a SE of dimensions 3x3x3
    StructuringElementType structuringElement = StructuringElementType::Box(elementRadius);

    // Create Binary dilate filter
    using BinaryDilateFilterType = itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType>;
    BinaryDilateFilterType::Pointer dilateFilter = BinaryDilateFilterType::New();
    dilateFilter->SetInput(itkData);
    dilateFilter->SetKernel(structuringElement);
    //dilateFilter->SetForegroundValue(2);
    //dilateFilter->SetBackgroundValue(2);
    dilateFilter->SetDilateValue(maxVal);
    dilateFilter->Update();
    itkData = dilateFilter->GetOutput();

    // Convert the output from dilation filter into vtk volume
    itkToVtk(itkData, voxData);
}

void ITKProcessing::generateNormals(vtkSmartPointer<vtkStructuredPoints> &voxData, QVector<vtkSmartPointer<vtkStructuredPoints> > &gradientVolumes)
{
    // Convert the vtk image to ITK image
    GradientImageType::ConstPointer itkData;
    floatVtkToItk(voxData, itkData);

    // Create 3D Sobel operators
    using SobelOperatorType = itk::SobelOperator<GradientPixelType, 3>;
    SobelOperatorType sobelOperatorX, sobelOperatorY, sobelOperatorZ;
    // Define kernel size -- 3x3x3
    itk::Size<3> radius;
    radius.Fill(1);
    // initialize sobel operators
    sobelOperatorX.SetDirection(0);
    sobelOperatorX.CreateToRadius(radius);
    sobelOperatorY.SetDirection(1);
    sobelOperatorY.CreateToRadius(radius);
    sobelOperatorZ.SetDirection(2);
    sobelOperatorZ.CreateToRadius(radius);

    // Filters for neighbourhood operations
    using NeighborhoodOperatorImageFilterType = itk::NeighborhoodOperatorImageFilter<GradientImageType, GradientImageType>;
    NeighborhoodOperatorImageFilterType::Pointer filterX = NeighborhoodOperatorImageFilterType::New();
    NeighborhoodOperatorImageFilterType::Pointer filterY = NeighborhoodOperatorImageFilterType::New();
    NeighborhoodOperatorImageFilterType::Pointer filterZ = NeighborhoodOperatorImageFilterType::New();

    // Filters for rescaling intensity values within [-1,1]
    using RescaleIntensityFiltertype = itk::RescaleIntensityImageFilter<GradientImageType, GradientImageType>;
    RescaleIntensityFiltertype::Pointer rescaleFilterX = RescaleIntensityFiltertype::New();
    RescaleIntensityFiltertype::Pointer rescaleFilterY = RescaleIntensityFiltertype::New();
    RescaleIntensityFiltertype::Pointer rescaleFilterZ = RescaleIntensityFiltertype::New();

    // apply sobel filter in X
    GradientImageType::ConstPointer gradientImageX;
    filterX->SetOperator(sobelOperatorX);
    filterX->SetInput(itkData);
    filterX->Update();
    // rescale the output gradient intensities
    rescaleFilterX->SetInput(filterX->GetOutput());
    rescaleFilterX->SetOutputMinimum(-1.0);
    rescaleFilterX->SetOutputMaximum(1.0);
    rescaleFilterX->Update();
    gradientImageX = rescaleFilterX->GetOutput();

    // Apply sobel filter in Y
    GradientImageType::ConstPointer gradientImageY;
    filterY->SetOperator(sobelOperatorY);
    filterY->SetInput(itkData);
    filterY->Update();
    // rescale the output gradients
    rescaleFilterY->SetInput(filterY->GetOutput());
    rescaleFilterY->SetOutputMinimum(-1.0);
    rescaleFilterY->SetOutputMaximum(1.0);
    rescaleFilterY->Update();
    gradientImageY = rescaleFilterY->GetOutput();

    // Apply sobel filter in Z
    GradientImageType::ConstPointer gradientImageZ;
    filterZ->SetOperator(sobelOperatorZ);
    filterZ->SetInput(itkData);
    filterZ->Update();
    // rescale output gradients
    rescaleFilterZ->SetInput(filterZ->GetOutput());
    rescaleFilterZ->SetOutputMinimum(-1.0);
    rescaleFilterZ->SetOutputMaximum(1.0);
    rescaleFilterZ->Update();
    gradientImageZ = rescaleFilterZ->GetOutput();

    // Add the gradient images into the output QVector after conversion into VTK image
    vtkSmartPointer<vtkStructuredPoints> gradientVTkImageX = vtkSmartPointer<vtkStructuredPoints>::New();
    floatItkToVtk(gradientImageX, gradientVTkImageX);
    gradientVolumes.append(gradientVTkImageX);

    vtkSmartPointer<vtkStructuredPoints> gradientVTkImageY = vtkSmartPointer<vtkStructuredPoints>::New();
    floatItkToVtk(gradientImageY, gradientVTkImageY);
    gradientVolumes.append(gradientVTkImageY);

    vtkSmartPointer<vtkStructuredPoints> gradientVTkImageZ = vtkSmartPointer<vtkStructuredPoints>::New();
    floatItkToVtk(gradientImageZ, gradientVTkImageZ);
    gradientVolumes.append(gradientVTkImageZ);
}

void ITKProcessing::getOccupiedExtent(vtkSmartPointer<vtkStructuredPoints> &voxData, int (&extentObj)[6])
{
    // get the extent of the volume
    int extent[6];
    voxData->GetExtent(extent);
    // initialize the extent of occupied voxels
    extentObj[0] = extent[1];
    extentObj[1] = extent[0];
    extentObj[2] = extent[3];
    extentObj[3] = extent[2];
    extentObj[4] = extent[5];
    extentObj[5] = extent[4];
    // iterate over the volume
    for(int z = 0; z <= extent[5]; z++)
    {
        for(int y = 0; y <= extent[3]; y++)
        {
            for(int x = 0; x <= extent[1]; x++)
            {
                // get the voxel
                unsigned int* val = static_cast<unsigned int*>(voxData->GetScalarPointer(x,y,z));
                // Compare coordinates of occupied voxels and establish boundaries around the vehicle
                if(*val != 0)
                {
                    // boundary in x-axis
                    extentObj[0] = x < extentObj[0] ? x : extentObj[0];
                    extentObj[1] = x > extentObj[1] ? x : extentObj[1];
                    // boundary in y-axis
                    extentObj[2] = y < extentObj[2] ? y : extentObj[2];
                    extentObj[3] = y > extentObj[3] ? y : extentObj[3];
                    // boundary in z-axis
                    extentObj[4] = z < extentObj[4] ? z : extentObj[4];
                    extentObj[5] = z > extentObj[5] ? z : extentObj[5];
                }
            }
        }
    }
}

void ITKProcessing::slicSegmentation(QVector<QVector3D> &boundaryVoxels, QVector<QVector3D> &boundaryNormals, vtkSmartPointer<vtkStructuredPoints> &voxData, int clusters, QVector<QVector3D> &segmentedVoxelCenters, QVector<QVector3D> &segmentedVoxelNormals)
{
    int numCamPos = boundaryVoxels.size();
    int clusterSize = numCamPos/clusters;
    /* Grid interval as defined in the paper -
     * Ideally for volume interval should computed as cube root of 'clusterSize' for a volume
     * But, we use the definition as in the paper considering that we are segmenting a surface
     * (boundary has thickness of only one voxel)
     */
    int S = sqrt(numCamPos/clusters);
    float residualThreshold = 10; // threshold for Residual error E in the paper
    float residualE = 50;         // Initialize the residual error
    float paramM = 10;            // Parameter m in the paper

    /* set of cluster centers
     * Holds normals and associated positions [x',y',z',x,y,z]^T
     */
    //QMap<QVector3D, QVector3D> clusterCenters;
    QVector<QVector<QVector3D>> clusterCenters;

    // Counter for number of centers adjusted
    int adjustedCount = 0;

    // iterate over the set of voxels and initialize cluster centers
    for(int i = 0; i < numCamPos; i++)
    {
        if(i % clusterSize == 0)
        {
            // calculate the magnitude of the current surface normal for neighborhood checking
            float currMag = boundaryNormals[i].length();
            // Vector to hold current cluster center for neighborhood checking
            int centerInd = i;
            // Iterate over the 3x3x3 neighbourhood
            for(int x = boundaryVoxels[i].x() - 1; x <= boundaryVoxels[i].x() + 1; x++)
            {
                for(int y = boundaryVoxels[i].y() - 1; y <= boundaryVoxels[i].y() + 1; y++)
                {
                    for(int z = boundaryVoxels[i].z() - 1; z <= boundaryVoxels[i].z() + 1; z++)
                    {
                        /* Check if the current voxel index exists in set of boundary voxels
                         * if voxel doesn't exist skip
                         */
                        int neighbour = boundaryVoxels.indexOf(QVector3D(x,y,z));
                        if(neighbour >= 0)
                        {
                            // calculate the magnitude of the normal of that neighbour
                            float neighbourMag = boundaryNormals[neighbour].length();
                            // if magnitude of the neighbour is less update the cluster center
                            if(neighbourMag < currMag)
                            {
                                currMag = neighbourMag;
                                centerInd = neighbour;
                            }
                        }

                    }
                }
            }
            // Update the adjusted centers count
            if(centerInd != i)
                adjustedCount++;
            //QVector<float> center = {boundaryNormals[i].x(),boundaryNormals[i].y(),boundaryNormals[i].z(),boundaryVoxels[i].x(),boundaryVoxels[i].y(),boundaryVoxels[i].z()};
            clusterCenters.append({boundaryNormals[centerInd], boundaryVoxels[centerInd]});
        }
    }

    qDebug() << "Number of clusters = " << clusterCenters.size();
    qDebug() << "Number of cluster centers adjusted = " << adjustedCount;
    
    // ------------------------------- SLIC loop --------------------------------
    int numClusters = clusterCenters.size();

    // ----Initialize the vector to hold distances----
    QVector<float> distances;
    distances.resize(numCamPos);
    /* The distances for all voxels are initialized with the maximum possible distances
     * surface normal distances are maximum between (-1,-1,-1) and (1,1,1)
     * while spatial distance is maximum at the extents of the voxel environment
     */
    int extent[6];
    voxData->GetExtent(extent);
    //getOccupiedExtent(voxData, extentObj);
    // calculate maximum possible distance between voxels
    //int nearIndex = boundaryVoxels.indexOf(QVector3D(extentObj[0],extentObj[2],extentObj[4]));
    //int farIndex = boundaryVoxels.indexOf(QVector3D(extentObj[1],extentObj[3],extentObj[5]));
    //QVector3D nearNormal = boundaryNormals[nearIndex];
    //QVector3D farNormal = boundaryNormals[farIndex];
    float maxNormDist = QVector3D(-1,-1,-1).distanceToPoint(QVector3D(1,1,1));
    float maxDist = QVector3D(extent[0],extent[2],extent[4]).distanceToPoint(QVector3D(extent[1],extent[3],extent[5]));
    // distance formula as mentioned in the paper
    float maxTotalDist = maxNormDist + ((paramM/S) * maxDist);
    for(int i = 0; i < numCamPos; i++)
    {
        // initialize all distances with max possible distance
        distances[i] = maxTotalDist;
    }
    // -----Initialize the vector of lables-------
    QVector<int> labels;
    labels.resize(numCamPos);
    for(int i = 0; i < numCamPos; i++)
    {
        // initialize the labels with -1 to keep track of voxels not associated with any cluster center
        labels[i] = -1;
    }

    // Counter for SLIC iterations
    int itCounter = 0;
    // SLIC iterator loop
    while(residualE >= residualThreshold)
    {
        // Loop over the cluster centers
        for(int i = 0; i < clusterCenters.size(); i++)
        {
            // Search in the neighbourhood of 2Sx2Sx2S of the cluster center
            QVector3D currCenter = clusterCenters[i][1];
            QVector3D currNormal = clusterCenters[i][0];
            /* Create a neighbourhood iterator and Check for image boundary conditions
             * -> ideally, image boundaries cannot be encountered in x,z directions as
             *    boundary voxels are very far from image boundaries
             * -> So, check for image boundaries only in y direction
             */
            // Check for bottom boundary
            int y1 = currCenter.y() - S;
            if(y1 < extent[2])
                y1 = extent[2];
            // Check for top boundary
            int y2 = currCenter.y() + S;
            if(y2 > extent[3])
                y2 = extent[3];

            for(int x = currCenter.x() - S; x <= currCenter.x() + S; x++)
            {
                for(int y = y1; y <= y2; y++)
                {
                    for(int z = currCenter.z() - S; z <= currCenter.z() + S; z++)
                    {
                        // Get the value at the neighbourhood voxel
                        unsigned int* val = static_cast<unsigned int*>(voxData->GetScalarPointer(x,y,z));
                        // process only the occupied voxels in the neighborhood
                        if(*val != 0)
                        {
                            // Get the index of this voxel location from the set of boundary voxels
                            int currIndex = boundaryVoxels.indexOf(QVector3D(x,y,z));
                            // Calculate distances according to the paper
                            float normalDistance = boundaryNormals[currIndex].distanceToPoint(currNormal);
                            float spatialDistance = boundaryVoxels[currIndex].distanceToPoint(currCenter);
                            float totalDistance = normalDistance + ((paramM/S) * spatialDistance);
                            // Compare with existing distance and update the label with that of closest cluster center
                            if(totalDistance < distances[currIndex])
                            {
                                // update the distance for later comparison
                                distances[currIndex] = totalDistance;
                                // update the label for the current voxel
                                labels[currIndex] = i;
                            }
                        }
                    }
                }
            }
        }
        /* After all cluster centers are processed for this iteration,
         * generate new centers as the average normal,x,y,z vector of all voxels belonging to that cluster
         */

        // initialize vector to hold sums of boundary normals and boundary voxels as per the labels
        QVector<QVector3D> normalSum;
        normalSum.resize(numClusters);
        QVector<QVector3D> voxelSum;
        voxelSum.resize(numClusters);
        // Initialize a vector to hold the count of associated voxels with each cluster center
        QVector<int> clusterSizes;
        clusterSizes.resize(numClusters);
        // initialize all the vectors with zeros
        for(int i = 0; i < numClusters; i++)
        {
            normalSum[i] = QVector3D(0,0,0);
            voxelSum[i] = QVector3D(0,0,0);
            clusterSizes[i] = 0;
        }

        // Iterate over the set of labels
        for(int i = 0; i < numCamPos; i++)
        {
            /* get the cluster index from the labels
             * and keep adding the normals and voxel positions corresponding to that cluster
             */
            int tempIndex = labels[i];
            normalSum[tempIndex] += boundaryNormals[i];
            voxelSum[tempIndex] += boundaryVoxels[i];
            clusterSizes[tempIndex] += 1;
        }

        // make a copy of clusterCenters for calculating residual error
        QVector<QVector3D> clusterCentersOld;
        clusterCentersOld.resize(numClusters);
        for(int i = 0; i < numClusters; i++)
        {
            clusterCentersOld[i] = clusterCenters[i][1];
        }

        // Update the cluster centers with new average values
        for(int i = 0; i < numClusters; i++)
        {
            // Number of voxels in that cluster
            int tempClusterSize = clusterSizes[i];
            // Average the normals at each label
            clusterCenters[i][0] = normalSum[i]/tempClusterSize;
            /* Average the voxel points
             * Round off the average value to the nearest integer
             * because voxel positions cannot take floating point values
             */
            clusterCenters[i][1].setX(round(voxelSum[i].x()/tempClusterSize));
            clusterCenters[i][1].setY(round(voxelSum[i].y()/tempClusterSize));
            clusterCenters[i][1].setZ(round(voxelSum[i].z()/tempClusterSize));
        }

        // Calculate the residual error
        double errorSum = 0;
        for(int i = 0; i < numClusters; i++)
        {
            QVector3D tempL1 = clusterCentersOld[i] - clusterCenters[i][1];
            errorSum += abs(tempL1.x()) + abs(tempL1.y()) + abs(tempL1.z());
        }
        // update the residual error as the average of all the l1 distances
        residualE = errorSum/numClusters;
        itCounter++;
        qDebug() << " Completed " << itCounter << " iterations.";
    }

    // setup output vectors
    for(int i = 0; i < clusterCenters.size(); i++)
    {
        segmentedVoxelCenters.append(clusterCenters[i][1]);
    }

}
