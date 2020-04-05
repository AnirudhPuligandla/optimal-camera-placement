/* Handles the simulation environment and all voxel processing and representation
 * voxWorld represents the complete environment and the values represent the following:
 * 0 -> Background points lying outisde the circumference of 12m
 * 1 -> Actually considered background voxels for optimization
 * 2 -> Represents vehicle boundaries / initial set of camera locations for optimization
 * 3 -> Represents vehicle (occupied voxels for object)
 * 4 -> Selected cluster centers (initial seeds)
 */

#include "vtkwidget.h"

// To avoid Error: no override found for 'vtkPolyDataMapper'
#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);

VTKWidget::VTKWidget(QWidget *parent) : QVTKOpenGLNativeWidget(parent)
{
    // Setup objects for initial scene
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkSphereSource> sphereSource;

    vtkNew<vtkPolyDataMapper> sphereMapper;
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

    vtkNew<vtkActor> sphereActor;
    sphereActor->SetMapper(sphereMapper);
    sphereActor->GetProperty()->SetColor(colors->GetColor4d("Tomato").GetData());

    // Call visualize function to Display initial scene
    visualizeVTK(sphereActor);

    // initialize the file
    fileNameVTK = "/media/anirudh/Data/Documents/PhD/Qt_projects/ocp_simulation/resources/bulldoser.vtk";
}

QSize VTKWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize VTKWidget::sizeHint() const
{
    return QSize(400, 400);
}

void VTKWidget::loadModelPressed()
{
    // Split the string check file extension
    std::vector<std::string> splitFileName;
    boost::split(splitFileName, fileNameVTK, [](char c){return c == '.';});

    // If extension is obj we need to voxelize it using binvox
    if(splitFileName[1].compare("obj") == 0)
    {
        // Setup Qprocess environment to voxelize using binvox
        //QString command = "xvfb :99 -screen 0 640x480x24 & setenv DISPLAY=:99 & xvfb-run -s \"-screen 0 640x480x24\" /home/anirudh/Downloads/binVox/binvox -d 128 -t vtk -pb " + QString::fromUtf8(fileNameOBJ.c_str());
        QString command = "/home/anirudh/Downloads/binVox/binvox -d 128 -t vtk " + QString::fromUtf8(fileNameVTK.c_str());
        //QString command = "/home/anirudh/opt/foxitsoftware/foxitreader";
        QProcess *binVoxProcess = new QProcess();
        qDebug().noquote()<<"Voxelizing using binvox with the command: " << command;
        binVoxProcess->setProcessChannelMode(QProcess::ForwardedChannels);

        // extract and display qprocess errors
        connect(binVoxProcess, &QProcess::errorOccurred, this, &VTKWidget::processError);

        binVoxProcess->start(command);
        // update filename with newly created vtk file
        fileNameVTK = splitFileName[0] + ".vtk";
    }

    //vtkNew<vtkOBJReader> reader;
    //auto reader = vtkSmartPointer<vtkOBJImporter>::New();

    // read data from vtk file as a structured grid
    vtkNew<vtkGenericDataObjectReader> reader;
    reader->SetFileName(fileNameVTK.c_str());
    reader->Update();

    // Check for type of data in the .vtk file
    if(reader->IsFileStructuredPoints())
    {

        // Read underlying as a 'structuredPoints' dataset
        vtkSmartPointer<vtkStructuredPoints> voxData = vtkSmartPointer<vtkStructuredPoints>::New();
        voxData = reader->GetStructuredPointsOutput();

        // get extent of volume
        int extent[6];
        voxData->GetExtent(extent);
        qDebug() << "Extent Old: " << "(" << extent[0] << "," << extent[1] << ")" << "(" << extent[2] << "," << extent[3] << ")" << "(" << extent[4] << "," << extent[5] << ")";

        // get data dimensions along each axis
        int dimsOld[3];
        voxData->GetDimensions(dimsOld);
        qDebug() << "Dimensions: " << "(" << dimsOld[0] << "," << dimsOld[1] << "," << dimsOld[2] << ")";

        // Containers to hold extent of object boundaries (in order (x1,x2,y1,y2,z1,z2)) 1 -> min. 2-> max.
        int extentObj[6];
        // initialize obejct extent
        extentObj[0] = extent[1];
        extentObj[1] = extent[0];
        extentObj[2] = extent[3];
        extentObj[3] = extent[2];
        extentObj[4] = extent[5];
        extentObj[5] = extent[4];
        // Iterate over the data and get boundaries of occupied voxels within the 128x128x128 grid
        for(int z = 0; z < dimsOld[2]-1; z++)
        {
            for(int y = 0; y < dimsOld[1]-1; y++)
            {
                for(int x = 0; x < dimsOld[0]-1; x++)
                {
                    // Get the value at that (x,y,z) location
                    bool* val = static_cast<bool*>(voxData->GetScalarPointer(x,y,z));
                    //qDebug () << *val;
                    // Compare coordinates of value 1 and establish boundaries around the vehicle
                    if(*val)
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
        qDebug() << "object boundaries: " << "(" << extentObj[0] << "," << extentObj[1] << ")" << "(" << extentObj[2] << "," << extentObj[3] << ")" << "(" << extentObj[4] << "," << extentObj[5] << ")";
        // get object center
        int objCenter[3];
        objCenter[0] = (extentObj[1] - extentObj[0])/2;
        objCenter[1] = (extentObj[3] - extentObj[2])/2;
        objCenter[2] = (extentObj[5] - extentObj[4])/2;

        /* *** Extend this for a general case (i.e, take inputs from user)
         * dimensions of bulldoser is ~ 8.1 x 4.5 x 4 m --- translates to ~ 63.3 mm / voxel on z-axis (length)
         * a distance of 12m roughly translates to 190 voxels, so, we create a 400x(maxBound_y)x400 grid with the vehicle at the center
         */
        int worldCenter = 200;      // 400/2
        // Translate the object in x-z such that the object center lies at 100 in x-z
        int transFactor[3];
        transFactor[0] = worldCenter - objCenter[0];
        transFactor[1] = 0;                         // No translation in y as anything below ground plane is not of interest
        transFactor[2] = worldCenter - objCenter[2];

        vtkNew<vtkImageTranslateExtent> transExtent;
        transExtent->SetTranslation(transFactor);
        transExtent->SetInputData(voxData);
        transExtent->Update();
        voxData->DeepCopy(transExtent->GetOutput());

        // update object center w.r.t. world
        objCenter[0] += transFactor[0];
        objCenter[1] += transFactor[1];
        objCenter[2] += transFactor[2];
        qDebug() << "Object center w.r.t. world: (" << objCenter[0] << "," << objCenter[1] << "," << objCenter[2] << ")";

        // Define and Set new extent for the simulation environment
        int extentNew[6];
        extentNew[0] = 0;
        extentNew[1] = 399;
        extentNew[2] = 0;
        extentNew[3] = extentObj[3] + 2;                // Space over the top of the vehicle is not of interest
        extentNew[4] = 0;
        extentNew[5] = 399;

        // Create a new volume and copy contents from translated model
        voxWorld = vtkSmartPointer<vtkStructuredPoints>::New();
        voxWorld->SetExtent(extentNew);
        voxWorld->AllocateScalars(VTK_UNSIGNED_INT, 1);
        // get extent of translated object
        voxData->GetExtent(extent);

        for(int z = extentNew[4]; z <= extentNew[5]; z++)
        {
            for(int y = extentNew[2]; y <= extentNew[3]; y++)
            {
                for(int x = extentNew[0]; x <= extentNew[1]; x++)
                {
                    unsigned int* valFinal = static_cast<unsigned int*>(voxWorld->GetScalarPointer(x,y,z));
                    *valFinal = 0;
                    // Check if the index falls within the dimensions of translated volume
                    if(x >= extent[0] && x <= extent[1] && y >= extent[2] && y <= extent[3] && z >= extent[4] && z <= extent[5])
                    {
                        // get the value at that location and update in the voxWorld volume if it is one
                        bool* valTemp = static_cast<bool*>(voxData->GetScalarPointer(x,y,z));
                        if(*valTemp)
                            *valFinal = 3;
                    }
                }
            }
        }

        voxWorld->GetExtent(extent);
        qDebug() << "Extent New: " << "(" << extent[0] << "," << extent[1] << ")" << "(" << extent[2] << "," << extent[3] << ")" << "(" << extent[4] << "," << extent[5] << ")";

        visualizeVolume(voxWorld);
    }
}

void VTKWidget::processError(QProcess::ProcessError error)
{
    qDebug() << "error enum val = " << error << endl;
}

void VTKWidget::dilatePressed()
{
    /* First extract the object boundary voxels, use one of the following:
     * 1) morphological dilation on the volume and get difference of volumes
     * 2) marching cubes (produces triangles or polygons)
     * 3) Dividing cubes (generates points representing the boundary)
     */

    // apply dilation on the volume using itkProcessing class
    vtkSmartPointer<vtkStructuredPoints> dilatedWorld = vtkSmartPointer<vtkStructuredPoints>::New();
    dilatedWorld->DeepCopy(voxWorld);
    int extent[6];
    voxWorld->GetExtent(extent);

    // Use ITK toolkit for morphological dilation
    volumeProcessor.volumeDilation(dilatedWorld);
    //qDebug() << "Number of cells before dilation: " << voxWorld->GetNumberOfCells();
    //qDebug() << "Number of cells after dilation: " << dilatedWorld->GetNumberOfCells();

    // Convert the image datatype to float to be able to process gradients
    vtkSmartPointer<vtkStructuredPoints> floatDilatedWorld = vtkSmartPointer<vtkStructuredPoints>::New();
    qDebug() << dilatedWorld->GetScalarType();
    vtkNew<vtkImageShiftScale> shiftFilter;
    shiftFilter->SetOutputScalarTypeToFloat();
    shiftFilter->SetInputData(dilatedWorld);
    shiftFilter->SetScale(1);
    shiftFilter->SetShift(0);
    shiftFilter->Update();
    floatDilatedWorld->DeepCopy(shiftFilter->GetOutput());
    //qDebug() << floatDilatedWorld->GetScalarType();
    //qDebug() << floatDilatedWorld->GetScalarRange()[0] << " " << floatDilatedWorld->GetScalarRange()[1];

    // Get Gradient images for the dilated volume
    gradientVolumes.clear();
    volumeProcessor.generateNormals(floatDilatedWorld, gradientVolumes);

    // volume image showing only the voxels added after dilation
    volBoundary = vtkSmartPointer<vtkStructuredPoints>::New();
    volBoundary->SetExtent(extent);
    volBoundary->AllocateScalars(VTK_UNSIGNED_INT, 1);

    // initialize the gradient magnitude image
    volGradientMagnitude = vtkSmartPointer<vtkStructuredPoints>::New();
    volGradientMagnitude->SetExtent(extent);
    volGradientMagnitude->AllocateScalars(VTK_FLOAT, 1);

    // ----------------------------Extract data Points and setup additional descriptive volumes-------------------------------------

    // initialize the containers that hold data points
    boundaryVoxels.clear();
    boundaryNormals.clear();
    backgroundVoxels.clear();

    qDebug() << gradientVolumes[0]->GetScalarRange()[0] << " " << gradientVolumes[0]->GetScalarRange()[1];

    int volDims[3];
    voxWorld->GetDimensions(volDims);

    // get world center
    int center[3] = {volDims[0]/2, volDims[1]/2, volDims[2]/2};

    int count = 0;
    int countBackground = 0;
    for(int z = 0; z <= volDims[2]-1; z++)
    {
        for(int y = 0; y <= volDims[1]-1; y++)
        {
            for(int x = 0; x <= volDims[0]-1; x++)
            {
                // Get the values at the location (x,y,z) from different voxel environments
                unsigned int* valWorld = static_cast<unsigned int*>(voxWorld->GetScalarPointer(x,y,z));
                unsigned int* valDilated = static_cast<unsigned int*>(dilatedWorld->GetScalarPointer(x,y,z));
                unsigned int* boundVal = static_cast<unsigned int*>(volBoundary->GetScalarPointer(x,y,z));
                float* normalX = static_cast<float*>(gradientVolumes[0]->GetScalarPointer(x,y,z));
                float* normalY = static_cast<float*>(gradientVolumes[1]->GetScalarPointer(x,y,z));
                float* normalZ = static_cast<float*>(gradientVolumes[2]->GetScalarPointer(x,y,z));
                float* gradientMagnitude = static_cast<float*>(volGradientMagnitude->GetScalarPointer(x,y,z));
                *gradientMagnitude = 0.0f;
                *boundVal = 0;

                /* Calculate the distance of the point from the center of voxWorld for collecting background voxels within 12m
                 * According to environment definition 12m ~= 190 voxels
                 */
                float distance = sqrt(pow(x - center[0], 2) + pow(y - center[1], 2) + pow(z - center[2], 2));

                /* Value 1 at (x,y,z) in dilatedWorld but not in voxWorld represents boundary voxels.
                 * So, change those values to '2' and leave the rest as they are
                 */
                if(*valDilated==3 && *valWorld==0)
                {
                    // Calculate magnitude of the gradient vector at each boundary voxel
                    *gradientMagnitude = sqrt(pow(*normalX,2) + pow(*normalY,2) + pow(*normalZ,2));
                    // Add dilated voxels into the global volume voxWorld
                    *valWorld = 2;
                    // Do not add points with normals parallel to y-axis { (0,1,0) and (0,-1,0) } cases (invert the values because sobel filter produces gradients facing into the object)
                    QVector3D voxNormal = QVector3D(-*normalX, -*normalY, -*normalZ);
                    if((voxNormal != QVector3D(0.0f, 1.0f, 0.0f)) && (voxNormal != QVector3D(0.0f, -1.0f, 0.0f)))
                    {
                        // Collect all the voxel coordinates lying on the vehicle boundary
                        boundaryVoxels.append(QVector3D(x,y,z));
                        // Also collect the corresponding surface normal
                        boundaryNormals.append(voxNormal);
                        // Highlight the boundary voxels with a value '2' in a boundary volume dataset
                        *boundVal = 2;
                        count++;
                    }
                }
                /* Collect all the boundary voxels (representing zero in both voxWorld and dilatedWorld)
                 * We also want to highlight these voxels with different opacity from the uncollected background voxels
                 */
                else if(*valDilated==0 && distance <= 190.0f)
                {
                    backgroundVoxels.append(QVector3D(x,y,z));
                    // highlight the voxels with value 1
                    *valWorld = 1;
                    countBackground++;
                }
            }
        }
    }

    // Call the function perform SLIC based segmentation
    superVoxels.clear();
    superNormals.clear();
    int clusters = 250;             // Number of required boundary voxels after segmentation
    volumeProcessor.slicSegmentation(boundaryVoxels,boundaryNormals,volBoundary,clusters,superVoxels,superNormals);

    // Update valWorld with initial cluster center seeds
    for(int i = 0; i < superVoxels.size(); i++)
    {
        unsigned int* valWorld = static_cast<unsigned int*>(voxWorld->GetScalarPointer(superVoxels[i].x(), superVoxels[i].y(), superVoxels[i].z()));
        // Set the cluster center seeds with value 4
        *valWorld = 100;
    }

    // Save all the images in .vtk files for verification on paraview
    std::vector<std::string> splitFileName;
    boost::split(splitFileName, fileNameVTK, [](char c){return c == '.';});

    vtkNew<vtkGenericDataObjectWriter> writer;
    vtkNew<vtkGenericDataObjectWriter> writerFloatDilated;
    vtkNew<vtkGenericDataObjectWriter> writerBoundary;
    vtkNew<vtkGenericDataObjectWriter> writerGradient;
    vtkNew<vtkGenericDataObjectWriter> writerGradientX;
    // Save the global volume
    writer->SetInputData(voxWorld);
    writer->SetFileName((splitFileName[0] + "_world.vtk").c_str());
    writer->Write();
    // Save the volume model with only boundary voxels
    writerBoundary->SetInputData(volBoundary);
    writerBoundary->SetFileName((splitFileName[0] + "_boundary.vtk").c_str());
    writerBoundary->Write();
    // Save the floating point dilated volume
    writerFloatDilated->SetInputData(floatDilatedWorld);
    writerFloatDilated->SetFileName((splitFileName[0] + "_floatDilated.vtk").c_str());
    writerFloatDilated->Write();
    // Save the gradient magnitude image
    writerGradient->SetInputData(volGradientMagnitude);
    writerGradient->SetFileName((splitFileName[0] + "_gradientMagnitude.vtk").c_str());
    writerGradient->Write();
    // Save gradients image in X direction
    writerGradientX->SetInputData(gradientVolumes[0]);
    writerGradientX->SetFileName((splitFileName[0] + "_gradientX.vtk").c_str());
    writerGradientX->Write();

    // Visualize only the boundary voxels after dilation
    qDebug() << "Number of boundary voxels: " << count;
    qDebug() << "Number of background points collected: " << countBackground;
    visualizeVolume(dilatedWorld);

    //-------------------------------------Perform Visibility Checks-------------------------------
    qDebug() << "Initial Data Extracted...Performing visibility checks...";

    // Create Visibility check object
    //VisibilityCheck visibilityCheck(boundaryVoxels, boundaryNormals, backgroundVoxels);
}

void VTKWidget::visualizeVTK(vtkNew<vtkActor> &actor)
{
    // Colors object to set background color
    vtkNew<vtkNamedColors> colors;
    vtkColor3d backgroundColor = colors->GetColor3d("SteelBlue");

    // Setup openGL render window
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;

#if VTK890
    setRenderWindow(renderWindow);
#else
    SetRenderWindow(renderWindow);
#endif

    // Setup renderer
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(actor);
    renderer->SetBackground(backgroundColor.GetData());
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Azimuth(30);
    renderer->GetActiveCamera()->Elevation(30);
    renderer->GetActiveCamera()->Dolly(1.5);
    renderer->ResetCameraClippingRange();

#if VTK890
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("VTK Widget");
#else
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("VTK Widget");
#endif

    renderWindow->Render();
}

void VTKWidget::visualizeVolume(vtkSmartPointer<vtkStructuredPoints> &voxelData)
{
    // setup data for visualization
    vtkNew<vtkSmartVolumeMapper> mapper;
    mapper->SetBlendModeToComposite();
    //mapper->SetInputConnection(reader->GetOutputPort());
    mapper->SetInputData(voxelData);
    mapper->SetRequestedRenderModeToRayCast();
    // VTK volume property
    vtkNew<vtkVolumeProperty> volumeProperty;
    volumeProperty->ShadeOff();
    volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

    // Define opacity for data
    vtkNew<vtkPiecewiseFunction> compositeOpacity;
    compositeOpacity->AddPoint(0.0, 0.001);
    compositeOpacity->AddPoint(1.0, 0.01);
    compositeOpacity->AddPoint(2.0, 0.6);                       // Give slightly different opacity for boundary voxels
    compositeOpacity->AddPoint(3.0, 1.0);
    //volumeProperty->SetScalarOpacity(compositeOpacity);       // composite opacity first
    volumeProperty->SetGradientOpacity(compositeOpacity);

    // Define colors
    vtkNew<vtkColorTransferFunction> color;
    color->AddRGBPoint(0.0, 1.0, 1.0, 1.0);                     // background points are white
    //color->AddRGBPoint(1.0, 1.0, 0.0, 0.0);                     // vehicle is in red
    color->AddRGBPoint(3.0, 1.0, 0.0, 0.0);                     // Boundary voxels in green
    volumeProperty->SetColor(color);

    // VTK volume (replaces actor)
    vtkNew<vtkVolume> volume;
    volume->SetMapper(mapper);
    volume->SetProperty(volumeProperty);

    // Colors object to set background color
    vtkNew<vtkNamedColors> colors;
    vtkColor3d backgroundColor = colors->GetColor3d("Grey");

    // Setup openGL render window
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;

#if VTK890
    setRenderWindow(renderWindow);
#else
    SetRenderWindow(renderWindow);
#endif

    // Setup renderer
    vtkNew<vtkRenderer> renderer;
    renderer->AddViewProp(volume);
    renderer->SetBackground(backgroundColor.GetData());
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Azimuth(30);
    renderer->GetActiveCamera()->Elevation(30);
    renderer->GetActiveCamera()->Dolly(1.5);
    renderer->ResetCameraClippingRange();

#if VTK890
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("VTK Widget");
#else
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("VTK Widget");
#endif

    renderWindow->Render();
}
