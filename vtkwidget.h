#ifndef VTKWIDGET_H
#define VTKWIDGET_H

#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkVersion.h>
#include <vtkOBJReader.h>
#include <vtkCamera.h>
//#include <vtkOBJImporter.h>
#include <vtkTexture.h>
#include <vtkPolyData.h>
//#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredPoints.h>
//#include <vtkStructuredPointsReader.h>
//#include <vtkImageDataGeometryFilter.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkImageTranslateExtent.h>
//#include <vtkImageDilateErode3D.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkImageShiftScale.h>

#include <string>
#include <math.h>

#include <QProcess>
#include <QDebug>
#include <boost/algorithm/string.hpp>
#include <QString>
#include <QVector>
#include <QVector3D>

#include "itkProcessing.h"
#include "visibilitycheck.h"

#if VTK_VERSION_NUMBER >= 89000000000ULL
#define VTK890 1
#endif

#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>

class VTKWidget : public QVTKOpenGLNativeWidget
{
    Q_OBJECT

public:
    VTKWidget(QWidget *parent = 0);
    ~VTKWidget() = default;

    // Functions to set vtk widget dimensions
    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;

    // Renderer function
    void visualizeVTK(vtkNew<vtkActor> &actor);
    // Volume renderer function
    void visualizeVolume(vtkSmartPointer<vtkStructuredPoints> &voxData);
    //

public slots:
    void loadModelPressed();            // Load model from .obj files using vtk
    void processError(QProcess::ProcessError error);                // slot to display QProcess errors
    void dilatePressed();               // Initialize data from model for optimization

private:
    vtkSmartPointer<vtkStructuredPoints> voxWorld;                      // private data memeber to store voxel data
    vtkSmartPointer<vtkStructuredPoints> volBoundary;                   // container for only the dilated boundary
    vtkSmartPointer<vtkStructuredPoints> volGradientMagnitude;          // volume showing gradient magnitudes in 3D
    QVector<vtkSmartPointer<vtkStructuredPoints>> gradientVolumes;      // vector holding gradient volumes after dilation in three dimensions
    //vtkSmartPointer<vtkStructuredPoints> boundaryVoxels;              // Data member to save only the voxels on vehicle boundaries
    std::string fileNameVTK;                                            // Global base file

    ITKProcessing volumeProcessor;                                      // Private itkProcessing object to acess image processing on the volume

    QVector<QVector3D> boundaryVoxels;                                  // container for holding the list of voxels on boundary of the vehicle (represented with x,y,z index)
    QVector<QVector3D> backgroundVoxels;                                // Container to hold list of background points that need to be covered (represented with x,y,z index)
    QVector<QVector3D> boundaryNormals;                                 // Holds the list of surface normals at boundary voxels
    QVector<QVector3D> superVoxels;                                     // Holds the set of super voxels after SLIC segmentation
    QVector<QVector3D> superNormals;                                    // Holds the set of normals to superVoxels after SLIC segmentation
};

#endif // VTKGUI_H
