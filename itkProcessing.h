#ifndef ITKPROCESSING_H
#define ITKPROCESSING_H

#include <vtkImageData.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>

#include <itkVTKImageToImageFilter.h>
#include <itkFlatStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkSobelOperator.h>
#include <itkNeighborhoodOperatorImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

#include <QVector>
#include <QVector3D>
#include <QDebug>

#include <math.h>

using PixelType = unsigned int;
using GradientPixelType = float;
//constexpr unsigned int Dimension = 3;
using ImageType = itk::Image<PixelType, 3>;
using GradientImageType = itk::Image<GradientPixelType, 3>;

class ITKProcessing
{
public:
    // VTK image to ITK image conversion
    void vtkToItk(vtkSmartPointer<vtkStructuredPoints> &vtkData, ImageType::ConstPointer &itkData);
    // ITK image to VTK image conversion
    void itkToVtk(ImageType::ConstPointer &itkData, vtkSmartPointer<vtkStructuredPoints> &vtkData);
    // VTK to ITK conversion helper function to handle images of type float (for gradients)
    void floatVtkToItk(vtkSmartPointer<vtkStructuredPoints> &vtkData, GradientImageType::ConstPointer &itkData);
    // ITK to VTK conversion for floating point images
    void floatItkToVtk(GradientImageType::ConstPointer &itkData, vtkSmartPointer<vtkStructuredPoints> &vtkData);
    // Function for morphological dilation
    void volumeDilation(vtkSmartPointer<vtkStructuredPoints> &voxData);
    /* Function to calculate surface normals on the boundary voxels -- uses sobel filters
     * voxData         -> input voxel volume
     * gradientVolumes -> output 3 volumes for x,y,z gradients
     */
    void generateNormals(vtkSmartPointer<vtkStructuredPoints> &voxData, QVector<vtkSmartPointer<vtkStructuredPoints>> &gradientVolumes);
    /* Function to do SLIC segmentation on the boundary voxels
     * Inputs : set of 'boundaryVoxels' representing the voxels on the boundary after dilation
     *          set of 'boundaryNormals' representing the voxel normals
     *          Number of required super voxels 'clusters'
     * Outputs : set of 'segmentedVoxelCenters' containing the centers of segmented super voxels
     *           set of 'segmentedVoxelNormals' representing the normal vectors of super voxels
     */
    void slicSegmentation(QVector<QVector3D> &boundaryVoxels,
                          QVector<QVector3D> &boundaryNormals,
                          vtkSmartPointer<vtkStructuredPoints> &voxData,
                          int clusters,
                          QVector<QVector3D> &segmentedVoxelCenters,
                          QVector<QVector3D> &segmentedVoxelNormals);

private:
    /* function to return the minimum and maximum extent of occupied voxels
     * input - vtk voxel volume
     * output: 'extentObj' 6 integers specifiying the extent of occupied voxels
     */
    void getOccupiedExtent(vtkSmartPointer<vtkStructuredPoints> &voxData, int (&extentObj)[6]);
};

#endif // ITKPROCESSING_H
