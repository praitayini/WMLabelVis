#include <vtkActor.h>
#include <vtkBalloonRepresentation.h>
#include <vtkBalloonWidget.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkRegularPolygonSource.h>
#include <vector>
#include <stack>
#include <string>
#include <iostream>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkDiffusionTensor3D.h>
#include <itkVector.h>
#include <itkVectorImage.h>
#include "vtkSmartPointer.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkImageProperty.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"
#include "vtkCellPicker.h"
#include "vtkPointData.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkStreamTracer.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkLineSource.h"
#include "vtkCellArray.h"
#include "vtkNamedColors.h"
#include "vtkNew.h"
#include <vector>
#include <stack>
#include <string>
#include <iostream>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>
#include <itkDiffusionTensor3D.h>
#include <itkVector.h>
#include <itkVectorImage.h>
#include "vtkSmartPointer.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"

typedef itk::Image < int,3 > AtlasImageType;


///nfs/masi/kanakap/open_srx
class CustomMouseCallback : public vtkCommand
{
public:
  // To use in the CustomTimerCallback
  //  double pickedVoxel ;
  //  double window_pickedVoxel ;
static CustomMouseCallback * New ()
  {
    CustomMouseCallback *callback = new CustomMouseCallback ;
    return callback ;
  }
  virtual void Execute (vtkObject *caller, unsigned long eventId, void *callData)
  {
   vtkSmartPointer < vtkRenderWindowInteractor > interactor = dynamic_cast < vtkRenderWindowInteractor * > ( caller ) ;
    vtkSmartPointer < vtkRenderWindow > window = interactor->GetRenderWindow () ;
    vtkSmartPointer < vtkRenderer > renderer = window->GetRenderers ()->GetFirstRenderer () ;
    vtkSmartPointer < vtkCellPicker > picker = dynamic_cast < vtkCellPicker * > ( interactor->GetPicker() ) ;    // Get the clicked location
     int clickLocation[2] ;
    interactor->GetEventPosition ( clickLocation ) ;
    // Detect where the click falls on the mesh and image
       int pickedId = picker->Pick ( clickLocation[0], clickLocation[1], 0, renderer ) ;
    int  pickedPos[3] ;
    double window_pickedPos[3];
    picker->GetPickPosition ( window_pickedPos ) ; // click on the mesh
    picker->GetPointIJK ( pickedPos ) ; // click on the image

    // Create SphereSource for seed point
        vtkSmartPointer<vtkSphereSource> sphere_source = vtkSmartPointer<vtkSphereSource>::New();
    sphere_source->SetCenter(window_pickedPos[0], window_pickedPos[1], window_pickedPos[2]);
    sphere_source->SetRadius(3.0);
    sphere_source->Update();
    vtkSmartPointer<vtkPolyDataMapper> sphere_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    sphere_mapper->SetInputData(sphere_source->GetOutput());
    vtkSmartPointer<vtkActor> sphere_actor = vtkSmartPointer<vtkActor>::New();
    sphere_actor->SetMapper(sphere_mapper);
    vtkNew<vtkNamedColors> colors;
    sphere_actor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
    renderer->AddActor(sphere_actor);
    vtkSmartPointer<vtkBalloonRepresentation> balloonRep = vtkSmartPointer <vtkBalloonRepresentation>::New();
    balloonRep->SetBalloonLayoutToImageRight();
    vtkSmartPointer<vtkBalloonWidget> balloonWidget = vtkSmartPointer<vtkBalloonWidget>::New();
    balloonWidget->SetInteractor(interactor);
    balloonWidget->SetRepresentation(balloonRep);
    AtlasImageType::IndexType pos;
    pos[0] = pickedPos[0];
    pos[1] = pickedPos[1];
    pos[2] = pickedPos[2];
    if (this->AtlasImage->GetPixel(pos) == 1)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Middle cerebellar peduncle",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 2)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Pontine crossing tract (a part of MCP)",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 3)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Genu of corpus callosum",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 4)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Body of corpus callosum",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 5)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Splenium of corpus callosum",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 6)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Fornix (column and body of fornix)",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 7)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Corticospinal tract R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 8)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Corticospinal tract L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 9)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Medial lemniscus R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 10)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Medial lemniscus L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 11)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Inferior cerebellar peduncle R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 12)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Inferior cerebellar peduncle L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 13)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior cerebellar peduncle R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 14)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior cerebellar peduncle L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 15)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cerebral peduncle R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 16)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cerebral peduncle L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 17)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Anterior limb of internal capsule R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 18)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Anterior limb of internal capsule L",NULL);
    }
            if (this->AtlasImage->GetPixel(pos) == 19)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior limb of internal capsule R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 20)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior limb of internal capsule L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 21)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Retrolenticular part of internal capsule R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 22)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Retrolenticular part of internal capsule L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 23)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Anterior corona radiata R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 24)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Anterior corona radiata L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 25)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior corona radiata R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 26)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior corona radiata L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 27)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior corona radiata R",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 28)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior corona radiata L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 29)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior thalamic radiation (include optic radiation) R",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 30)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Posterior thalamic radiation (include optic radiation) L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 31)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Sagittal stratum R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 32)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Sagittal stratum L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 33)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "External capsule R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 34)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "External capsule L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 35)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cingulum (cingulate gyrus) R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 36)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cingulum (cingulate gyrus) L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 37)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cingulum (hippocampus) R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 38)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Cingulum (hippocampus) L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 39)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Fornix (cres) / Stria terminalis R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 40)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Fornix (cres) / Stria terminalis L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 41)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior longitudinal fasciculus R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 42)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior longitudinal fasciculus L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 43)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior fronto-occipital fasciculus R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 44)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Superior fronto-occipital fasciculus L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 45)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Inferior fronto-occipital fasciculus R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 46)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Inferior fronto-occipital fasciculus L",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 47)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Uncinate fasciculus R",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 48)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Uncinate fasciculus L",NULL);
    }
    if (this->AtlasImage->GetPixel(pos) == 49)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Tapetum R",NULL);
    }
        if (this->AtlasImage->GetPixel(pos) == 49)
    {
    	balloonWidget->AddBalloon(sphere_actor,
                            "Tapetum L",NULL);
    }
    
    balloonWidget->EnabledOn();
    window->Render();

/*    int iter_no = 0 ;
    this->pickedVoxel[0] = pickedPos[0] ; this->pickedVoxel[1] = pickedPos[1] ; this->pickedVoxel[2] = pickedPos[2];
    this->window_pickedVoxel[0] = window_pickedPos[0]; this->window_pickedVoxel[1] = window_pickedPos[1]; this->window_pickedVoxel[2] = window_pickedPos[2] ;*/
  }
void SetParam (  AtlasImageType::Pointer AtlasImagePassed)
  {
    this->AtlasImage = AtlasImagePassed ;
  }
  
  private:
  	AtlasImageType::Pointer AtlasImage;
  
} ;

int main(int argc, char* argv[])
{

	    //
    // PART 1: Set up typedefs and read image(s) in 
    //

    if( argc < 1 )
    {
      std::cerr << "Usage: " << std::endl ;
      std::cerr << argv[0] << " inputImageFile" << std::endl ; 
      return -1 ;
    }

    // Make data type and image type for the input diffusion image
    typedef itk::Image < double,3 > DiffusionImageType;
    
    typedef itk::Image < int, 3 > AtlasImageType ;

	
    // Make file reader for input diffusion .nrrd image
    typedef itk::ImageFileReader<DiffusionImageType> DiffusionFileReaderType;
    DiffusionFileReaderType::Pointer DiffReader = DiffusionFileReaderType::New();
    DiffReader->SetFileName(argv[1]);
    DiffReader->Update();

    DiffusionImageType::Pointer DiffusionImage = DiffReader->GetOutput() ; 
    
    typedef itk::ImageFileReader<AtlasImageType> AtlasFileReaderType;
    AtlasFileReaderType::Pointer AtlasReader = AtlasFileReaderType::New();
    AtlasReader->SetFileName(argv[2]);
    AtlasReader->Update();
    
    AtlasImageType::Pointer AtlasImage = AtlasReader->GetOutput();
  /*
  // Sphere
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetCenter(-4.0, 0.0, 0.0);
  sphereSource->SetRadius(4.0);
  vtkSmartPointer<vtkPolyDataMapper> sphereMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
  vtkSmartPointer<vtkActor> sphereActor =
    vtkSmartPointer<vtkActor>::New();
  sphereActor->SetMapper(sphereMapper);
  // Regular Polygon
  vtkSmartPointer<vtkRegularPolygonSource> regularPolygonSource =
    vtkSmartPointer<vtkRegularPolygonSource>::New();
  regularPolygonSource->SetCenter(4.0, 0.0, 0.0);
  regularPolygonSource->SetRadius(4.0);
  vtkSmartPointer<vtkPolyDataMapper> regularPolygonMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  regularPolygonMapper->SetInputConnection(regularPolygonSource->GetOutputPort());
  vtkSmartPointer<vtkActor> regularPolygonActor =
    vtkSmartPointer<vtkActor>::New();
  regularPolygonActor->SetMapper(regularPolygonMapper);
  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  // Create the widget
  vtkSmartPointer<vtkBalloonRepresentation> balloonRep =
    vtkSmartPointer<vtkBalloonRepresentation>::New();
  balloonRep->SetBalloonLayoutToImageRight();
  vtkSmartPointer<vtkBalloonWidget> balloonWidget =
    vtkSmartPointer<vtkBalloonWidget>::New();
  balloonWidget->SetInteractor(renderWindowInteractor);
  balloonWidget->SetRepresentation(balloonRep);
  balloonWidget->AddBalloon(sphereActor,
                            "This is a sphere",NULL);
  balloonWidget->AddBalloon(regularPolygonActor,
                            "This is a regular polygon",NULL);
  // Add the actors to the scene
  renderer->AddActor(sphereActor);
  renderer->AddActor(regularPolygonActor);
  // Render an image (lights and cameras are created automatically)
  renderWindow->Render();
  balloonWidget->EnabledOn();
  // Begin mouse interaction
  renderWindowInteractor->Start();
	*/
// Connect ITK portion to the VTK portion
  typedef itk::ImageToVTKImageFilter < DiffusionImageType > ITKToVTKFilterType ;
  ITKToVTKFilterType::Pointer itkToVTKfilter = ITKToVTKFilterType::New() ;
  itkToVTKfilter->SetInput ( DiffusionImage ) ;
  itkToVTKfilter->Update() ;

  // VTK Portion of the code - visualization pipeline
  // mapper
  vtkSmartPointer < vtkImageSliceMapper > imageMapper = vtkSmartPointer < vtkImageSliceMapper > ::New() ;
  imageMapper->SetInputData ( itkToVTKfilter->GetOutput() ) ;
  imageMapper->SetOrientationToX () ;
  imageMapper->SetSliceNumber ( 55 ) ;
  std::cout << "default for atfocalpoint: " << imageMapper->GetSliceAtFocalPoint () << std::endl ;
  std::cout << "default for faces camera: " << imageMapper->GetSliceFacesCamera () << std::endl ;
  imageMapper->SliceAtFocalPointOn () ;
  imageMapper->SliceFacesCameraOn () ;

  // actor 
  vtkSmartPointer < vtkImageActor > imageActor = vtkSmartPointer < vtkImageActor > ::New() ;
  imageActor->SetMapper ( imageMapper ) ;

  // set up the scene, window, interactor
  vtkSmartPointer < vtkRenderer > renderer = vtkSmartPointer < vtkRenderer >::New() ;
  renderer->AddActor ( imageActor ) ;
  
/*
 * vtkPoints x
 *
 * x->insertNextPoint(double x[3])
 *
 * vtkPointsSet y
 *
 * y->setPoints(x)
 *
 * vtkActor actors
 *
 */ 
 
     vtkPoints* newp= vtkPoints::New();
     vtkCellArray* newv= vtkCellArray::New();
     vtkPolyData* pd= vtkPolyData::New();
     pd->SetPoints( newp );
     pd->SetVerts( newv );
 
    //AtlasImageType::RegionType wholeImage;
    AtlasImageType::RegionType region = DiffusionImage->GetLargestPossibleRegion();
    AtlasImageType::SizeType size = region.GetSize();
    typedef itk::ImageRegionIterator < AtlasImageType > AtlasImageIteratorType ;
    AtlasImageIteratorType AtlasIterator ( AtlasImage, region );
    while (! AtlasIterator.IsAtEnd())
    {
     // AtlasIterator.GetIndex() ; 
     //double a = AtlasImage->GetPixel() ; 
     
     if (AtlasImage->GetPixel(AtlasIterator.GetIndex()) == 3.0)
     {
        vtkIdType id= newp->InsertNextPoint(AtlasIterator.GetIndex()[0], AtlasIterator.GetIndex()[1], AtlasIterator.GetIndex()[2]);
        newv->InsertNextCell( 1, &id );
        std::cout << "current index " << AtlasIterator.GetIndex() << std::endl ;
     	std::cout << "current pixel " << AtlasImage->GetPixel(AtlasIterator.GetIndex()) << std::endl ;
        
     }
     ++AtlasIterator ; 
    }

   vtkPolyDataMapper* my_mapper= vtkPolyDataMapper::New();
   vtkActor* my_actor= vtkActor::New();
   my_mapper->SetInputData( pd );
   my_actor->SetMapper( my_mapper );

    /*AtlasImageType::IndexType corner;
    corner[0] = 0;
    corner[1] = 0;
    corner[2] = 0;
    //wholeImage.SetSize(size);
   // wholeImage.SetIndex(corner);
  
   // make sure that the PrinciplEigeniamge has the same parameters as the DiffusionImage!!!
    ParcelImage->SetOrigin(AtlasImage->GetOrigin());
    ParcelImage->SetDirection(AtlasImage->GetDirection());
    ParcelImage->SetSpacing(AtlasImage->GetSpacing());
    ParcelImage->SetRegions(wholeImage);
    ParcelEigenImage->Allocate();
    
    typedef itk::ImageRegionIterator <AtlasImageType> InIteratorType;
    typedef itk::ImageRegionIterator <AtlasImageType> OutIteratorType;
    InIteratorType inIterator (AtlasImage, wholeImage);
    OutIteratorType outIterator(ParcelImage, wholeImage);
    inIterator.GoToBegin();
    outIterator.GoToBegin();
    DiffusionTensorType currentTensor;
    VectorType currentVector;
    DiffusionTensorType::EigenValuesArrayType eigenArray;
    DiffusionTensorType::EigenVectorsMatrixType eigenMatrix;
    while (! outIterator.IsAtEnd())
    {
        currentTensor = inIterator.Value();
        
        outIterator.Set(0);
        ++outIterator;
        ++inIterator;
    }
  
  
  */
  
  /*
    // vtk data set to fill
	//   vtkNew<vtkPolyData> myPolyData;
	vtkSmartPointer<vtkPoints> myvtkPoints = vtkSmartPointer <vtkPoints>::New();
	vtkSmartPointer < vtkIntArray > intArray = vtkSmartPointer < vtkIntArray >::New() ;
  parcelArray->SetNumberOfTuples ( 50 ) ;
  for ( unsigned int i = 0 ; i < 50 ; i++ )
    {
      parcelArray->SetValue ( i, 1 ) ;
    }
	
	myvtkPoints->SetData(parcelArray);
   // compute a vtkIdList that contains ids of each points
  vtkIdType numPoints = myvtkPoints->GetNumberOfPoints();
  vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
  pointIds->SetNumberOfIds(numPoints);
  for (vtkIdType i = 0; i < numPoints; ++i)
  {
    pointIds->SetId(i, i);
  }
   // create a vtkCellArray from this list
  vtkSmartPointer<vtkCellArray> polyPoint = vtkSmartPointer<vtkCellArray>::New();
  polyPoint->InsertNextCell(pointIds);
   // give the points and cells to the final poly data
  myPolyData->SetPoints(myVtkPoints);
  myPolyData->SetVerts(polyPoint);
  
  */
  
    vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetCenter(55, 229, 11);
  sphereSource->SetRadius(10.0);


  vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

  vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
  sphereActor->SetMapper(sphereMapper);
  sphereActor->GetProperty()->SetOpacity(0);

	renderer->AddActor(sphereActor);


    
  // get the camera so we can position it better
  vtkSmartPointer < vtkCamera > camera = renderer->GetActiveCamera() ;

  double position[3],  imageCenter[3] ;
  itkToVTKfilter->GetOutput()->GetCenter ( imageCenter ) ;
  position[0] = imageCenter[0] ;
  position[1] = imageCenter[1] ;
  position[2] = -160 ;
  std::cout << "Image center: " << imageCenter[0] << " " << imageCenter[1] << " " << imageCenter[2] << std::endl ;
  double spacing[3] ;
  int imageDims[3] ;
  itkToVTKfilter->GetOutput()->GetSpacing ( spacing ) ;
  itkToVTKfilter->GetOutput()->GetDimensions ( imageDims ) ;
  double imagePhysicalSize[3] ;
  for ( unsigned int d = 0 ; d < 3 ; d++ )
    {
      imagePhysicalSize[d] = spacing[d] * imageDims[d] ;
    }


  camera->ParallelProjectionOn () ; 
  camera->SetFocalPoint ( imageCenter ) ;
  camera->SetPosition ( position ) ;
  //  camera->SetDistance ( imagePhysicalSize[2] * -1 ) ;
  std::cout << "Parallel scale: " << camera->GetParallelScale() << std::endl ;
  std::cout << imageDims[0] << " " << imageDims[1] << " " << imageDims[2] << std::endl ;
  camera->SetParallelScale ( imageDims[2] * 1.2 ) ;

  // set up window
  vtkSmartPointer < vtkRenderWindow > window = vtkSmartPointer < vtkRenderWindow >::New() ;
  window->AddRenderer ( renderer ) ;
  window->SetSize ( 500, 500 ) ;

  vtkSmartPointer < vtkRenderWindowInteractor > interactor = vtkSmartPointer < vtkRenderWindowInteractor >::New() ;
  interactor->SetRenderWindow ( window ) ;

  vtkSmartPointer < vtkInteractorStyleImage > style = vtkSmartPointer < vtkInteractorStyleImage >::New() ;
  //style->SetInteractionModeToImage3D() ;
  style->SetInteractionModeToImageSlicing() ;

  interactor->SetInteractorStyle ( style ) ;
  interactor->Initialize() ;
  interactor->CreateRepeatingTimer ( 300 ) ;
   
  // pk mouse stuff
  vtkSmartPointer < CustomMouseCallback > myMouseCallback = vtkSmartPointer < CustomMouseCallback >::New() ;
  interactor->AddObserver ( vtkCommand::LeftButtonPressEvent , myMouseCallback, 0 ) ;
  myMouseCallback->SetParam ( AtlasImage );
  vtkSmartPointer < vtkCellPicker > picker = vtkSmartPointer < vtkCellPicker >::New() ;
  interactor->SetPicker ( picker ) ;
  //vtkSmartPointer < CustomTimerCallback > myCallback = vtkSmartPointer < CustomTimerCallback >::New() ;
  //int timerId = interactor->AddObserver ( vtkCommand::TimerEvent, myCallback, 0 ) ;
  //myCallback->SetTimerId ( timerId ) ;
  //  interactor->DestroyTimer ( timerId ) ;
  // run!
  // Create the widget
  /*renderer->AddActor(my_actor);
  vtkSmartPointer<vtkBalloonRepresentation> balloonRep = vtkSmartPointer <vtkBalloonRepresentation>::New();
  balloonRep->SetBalloonLayoutToImageRight();
  vtkSmartPointer<vtkBalloonWidget> balloonWidget = vtkSmartPointer<vtkBalloonWidget>::New();
  balloonWidget->SetInteractor(interactor);
  balloonWidget->SetRepresentation(balloonRep);
  balloonWidget->AddBalloon(my_actor,
                            "This is a intensity 20",NULL);
  balloonWidget->EnabledOn(); // must come before interactor->start and interactor->render*/
  interactor->Start() ;
	interactor->Render();

  // Done.
  return 0 ;
}
