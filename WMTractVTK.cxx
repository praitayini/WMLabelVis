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

///nfs/masi/kanakap/open_srx

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
    typedef itk::Image<double,3> DiffusionImageType;
    
    typedef itk::Image < double, 3 > AtlasImageType ;

	
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
  
  
  
    //AtlasImageType::RegionType wholeImage;
    AtlasImageType::RegionType region = DiffusionImage->GetLargestPossibleRegion();
    AtlasImageType::SizeType size = region.GetSize();
    typedef itk::ImageRegionIterator < AtlasImageType > AtlasImageIteratorType ;
    AtlasImageIteratorType AtlasIterator ( AtlasImage, region );
    while (! AtlasIterator.IsAtEnd())
    {
     // AtlasIterator.GetIndex() ; 
     //double a = AtlasImage->GetPixel() ; 
     
     if (AtlasImage->GetPixel(AtlasIterator.GetIndex()) == 20.0)
     {
        std::cout << "current index " << AtlasIterator.GetIndex() << std::endl ;
     	std::cout << "current pixel " << AtlasImage->GetPixel(AtlasIterator.GetIndex()) << std::endl ;
        
     }
     ++AtlasIterator ; 
    }

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
  sphereSource->SetRadius(4.0);


  vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

  vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
  sphereActor->SetMapper(sphereMapper);

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

  //vtkSmartPointer < CustomTimerCallback > myCallback = vtkSmartPointer < CustomTimerCallback >::New() ;
  //int timerId = interactor->AddObserver ( vtkCommand::TimerEvent, myCallback, 0 ) ;
  //myCallback->SetTimerId ( timerId ) ;
  //  interactor->DestroyTimer ( timerId ) ;
  // run!
  // Create the widget
  vtkSmartPointer<vtkBalloonRepresentation> balloonRep = vtkSmartPointer <vtkBalloonRepresentation>::New();
  balloonRep->SetBalloonLayoutToImageRight();

  vtkSmartPointer<vtkBalloonWidget> balloonWidget = vtkSmartPointer<vtkBalloonWidget>::New();
  balloonWidget->SetInteractor(interactor);
  balloonWidget->SetRepresentation(balloonRep);
  balloonWidget->AddBalloon(sphereActor,
                            "This is a sphere",NULL);
  balloonWidget->EnabledOn(); // must come before interactor->start and interactor->render
  interactor->Start() ;
	interactor->Render();

  // Done.
  return 0 ;
}
	
	

