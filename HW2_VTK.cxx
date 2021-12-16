#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "vtkSmartPointer.h"
#include "vtkImageSliceMapper.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkImageProperty.h"
#include "vtkInteractorStyleImage.h"
#include "vtkCommand.h"
#include "vtkCellData.h"
#include "vtkPolyDataMapper.h"
#include "vtkUnsignedCharArray.h"
#include "vtkRendererCollection.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"
#include "vtkCellPicker.h"
#include "vtkPointData.h"
#include "vtkDecimatePro.h"
#include "vtkSphereSource.h"
#include "vtkProperty.h"
#include "vtkVector.h"
#include "vtkPolyDataMapper.h"
#include "vtkLine.h"
#include "vtkNew.h"
#include "vtkLineSource.h"
#include "vtkNamedColors.h"
#include <vtkBalloonRepresentation.h>
#include <vtkBalloonWidget.h>

const unsigned int nDims = 3 ;
typedef itk::Vector < double, nDims > vectorType;
typedef itk::Image < vectorType, nDims > pevImageType ;
typedef itk::Image < double, nDims > faImageType ;
typedef itk::Image < double, nDims > wmtractImageType ;

class CustomTimerCallback : public vtkCommand
{
public:

    static CustomTimerCallback * New () 
  {
    CustomTimerCallback *callback = new CustomTimerCallback ;
    return callback ;
  }

  virtual void Execute (vtkObject *caller, unsigned long eventId, void *callData)
  {
    vtkSmartPointer < vtkRenderWindowInteractor > interactor = dynamic_cast < vtkRenderWindowInteractor * > ( caller ) ;
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer < vtkRenderWindow > window = interactor->GetRenderWindow () ;
    vtkSmartPointer < vtkRenderer > renderer = window->GetRenderers ()->GetFirstRenderer () ;
    if ( this->m_counter > 100000 ) 
      {
       interactor->DestroyTimer ( this->m_timerId ) ;
	std::cout << "destroyed timer" << std::endl ;
      }
    else
      {
        if (this->voxelList.size() > 0 ) {
  	   faImageType::IndexType ret[3] ;
           int ret_segtract= segmentTract(this->voxelList.front(), wmtractImage, faImage, pevImage, &iteration_no, ret, this->window_voxelList.front());
           this->voxelList.pop_front() ;
           this->window_voxelList.pop_front() ;
           if (ret_segtract == 0) {
               plotLineSource( ret[0], ret[1] , mapper, actor, renderer, window ) ;
               plotLineSource( ret[0], ret[2] , mapper, actor, renderer, window ) ;        
           }
        }
      }
    this->m_counter++ ;
  }
 
// Function to plot a line between current point and next point on tract
 int plotLineSource( faImageType::IndexType  ret0 , faImageType::IndexType ret1 , vtkSmartPointer<vtkPolyDataMapper> mapper, vtkSmartPointer<vtkActor> actor,  vtkSmartPointer < vtkRenderer > renderer, vtkSmartPointer < vtkRenderWindow > window)
 {
  vtkNew<vtkLineSource> lineSource;
  vtkNew<vtkNamedColors> colors;
  double p0[3] = {ret0[0], ret0[1], ret0[2]}; 
  double p1[3] = {ret1[0], ret1[1], ret1[2]}; 
  lineSource->SetPoint1(p0);
  lineSource->SetPoint2(p1);
  mapper->SetInputConnection(lineSource->GetOutputPort());
  actor->SetMapper(mapper);
  actor->GetProperty()->SetLineWidth(4);
  actor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
  renderer->AddActor(actor);
  window->Render( );
}

// Function to trace the tracts
 int segmentTract( faImageType::IndexType cVoxel, wmtractImageType::Pointer  wmtractImage, faImageType::Pointer faImage, pevImageType::Pointer pevImage, int *iteration_no, faImageType::IndexType *ret, faImageType::IndexType window_pickedVoxel )
{
   // Termination criteria 
   if (*iteration_no == 100000)
   {  
return 1;}
   else if (!wmtractImage->GetLargestPossibleRegion().IsInside(cVoxel)) 
   {
  return 1;} 
   else if ( faImage->GetPixel( cVoxel ) < 0.2)
   {
 return 1;}
   else if ( wmtractImage->GetPixel( cVoxel ) == 1.0 ) 
   {
 return 1;}
   else
   {

   // At the voxel set the voxel to iteration number
   wmtractImage->SetPixel( cVoxel, 1.0 );

   // Eigen vector at the voxel d[x]
   vectorType dx = pevImage->GetPixel( cVoxel ); 

   // Compute next voxel to visit
   faImageType::IndexType AB ;
   faImageType::IndexType BA ;
  for (int j = 0; j < 3; j++) 
   {
    AB[j] = round(cVoxel[j] + dx[j] * 0.8)  ;
    BA[j] = round(cVoxel[j] - dx[j] * 0.8) ;
   }
   *iteration_no += 1; 
    
   // Compute next window position to visit 
   faImageType::IndexType window_AB ;
   faImageType::IndexType window_BA ;
   for (int k = 0; k < 3; k++)
   {
    window_AB[k] = window_pickedVoxel[k] + dx[k] * 2; 
    window_BA[k] = window_pickedVoxel[k] - dx[k] * 2; 
   }

   // Add backward and forward voxel to the list
   this->voxelList.push_back(AB);
   this->voxelList.push_back(BA);

   // Keep track of the window (mesh) postion
   this->window_voxelList.push_back(window_AB);
   this->window_voxelList.push_back(window_BA);

   // Return postion on the window
   ret[0] = window_pickedVoxel;
   ret[1] = window_AB;
   ret[2] = window_BA;
   
   return 0;
   }
}
 
  void SetTimerId ( int id ) 
  {
    this->m_timerId = id ;
    this->m_counter = 0 ;
  }

  void SetSegmentTractParams (  faImageType::IndexType pickedVoxel, wmtractImageType::Pointer  wmtractImage, faImageType::Pointer faImage, pevImageType::Pointer pevImage, int iteration_no , faImageType::IndexType window_pickedVoxel)
  {
    this->voxelList.push_back(pickedVoxel) ;
    this->window_voxelList.push_back(window_pickedVoxel) ;
    this->wmtractImage = wmtractImage ;
    this->faImage = faImage ;
    this->pevImage = pevImage ;
    this->iteration_no = iteration_no;
  }
  
private: 
  int m_timerId ;
  int m_counter ;
  std::list<faImageType::IndexType> voxelList;
  std::list<faImageType::IndexType> window_voxelList;
  wmtractImageType::Pointer wmtractImage ;
  faImageType::Pointer faImage ;
  pevImageType::Pointer pevImage ;
  int iteration_no ;
} ;

class CustomMouseCallback : public vtkCommand
{
public:
  // To use in the CustomTimerCallback
  faImageType::IndexType pickedVoxel ;
  faImageType::IndexType window_pickedVoxel ;

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
    vtkSmartPointer < vtkCellPicker > picker = dynamic_cast < vtkCellPicker * > ( interactor->GetPicker() ) ;
    // Get the clicked location
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
    balloonWidget->AddBalloon(sphere_actor,
                            "This is a intensity 20",NULL);
    balloonWidget->EnabledOn();
    window->Render();

    int iter_no = 0 ;
    this->pickedVoxel[0] = pickedPos[0] ; this->pickedVoxel[1] = pickedPos[1] ; this->pickedVoxel[2] = pickedPos[2];
    //this->pickedVoxel[0] = 69 ; this->pickedVoxel[1] = 47 ; this->pickedVoxel[2] = 45 ; //testing
    this->window_pickedVoxel[0] = window_pickedPos[0]; this->window_pickedVoxel[1] = window_pickedPos[1]; this->window_pickedVoxel[2] = window_pickedPos[2] ; 
  }

} ;


int main ( int argc, char * argv[] )
{
 // Verify command line arguments
 if( argc < 2 )
    {
      std::cerr << "Usage: " << std::endl ;
      std::cerr << argv[0] << " inputPEVImageFile inputFAImageFile" << std::endl ; 
      return -1 ;
    }
 

 // ITK Portion of the code - reading the file in
 // Setup types
 typedef itk::ImageFileReader < pevImageType > ImageReaderType ;
 ImageReaderType::Pointer myReader = ImageReaderType::New() ;
 myReader->SetFileName ( argv[1] ) ;
 pevImageType::Pointer myITKImage = myReader->GetOutput() ;

 // Connect ITK portion to VTK portion
 typedef itk::ImageToVTKImageFilter < pevImageType > ITKToVTKFilterType ;
 ITKToVTKFilterType::Pointer itkToVTKfilter = ITKToVTKFilterType::New() ;
 itkToVTKfilter->SetInput ( myITKImage ) ;
 itkToVTKfilter->Update() ;
  
  // Read the FA image in too
  typedef itk::ImageFileReader < faImageType > faImageReaderType ;
  faImageReaderType::Pointer myfaReader = faImageReaderType::New() ;
  myfaReader->SetFileName ( argv[2] ) ;
  myfaReader->Update() ;
  faImageType::Pointer myfaITKImage = myfaReader->GetOutput() ;
  faImageType::RegionType myRegion = myfaITKImage->GetLargestPossibleRegion() ;
  faImageType::SizeType mySize = myRegion.GetSize() ;
  faImageType::IndexType myStart ;
  myStart[0] =   0;  // first index on X
  myStart[1] =   0;  // first index on Y
  myStart[2] =   0;  // first index on Z
  faImageType::SpacingType mySpacing = myfaITKImage->GetSpacing() ;
  faImageType::PointType myOrigin = myfaITKImage->GetOrigin() ;
  faImageType::DirectionType myDirection = myfaITKImage->GetDirection() ;

  typedef itk::ImageToVTKImageFilter < faImageType > faITKToVTKFilterType ;
  faITKToVTKFilterType::Pointer faitkToVTKfilter = faITKToVTKFilterType::New() ;
  faitkToVTKfilter->SetInput ( myfaITKImage ) ;
  faitkToVTKfilter->Update() ;

  // WM tract image
  wmtractImageType::Pointer wmtractImage = wmtractImageType::New() ;
  wmtractImage->SetRegions( myRegion ) ;
  wmtractImage->SetOrigin( myOrigin ) ;
  wmtractImage->SetSpacing( mySpacing );
  wmtractImage->SetDirection( myDirection );
  wmtractImage->Allocate() ;

  // VTK Portion of the code - visualization pipeline
  // Mapper
  vtkSmartPointer < vtkImageSliceMapper > imageMapper = vtkSmartPointer < vtkImageSliceMapper > ::New() ;
  imageMapper->SetInputData ( itkToVTKfilter->GetOutput() ) ;
  imageMapper->SetOrientationToX () ;
  imageMapper->SetSliceNumber ( 69 ) ;
  imageMapper->SliceAtFocalPointOn () ;
  imageMapper->SliceFacesCameraOn () ;
   
  // Image property 
  vtkSmartPointer < vtkImageProperty > image_property = vtkSmartPointer <vtkImageProperty>::New() ;
  image_property->SetColorWindow(1.0) ;
  image_property->SetColorLevel(0.5) ;

  // Actor 
  vtkSmartPointer < vtkImageActor > imageActor = vtkSmartPointer < vtkImageActor > ::New() ;
  imageActor->SetMapper ( imageMapper ) ;
  imageActor->SetProperty( image_property ) ; 

  // Set up the scene, window, interactor
  vtkSmartPointer < vtkRenderer > renderer = vtkSmartPointer < vtkRenderer >::New() ;
  renderer->AddActor ( imageActor ) ;

  // Get the camera so we can position it better
  vtkSmartPointer < vtkCamera > camera = renderer->GetActiveCamera() ; 
  double position[3],  imageCenter[3] ;
  itkToVTKfilter->GetOutput()->GetCenter ( imageCenter ) ;
  position[0] = imageCenter[0] ;
  position[1] = imageCenter[1] ;
  position[2] = -160 ;
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
  camera->SetParallelScale ( imageDims[2] * 1.3) ; // custom to this PEV image

  // Set up window
  vtkSmartPointer < vtkRenderWindow > window = vtkSmartPointer < vtkRenderWindow >::New() ;
  window->AddRenderer ( renderer ) ;
  window->SetSize ( 500, 500 ) ;
  
  // Create the interactor
  vtkSmartPointer < vtkRenderWindowInteractor > interactor = vtkSmartPointer < vtkRenderWindowInteractor >::New() ;
  interactor->SetRenderWindow ( window ) ;  

  vtkSmartPointer < vtkInteractorStyleImage > style = vtkSmartPointer < vtkInteractorStyleImage >::New() ;
  //style->SetInteractionModeToImage3D() ; // for 3D 
  style->SetInteractionModeToImageSlicing() ;
  
  interactor->SetInteractorStyle ( style ) ;
  interactor->Initialize() ;

  interactor->CreateRepeatingTimer( 0.1 ); 

  // create callback and add it to the interactor
  vtkSmartPointer < CustomMouseCallback > myMouseCallback = vtkSmartPointer < CustomMouseCallback >::New() ;
  interactor->AddObserver ( vtkCommand::LeftButtonPressEvent , myMouseCallback, 0 ) ;

  vtkSmartPointer < vtkCellPicker > picker = vtkSmartPointer < vtkCellPicker >::New() ;
  interactor->SetPicker ( picker ) ;
  
  // Run
  interactor->Start() ;
  std::cout << "Pick a starting seed point" <<  std::endl ;

  //pkvtkSmartPointer < CustomTimerCallback > myCallback = vtkSmartPointer < CustomTimerCallback >::New() ;
  //pkint timerId = interactor->AddObserver ( vtkCommand::TimerEvent, myCallback, 0) ;
  //pkmyCallback->SetTimerId ( timerId ) ;
  //pkint iter_no = 0 ;
  //pkmyCallback->SetSegmentTractParams ( myMouseCallback->pickedVoxel, wmtractImage, myfaITKImage, myITKImage, iter_no, myMouseCallback->window_pickedVoxel ) ;
  
  // Run
  //pkstd::cout << "Generating WM tracts from the selected seed point" <<  std::endl ;
  //pkinteractor->Start() ;
  //interactor->ReInitialize();

  // Done
 return 0;
}
