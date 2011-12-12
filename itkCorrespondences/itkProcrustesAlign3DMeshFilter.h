#ifndef _itkProcrustesAlign3DMeshFilter_h
#define _itkProcrustesAlign3DMeshFilter_h

#include "itkProcessObject.h"
#include "itkScalableAffineTransform.h"
#include "itkTransformMeshFilter.h"


namespace itk
{
  /** \class ProcrustesAlign3DMeshFilter
  *   \brief Class for groupwise aligning a set of 3D meshes.
  * 
  * All input meshes must have the same number of points (with corresponding 
  * indices).
  * Default mode: Iterative Generalized Procrustes alignment with initialization
  * from first object
  * Option: UseInitialAverageOn/Off() Use  average structure as initialization (off by default)
  *         Only appropriate if objects are already pre-aligned 
  * Option: UseSingleIterationOn/Off() Only run one iteration (off by default) 
  * 
  * All output meshes will be centered at the origin and scaled to 
  * match a mean shape with norm 1.
  * Option: UseScalingOn/Off() enables/disables scaling normalization (on by default)
  *
  * GetTransform() can be used to query the employed transformation for each 
  * input mesh.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  * change  Martin Styner, UNC support for single template and initialization with average
  */
  template <class TInputMesh, class TOutputMesh>
  class ProcrustesAlign3DMeshFilter : public ProcessObject
  {

  public:

    /** Standard typedefs. */
    typedef ProcrustesAlign3DMeshFilter   Self;
    typedef ProcessObject                 Superclass;
    typedef SmartPointer<Self>            Pointer;
    typedef SmartPointer<const Self>      ConstPointer;

    /** Convenient typedefs. */
    typedef TInputMesh                                  InputMeshType;
    typedef TOutputMesh                                 OutputMeshType;
    typedef typename InputMeshType::Pointer             InputMeshPointer;
    typedef typename OutputMeshType::Pointer            OutputMeshPointer;
    typedef typename InputMeshType::PointType           InputPointType;
    typedef typename OutputMeshType::PointType          OutputPointType;
    typedef typename OutputMeshType::PointsContainer    OutputPointsContainer;
    typedef DataObject::Pointer                         DataObjectPointer;
    typedef typename OutputMeshType::CoordRepType       CoordRepType;
    typedef vnl_matrix<CoordRepType>                    MatrixType;
    typedef ScalableAffineTransform<CoordRepType, 3>    TransformType;
    typedef typename TransformType::Pointer             TransformPointer;
    typedef TransformMeshFilter
      <OutputMeshType, OutputMeshType, TransformType>   TransformMeshType;
    typedef typename TransformMeshType::Pointer         TransformMeshPointer;
    typedef std::vector<TransformMeshPointer>           TransformMeshArray;
    typedef typename TransformType::OffsetType          TranslationType;
    typedef typename TransformType::MatrixType          RotationType;        
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ProcrustesAlign3DMeshFilter, ProcessObject);

    /** Sets the number of input meshes that have to be aligned.
    * Call this before any other methods.
    */
    void SetNumberOfInputs( unsigned int num );

    /** Sets one input mesh (starting with idx=0). */
    void SetInput( unsigned int idx, InputMeshPointer mesh );

    /** Returns one of the input meshes (starting with idx=0). */
    InputMeshType* GetInput( unsigned int idx )
    {
      return static_cast<InputMeshType*>(this->ProcessObject::GetInput( idx ));
    }

    /** Gets one transformed output mesh (starting with idx=0). */
    OutputMeshType* GetOutput(unsigned int idx);

    /** Gets the transformed mean of all aligned meshes. */
    itkGetConstObjectMacro( Mean, OutputMeshType );

    /** Returns the transform used to align a mesh (starting with idx=0). */
    TransformType* GetTransform( unsigned int idx )
    {
      return m_MeshTransform[idx]->GetTransform();
    } 

    /** Get/Set the convergence value that determines when the iterative
    * calculation of alignment is stopped. The smaller the value, the higher
    * the accuracy (the default should be sufficient for all applications).
    */
    itkGetConstMacro( Convergence, double );
    itkSetMacro (Convergence, double );

    /** Scaling normalization on/off.
    * Default is on.
    */
    itkGetConstMacro( UseScaling, bool );
    itkSetMacro( UseScaling, bool );
    itkBooleanMacro( UseScaling );

    /** Determines if an average of all shapes is used for initialization
    * or just the first shape.
    * Default is off (use just the first shape).
    */
    itkGetConstMacro( UseInitialAverage, bool );
    itkSetMacro( UseInitialAverage, bool );
    itkBooleanMacro( UseInitialAverage );

    /** Only do one iteration.
    * Default is off.
    */
    itkGetConstMacro( UseSingleIteration, bool );
    itkSetMacro( UseSingleIteration, bool );
    itkBooleanMacro( UseSingleIteration );

    /** Determines if meshes are translated individually for best alignment. 
    * Default is on.
    */
    itkGetConstMacro( AlignTranslation, bool );
    itkSetMacro( AlignTranslation, bool );
    itkBooleanMacro( AlignTranslation );

    /** Determines if meshes are scaled individually for best alignment. 
    * Default is on.
    */
    itkGetConstMacro( AlignScale, bool );
    itkSetMacro( AlignScale, bool );
    itkBooleanMacro( AlignScale );

    /** Determines if meshes are rotated for best alignment. 
    * Default is on.
    */
    itkGetConstMacro( AlignRotation, bool );
    itkSetMacro( AlignRotation, bool );
    itkBooleanMacro( AlignRotation );

    /** Gets the mean scale factor that has been used to transform the meshes. */
    itkGetConstMacro( MeanScale, double );

    /** Creates an ouput object. */
    virtual DataObjectPointer MakeOutput(unsigned int idx);


  protected:

    /** Standard constructor. */
    ProcrustesAlign3DMeshFilter();

    /** Standard destructor. */
    ~ProcrustesAlign3DMeshFilter();

    /** Performs the alignment. */
    virtual void GenerateData();

    /** Calculates the center coordinates of the specified mesh. */
    TranslationType GetMeshCenter( InputMeshPointer mesh );

    /** Uses the current transformations to calculate a new mean.*/
    void CalculateMean();

    /** Returns the best transform to fit input mesh idx, translated to  
     * zero origin by using the values in m_Center, to the given targetMesh. 
     * targetMesh has to be centered at zero origin as well!
     */
    TransformPointer GetProcrustesMatch( unsigned int idx, OutputMeshPointer targetMesh );

    
  private:

    /** When the difference between two consecutive mean calculations gets
    * samller than m_Convergence, the alignment is terminated. */
    double                        m_Convergence;
    /** Holds the transforms for all input meshes.*/
    TransformMeshArray            m_MeshTransform;
    /** Holds the center coordinates for all input meshes.*/
    std::vector<TranslationType>  m_Center;
    /** The mean of all transformed meshes. */
    OutputMeshPointer             m_Mean;
    /** The mean of all transformed meshes from the preceding iteration. */
    OutputMeshPointer             m_OldMean;
    /** Determines if different meshes are translated for best fit. */
    bool                          m_AlignTranslation;
    /** Determines if different meshes are resized for best fit. */
    bool                          m_AlignScale;
    /** Determines if different meshes are rotated for best fit. */
    bool                          m_AlignRotation;
    /** Scaling normalization on/off */
    bool                          m_UseScaling;
    /** Determines if the first iteration starts with an average or with a single shape. */
    bool                          m_UseInitialAverage;
    /** Only run one single iteration on/off */
    bool                          m_UseSingleIteration;
    /** The mean scale factor for all meshes in the final alignment. */
    double                        m_MeanScale;
    /** The mean center of all meshes. */
    TranslationType               m_MeanCenter;
  };

  
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProcrustesAlign3DMeshFilter.txx"
#endif

#endif
