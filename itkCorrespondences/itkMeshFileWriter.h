#ifndef MESHFILEWRITER_H_HEADER_INCLUDED
#define MESHFILEWRITER_H_HEADER_INCLUDED

#include <itkProcessObject.h>

namespace itk {


/** \class MeshFileWriter
  *   \brief Base class for mesh writers in different formats.
  *
  * The class offers a Write() method which calls GenerateData() to write the 
  * mesh to disk. In addition, there are several methods to get and set 
  * FilePrefix, Filename and FilePattern. The subclasses can decide which of 
  * these naming methods they want to use.
  *
  * \author Tobias Heimann. Division Medical and Biological Informatics, 
  *         German Cancer Research Center, Heidelberg, Germany.
  */
template <class TMeshType>
class MeshFileWriter : public ProcessObject
{
  public:

    typedef MeshFileWriter<TMeshType> Self;
    typedef ProcessObject Superclass;
    typedef SmartPointer < Self > Pointer;
    typedef SmartPointer < const Self > ConstPointer;

    itkTypeMacro( MeshFileWriter<TMeshType>, ProcessObject );


        //##Documentation
    //## @brief Get the specified the file to write
    //## 
    //## Either the FileName or FilePrefix plus FilePattern are used to write.
    virtual const char* GetFileName() const { return m_fileName; }

    //##Documentation
    //## @brief Specify the file to write.
    //## 
    //## Either the FileName or FilePrefix plus FilePattern are used to write.
    virtual void SetFileName(const char* aFileName) { strcpy( m_fileName, aFileName ); }

    //##Documentation
    //## @brief Get the specified file prefix for the file(s) to write. 
    //## 
    //## You should specify either a FileName or FilePrefix. Use FilePrefix if
    //## the data is stored in multiple files.
    virtual const char* GetFilePrefix() const { return m_filePrefix; }

    //##Documentation
    //## @brief Specify file prefix for the file(s) to write.
    //## 
    //## You should specify either a FileName or FilePrefix. Use FilePrefix if
    //## the data is stored in multiple files.
    virtual void SetFilePrefix(const char* aFilePrefix) { strcpy( m_filePrefix, aFilePrefix ); };

    //##Documentation
    //## @brief Get the specified file pattern for the file(s) to write. The
    //## sprintf format used to build filename from FilePrefix and number.
    //## 
    //## You should specify either a FileName or FilePrefix. Use FilePrefix if
    //## the data is stored in multiple files.
    virtual const char* GetFilePattern() const { return m_filePattern; }

    //##Documentation
    //## @brief Specified file pattern for the file(s) to write. The sprintf
    //## format used to build filename from FilePrefix and number.
    //## 
    //## You should specify either a FileName or FilePrefix. Use FilePrefix if
    //## the data is stored in multiple files.
    virtual void SetFilePattern(const char* aFilePattern) { strcpy( m_filePattern, aFilePattern ); };


    virtual DataObject* GetInput() 
    {
      return ProcessObject::GetInput( 0 );
    }

    virtual void SetInput( DataObject *input )
    {
      SetNthInput( 0, input );
    }

    virtual void Write()
    {
      if ( this->GetInput() == NULL ) {
        itkExceptionMacro(<<"Write:Please specify an input!");
        return;
      }
      this->UpdateOutputInformation();
      (*this->GetInputs().begin())->SetRequestedRegionToLargestPossibleRegion();
      this->PropagateRequestedRegion(NULL);
      this->UpdateOutputData(NULL);
    }

    virtual void Update()
    {
      Write();
    }
    
    virtual void GenerateData() = 0;

  protected:

    MeshFileWriter() 
    {
      m_fileName[0] = 0;
      m_filePrefix[0] = 0; 
      m_filePattern[0] = 0; 
    };

    virtual ~MeshFileWriter() {};

    char  m_fileName[1024];       // filename for mesh data
    char  m_filePrefix[1024];     // fileprefix for mesh data
    char  m_filePattern[1024];    // filepattern for mesh data

};
                                                   

} // namespace itk
#endif /* MESHFILEWRITER_H_HEADER_INCLUDED */
