#ifndef READWRITEMS_H
#define READWRITEMS_H

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <istream>
#include <vector>
#include <string>
#include <array>

class Alignment;
class BaseFormatHandler;

/**
 \brief Class to handle \link BaseFormatHandler Format Handlers \endlink.\n
    It serves as a proxy to the handlers, so the code outside the FormatManager
        is format-agnostic.
    \note The Format Handlers are added automatically by CMake, on CMake configuration.\n
    This is achieved by creating the file include/FormatHandling/formats_header.h ,
        which implements the FormatManager constructor.
 */
class FormatManager
{
public:
    FormatManager();
    ~FormatManager();
    
private:
    /**
     \brief Vector that contains the available formats to load/save from.\n
            They are loaded into the format in the constructor function.
     \note The formats (\link BaseFormatHandler \endlink) are hard-coded loaded, 
        and thus so when a new format is implemented, it has to be added 'manually'.
        There is a workaround for this, done on CMake, which creates the constructor of the class.
        This on-the-fly constructor implements all states found, using some rules.\n
<br>
\note 
This file includes all States found, and also defines the ReadWriteMS constructor.\n
To be able to be automatically recognized, the new state should:\n
       -# Have the same Class Name as File Name (without the extension)\n
       -# Name must end with '_state'\n
       -# Be placed on ReadWriteMS folder\n
     */
    std::vector<BaseFormatHandler*> available_states;
    /**
     \brief Function that adds a newState to the
             \link FormatManager::available_states \endlink vector.\n
            This should be called on the constructor function foreach format existent.
     \param newState Pointer to the newState we want to instantiate.
     */
    void addState(BaseFormatHandler* newState);
    
public:
    
    /** \brief Tag to know if the machine
        has an output file or it has to output to console. */
    bool hasOutputFile  = true;
    /** \brief Tag to know if the machine should keep original headers.*/
    bool keepHeader     = false;
    /** \brief Tag to know if sequences should be reversed before saving them.*/
    bool reverse        = false;
    
    // LEGACY PARAMETERS
    /** \brief Tag to know if the machine should output
        the format information about the alignment.*/
    bool format         = false;
    /** \brief Tag to know if the machine should output
        the type of the alignment.*/
    bool type           = false;
    /** \brief Tag to know if the machine should output
        the information of the alignment.*/
    bool info           = false;
    
    /**
    \brief Function that loads an alignment given a file path.
        It automatically detects the format of the file.
    \param inFile File path of the alignment to load.
    \return <b>Pointer to the alignment</b> if it could be loaded.\n
            <b>Null</b> if not.
        */
    Alignment* loadAlignment(std::string inFile);

    /**
    \brief Function to save an alignment to a file.
        It searches among the available_states one that can write the alignment
         in the specified format.
    \param outPattern File path to save the alignment.
    \param outFormats Format in which save the alignment.
    \param alignment Alignment
            */
    bool saveAlignment(std::string outPattern, std::vector< std::string >* outFormats, Alignment* alignment);
    
    /**
     \brief Function that takes multiple files,
             loads them and saves in a cumulus of formats, using an outPattern.
     \param inFile Vector of files to load, reformat and save.
     \param outPattern Path and name of the new files.\n
               The function changes some optional tokens on the original string to obtain multiple versions: \n
                - <b> [in] </b>        Token that is changed with the original filename without extension.\n
                - <b> [format] </b>    Token that is changed with the new format name.\n
                - <b> [extension] </b> Token that is changed with the format file extensions.
     \note This method should be used in combination of the output tag system, which allows to reuse the name of the original alignment, 
        in the new output, along the new format and extension.
        Otherwise, the system would overwrite the same file over and over.
     \param outFormats Output formats that original files should reformat to. 
     */
    void loadAndSaveMultipleAlignments(std::vector< std::string >* inFile, std::string* outPattern, std::vector< std::string >* outFormats);
    
    /**
     \brief Function to obtain the format name of a given file.
     \param inFile File path of the file which we want to obtain it's format.
     */
    std::string getFileFormatName(std::string inFile);
    
    /**
     \brief Function to obtain all format names available
      by this object that can \b load an alignment.
     */
    std::string getInputFormatsAvailable();
    
    /**
     \brief Function to obtain all format names available
      by this object that can \b save an alignment.
     */
    std::string getOutputFormatsAvailable();
    
    /**
     \brief Function to divide an alignment into different alignments,
             each one with a sequence from the original.\n
            This function does a deep copy of each sequence,
             so the original alignment can be deleted after being splitted*/
    std::vector<Alignment*> splitAlignmentKeeping(Alignment& alignment);
};


#endif // READWRITEMS_H