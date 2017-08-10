#include "../../include/ReadWriteMS/mega_sequential_state.h"
#include "../../include/ReadWriteMS/ReadWriteMachineState.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "../../include/defines.h"
#include "../../include/values.h"

using namespace std;

int MegaSequentialState::CheckAlignment(istream* origin)
{
    origin->seekg(0);
    origin->clear();
 
    char c, *firstWord = NULL, *line = NULL;
    int format = 0, blocks = 0;
    string nline;

    /* Read first valid line in a safer way */
    do {
        line = utils::readLine(*origin);
    } while ((line == NULL) && (!origin->eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (origin->eof())
        return 0;

    /* Otherwise, split line */
    firstWord = strtok(line, OTHDELIMITERS);

        /* Mega Format */
    if((!strcmp(firstWord, "#MEGA")) || (!strcmp(firstWord, "#mega"))) {

        /* Determine specific mega format: sequential or interleaved.
         * Counting the number of blocks (set of lines starting by "#") in
         * the input file. */
        blocks = 0;
        do {
            origin->read(&c, 1);
        } while((c != '#') && (!origin->eof()));

        do {
            while((c != '\n') && (!origin->eof()))
                origin->read(&c, 1);
            origin->read(&c, 1);
            if(c == '#')
                blocks++;
        } while((c != '\n') && (!origin->eof()));

        /* MEGA Sequential (22) or Interleaved (21) */
        return (!blocks) ? 1 : 0;
    }
    return 0;
}

newAlignment* MegaSequentialState::LoadAlignment(std::__cxx11::string filename)
{
    newAlignment * _alignment = new newAlignment();
   /* MEGA sequential file format parser */

    char *frag = NULL, *str = NULL, *line = NULL;
    ifstream file;
    int i;

    /* Check the file and its content */
    file.open(filename, ifstream::in);
    if(!utils::checkFile(file))
        return nullptr;

    /* Filename is stored as a title for MEGA input alignment.
     * If it is detected later a label "TITLE" in input file, this information
     * will be replaced for that one */
    _alignment->filename.append("!Title ");
    _alignment->filename.append(filename);
    _alignment->filename.append(";");

    /* Skip first valid line */
    do {
        line = utils::readLine(file);
    } while ((line == NULL) && (!file.eof()));

    /* If the file end is reached without a valid line, warn about it */
    if (file.eof())
        return nullptr;

    /* Try to get input alignment information */
    while(!file.eof()) {

        /* Destroy previously allocated memory */
        if (line != NULL)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* If a sequence name flag is found, go out from getting information loop */
        if(!strncmp(line, "#", 1))
            break;

        /* Destroy previously allocated memory */
        if (frag != NULL)
            delete [] frag;

        /* Create a local copy from input line */
        frag = new char[strlen(line) + 1];
        strcpy(frag, line);

        /* Split input line copy into pieces and analize it
         * looking for specific labels */
        str = strtok(frag, "!: ");
        for(i = 0; i < (int) strlen(str); i++)
            str[i] = toupper(str[i]);

        /* If TITLE label is found, replace previously stored information with
         * this info */
        if(!strcmp(str, "TITLE")) {
             _alignment->filename.clear();
            if(strncmp(line, "!", 1))
                 _alignment->filename += "!";
             _alignment->filename += line;
        }

            /* If FORMAT label is found, try to get some details from input file */
        else if(!strcmp(str, "FORMAT"))
            _alignment -> aligInfo.append(line, strlen(line));
    }

    /* Deallocate local memory */
    if (frag != NULL)
        delete [] frag;

    /* Count how many sequences are in input alignment file */
    do {

        /* Check whether input line is valid or not */
        if (line == NULL) {
            line = utils::readLine(file);
            continue;
        }

        /* If current line starts by a # means that it is a sequence name */
        if (!strncmp(line, "#", 1))
            _alignment->sequenNumber++;

        /* Destroy previously allocated memory */
        if (line != NULL)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);

    } while(!file.eof());

    /* Move file pointer to the beginner */
    file.clear();
    file.seekg(0);

    /* Allocate memory */
    _alignment->seqsName  = new string[_alignment->sequenNumber];
    _alignment->sequences = new string[_alignment->sequenNumber];

    /* Skip first line */
    line = utils::readLine(file);

    /* Skip lines until first sequence name is found */
    while(!file.eof()) {

        /* Destroy previously allocated memory */
        if (line != NULL)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
        if (line == NULL)
            continue;

        /* If sequence name label is found, go out from loop */
        if (!strncmp(line, "#", 1))
            break;
    }

    /* First sequence is already detected so its name should be stored */
    i = -1;

    /* This loop is a bit tricky because first sequence name has been already
     * detected, so it is necessary to process it before moving to next line.
     * That implies that lines are read at loop ends */
    while(!file.eof()) {

        /* Skip blank lines */
        if (line == NULL) {
            line = utils::readLine(file);
            continue;
        }

        /* Skip lines with comments */
        if (!strncmp(line, "!", 1)) {
            /* Deallocate memory and read a new line */
            delete [] line;
            line = utils::readLine(file);
            continue;
        }

        /* Remove comments inside of sequences and split input line */
        frag = utils::trimLine(line);

        /* Skip lines with only comments */
        if (frag == NULL) {

            /* Deallocate memory and read a new line */
            if (line != NULL)
                delete [] line;
            line = utils::readLine(file);
            continue;
        }

        /* Otherwise, split it into fragments */
        str = strtok(frag, " #\n");

        /* Sequence Name */
        if (!strncmp(line, "#", 1)) {
            i += 1;
            _alignment->seqsName[i].append(str, strlen(str));
            str = strtok(NULL, " #\n");
        }

        /* Sequence itself */
        while(str != NULL) {
            _alignment->sequences[i].append(str, strlen(str));
            str = strtok(NULL, " \n");
        }

        /* Deallocate dynamic memory */
        if (frag != NULL)
            delete [] frag;

        if (line != NULL)
            delete [] line;

        /* Read a new line in a safe way */
        line = utils::readLine(file);
    }

    /* Close input file */
    file.close();

    /* Deallocate dynamic memory */
    if (line != NULL)
        delete [] line;

    /* Check the matrix's content */
    _alignment->fillMatrices(true);
    return _alignment;
}

bool MegaSequentialState::SaveAlignment(newAlignment* alignment, std::ostream* output, std::string* FileName)
{
    /* Generate output alignment in MEGA format */

    int i, j, k;
    string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment->isAligned) {
        cerr << endl << "ERROR: Sequences are not aligned. Format (Mega) "
             << "not compatible with unaligned sequences." << endl << endl;
        return false;
    }

    /* Allocate local memory for generating output alignment */
    tmpMatrix = new string[alignment->sequenNumber];

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    for(i = 0; i < alignment->sequenNumber; i++)
        tmpMatrix[i] = (!Machine->reverse) ?
                       alignment->sequences[i] :
                       utils::getReverse(alignment->sequences[i]);

    /* Compute output file datatype */
    alignment->getAlignmentType();

    /* Print output alignment header */
    *output << "#MEGA" << endl << alignment->filename << endl;

    /* Print alignment datatype */
    if (alignment->dataType & SequenceTypes::DNA)
        *output << "!Format DataType=DNA ";
    else if (alignment->dataType & SequenceTypes::RNA)
        *output << "!Format DataType=RNA ";
    else if (alignment->dataType & SequenceTypes::AA)
        *output << "!Format DataType=protein ";

    /* Print number of sequences and alignment length */
    *output << "NSeqs=" << alignment->sequenNumber << " Nsites=" << alignment->residNumber
         << " indel=- CodeTable=Standard;" << endl << endl;

    /* Print sequences name and sequences divided into blocks of 50 residues */
    for(i = 0; i < alignment->sequenNumber; i++) {
        *output << "#" << alignment->seqsName[i] << endl;
        for(j = 0; j < alignment->residuesNumber[i]; j += 50) {
            for(k = j; ((k < alignment->residNumber) && (k < j + 50)); k += 10)
                *output << tmpMatrix[i].substr(k, 10) << " ";
            *output << endl;
        }
        *output << endl;
    }

    /* Deallocate local memory */
    delete [] tmpMatrix;
    
    return true;
}

bool MegaSequentialState::RecognizeOutputFormat(std::string FormatName)
{
    if (ReadWriteBaseState::RecognizeOutputFormat(FormatName)) return true;
    if (FormatName == "mega") return true;
    return false;
}