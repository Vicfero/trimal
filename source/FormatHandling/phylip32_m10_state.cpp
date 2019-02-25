#include "FormatHandling/phylip32_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "utils.h"

namespace FormatHandling {

int phylip32_m10_state::CheckAlignment([[maybe_unused]] std::istream* origin)
{
    return 0;
}

Alignment* phylip32_m10_state::LoadAlignment([[maybe_unused]] const std::string &filename)
{
    return nullptr;
}

bool phylip32_m10_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
    /* Generate output alignment in PHYLIP 3.2 format (interleaved) */

    int i, j, k, maxLongName;
    std::string *tmpMatrix;

    /* Check whether sequences in the alignment are aligned or not.
     * Warn about it if there are not aligned. */
    if (!alignment.isAligned) {
        debug.report(ErrorCode::UnalignedAlignmentToAlignedFormat, new std::string[1] { this->name });
        return false;
    }

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    if (Machine->reverse)
    {
        /* Allocate local memory for generating output alignment */
        tmpMatrix = new std::string[alignment.originalNumberOfSequences];
        for(i = 0; i < alignment.originalNumberOfSequences; i++)
            tmpMatrix[i] = utils::getReverse(alignment.sequences[i]);
    }
    else tmpMatrix = alignment.sequences;

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length */
    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment.originalNumberOfSequences); i++)
        if  (alignment.saveSequences[i] != -1)
            maxLongName = utils::max(maxLongName, (int)alignment.seqsName[i].size());

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        debug.report(WarningCode::HeaderWillBeCut, new std::string[1]{std::string(name)});
    }

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    (*output) << " " << alignment.numberOfSequences << " " << alignment.numberOfResidues;

    /* Alignment */
    /* For each sequence, print its identifier and then the sequence itself in
     * blocks of 50 residues */
    for(i = 0; i < alignment.originalNumberOfSequences; i++) {
        /* Sequence Name */
        if (alignment.saveSequences[i] == -1) continue;
        (*output) << "\n" << std::setw(maxLongName + 3) << std::left << alignment.seqsName[i].substr(0, maxLongName);
        /* Sequence. Each line contains a block of 5 times 10 residues. */
        
        for (j = 0, k = 0; j < alignment.originalNumberOfResidues; j++)
        {
            if (alignment.saveResidues[j] == -1) continue;
            if (k == 50)
            {
                *output << "\n" << std::setw(maxLongName + 3) << std::left << " " ;
                k = 0;
            }
            *output << alignment.sequences[i][j];
            k++;
            if (k % 10 == 0) 
                *output << " ";
        }
        if (k % 10 != 0)
            *output << " ";
        
        /* Print a blank line to mark sequences separation */
//         if (k % 50 != 0)
        (*output) << "\n";
    }
    *output << "\n";

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool phylip32_m10_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylip32_m10";
}

}