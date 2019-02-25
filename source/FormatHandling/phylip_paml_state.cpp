#include "FormatHandling/phylip_paml_state.h"

#include "FormatHandling/FormatManager.h"
#include "defines.h"
#include "utils.h"

namespace FormatHandling {

int phylip_paml_state::CheckAlignment([[maybe_unused]] std::istream* origin)
{
    return 0;
}

Alignment* phylip_paml_state::LoadAlignment([[maybe_unused]] const std::string &filename)
{
    return nullptr;
}

bool phylip_paml_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
    /* Generate output alignment in PHYLIP format compatible with PAML program */

    int i, maxLongName;
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

    maxLongName = PHYLIPDISTANCE;
    for(i = 0; (i < alignment.numberOfSequences); i++)
        maxLongName = utils::max(maxLongName, (int)alignment.seqsName[i].size());

    /* Generating output alignment */
    /* First Line: Sequences Number & Residued Number */
    *output << " " << alignment.numberOfSequences << " " << alignment.numberOfResidues << "\n";

    /* Print alignment */
    /* Print sequences name follow by the sequence itself in the same line */
    for(i = 0; i < alignment.numberOfSequences; i++)
        *output << std::setw(maxLongName + 3) << std::left << alignment.seqsName[i].substr(0, maxLongName)
             << alignment.sequences[i] << "\n";
    *output << "\n";

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool phylip_paml_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "phylippaml";
}

}