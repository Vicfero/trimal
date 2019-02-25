#include "FormatHandling/fasta_m10_state.h"

#include "FormatHandling/FormatManager.h"
#include "utils.h"

namespace FormatHandling {
int fasta_m10_state::CheckAlignment([[maybe_unused]] std::istream* origin)
{
    return 0;
}

Alignment* fasta_m10_state::LoadAlignment([[maybe_unused]] const std::string &filename)
{
    return nullptr;
}

bool fasta_m10_state::SaveAlignment(const Alignment &alignment, std::ostream *output)
{
    /* Generate output alignment in FASTA format. Sequences can be unaligned. */

    ulong i, j, k, maxLongName;
    std::string *tmpMatrix;
    bool lastcharIsnewline = true;

    /* Depending on alignment orientation: forward or reverse. Copy directly
     * sequence information or get firstly the reversed sequences and then
     * copy it into local memory */
    if (Machine->reverse)
    {
        /* Allocate local memory for generating output alignment */
        tmpMatrix = new std::string[alignment.originalNumberOfSequences];
        for(i = 0; i < (ulong)alignment.originalNumberOfSequences; i++)
            tmpMatrix[i] = utils::getReverse(alignment.sequences[i]);
    }
    else tmpMatrix = alignment.sequences;

    /* Depending on if short name flag is activated (limits sequence name up to
     * 10 characters) or not, get maximum sequence name length. Consider those
     * cases when the user has asked to keep original sequence header */
    maxLongName = 0;
    for(i = 0; i < (unsigned)alignment.originalNumberOfSequences; i++)
    {
        if (alignment.saveSequences && alignment.saveSequences[i] == -1) continue;
        if (!Machine->keepHeader)
            maxLongName = utils::max(maxLongName, alignment.seqsName[i].size());
        else if (alignment.seqsInfo != nullptr)
            maxLongName = utils::max(maxLongName, alignment.seqsInfo[i].size());
    }

    if (maxLongName > PHYLIPDISTANCE) {
        maxLongName = PHYLIPDISTANCE;
        debug.report(WarningCode::HeaderWillBeCut, new std::string[1]{std::string(name)});
    }
    /* Print alignment. First, sequences name id and then the sequences itself */
    for(i = 0; i < (unsigned)alignment.originalNumberOfSequences; i++) {
        
        if (alignment.saveSequences != nullptr && alignment.saveSequences[i] == -1) continue;
        
        if (!Machine->keepHeader)
            (*output) << ">" << alignment.seqsName[i].substr(0, (unsigned)maxLongName) << "\n";
        
        else if (alignment.seqsInfo != nullptr)
            (*output) << ">" << alignment.seqsInfo[i].substr(0, (unsigned)maxLongName) << "\n";
        
        
        for (j = 0, k = 0; (ulong)j < (ulong)alignment.sequences[i].length(); j++)
        {
            if (alignment.saveResidues != nullptr && alignment.saveResidues[j] == -1) 
            {
                if (!lastcharIsnewline && (unsigned)j == alignment.sequences[i].length() -1 )
                {
                    (*output) << "\n";
                    lastcharIsnewline = true;
                }
            }
            else
            {
                (*output) << tmpMatrix[i][j];
                k++;
                lastcharIsnewline = false;
                if ((k % 60 == 0 || j == alignment.sequences[i].length() - 1))
                {
                    (*output) << "\n";
                    lastcharIsnewline = true;
                }
            }
        }
    }

    /* Deallocate local memory */
    if (Machine->reverse)
        delete [] tmpMatrix;
    
    return true;
}

bool fasta_m10_state::RecognizeOutputFormat(const std::string &FormatName)
{
    if (BaseFormatHandler::RecognizeOutputFormat(FormatName)) return true;
    return FormatName == "fasta_m10";
}
}