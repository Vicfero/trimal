//
// Created by vfernandez on 25/01/19.
//

#include <Alignment/Alignment.h>
#include <FormatHandling/FormatManager.h>
#include <utils.h>
#include <reportsystem.h>

#include <iomanip>

#define DEBUG !NDEBUG

enum returnCodes
{
    correct = 0,
    notPerfectMatch = 1,
    notSameSize = 2,
    referenceNotLoaded = 4,
    compareAlignmentNotLoaded = 8,
    scoreNotRecognized = 16
};

bool compareAlig(
        const Alignment * const ref,
        const Alignment * const cmp,
        const float score,
        float * scoreVector,
        int * counterVector,
        bool verbose = true)
{

    int
        // Reference alignment sequence column index
        it0 = 0,
        // Compare alignment sequence column index
        it1 = 0;
    int SequenceIndex;

#if DEBUG
    // Perfect match between alignments?
    // Used as SanityCheck when comparing two identical alignments
    bool perfectMatch = true;
#endif

    if (verbose) std::cout << cmp->filename << "\n\t";

    // Start comparing by the current alignment column
    for(; it1 < cmp->originalNumberOfResidues; it1++)
    {
        // Start comparing, column by column, of the original alignment
        for(; it0 < ref->originalNumberOfResidues; it0++)
        {
            // Go residue by residue (iterate over all sequences al fixed residue position)
            for (SequenceIndex = 0; SequenceIndex < ref->originalNumberOfSequences; SequenceIndex++)
            {
                // If one residue is not the same, skip the loop.
                if (ref->sequences[SequenceIndex][it0] != cmp->sequences[SequenceIndex][it1])
                    break;
            }

            // If the loop hasn't been skipped, we found a perfect match
            if (SequenceIndex == ref->originalNumberOfSequences) {
                // Increase the N for the column
                counterVector[it0]++;
                // Mark as kept
                scoreVector[it0] += score;
                // Go over next column in the cmp alignment
                it1++;

                // If verbose mode, print the score for the column
                if (verbose)
                    std::cout << ' ' << score ;

            // If we haven't found a perfect match
            } else {
#if DEBUG
                // Mark as removed
                perfectMatch = false;
#endif
                // If verbose mode, print the score for the column ~ 0.0F
                if (verbose)
                    std::cout << ' ' << 0.0F ;
            }

            // Go over next column in the ref alignment is done within the loop
        }
    }

    if (verbose) std::cout << '\n';
#if DEBUG
    return perfectMatch;
#else
    return true;
#endif
}

void printUsage()
{
    std::cout << "Usage:"
                 "\n\tmultipleNumCol [-verbose]"
                 "\n\treferenceAlignment referenceScore"
                 "\n\t[compare_alignment1 alignment_score1]"
                 "\n\t[compare_alignment2 ...]" << std::endl;
    exit(0);
}

int main(int argc, char *argv[]) {

#if DEBUG
    std::cout << "RUNNING IN DEBUG - "
                 "Please, run \"cmake -DCMAKE_BUILD_TYPE=RELEASE . ; make utils\" "
                 "if this wasn't intentional\n\n";
#endif

    // We're going to use Debug Messages
    debug.IsDebug = true;

    // argv - iterator
    int argv_i = 1;

    // Check if enough arguments have been provided
    if ((argc - argv_i) < 2)
    {
        printUsage();
    }

    // Verbose Option
    bool verbose = false;
    if (argc > 2 && !strcmp(argv[argv_i], "-verbose"))
    {
        verbose = true;
        argv_i++;
    }

    // Check there are paired arguments ~ Alignment + Score
    if ((argc - argv_i) % 2 != 0)
    {
        printUsage();
    }

    // Have a return value to provide meaningful error codes
    int returnValue = returnCodes::correct;

    // Configure cout to print float using the same space
    std::cout << std::fixed << std::setw(2) << std::setprecision(5);

    // Use the Format Manager to load alignments
    FormatHandling::FormatManager * formatManager = new FormatHandling::FormatManager();

    // Load the reference alignment argc-i = i + 1
    Alignment * ref = formatManager->loadAlignment(argv[argv_i++]), * alig;

    // Dont continue if reference alignment couldn't be loaded
    if (ref == nullptr)
        exit(ReferenceFileNotLoaded);

    // Get the score for the reference alignment
    float refScore;

    // Check is a number
    if (utils::isNumber(argv[argv_i]))
    {
        // Store the score
        refScore = (float) atof(argv[argv_i++]);
    } else {
        // Report the error
        debug.log(ERROR) << alig->filename
                         << "\n\t Score for reference not recognized: \""
                         << argv[argv_i] << "\"\n";
        exit(scoreNotRecognized);
    }

    // Keep a count of the number of sequences
    int seqCount = ref->originalNumberOfSequences;

    // Initialize the values array
    float * scoreVector = new float[ref->originalNumberOfResidues];

    // Initialize the counter array
    int * counterVector = new int[ref->originalNumberOfResidues];

    // Initialize the values array values to 0.
    for (int x = 0; x < ref->originalNumberOfResidues; ++x) {
        scoreVector[x] = 0;
        counterVector[x] = 0;
    }

#if DEBUG
    // Compare ref vs ref as a sanity check
    if (!compareAlig(
            // Alignments: Reference - Compare
            ref, ref,
            // Score of Compare alignment
            refScore,
            // Score and Counter vectors
            scoreVector, counterVector,
            // Verbose option
            verbose))
    {
        returnValue |= returnCodes::notPerfectMatch;
        debug.log(ERROR) << "This shouldn't happen, contact with developer.";
    } else if (verbose) {
        std::cout << '\n';
    }
#else
    compareAlig(
            // Alignments: Reference - Compare
            ref, ref,
            // Score of Compare alignment
            refScore,
            // Score and Counter vectors
            scoreVector, counterVector,
            // Verbose option
            verbose);
#endif

    // Iterate over all other alignments -> argc-i = 2+
    for (; argv_i < argc; argv_i++)
    {
        // Load the alignment to compare
        alig = formatManager->loadAlignment(argv[argv_i++]);

        if (alig == nullptr) {
            returnValue |= returnCodes::referenceNotLoaded;
            continue;
        }

        // SanityCheck: Compare all the alignments sequence count
        if (alig->originalNumberOfSequences != seqCount)
        {
            returnValue |= returnCodes::notSameSize;
            debug.log(ERROR) << "Skipping file " << alig->filename
                             << " as it doesn't have the same number of sequences as "
                                << ref->filename << '\n';
            delete alig;
            continue;
        }

        if (utils::isNumber(argv[argv_i]))
        {
            // Compare the reference vs the alignment
            compareAlig(
                    // Alignments: Reference - Compare
                    ref, alig,
                    // Score of Compare alignment
                    (const float) atof(argv[argv_i]),
                    // Score and Counter vectors
                    scoreVector, counterVector,
                    // Verbose Option
                    verbose);
        } else {
            returnValue |= scoreNotRecognized;
            debug.log(ERROR) << alig->filename
                             << "\n[ERROR] Score \"" << argv[argv_i]
                             << "\"for alignment \"" << alig->filename
                             << "\" not recognized.\n";
            delete alig;
            continue;
        }

        // Delete the cmp alignment
        delete alig;
    }

    if (verbose)
        std::cout << "\nFinal score Vector: \n\t";

    // Iterate over all scores, dividing their value by
    //      the number of elements that gave score
    for (int x = 0; x < ref->originalNumberOfResidues; ++x) {
        scoreVector[x] /= counterVector[x];
        std::cout << ' ' << scoreVector[x];
    }

    if (verbose)
    {
        std::cout << "\nCounter Vector: \n\t" << std::setprecision(0);
        // Iterate over all scores, dividing their value by
        //      the number of elements that gave score
        for (int x = 0; x < ref->originalNumberOfResidues; ++x) {
            std::cout << ' ' << counterVector[x] << std::setw(7);
        }
    }

    // Delete the format manager
    delete formatManager;

    // Delete the reference alignment
    delete ref;

    // Delete the temporal arrays
    delete [] scoreVector;
    delete [] counterVector;

    // Exit, returning the 'return value'
    exit(returnValue);
}
