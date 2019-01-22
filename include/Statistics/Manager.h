//
// Created by bioinfo on 5/06/17.
//

#ifndef TRIMAL_STATISTICSMANAGER_H
#define TRIMAL_STATISTICSMANAGER_H
#include "similarityMatrix.h"

// Forward declarations
class Alignment;



/// Namespace containing all classes related to statistics handling
namespace statistics {
// Forward declarations
    class Gaps;

    class Similarity;

    class Consistency;

    class Entropy;

    /**
     * \brief Class to handle the interaction with Alignment and statistics objects.\n
     * It serves as a wrapper or intermediate between the alignment and each specific stat.\n
     * It also encapsulates the similarityMatrix.
     */
    class Manager {
    public:

        /**
         * \brief Gaps submodule
         * */
        Gaps *gaps = nullptr;
        /**
         * \brief Similarity submodule
         * */
        Similarity *similarity = nullptr;
        /**
         * \brief Consistency submodule
         * */
        Consistency *consistency = nullptr;

        /**
         * \brief Entropy submodule
         * */
        Entropy *entropy = nullptr;

        /**
         * \brief SimilarityMatrix used by Similarity
         * */
        similarityMatrix *_similarityMatrix = nullptr;

        /**
         * \brief Gap window
         * */
        int ghWindow;
        /**
         * \brief Similarity window
         * */
        int shWindow;

        /**
         * \brief Entropy window
         * */
        int ehWindow;

        /**
         * \brief Method to set a similarity matrix
         * */
        bool setSimilarityMatrix(similarityMatrix *sm);

        /**
         * \brief Method to handle gap stat calculation\n
         * It checks if the #gaps submodule has been created, otherwise, creates it
         * */
        bool calculateGapStats();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsColumns()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsGapsColumns();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsAcl()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsGapsTotal();

        /**
         * \brief Method to handle gap stat calculation\n
         * It checks if the #gaps submodule has been created, otherwise, creates it
         * */
        bool calculateEntropyStats();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsColumns()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsEntropyColumns();

        /**
         * \brief Wrapper to Statistics::Gaps::printGapsAcl()\n
         * It calls to calculateGapStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsEntropyTotal();

        /**
         * \brief Method to handle similarity stat calculation\n
         * It checks if the #similarity submodule has been created, otherwise, creates it
         * */
        bool calculateConservationStats();

        /**
         * \brief Wrapper to Statistics::Similarity::printConservationAcl()\n
         * It calls to calculateConservationStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsConservationColumns();

        /**
         * \brief Wrapper to Statistics::Similarity::printGapsTotal()\n
         * It calls to calculateConservationStats() to make sure the information is available before reporting the requested values
         * */
        void printStatisticsConservationTotal();

        /**
         * \brief Method to print the vector containing the keep/reject (Alignment::saveResidues) values of the associated Alignment\n
         * If \b -1 , the residue would be rejected, otherwise, it prints the position on the alignment of the kept residue
         * */
        void printCorrespondence();

    private:
        /// Making Alignment a friend class allows us to have private constructor and destructors.\n
        /// This gives us more control over the statistics classes, hinting that they shouldn't be created outside the scope of a newAlignmet.
        friend class ::Alignment;

        Alignment *alig;

        explicit Manager(Alignment *parent);

        Manager(Alignment *parent, Manager *mold);

        ~Manager();
    };

}
#endif //TRIMAL_STATISTICSMANAGER_H
