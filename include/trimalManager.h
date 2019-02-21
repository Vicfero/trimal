//
// Created by bioinfo on 8/06/17.
//

#ifndef TRIMAL_TRIMALARGUMENTPARSER_H
#define TRIMAL_TRIMALARGUMENTPARSER_H


#include "FormatHandling/FormatManager.h"
#include "Alignment/Alignment.h"

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <iosfwd>
#include <string>


// Forward declatarion
class similarityMatrix;
namespace statistics
{
    class Consistency;
    class similarityMatrix;
}

class trimAlManager
{
public:
    std::vector<std::string> * vcfs = nullptr;
    
    
    bool 
        appearErrors        = false,

        getComplementary    = false, 
        columnNumbering     = false,
        nogaps              = false, 
        noallgaps           = false,
        gappyout            = false, 
        strict              = false,
        strictplus          = false, 
        automated1          = false,
        patternTrim         = false,
        complexPatternTrim  = false,
        sgc                 = false,
        sgt                 = false, 
        ssc                 = false, 
        sst                 = false,
        sfc                 = false, 
        sft                 = false, 
        sident              = false, 
        soverlap            = false, 
        selectSeqs          = false,
        selectCols          = false, 
        splitByStopCodon    = false,
        terminalOnly        = false, 
        keepSeqs            = false,
        ignoreStopCodon     = false,
        ignoreFilter        = false;

    float 
        conservationThreshold   = -1,
        gapThreshold            = -1,
        similarityThreshold     = -1,
        consistencyThreshold    = -1,
        residuesOverlap         = -1,
        sequenceOverlap         = -1,
        maxIdentity             = -1,
        minCoverage             = -1,
        minQuality              = -1;
    
    int
        i                       = 1,
        stats                   = 0,
        windowSize              = -1, 
        gapWindow               = -1, 
        similarityWindow        = -1,
        consistencyWindow       = -1, 
        blockSize               = -1, 
        clusters                = -1
        ,
        automatedMethodCount    = -1,
        alternative_matrix      = -1,
        
        *delColumns         = nullptr,
        *delSequences       = nullptr,
        *sequencesLengths   = nullptr;
    
    size_t
            argumentLength = size_t(-1);

    std::string 
        *sequencesNames     = nullptr;
    

    /* Others variables */
    std::ifstream compare;
    statistics::similarityMatrix *similMatrix   = nullptr;
    
    Alignment
        *origAlig                   = nullptr,
        *singleAlig                 = nullptr,
        *tempAlig                   = nullptr,
        *backtranslationAlig        = nullptr,
        **compareAlignmentsArray    = nullptr;

    char 
        c, 
        *forceFile          = nullptr,
        
        *infile             = nullptr,
        
        *backtransFile      = nullptr,
        
        *outfile            = nullptr,
        *htmlOutFile        = nullptr,
        *svgOutFile         = nullptr,
        *svgStatsOutFile    = nullptr,

        *compareset         = nullptr,

        
        *matrixFile         = nullptr,
                
        **filesToCompare    = nullptr,
        line[256];
             
    std::vector<std::string> oformats;

    /// Consistency Manager. 
    ///     We have to save the reference, 
    ///     while it is still not part of an alignment
    statistics::Consistency *CS = nullptr;

public:
    trimAlManager();
    ~trimAlManager();
    void parseArguments(int argc, char *argv[]);
private: // Parse Arguments Methods
        void verbosity_argument             (const int* argc, char* argv[]);
        
        void help_arguments                 (const int *argc, char* argv[], int *currentArg);

        bool in_argument                    (const int* argc, char* argv[], int* currentArg);
        bool vcf_argument                   (const int* argc, char* argv[], int* currentArg);
        bool out_argument                   (const int* argc, char* argv[], int* currentArg);
        bool html_out_argument              (const int* argc, char* argv[], int* currentArg);
        bool timetracker_out_argument       (const int* argc, char* argv[], int* currentArg);
        bool svg_out_argument               (const int* argc, char* argv[], int* currentArg);
        bool svg_stats_argument             (const int* argc, char* argv[], int* currentArg);
        bool out_format_arguments           (const int* argc, char* argv[], int* currentArg);
        bool matrix_argument                (const int* argc, char* argv[], int* currentArg);
        bool compareset_argument            (const int* argc, char* argv[], int* currentArg);
        bool force_select_argument          (const int* argc, char* argv[], int* currentArg);
        bool back_trans_argument            (const int* argc, char* argv[], int* currentArg);
        bool gap_threshold_argument         (const int* argc, char* argv[], int* currentArg);
        bool similarity_threshold_argument  (const int* argc, char* argv[], int* currentArg);
        bool consistency_threshold_argument (const int* argc, char* argv[], int* currentArg);
        bool conservation_threshold_argument(const int* argc, char* argv[], int* currentArg);
        bool select_cols_argument           (const int* argc, char* argv[], int* currentArg);
        bool no_gaps_argument               (const int* argc, char* argv[], int* currentArg);
        bool no_all_gaps_argument           (const int* argc, char* argv[], int* currentArg);
        bool keep_seqs_argument             (const int* argc, char* argv[], int* currentArg);
        bool keep_header_argument           (const int* argc, char* argv[], int* currentArg);
        bool gappy_out_argument             (const int* argc, char* argv[], int* currentArg);
        bool strict_argument                (const int* argc, char* argv[], int* currentArg);
        bool strict_plus_argument           (const int* argc, char* argv[], int* currentArg);
        bool automated1_argument            (const int* argc, char* argv[], int* currentArg);
        bool residue_overlap_argument       (const int* argc, char* argv[], int* currentArg);
        bool sequence_overlap_argument      (const int* argc, char* argv[], int* currentArg);
        bool seqs_select_argument           (const int* argc, char* argv[], int* currentArg);
        bool max_identity_argument          (const int* argc, char* argv[], int* currentArg);
        bool clusters_argument              (const int* argc, char* argv[], int* currentArg);
        bool terminal_only_argument         (const int* argc, char* argv[], int* currentArg);
        bool window_argument                (const int* argc, char* argv[], int* currentArg);
        bool gap_window_argument            (const int* argc, char* argv[], int* currentArg);
        bool similarity_window_argument     (const int* argc, char* argv[], int* currentArg);
        bool consistency_window_argument    (const int* argc, char* argv[], int* currentArg);
        bool block_argument                 (const int* argc, char* argv[], int* currentArg);
        bool stats_arguments                (const int* argc, char* argv[], int* currentArg);
        bool complementary_argument         (const int* argc, char* argv[], int* currentArg);
        bool col_numbering_argument         (const int* argc, char* argv[], int* currentArg);
        bool split_by_stop_codon_argument   (const int* argc, char* argv[], int* currentArg);
        bool ignore_stop_codon_argument     (const int* argc, char* argv[], int* currentArg);
        bool ignore_filter_argument         (const int* argc, char* argv[], int* currentArg);
        bool min_quality_argument           (const int* argc, char* argv[], int* currentArg);
        bool min_coverage_argument          (const int* argc, char* argv[], int* currentArg);
        bool patternCleaning                (const int* argc, char* argv[], int* currentArg);
        bool complexPatternCleaning         (const int* argc, char* argv[], int* currentArg);

public:
    bool processArguments(char* argv[]);
private: // Process Arguments Methods
    bool check_arguments_incompatibilities();
        bool check_inFile_incompatibilities();
        bool check_select_cols_and_seqs_incompatibilities();
        bool check_thresholds_incompatibilities();
        bool check_automated_methods_incompatibilities();
        bool check_max_identity_incompatibilities();
        bool check_clusters_incompatibilities();
        bool check_windows_incompatibilities();
        bool check_stats_incompatibilities();
        bool check_codon_behaviour_incompatibility();
        bool check_combinations_among_thresholds_incompatibility();
        bool check_automated_manual_incompatibilities();
        bool check_vcf_incompatibility();

    bool check_arguments_needs(char* argv[]);

        bool check_force_selection();
        bool check_input_file_with_coding_sequences_argument();
        bool check_file_aligned();
        bool check_similarity_matrix();
        bool check_outputs_coincidence();
        bool check_col_numbering();
        bool check_residue_and_sequence_overlap();
        bool check_output_relevance();
        bool check_output_file_with_statistics();
        bool check_block_size();
        bool check_backtranslations();
        bool check_coding_sequences_type();
        bool check_and_prepare_coding_sequence();
        bool check_backtranslation_infile_names_correspondence();
        void check_compareset_window_argument();
        void check_output_format();

public:
    int perform();
    void delete_variables();

private: // General, private, methods
    int innerPerform();
    bool performCompareset();

    void output_reports();
    void save_alignment();

    void svg_stats_out();
    void print_statistics();
    bool create_or_use_similarity_matrix();
    void clean_alignment();
    void postprocess_alignment();

        void set_window_size();

    // NON COMPLEX OPTIONS
    void menu();
    void examples();


    void CleanSequences();
    void CleanResiduesAuto();
    void CleanResiduesNonAuto();


private:
    FormatHandling::FormatManager formatManager;
};



#endif //TRIMAL_TRIMALARGUMENTPARSER_H
