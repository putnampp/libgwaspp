/*
* Copyright (c) 2012, Patrick Putnam
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met: 
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer. 
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution. 
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* The views and conclusions contained in the software and documentation are those
* of the authors and should not be interpreted as representing official policies, 
* either expressed or implied, of the FreeBSD Project.
*/
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <cstring>
#include <memory>
#include <cmath>
#include <fstream>

#include <boost/program_options.hpp>

#include "libgwaspp.h"
#include "genetics/genetic_data.h"
#include "genetics/genetic_data_file.h"
#include "genetics/genotype/genotype_record.h"

#include "genetics/individual/individual_genotype_file.h"
#include "genetics/individual/tped_genotype_file.h"

#include "genetics/individual/individual_phenotype_file.h"
#include "genetics/individual/tfam_phenotype_file.h"

#include "genetics/individual/illumina_annotation_file.h"
#include "genetics/individual/tfam_annotation_file.h"

#include "genetics/genotype/geno_table.h"

#include "genetics/analyzable/case_control_set.h"

#include "algorithms/computation_engine.h"
#include "algorithms/epistasis_func.h"
#include "algorithms/maf_func.h"

#include "validation_tests/validate_func.h"

#include "util/time/timing.h"

#define MATHLIB_STANDALONE
#include "Rmath.h"

#define AA 0
#define AB 1
#define BB 2

using namespace std;
using namespace libgwaspp::genetics;
using namespace libgwaspp::algorithms;
using namespace util;

namespace po = boost::program_options;

const int QUIT = -1;
const int BUILD_CONTINGENCY_TABLES = 0;
const int CASE_CONTROL_EPISTASIS = 1;

int selectAnalysis();

const string HELP_KEY = "help";
const string VERSION_KEY = "version";

const string TPLINK_KEY = "tplink";
const string ILLUMINA_KEY = "illu";

const string GENOTYPE_FILE_KEY = "geno";
const string PHENOTYPE_FILE_KEY = "pheno";
const string CASE_CONTROL_ANNOTATION_FILE = "annot";

const string OUTPUT_FILE_KEY = "output";

const string TEST_CONTINGENGY_PERFORMANCE_KEY = "contin-perform";
const string TEST_CONTINGENCY_DEBUG_KEY = "contin-debug";
const string TEST_EPISTASIS_PERFORMANCE_KEY = "epi-perform";
const string TEST_EPISTASIS_DEBUG_KEY = "epi-debug";

const string TEST_DISTRIBUTION_PERFORMANCE_KEY = "dist-perform";
const string TEST_DISTRIBUTION_DEBUG_KEY = "dist-debug";

const string TEST_CC_SELECT_DIST_PERFORMANCE_KEY = "select-cc-maf";
const string TEST_CC_INLINE_DIST_PERFORMANCE_KEY = "inline-cc-maf";

const string TEST_INLINE_DIST_DEBUG_KEY = "test-inline-maf";
const string TEST_SELECT_DIST_DEBUG_KEY = "test-select-maf";

const string TEST_BOOST_KEY = "test-boost-epi";

const string VALIDATE_CALL_KEY = "valid-calls";
const string VALIDATE_GENO_KEY = "valid-geno";

const string COMPRESSION_LEVEL_KEY = "comp-level";

bool parseArguments( int argc, char **argv, po::variables_map &vm );

enum FileType { UNK, TPLINK, ILLUMINA };

FileType g_ft = UNK;

int main( int argc, char **argv ) {
    po::variables_map vm;

    if( !parseArguments(argc, argv, vm) ) {
        return 1;
    }

    string pheno_file = vm[ PHENOTYPE_FILE_KEY.c_str() ].as<string>();
    string geno_file = vm[ GENOTYPE_FILE_KEY.c_str() ].as<string>();
    string annot_file = vm[ CASE_CONTROL_ANNOTATION_FILE.c_str() ].as<string>();
    string out_file = vm[ OUTPUT_FILE_KEY.c_str() ].as<string>();

    auto_ptr<GeneticData> gd( new GeneticData( (eCompressionLevel) vm[ COMPRESSION_LEVEL_KEY ].as< int >() ) );
    auto_ptr<GeneticDataFile> ipf;  // phenotype file parser
    auto_ptr<GeneticDataFile> igf;  // genotype file parser
    auto_ptr<GeneticDataFile> iaf;  // annotation file parser

    auto_ptr< set< string > > marker_ids( new set<string>() ), individual_ids(new set<string>());

    char delim = ' ';
    if( g_ft == TPLINK ) {
        ipf = auto_ptr<GeneticDataFile>(new TfamPhenotypeFile());
        igf = auto_ptr<GeneticDataFile>(new TpedGenotypeFile());
        iaf = auto_ptr<GeneticDataFile>(new TFamAnnotationFile());
        annot_file = pheno_file;
    } else if( g_ft == ILLUMINA ) {
        ipf = auto_ptr<GeneticDataFile>(new IndividualPhenotypeFile());
        igf = auto_ptr<GeneticDataFile>(new IndividualGenotypeFile());
        iaf = auto_ptr<GeneticDataFile>(new IlluminaAnnotationFile());
    } else {
        cout << "Unknown Phenotype file" << endl;
        return 1;
    }

    ostream *out;
    if( out_file.empty() ) {
        out = &cout;
    } else {
        ofstream *tmp = new ofstream( out_file.c_str() );
        if( !tmp->is_open() ) {
            cout << "ERROR: Could not open the output file " << out_file << endl;
            return 1;
        }
        out = (ostream *) tmp;
    }

    cout << "Start to populate data" << endl;

    INIT_LAPSE_TIME;
    RECORD_START;

    ipf->populateGeneticData( pheno_file, &*gd, delim);

    RECORD_STOP;
    PRINT_LAPSE(cout, "Time to populate Phenotype Data: " );
    cout << endl;

    cout << "Found " << gd->getIndividualCount() << " individuals." << endl;
    cout << "Found " << gd->getPhenotypeCount() << " phenotyped traits." << endl;

    cout << "Starting to populate genotype data" << endl;
    igf->populateGeneticData( geno_file, &*gd, delim);

    cout << "Genotyped Individual Count: " << gd->getGenotypedIndividualsCount() << endl;
    cout << "Found " << gd->getMarkerCount() << " markers" << endl;
    iaf->populateGeneticData( annot_file, &*gd, delim );

    BasicInput inp_all( gd.get(), marker_ids.get(), individual_ids.get() );

    if( vm.count( VALIDATE_CALL_KEY ) || vm.count( VALIDATE_GENO_KEY ) ) {
        compute( printAllCalls, (void *) &inp_all, (void *) out );
    }

    if( vm.count(TEST_DISTRIBUTION_PERFORMANCE_KEY)) {
        compute(genotype_dist_performance, &*gd, out );
    }

    if( vm.count( TEST_CONTINGENGY_PERFORMANCE_KEY ) ) {
        compute( ContingencyPerformance, ( void * ) &inp_all, ( void * ) out );
    }

    if( vm.count( TEST_CONTINGENCY_DEBUG_KEY ) ) {
        compute( ContingencyDebug,  ( void * ) &inp_all, ( void * ) out );
    }

    if( vm.count( TEST_EPISTASIS_DEBUG_KEY ) ) {
        compute ( EpistasisDebug,  ( void * ) &inp_all, ( void * ) out );
    }

    if( vm.count( TEST_EPISTASIS_PERFORMANCE_KEY ) ) {
        compute ( EpistasisPerformance,  ( void * ) &inp_all, ( void * ) out );
    }

    if( vm.count( TEST_CC_INLINE_DIST_PERFORMANCE_KEY ) ) {
        compute ( inline_cc_maf, &*gd, out );
    }

    if( vm.count( TEST_CC_SELECT_DIST_PERFORMANCE_KEY ) ) {
        compute ( select_cc_maf, &*gd, out );
    }

    if( vm.count( TEST_INLINE_DIST_DEBUG_KEY ) ) {
        compute ( inline_maf_print, &*gd, out );
    }

    if( vm.count( TEST_BOOST_KEY ) ) {
        compute( computeBoost, &*gd, out );
    }

    marker_ids->clear();
    individual_ids->clear();

    cout << "DONE" << endl;
    return 0;
}

bool parseArguments( int argc, char **argv, po::variables_map &vm ) {
    po::options_description general( "General Options" );
    general.add_options()
    (( HELP_KEY + ",h" ).c_str(), "Help options" )
    (( VERSION_KEY + ",v" ).c_str(), "Version" )
    ;

    po::options_description data_opt( "Data Files" );
    data_opt.add_options()
    (( GENOTYPE_FILE_KEY + ",g" ).c_str(), po::value< string >()->default_value( "" ), "Genotype File" )
    (( PHENOTYPE_FILE_KEY + ",p").c_str(), po::value< string >()->default_value( "" ), "Phenotype File" )
    (( CASE_CONTROL_ANNOTATION_FILE + ",a").c_str(), po::value<string>()->default_value(""), "Case/Control set annotation file")
    (( OUTPUT_FILE_KEY + ",o").c_str(), po::value< string >()->default_value( "" ), "Results file")
    ;

    po::options_description annotations("Annotation Types (optional)");
    annotations.add_options()
    (( TPLINK_KEY + ",t").c_str(), "Providing genotype file in TPED format, and phenotype file in TFAM format")
    (( ILLUMINA_KEY + ",i").c_str(), "Providing genotype file in ILLU format, and a table of phenotype values")
    (( COMPRESSION_LEVEL_KEY ).c_str(), po::value< int >()->default_value( e2BitBlockCompression ), "Specify which compression method to use default is 2-bit block compression")
    ;

    po::options_description tests("Tests");
    tests.add_options()
    ((TEST_CONTINGENGY_PERFORMANCE_KEY ).c_str(), "Performance test of Contingency tables")
    ((TEST_CONTINGENCY_DEBUG_KEY).c_str(), "Build and Print Contingency tables")
    ((TEST_EPISTASIS_PERFORMANCE_KEY).c_str(), "Performance test of Epistasis algorithm")
    ((TEST_EPISTASIS_DEBUG_KEY).c_str(), "Build and Print Case/Control Contingency Tables for Epistasis analysis" )
    ((TEST_DISTRIBUTION_PERFORMANCE_KEY).c_str(), "Performance test of Genotype Distribution table")
    ((TEST_CC_SELECT_DIST_PERFORMANCE_KEY).c_str(), "Performance test of Case/Control Genotype Distributions; Select C/C first, then compute MAF")
    ((TEST_CC_INLINE_DIST_PERFORMANCE_KEY).c_str(), "Performance test of Case/Control Genotype Distributions; Compute MAF using C/C inline")
    ((TEST_INLINE_DIST_DEBUG_KEY).c_str(), "Print all maf distributions for the input set, using an inline distribution counting method")
    ((TEST_BOOST_KEY).c_str(), "Perform Epistasis analysis using an optimized BOOST algorithm")
    ;

    po::options_description validate( "Validations" );
    validate.add_options()
    ((VALIDATE_CALL_KEY).c_str(), "Print ALL input calls")
    ((VALIDATE_GENO_KEY).c_str(), "Print ALL input genotypes")
    ;

    po::options_description cmdline;
    cmdline.add( general ).add(data_opt).add( annotations ).add( tests ).add(validate);

    po::positional_options_description p;
    p.add(( GENOTYPE_FILE_KEY).c_str(), 1).add(( PHENOTYPE_FILE_KEY).c_str(), 1).add(( CASE_CONTROL_ANNOTATION_FILE).c_str(), 1);

    po::store( po::command_line_parser( argc, argv ).options( cmdline ).positional( p ).run(), vm );
    po::notify( vm );

    if( vm.count( HELP_KEY.c_str() ) ) {
        cout << cmdline << "\n";
        return false;
    }

    if( vm[ GENOTYPE_FILE_KEY ].as<string>().empty() ) {
        cout << "ERROR: A genotype file must be specified." << endl;
        return false;
    } else if( vm[ PHENOTYPE_FILE_KEY ].as<string>().empty()) {
        cout << "ERROR: A phenotype file must be specified." << endl;
        return false;
    }

    if( vm.count( TPLINK_KEY ) ) {
        cout << "Transposed PLINK format" << endl;
        g_ft = TPLINK;
    } else if( vm.count( ILLUMINA_KEY ) ) {
        cout << "Illumina genotype table file format" << endl;
        g_ft = ILLUMINA;
    } else {
        cout << "Guessing file format from extensions" << endl;
    }

    switch( (eCompressionLevel) vm[ COMPRESSION_LEVEL_KEY ].as< int >() ) {
    case eBasicCompression:
    case eByteCompression:
    case eHalfByteCompression:
    case e2BitBlockCompression:
    case e3BitStream:
    case e2BitStream:
        break;
    default:
        cout << "Invalid Compression Level specified.";
        return false;
    }

    return true;
}
