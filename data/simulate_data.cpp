#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <time.h>
#include <cstdio>
#include <cstdlib>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

using namespace std;

struct percent_count {
    double first;
    double second;
    percent_count( double p, double c ) : first(p), second(c) {}
};

typedef vector< percent_count * > Distro;

struct bin {
    string name;
    Distro * distribution;
    double total_size;

    bin( const string & n ) : name( n ), distribution( new Distro() ), total_size(0.0) {}
    ~bin( ) {
        while( !distribution->empty() ) {
            percent_count * p = distribution->back();
            distribution->pop_back();
            delete p;
        }
        delete distribution;
    }
};

typedef vector< bin * > bins;

bool parseBins( const string & file, bins * b );
ostream & operator<<( ostream & o, const bin & b );
void buildDataSet( bin * b, unsigned int nSnps, unsigned int nSamples, const string &  );
bool printTFam( const string & , unsigned int samples );

void print( ostream &, unsigned int, const vector< char > * alleles );

int main( int argc, char ** argv ) {

    if( argc != 6 ) {
        cout << "Expected usage: simgeno <dist_table> <dist_idx> <nSnps> <nSamples> <outfile_prefix>" << endl;
        return 1;
    }

    srand( time(NULL) );

    string dist_file( argv[1] );
    unsigned int dist_idx = boost::lexical_cast< unsigned int >( argv[2] );
    unsigned int nSnps = boost::lexical_cast< unsigned int >( argv[3] );
    unsigned int nSamples = boost::lexical_cast< unsigned int >( argv[4] );

    string out_prefix( argv[5] );
    
    bins b;

    cout << "Parsing distributions:" << endl;
    if( !parseBins( dist_file, &b ) ) {
        cout << "Failed to parse bins" << endl;
        return 1;
    }

    buildDataSet( b[dist_idx], nSnps, nSamples, out_prefix );

    while(! b.empty() ) {
        bin * del = b.back();
        b.pop_back();
        delete del;
    }

    return 0;
}

ostream & operator<<( ostream & o, const bin & b ) {
    o << "{ " << b.name << ", " << b.total_size << " }";
    return o;
}

bool parseBins( const string & file, bins * b ) {
    ifstream iFile( file.c_str() );

    if( !iFile.is_open() ) {
        cout << "Failed to open file: " << file << endl;
        return false;
    }

    string l = "";
    // parse column headers
    std::getline(iFile, l );
    if( iFile.bad() ) {
        cout << "Failed to read column headers" << endl;
        iFile.close();
        return false;
    }

    if( l.back() == '\r' ) {
        l.pop_back();
    }

    typedef boost::char_separator< char > seperator;
    typedef boost::tokenizer< seperator > tokenizer;
    seperator sep("\t");
    tokenizer tokens( l, sep );

    for( tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it ) {
        if(! it->empty() ) {
            b->push_back( new bin(*it) );
        }
    }

    while( !iFile.eof() ) {
        getline( iFile, l );
        if( l.empty() ) {
            continue;
        } else if( l.back() == '\r') {
            l.pop_back();
        }

        tokenizer t( l, sep );
        tokenizer::iterator it = t.begin();
        double percent = boost::lexical_cast< double >( (*it++) ) / 100.0;

        int i = 0;
        while( it != t.end() ) {
            if(! it->empty() ) {
                double count = boost::lexical_cast< double >( (*it) );

                (*b)[i]->total_size += count;
                (*b)[i++]->distribution->push_back( new percent_count( percent, count) );
                it++;
            }
        }
    }

    iFile.close();
    return true;
}

void buildDataSet( bin * b, unsigned int nSnps, unsigned int nSamples, const string & outprefix ) {
    cout << "Building data set for : " << *b << ", " << nSnps << ", " << nSamples << endl;
    vector< char > alleles;

    unsigned int nAlleles = 2 * nSamples;
    alleles.reserve( nAlleles );

    string tmp_filename = outprefix + "." + boost::lexical_cast<string>(nSnps) + "." + boost::lexical_cast<string>(nSamples) + ".tfam";
    
    if( !printTFam( tmp_filename, nSamples) ) {
        return;
    }

    tmp_filename = outprefix + "." + boost::lexical_cast<string>(nSnps) + "." + boost::lexical_cast<string>(nSamples) + ".tped";

    ofstream out( tmp_filename.c_str() );

    if( !out.is_open() ) {
        cout << "Unable to open: " << tmp_filename << endl;
        return;
    }

    vector< pair< double, unsigned int > > snp_sets;

    /*
    * pre-compute the expected number of snps from each bin
    */
    Distro::iterator it = b->distribution->begin();
    unsigned int t = 0;
    while( it != b->distribution->end() ) {
        percent_count * lb = (*it++);

        unsigned int m = ceil((lb->second / b->total_size) * (double) nSnps);
        snp_sets.push_back( make_pair( lb->first, m ) );
    }

    for( unsigned int t = 0; t < nSnps; ) {
        random_shuffle( snp_sets.begin(), snp_sets.end() );
        for( vector< pair< double, unsigned int > >::iterator sset = snp_sets.begin(); sset != snp_sets.end(); sset++ ) {
            if( sset->second > 0 ) {
                alleles.clear();
                int r = rand() % 1000; //  [0, 1000)
                double p = sset->first + ((double) r / 100000.0); // percent + [0.00000, 0.00999]
                unsigned int s = (unsigned int) floor( p * 2.0 * (double) nSamples);

                for( unsigned int j = 0; j < s; ++j ) {
                    alleles.push_back( 'B' );
                }

                for( ; s < nAlleles; ++s ) {
                    alleles.push_back( 'A' );
                }

                random_shuffle( alleles.begin(), alleles.end() );

                print(out, t++, &alleles );
                sset->second--;
            }
        }
    }
}

bool printTFam( const string & file, unsigned int samples ) {
    ofstream out( file.c_str() );
    
    if( !out.is_open() ) {
        cout << "Unable to open: " << file << endl;
        return false;
    }
    unsigned int rnd = rand();
    int j = 0;
    for( unsigned int i = 0; i < samples; ++i ) {
        out << "Fam_" << i 
            << "\tInd_" << i 
            << "\tPat_" << i 
            << "\tMat_" << i
            << "\tx";

        out << "\t" << (int)(rnd & 1) << "\n";

        rnd >>= 1;
        if( ++j >= 32 ) {
            rnd = rand();
            j = 0;
        }
    }

    return true;
}

void print( ostream & o, unsigned int snp_idx, const vector< char > * alleles ) {
    if( alleles->empty() ) return;

    vector< char >::const_iterator it = alleles->begin();
    o << "0\trs" << snp_idx << "\t0\t0";
    while( it != alleles->end() ) {
        o << "\t" << (*it++);
    }
    o << "\n";
}
