#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

using namespace std;
// declare functions

void check_args( int argc, char** argv);

class Parse_args
{
public:
  int k, m;
  string fa_file;
  void set_args( int argc, char** argv );
};

void Parse_args::set_args( int argc, char** argv )
{
  check_args( argc, argv );
  fa_file = argv[1];
  k = stoi(argv[2]);
  m = stoi(argv[3]);

  cout << "fasta=" << fa_file << "\tk=" << k << "\tm=" << m << endl;
}

void check_args( int argc, char** argv)
{
  string docstring =  "A simple mans C++ kmeriser tool. ArthurVM '20. ";
  string usage =      "Usage: kmerise <seq.fa> <k> <min count> ";
  string kstring =    "       <k>         The kmer size to assess. ";
  string mcstring =   "       <m>         The maximum count to return kmers. ";
  string helpstring = "       -h          Show this message and exit. ";

  string __doc__ =
  docstring + "\n\n" +
  usage + "\n" +
  kstring + "\n" +
  mcstring + "\n" +
  helpstring;

  if (argc < 3)
  {
    cerr << __doc__ << endl;
    exit (EXIT_FAILURE);
  }
}

void printNtimes( char c, int n )
{
  // character c will be printed n times
  cout << string(n, c);
}
