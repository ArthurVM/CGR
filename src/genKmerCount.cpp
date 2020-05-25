#include "genKmerCount.hpp"

using namespace std;

// declare functions
map<string, int> read_fasta( string fa_file, int k );
void view_kmer_array( map<string, int>& kmer_map, int max_count );
Parse_args args;

int main( int argc, char** argv )
{
  // check and parse arguments
  args.set_args(argc, argv);
  // call functions
  map<string, int> kmer_map = read_fasta( args.fa_file, args.k );
  view_kmer_array( kmer_map, args.m );
}

map<string, int> read_fasta( string fa_file, int k )
{

  // Reads in a fasta file and geerates a counted kmer map structure

  void kmerise( string seq, int k, map<string, int>& kmer_map );

  map<string, int> kmer_map;

  std::ifstream input(fa_file);

    // check fa file can be opened
    if (!input.good())
    {
      cerr << "Error opening: " << fa_file << endl;
      return kmer_map;
    }

    std::string line, id, DNA_sequence;

    while (std::getline(input, line))
    {

      if(line.empty())
      {
        continue;
      }

      if (line[0] == '>')
      {
        // output previous line before overwriting id if it is non-empty
        if(!id.empty())
        {
          kmerise( DNA_sequence, k, kmer_map );
        }

        id = line.substr(1);
        DNA_sequence.clear();
      }

      else
      {
        DNA_sequence += line;
      }
    }

    // output final entry if id is non-empty
    if(!id.empty())
    {
      kmerise( DNA_sequence, k, kmer_map );
    }

  return kmer_map;

}

void view_kmer_array( map<string, int>& kmer_map, int max_count )
{
  // Take the kmer array map structure and print to stdout as:
  // '  args={argc}                   '
  // '  fasta=/path/to/fasta   k={k}  '
  // '  kmer count                    '
  // '  ...  ...                      '

  // map<string, int>::iterator it = kmer_map.begin();

  for( auto i = kmer_map.begin(); i != kmer_map.end(); i++)
  {
  // while(it != kmer_map.end()){
    string kmer = i -> first;
    float count = i -> second;

    // check max count to return kmers
    if(max_count != 0)
    {
      if(count <= max_count)
      {
        cout << kmer << "\t" << count << endl;
      }
    }
    // else return all kmers
    else
    {
      cout << kmer << "\t" << count << endl;
    }
  }
}

void kmerise( string seq, int k, map<string, int>& kmer_map )
{
  // takes a sequence and a value of k and generates a map of kmer counts

  for(int i = 0; i < seq.size()-k+1; i++)
  {
    string kmer = seq.substr(i, k);

    // transforms all kmers to uppercase
    // apart from lowercase DNA seq being gross to look at, it may be interpreted as soft masking by some software. This is bad.
    std::transform(kmer.begin(), kmer.end(), kmer.begin(), std::ptr_fun<int, int>(std::toupper));

    // counts kmers which do not contain "N" characters
    // IDEA :: extend to non "ACTG" characters
    if(kmer.find("N") == string::npos)
    {
      kmer_map[kmer]++;
    }
  }
}
