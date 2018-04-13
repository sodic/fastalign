#include <iostream>
#include <winnowing/winnowing.hpp>
#include "bioparser/bioparser.hpp"
#include <vector>
#include <unordered_map>
#include "FASTARead/FASTARead.h"

#define K 15
#define W 5

using namespace std;

int main(int argc, char const *argv[]) {

    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTARead>("../examples/escherichia_coli_reference.fasta");
    fasta_parser->parse_objects(fasta_reads, static_cast<uint64_t>(-1));

    const char* reference = fasta_reads[0]->get_data();
    const uint32_t length = fasta_reads[0]->get_data_length();

    vector<winnowing::minimizer> minimizers;
    unordered_map<winnowing::minhash_t, vector<int32_t>> lookup_table;
    winnowing::index_sequence(reference, length, W,K ,minimizers, lookup_table);

    cout << "Reference sequence successfully indexed." << endl;

    return 0;
}