#include <iostream>
#include <winnowing/winnowing.hpp>
#include "bioparser/bioparser.hpp"
#include "FASTARead/FASTARead.h"

using namespace std;

int main(int argc, char const *argv[]) {

    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTARead>(
            "../examples/escherichia_coli_reference.fasta");
    fasta_parser->parse_objects(fasta_reads, static_cast<uint64_t>(-1));

    const char *reference = fasta_reads[0]->get_data();
    const uint32_t length = fasta_reads[0]->get_data_length();

    vector<winnowing::minimizer> minimizers;
    unordered_map<winnowing::minhash_t, vector<uint32_t>> lookup_table;
    winnowing::index_sequence(reference, length, winnowing::DEFAULT_W, winnowing::DEFAULT_K, minimizers, lookup_table);

    cout << "Reference sequence successfully indexed." << endl;

    return 0;
}