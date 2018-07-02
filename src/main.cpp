#include <iostream>
#include <sstream>
#include <winnowing_dac/winnowing.hpp>
#include <config/config.hpp>
#include <mapper/mapper.hpp>
#include "mapper/mapper.hpp"
#include "bioparser/bioparser.hpp"
#include "statistics/statistics.hpp"
#include "FASTARead/FASTARead.h"

using namespace std;

const char *query_name;
const char *ref_name;
uint32_t query_length;
uint32_t ref_length;


string join_file(vector<unique_ptr<FASTARead>> &reads) {
    stringstream ss;
    for (unique_ptr<FASTARead> &read : reads) {
        ss << read->get_data();
    }
    return ss.str();
}


void print_mappings(const vector<mapper::Mapping> &mappings, FILE *output) {
    for (const mapper::Mapping &mapping: mappings) {
        fprintf(output, "%s %d %d %d %c %s %d %d %d %lf\n",
                query_name,
                query_length,
                mapping.query_start,
                mapping.query_end,
                mapping.strand ? '+' : '-',
                ref_name,
                ref_length,
                mapping.ref_start,
                mapping.ref_end,
                mapping.identity_estimate
        );
    }
}

int main(int argc, char const *argv[]) {

    vector<unique_ptr<FASTARead>> query_reads;
    auto query_parser = bioparser::createParser<bioparser::FastaParser, FASTARead>(argv[1]);
    query_parser->parse_objects(query_reads, static_cast<uint64_t>(-1));
    string q_str = join_file(query_reads);
    query_length = q_str.size();
    query_name = query_reads[0]->get_name();
    const char *query = q_str.c_str();

    vector<unique_ptr<FASTARead>> ref_reads;
    auto reference_parser = bioparser::createParser<bioparser::FastaParser, FASTARead>(argv[2]);
    reference_parser->parse_objects(ref_reads, static_cast<uint64_t>(-1));
    string ref_str = join_file(ref_reads);
    ref_name = ref_reads[0]->get_name();
    ref_length = ref_str.size();
    const char *reference = ref_str.c_str();


    vector<winnowing::minimizer> ref_minimizers;
    unordered_map<winnowing::minhash_t, vector<uint32_t>> lookup_table;

    std::vector<mapper::Mapping> mappings;
    mapper::compute_mappings(reference, ref_length, query, query_length, mappings);
    FILE *output = fopen("fastalign.out", "w");
    print_mappings(mappings, output);
    return 0;
}