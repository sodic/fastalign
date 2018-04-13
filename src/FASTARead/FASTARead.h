//
// Created by filip on 03.04.18..
//

#ifndef ZAVRSNI_FASTAREAD_H
#define ZAVRSNI_FASTAREAD_H

#include <string>

/**
 * The class FASTARead encapsulates information about a read taken from a FASTA file.
 * The class with this format was required by the library bioparser (https://github.com/rvaser/bioparser) used for reading data from FASTA files.
 */
class FASTARead {

    std::string name_{};

    std::string data_{};

public:
    FASTARead(const char *name,
              uint32_t name_length,
              const char *data,
              uint32_t data_length) : name_(name, name_length), data_(data, data_length) {
    }

    const char*  get_data();

    const char*  get_name();

    uint32_t get_data_length();

    uint32_t get_name_length();

};


#endif //ZAVRSNI_FASTAREAD_H
