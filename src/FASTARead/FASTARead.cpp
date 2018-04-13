//
// Created by filip on 03.04.18..
//

#include <cstdint>
#include "FASTARead.h"

const char* FASTARead::get_name(){
    return name_.c_str();
};

const char* FASTARead::get_data(){
    return data_.c_str();
};

uint32_t FASTARead::get_data_length(){
    return static_cast<uint32_t>(data_.size());
}

uint32_t FASTARead::get_name_length() {
    return static_cast<uint32_t>(name_.size());
}