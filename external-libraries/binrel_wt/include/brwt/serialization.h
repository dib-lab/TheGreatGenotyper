#ifndef SERIALIZATION_BRWT_H
#define SERIALIZATION_BRWT_H

#include <vector>
#include <iostream>

namespace brwt {

void serialize_number(std::ostream &out, uint64_t number);

uint64_t load_number(std::istream &in);

void serialize_number_vector(std::ostream &out,
                             const std::vector<std::uint_fast64_t> &vector);

std::vector<std::uint_fast64_t> load_number_vector(std::istream &in);

uint64_t get_number_vector_size(std::istream &in);

} // namespace brwt

#endif // SERIALIZATION_BRWT_H
