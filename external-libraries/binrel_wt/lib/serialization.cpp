#include <brwt/serialization.h>

#include <sdsl/int_vector.hpp>

namespace brwt {

void serialize_number(std::ostream &out, uint64_t number) {
    if (!out.good()) {
        throw std::ostream::failure("Bad stream");
    }
    out.write(reinterpret_cast<const char *>(&number), sizeof(number));
}

uint64_t load_number(std::istream &in) {
    if (!in.good()) {
        throw std::istream::failure("Bad stream");
    }
    uint64_t value;
    in.read(reinterpret_cast<char *>(&value), sizeof(value));
    return value;
}

void serialize_number_vector(std::ostream &out,
                             const std::vector<std::uint_fast64_t> &vector) {
    if (!out.good()) {
        throw std::ostream::failure("Bad stream");
    }
    sdsl::int_vector<> int_vector(vector.size(), 0, sizeof(std::uint_fast64_t) * 8);
    for (size_t i = 0; i < vector.size(); ++i) {
        int_vector[i] = vector[i];
    }
    int_vector.serialize(out);
}

uint64_t get_number_vector_size(std::istream &in) {
    if (!in.good())
        throw std::istream::failure("Bad stream");

    // save the position
    auto position = in.tellg();

    typename sdsl::int_vector<>::size_type size;
    typename sdsl::int_vector<>::int_width_type int_width;
    sdsl::int_vector<>::read_header(size, int_width, in);
    if (!int_width || size % int_width)
        throw std::istream::failure("Error when loading size and int_width in vector");

    // restore the initial position in stream
    in.seekg(position, in.beg);

    return size / int_width;
}

std::vector<std::uint_fast64_t> load_number_vector(std::istream &in) {
    if (!in.good()) {
        throw std::istream::failure("Bad stream");
    }
    try {
        sdsl::int_vector<> int_vector;
        int_vector.load(in);
        return std::vector<std::uint_fast64_t>(int_vector.begin(), int_vector.end());
    } catch (...) {
        throw std::istream::failure("Bad stream");
    }
}

} // namespace brwt
