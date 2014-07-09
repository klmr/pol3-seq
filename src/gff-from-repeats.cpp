#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <sstream>
#include <tuple>
#include <vector>

#include "util/lang/range.hpp"
using util::lang::range;

// RepeatMasker file format:
//
//    SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
// score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID
// 
// 12937  10.6  1.3  0.2  chr1      3000001 3000097 (192471874) C  L1_Mus3        LINE/L1             (3055) 3592   3487      1
//    27   0.0  0.0  0.0  chr1      3000098 3000123 (192471848) +  (T)n           Simple_repeat            1   26    (0)      2
// 12937  10.6  1.3  0.2  chr1      3000124 3002128 (192469843) C  L1_Mus3        LINE/L1             (3161) 3486   1467      1

enum class strand : char {
    positive = '+',
    negative = '-',
    unknown = 'C'
};

std::istream& operator >>(std::istream& in, strand& value) {
    auto val = char{};
    if (in >> val)
        value = static_cast<strand>(val);
    return in;
}

std::ostream& operator <<(std::ostream& out, strand value) {
    return out << static_cast<char>(value);
}

struct repeat {
    using pos_t = unsigned int;
    int score;
    // float perc_div, perc_del, perc_ins
    std::string reference;
    pos_t begin;
    pos_t end;
    // pos_t left;
    strand strand;
    std::string matching_repeat;
    std::string family;
    // pos_t begin, end, left; // position in repeat
    int id;

    repeat() = default;
    repeat(repeat const&) = default;

    static int const format_flag_i;
};

int const repeat::format_flag_i = std::ios_base::xalloc();

struct repeat_annotation {
    std::vector<repeat> repeats;
};

struct repeat_input_helper {
    std::unordered_map<std::string, int>& ids;
    repeat& value;

    repeat_input_helper(std::unordered_map<std::string, int>& ids, repeat& value)
        : ids{ids}, value{value} {}

    int get_next(std::string const& type) {
        return ++ids[type];
    }
};

std::istream& operator >>(std::istream& in, repeat& value);

std::istream& operator >>(std::istream& in, repeat_input_helper&& value) {
    if (in >> value.value)
        value.value.id = value.get_next(value.value.family);

    return in;
}

std::istream& operator >>(std::istream& in, repeat_annotation& value) {
    auto ids = std::unordered_map<std::string, int>{};
    auto r = repeat{};

    value.repeats.clear();

    while (in >> repeat_input_helper{ids, r})
        value.repeats.push_back(std::move(r));

    return in;
}

std::ostream& operator <<(std::ostream& out, repeat_annotation const& value) {
    std::copy(begin(value.repeats), end(value.repeats),
        std::ostream_iterator<repeat>(out, "\n"));
    return out;
}

enum class file_formats {
    plain, gff, bed
};

std::istream& operator >>(std::istream& in, repeat& value) {
    auto line = std::string{};
    if (not getline(in, line))
        return in;

    std::istringstream istr{line};
    // Fields of no interest to skip
    auto ftmp = float{};
    auto itmp = int{};
    auto ctmp = char{};
    return istr >> value.score >> ftmp >> ftmp >> ftmp >>
        value.reference >> value.begin >> value.end >> ctmp >> itmp >>
        ctmp >> value.strand >> value.matching_repeat >> value.family;
}

void print_gff(std::ostream& out, repeat const& value) {
    auto gff_strand = [](strand s) {
        return s == strand::unknown ? '.' : static_cast<char>(s);
    };

    auto prefix = value.family.substr(0, value.family.find("/"));

    out << value.reference << '\t' << "RepeatMasker" <<
        '\t' << prefix <<
        '\t' << value.begin << '\t' << value.end <<
        '\t' << value.score <<
        '\t' << gff_strand(value.strand) << "\t." <<
        '\t' << prefix << "_ID=\"" << value.family << '-' << value.id << "\"" <<
        "; " << "type=\"" << value.family << "\"" <<
        "; " << "source=\"" << value.matching_repeat << "\"";
}

std::ostream& operator <<(std::ostream& out, repeat const& value) {
    switch (static_cast<file_formats>(out.iword(repeat::format_flag_i))) {
        case file_formats::plain:
            return out << value.reference << ':' << value.begin << '-' << value.end <<
                " (" << value.strand << ')' <<
                '\t' << value.family <<
                " (" << value.matching_repeat << ')' <<
                "\t[" << value.score << ']';
        case file_formats::gff:
            print_gff(out, value);
            return out;
        case file_formats::bed:
            throw "Not implemented";
    }
}

std::ostream& format_gff(std::ostream& out) {
    out.iword(repeat::format_flag_i) = static_cast<int>(file_formats::gff);
    return out;
}

std::ostream& format_plain(std::ostream& out) {
    out.iword(repeat::format_flag_i) = static_cast<int>(file_formats::plain);
    return out;
}

std::istream& skip_header(std::istream& in, unsigned int length) {
    auto constexpr infinity = std::numeric_limits<std::streamsize>::max();

    for (auto i : range(0u, length))
        in.ignore(infinity, '\n'), i;

    return in;
}

template <typename Ex> void fail(Ex ex) { throw ex; }

class invalid_input : std::exception { };

int main(int argc, char** argv)  try {
    (void) argc;
    (void) argv;
    auto repeats = repeat_annotation{};
    std::istream& infile = std::cin;

    if (not skip_header(infile, 3))
        fail(invalid_input{});

    infile >> repeats;
    std::cout << format_gff << repeats;
}
catch (std::exception const& ex) {
    std::cerr << ex.what() << '\n';
    return 1;
}
catch (...) {
    return 1;
}
