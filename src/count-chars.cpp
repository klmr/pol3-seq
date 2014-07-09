#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>

using map_type = std::unordered_map<int, int>;

map_type count_chars_per_line(std::istream& input, char character) {
    auto map = map_type{};

    auto line = std::string{};
    while (getline(input, line)) {
        auto count = std::count(begin(line), end(line), character);
        ++map[count];
    }

    return map;
}

void print_histogram(map_type const& map) {
    for (auto&& i : map)
        std::cout << i.first << '\t' << i.second << '\n';
}

int main(int argc, char** argv) {
    std::ios_base::sync_with_stdio(false);

    if (not (argc > 1 and std::strlen(argv[1]) > 0)) {
        std::cerr << "Invalid arguments\n";
        return EXIT_FAILURE;
    }

    auto character = argv[1][0];
    auto map = count_chars_per_line(std::cin, character);

    print_histogram(map);
}
