cmake_minimum_required(VERSION 3.10)

# Set the project name
project(HMM_NRP_Matcher)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP)

# Add the executable
add_executable(hmm_nrp_matcher
        src/cpp_matching/main.cpp
        src/cpp_matching/parsing/parse_hmms.cpp
        src/cpp_matching/parsing/parse_nrps.cpp
        src/cpp_matching/matching/matcher.cpp
        src/cpp_matching/matching/viterbi_algorithm.cpp
        src/cpp_matching/write_matches.cpp
        )

# Include directories for header files
target_include_directories(hmm_nrp_matcher PRIVATE src/cpp_matching)

if (OPENMP_FOUND)
    target_link_libraries(hmm_nrp_matcher PUBLIC OpenMP::OpenMP_CXX)
else ()
    if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        message(FATAL_ERROR "Nerpa in the fast mode (C++-based) requires OpenMP to be available")
    endif()
endif()
