cmake_minimum_required(VERSION 3.10)
include(FetchContent)   # <--- REQUIRED for FetchContent_Declare!

# Set extension name here
set(TARGET_NAME miint)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# from stobealign, seem needed for kseq++
# note: bumped from c++17 to 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# DuckDB's extension distribution supports vcpkg. As such, dependencies can be added in ./vcpkg.json and then
# used in cmake with find_package. Feel free to remove or replace with other dependencies.
# Note that it should also be removed from vcpkg.json to prevent needlessly installing it..
find_package(ZLIB)
find_package(Threads REQUIRED)

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})

# adapted kseq++ include use from strobealign
include_directories(src/include ext)
set(EXTENSION_SOURCES
    src/SequenceReader.cpp
    src/QualScore.cpp
    src/read_fastx.cpp
    src/miint_extension.cpp
    src/QualScore.cpp)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})


target_link_libraries(${EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB)
target_link_libraries(${LOADABLE_EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB)

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")

FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(catch2)
set(TEST_SOURCES
    src/SequenceReader.cpp
    test/cpp/test_SequenceReader.cpp
    src/QualScore.cpp
    test/cpp/test_QualScore.cpp
)

add_executable(tests ${TEST_SOURCES})

target_link_libraries(tests
  Catch2::Catch2WithMain
  Threads::Threads
  ZLIB::ZLIB
)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
