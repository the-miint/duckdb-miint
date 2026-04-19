# Bundle sortmerna's transitive OBJECT library outputs into a single static
# archive. Upstream's smr_api.a only contains smr_objs' objects because CMake
# doesn't propagate OBJECT library objects transitively through STATIC link.
# We pick up cmph/alp/build_version .o files from the upstream build tree and
# archive them together with libsmr_api.a into libsortmerna_bundle.a.
#
# Required input variables:
#   SORTMERNA_BIN_DIR      — sortmerna's CMake binary dir
#   SORTMERNA_API_LIB      — path to libsmr_api.a produced upstream
#   SORTMERNA_BUNDLE_LIB   — output path for the bundled archive
#   CMAKE_AR               — ar binary

file(GLOB CMPH_OBJS "${SORTMERNA_BIN_DIR}/3rdparty/cmph/CMakeFiles/cmph.dir/*.o")
file(GLOB ALP_OBJS "${SORTMERNA_BIN_DIR}/3rdparty/alp/CMakeFiles/alp.dir/*.o")
file(GLOB BV_OBJS "${SORTMERNA_BIN_DIR}/CMakeFiles/build_version.dir/*.cpp.o")

if(NOT CMPH_OBJS)
    message(FATAL_ERROR "No cmph object files found under ${SORTMERNA_BIN_DIR}/3rdparty/cmph")
endif()
if(NOT ALP_OBJS)
    message(FATAL_ERROR "No alp object files found under ${SORTMERNA_BIN_DIR}/3rdparty/alp")
endif()
if(NOT BV_OBJS)
    message(FATAL_ERROR "No build_version object files found under ${SORTMERNA_BIN_DIR}/CMakeFiles/build_version.dir")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SORTMERNA_API_LIB} ${SORTMERNA_BUNDLE_LIB}
    RESULT_VARIABLE rc
)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "Failed to seed ${SORTMERNA_BUNDLE_LIB} from ${SORTMERNA_API_LIB}")
endif()

execute_process(
    COMMAND ${CMAKE_AR} rs ${SORTMERNA_BUNDLE_LIB} ${CMPH_OBJS} ${ALP_OBJS} ${BV_OBJS}
    RESULT_VARIABLE rc
)
if(NOT rc EQUAL 0)
    message(FATAL_ERROR "ar failed to append OBJECT library objects to ${SORTMERNA_BUNDLE_LIB}")
endif()

message(STATUS "sortmerna bundle ready: ${SORTMERNA_BUNDLE_LIB}")
