include(FetchContent)

cmake_minimum_required(VERSION 3.0)

project(cds-a5)


include_directories("include")

#----------------------------------------
# Add sdsl
#----------------------------------------
FetchContent_Declare(
  sdsl
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")

  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
# compile with address sanitizer
# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)

add_executable(build src/build.cpp)
target_link_libraries(build sdsl)
target_compile_options(build PUBLIC -O2)

add_executable(load src/deserialize.cpp)
target_link_libraries(load sdsl)
target_compile_options(load PUBLIC -O2)

add_executable(movelet src/lf_movelet.cpp)
target_link_libraries(movelet sdsl)
target_compile_options(movelet PUBLIC -O2)

add_executable(lf src/lf.cpp)
target_link_libraries(lf sdsl)
target_compile_options(lf PUBLIC -O2)

add_executable(BWT_inverse_lf src/BWT_inverse_lf.cpp)
target_link_libraries(BWT_inverse_lf sdsl)
target_compile_options(BWT_inverse_lf PUBLIC -O2)