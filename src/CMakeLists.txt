set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/fht)

set(FILE_LIST
   fast_hankel_transform.cc
)

pybind11_add_module(libcxx ${FILE_LIST})

target_link_libraries(libcxx
   PRIVATE
   ${Boost_LIBRARIES}
   fmt-header-only
   ${FFTW_LIBRARIES}
)

target_include_directories(libcxx
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)
