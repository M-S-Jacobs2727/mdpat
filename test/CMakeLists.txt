find_package(MPI REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_definitions(-DBOOST_TEST_DYN_LINK)

find_file(MPI_CXX NAMES mpiCC mpicxx)
set(CMAKE_CXX_COMPILER ${MPI_CXX})

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")
#message("Include directory: ${INCLUDE_DIR}")


file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(testSrc ${TEST_SRCS})
    get_filename_component(testName ${testSrc} NAME_WE)
    add_executable(${testName}
        ${testSrc}
    )
    target_include_directories(${testName}
      PRIVATE
        ${INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
    )
    target_link_libraries(${testName}
      PRIVATE
        mdpatlib
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )
    add_test(NAME ${testName}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND mpirun -n 2 ${CMAKE_CURRENT_BINARY_DIR}/${testName} )
  
endforeach(testSrc)
