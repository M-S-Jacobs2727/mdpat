find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_definitions(-DBOOST_TEST_DYN_LINK)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../src")
#message("Include directory: ${INCLUDE_DIR}")

add_executable(mdpat-test
    # testReadTrajectories.cpp
    testPermuteDims.cpp
)
target_include_directories(mdpat-test
  PRIVATE
    ${INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(mdpat-test
  PRIVATE
    mdpatlib
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
