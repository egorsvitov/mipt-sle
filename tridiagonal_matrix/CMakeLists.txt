cmake_minimum_required(VERSION 3.14)

project(tester)

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)

FetchContent_MakeAvailable(googletest)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(tester tester.cpp)
target_link_libraries(tester gtest_main ${CMAKE_CURRENT_SOURCE_DIR}/tridiag/build/libtridiag.a)
add_test(NAME simple_test COMMAND tester)