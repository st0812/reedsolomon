cmake_minimum_required(VERSION 3.16)
project(rs C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")

# CppUTestを追加
find_package(CppUTest REQUIRED)
find_library(GF_LIB gf PATHS lib)
find_library(POLY_LIB poly PATHS lib)
find_library(MATRIX_LIB matrix PATHS lib)


# テスト対象
add_library(rs src/RS.c)

# テスト実行ファイル
add_executable(test_runner
	test/TestRS.cpp
    test/TestRunner.cpp
)

target_include_directories(rs PUBLIC include)
target_include_directories(rs PUBLIC lib)

## リンク
target_link_libraries(rs PRIVATE
    ${GF_LIB}
    ${POLY_LIB}
    ${MATRIX_LIB}
)

target_link_libraries(test_runner PRIVATE
    rs 
    CppUTest
    CppUTestExt
)

# テスト登録
enable_testing()
add_test(NAME rs COMMAND test_runner)
