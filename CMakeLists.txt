cmake_minimum_required(VERSION 3.5)
project(qfib)

file(GLOB_RECURSE source_files sources/*.cpp)
add_executable(qfib ${source_files})


set_property(TARGET qfib PROPERTY CXX_STANDARD 11)

target_include_directories(qfib	PUBLIC eigen3)

if (USE_OPENMP)
	find_package(OpenMP)
	target_link_libraries(qfib PUBLIC OpenMP::OpenMP_CXX)
	add_compile_definitions(USE_OPENMP)
	message("Using OpenMP")
else()
	find_package(Threads)
	target_link_libraries (qfib ${CMAKE_THREAD_LIBS_INIT})
endif()
