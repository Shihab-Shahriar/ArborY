cmake_minimum_required(VERSION 3.23)

# Set your project name
project(distance)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


include(FetchContent)
FetchContent_Declare(
  kokkos
  GIT_REPOSITORY https://github.com/kokkos/kokkos.git
  GIT_TAG        master 
)

# Make Kokkos available for `add_subdirectory`
FetchContent_MakeAvailable(kokkos)

add_executable(app program.cpp)
target_link_libraries(app Kokkos::kokkos)
target_include_directories(app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# For testing purpose
find_package(CGAL REQUIRED COMPONENTS Core)
add_executable(cgal cgall.cpp)
target_link_libraries(cgal CGAL::CGAL CGAL::CGAL_Core)
