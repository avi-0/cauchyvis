cmake_minimum_required(VERSION 3.10)
project(cauchyvis CXX)

include(FetchContent) 
FetchContent_Declare(mahi-gui GIT_REPOSITORY https://github.com/mahilab/mahi-gui.git) 
FetchContent_MakeAvailable(mahi-gui)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(cauchyvis cauchyvis.cpp)
target_link_libraries(cauchyvis mahi::gui)
target_link_libraries(cauchyvis Eigen3::Eigen)
