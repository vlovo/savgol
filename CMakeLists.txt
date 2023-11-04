# main CMakeList.txt
# (c) vlovo for savgol

cmake_minimum_required(VERSION 3.15)

set(ProjectName  savgol)

project(${ProjectName})

set(SRC  
main.cpp
)

set(HDR  
)

add_library(savgol INTERFACE) 
target_include_directories(savgol INTERFACE include/) 
target_compile_features(savgol INTERFACE cxx_std_17)

add_executable(main ${SRC} ${HDR})
target_link_libraries(main savgol) 
target_compile_features(main PRIVATE cxx_std_17) 
 