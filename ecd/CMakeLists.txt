cmake_minimum_required(VERSION 2.6)


#INCLUDE_DIRECTORIES(${openbabel_SOURCE_DIR}/include/ga)
#LINK_DIRECTORIES(${openbabel_SOURCE_DIR}/include/ga)
#TARGET_LINK_LIBRARIES(ecdconf ga)

add_executable(ecdconf ecdconf.cpp)
target_link_libraries(ecdconf openbabel)
target_link_libraries(ecdconf /home/liu/soft/babel/lib/libopenbabel.so)
include_directories(/home/liu/soft/babel/include/openbabel-2.0)
target_link_libraries(ecdconf ga)
target_link_libraries(ecdconf /home/liu/projects/conformation/GA/galib247/ga/libga.so)
include_directories(/home/liu/projects/conformation/GA/galib247)
