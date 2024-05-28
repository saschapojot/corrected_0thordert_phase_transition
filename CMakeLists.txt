cmake_minimum_required(VERSION 3.27)
project(t_phase_transition)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories("/home/polya/cppLocal/armadillo/include")
include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")
set(Python3_INCLUDE_DIRS "/home/polya/cppLocal/python/include/python3.12")
set(Python3_LIBRARIES "/home/polya/cppLocal/python/lib/libpython3.12.a")
include_directories(${Python3_INCLUDE_DIRS})



add_executable(run1d run1d.cpp 1d/1d.cpp)
add_executable(genA genA.cpp version1/LJPotPBC/version1LJPotPBC2Atom.cpp)




add_executable(runV1LJ2AtomPBC runV1LJ2AtomPBC.cpp version1/LJPotPBC/version1LJPotPBC2Atom.cpp)
add_executable(runV1LJ2AtomParseXMLPBC runV1LJ2AtomParseXMLPBC.cpp version1/LJPotPBC/parseXML.cpp version1/LJPotPBC/version1LJPotPBC2Atom.cpp)

set_source_files_properties(
        1d/1d.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)





set_source_files_properties(
        version1/LJPotPBC/version1LJPotPBC2Atom.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set_source_files_properties(
        version1/LJPotPBC/parseXML.cpp
        PROPERTIES

        COMPILE_OPTIONS "-O3"
)

set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

#find_package(Boost  REQUIRED COMPONENTS ALL)
find_package(Boost  REQUIRED serialization filesystem json python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()


file(GLOB LIBS_SO "/home/polya/cppLocal/armadillo/lib/*.so")
target_link_libraries(run1d ${LIBS_SO})
target_link_libraries(genA ${LIBS_SO})
target_link_libraries(runV1LJ2AtomPBC ${LIBS_SO})
target_link_libraries(runV1LJ2AtomParseXMLPBC ${LIBS_SO})

target_link_libraries(run1d Boost::serialization)
target_link_libraries(run1d Boost::filesystem)
target_link_libraries(run1d Boost::json)


target_link_libraries(genA Boost::serialization)
target_link_libraries(genA Boost::filesystem)
target_link_libraries(genA Boost::json)
target_link_libraries(genA Boost::python)
target_link_libraries(genA ${Python3_LIBRARIES})




target_link_libraries(runV1LJ2AtomPBC Boost::serialization)
target_link_libraries(runV1LJ2AtomPBC Boost::filesystem)
target_link_libraries(runV1LJ2AtomPBC Boost::json)
target_link_libraries(runV1LJ2AtomPBC Boost::python)
target_link_libraries(runV1LJ2AtomPBC ${Python3_LIBRARIES})

target_link_libraries(runV1LJ2AtomParseXMLPBC Boost::serialization)
target_link_libraries(runV1LJ2AtomParseXMLPBC Boost::filesystem)
target_link_libraries(runV1LJ2AtomParseXMLPBC Boost::json)
target_link_libraries(runV1LJ2AtomParseXMLPBC Boost::python)
target_link_libraries(runV1LJ2AtomParseXMLPBC ${Python3_LIBRARIES})
