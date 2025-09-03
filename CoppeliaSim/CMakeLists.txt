cmake_minimum_required(VERSION 3.5)
project(zmqRemoteApi-client-cpp3)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 17)

if(NOT COPPELIASIM_INCLUDE_DIR)
    if(DEFINED ENV{COPPELIASIM_ROOT_DIR})
        set(COPPELIASIM_INCLUDE_DIR $ENV{COPPELIASIM_ROOT_DIR}/programming/include)
    else()
        message(FATAL_ERROR "Environment variable COPPELIASIM_ROOT_DIR is not set")
    endif()
endif()

list(APPEND CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules
    ${COPPELIASIM_INCLUDE_DIR}/cmake
)
find_package(CoppeliaSim 4.1.0.0 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_package(OpenCV COMPONENTS core imgproc highgui calib3d REQUIRED)
message(STATUS "OpenCV libs: ${OpenCV_LIBS}")
find_package(Qt5 COMPONENTS Core Widgets Gui)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Inclui pasta de headers da biblioteca apriltags
include_directories(/home/impossivel/apriltags/AprilTags)

set(GENERATE OFF CACHE BOOL "Generate wrappers for objects and methods (requires CoppeliaSim to be running).")
set(GENERATE_INCLUDE_OBJECTS "sim,simIK" CACHE STRING "")
set(GENERATE_EXCLUDE_OBJECTS "" CACHE STRING "")
set(GENERATE_EXCLUDE_METHODS "sim.test,sim.auxFunc,sim.getStringSignal,sim.getInt32Signal,sim.getFloatSignal" CACHE STRING "")

include(FetchContent)

FetchContent_Declare(jsoncons
    GIT_REPOSITORY https://github.com/danielaparker/jsoncons
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/jsoncons
)
FetchContent_GetProperties(jsoncons)
if(NOT jsoncons_POPULATED)
    FetchContent_Populate(jsoncons)
endif()

set(CPPZMQ_BUILD_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(cppzmq
    GIT_REPOSITORY https://github.com/zeromq/cppzmq
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/cppzmq
)
FetchContent_GetProperties(cppzmq)
if(NOT cppzmq_POPULATED)
    FetchContent_Populate(cppzmq)
    add_subdirectory(${cppzmq_SOURCE_DIR} ${cppzmq_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/generated)

# Geração de JSONs (calltips, constants)
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/calltips.json"
    COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH="${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/../../tools/get_raw_calltips.py"
        "${CMAKE_CURRENT_BINARY_DIR}/generated/calltips.json"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../../tools/get_raw_calltips.py"
)
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/constants.json"
    COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH="${CMAKE_CURRENT_SOURCE_DIR}/../python/src"
        ${Python3_EXECUTABLE}
        "${CMAKE_CURRENT_SOURCE_DIR}/../../tools/get_constants.py"
        "${CMAKE_CURRENT_BINARY_DIR}/generated/constants.json"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../../tools/get_constants.py"
)

set(generatedFiles)
file(GLOB templateFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/templates/ ${CMAKE_CURRENT_SOURCE_DIR}/templates/*)
foreach(templateFile ${templateFiles})
    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/${templateFile}"
        COMMAND ${Python3_EXECUTABLE}
            "${COPPELIASIM_INCLUDE_DIR}/simStubsGen/external/pycpp/pycpp.py"
            -p "calltips_json=${CMAKE_CURRENT_BINARY_DIR}/generated/calltips.json"
            -p "constants_json=${CMAKE_CURRENT_BINARY_DIR}/generated/constants.json"
            -p "include_objects=${GENERATE_INCLUDE_OBJECTS}"
            -p "exclude_objects=${GENERATE_EXCLUDE_OBJECTS}"
            -p "exclude_methods=${GENERATE_EXCLUDE_METHODS}"
            -i "${CMAKE_CURRENT_SOURCE_DIR}/templates/${templateFile}"
            -o "${CMAKE_CURRENT_BINARY_DIR}/generated/${templateFile}"
            -P "${COPPELIASIM_PYTHON_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/../../tools"
            -P "${CMAKE_CURRENT_SOURCE_DIR}"
        DEPENDS
            "${COPPELIASIM_INCLUDE_DIR}/simStubsGen/external/pycpp/pycpp.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/templates/${templateFile}"
            "${CMAKE_CURRENT_BINARY_DIR}/generated/calltips.json"
            "${CMAKE_CURRENT_BINARY_DIR}/generated/constants.json"
            "${CMAKE_CURRENT_SOURCE_DIR}/cpp_utils.py"
            "${COPPELIASIM_PYTHON_DIR}/calltip.py"
            "${COPPELIASIM_PYTHON_DIR}/calltip.lark"
    )
    list(APPEND generatedFiles "${CMAKE_CURRENT_BINARY_DIR}/generated/${templateFile}")
endforeach()

add_custom_target(generate_code DEPENDS ${generatedFiles})

add_library(RemoteAPIClient STATIC RemoteAPIClient.cpp)
if(GENERATE)
    add_dependencies(RemoteAPIClient generate_code)
    target_include_directories(RemoteAPIClient BEFORE PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/generated)
else()
    set_source_files_properties(RemoteAPIClient.h OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RemoteAPIObjects.h)
    set_source_files_properties(RemoteAPIClient.cpp OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/RemoteAPIObjects.cpp)
endif()

target_compile_definitions(RemoteAPIClient PUBLIC -DSIM_REMOTEAPICLIENT_OBJECTS)
target_include_directories(RemoteAPIClient PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/jsoncons/include)
target_link_libraries(RemoteAPIClient PUBLIC cppzmq)

# Executáveis e linkagem com apriltags
set(APRILTAG_LIB /home/impossivel/apriltags/build/lib/libapriltags.a)

add_executable(example example.cpp)
target_include_directories(example PRIVATE /home/impossivel/apriltags/AprilTags)
target_link_libraries(example PRIVATE
    RemoteAPIClient
    Eigen3::Eigen
    v4l2
    ${APRILTAG_LIB}
    ${OpenCV_LIBS}   
)

add_executable(synchronousImageTransmission synchronousImageTransmission.cpp)
target_include_directories(synchronousImageTransmission PRIVATE /home/impossivel/apriltags/AprilTags)
target_link_libraries(synchronousImageTransmission PRIVATE RemoteAPIClient ${APRILTAG_LIB} Eigen3::Eigen)

add_executable(bubbleRobClient bubbleRobClient.cpp)
target_include_directories(bubbleRobClient PRIVATE /home/impossivel/apriltags/AprilTags)
target_link_libraries(bubbleRobClient PRIVATE RemoteAPIClient ${APRILTAG_LIB} Eigen3::Eigen)

if(Qt5_FOUND)
    add_executable(qt qt.cpp)
    target_include_directories(qt PRIVATE /home/impossivel/apriltags/AprilTags)
    target_link_libraries(qt PRIVATE RemoteAPIClient Qt5::Widgets ${APRILTAG_LIB})
endif()

