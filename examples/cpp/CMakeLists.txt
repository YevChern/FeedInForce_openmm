cmake_minimum_required(VERSION 2.8)

# We need to know where OpenMM is installed so we can access the headers and libraries.
SET(OPENMM_DIR "/home/yevhen/programs/openmm" CACHE PATH "Where OpenMM is installed")
INCLUDE_DIRECTORIES("${OPENMM_DIR}/include")
LINK_DIRECTORIES("${OPENMM_DIR}/lib" "${OPENMM_DIR}/lib/plugins")

project(template_openmm)
add_executable(${PROJECT_NAME} "main.cpp")

target_link_libraries(template_openmm OpenMM FeedInPlugin)
