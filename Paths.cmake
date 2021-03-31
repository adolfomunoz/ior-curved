####################################################################
# PATHS
####################################################################
if (NOT EXECUTABLE_OUTPUT_PATH)
	set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
endif()
if (NOT LIBRARY_OUTPUT_DIRECTORY)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
endif()


if (NOT IS_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
	file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
endif()
if (NOT IS_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
	file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
endif()
if (NOT IS_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

