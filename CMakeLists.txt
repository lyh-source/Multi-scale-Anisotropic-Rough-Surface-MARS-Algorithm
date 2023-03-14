if(NOT PROJECT_NAME)
    cmake_minimum_required(VERSION 3.1...3.15)
    project(MARS CXX C)
endif()

set(CMAKE_C_STANDARD 14)
set(CMAKE_CXX_STANDARD 14)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

SET(MARS_INCLUDE_DIR include)
include_directories( ${MARS_INCLUDE_DIR} )

file (GLOB MAIN_FILE
    src/main.cpp
)

file (GLOB MARS_FILE
    src/mars/*
)

add_executable(${PROJECT_NAME}
    	${MAIN_FILE}
	${MARS_FILE}
)



if (CMAKE_SYSTEM_NAME MATCHES "Windows")
	find_package(OpenCV REQUIRED)
	SET(OpenCV_LIBRARIES ${OpenCV_LIBS})
	INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
	target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
#elseif(CMAKE_SYSTEM_NAME MATCHES "Emscripten")
endif()







