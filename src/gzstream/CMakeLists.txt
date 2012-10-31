include_directories(${PROJECT_SOURCE_DIR}/src/gzstream)

SET(SRCS gzstream.cpp)

#FOREACH(SRC ${SRCS})
#	GET_FILENAME_COMPONENT(SRC_NAME ${SRC} NAME_WE)
#	add_library(${SRC_NAME} OBJECT ${SRC})
#ENDFOREACH(SRC)

ADD_LIBRARY(gzstream ${SRCS})
TARGET_LINK_LIBRARIES( gzstream ${ZLIB_LIBRARY} )
