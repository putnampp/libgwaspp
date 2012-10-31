include_directories(${BIOLIB2_SOURCE_DIR}/src/biolib2)

SET(SRCS lookup_table.cpp)

#FOREACH(SRC ${SRCS})
#	GET_FILENAME_COMPONENT(SRC_NAME ${SRC} NAME_WE)
#	add_library(${SRC_NAME} OBJECT ${SRC})
#ENDFOREACH(SRC)

ADD_LIBRARY(lookup_table OBJECT ${SRCS})
