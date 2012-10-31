include_directories(${BIOLIB2_SOURCE_DIR}/src/biolib2)

SET(SRCS qualitative_trait.cpp 
	quantitative_trait.cpp 
	trait.cpp)

#FOREACH(SRC ${SRCS})
#	GET_FILENAME_COMPONENT(SRC_NAME ${SRC} NAME_WE)
#	#add_library(${SRC_NAME} OBJECT ${SRC})
#ENDFOREACH(SRC)

ADD_LIBRARY(trait OBJECT ${SRCS})
