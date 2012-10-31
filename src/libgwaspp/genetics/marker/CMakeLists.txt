include_directories(${BIOLIB2_SOURCE_DIR}/src/biolib2)

SET(SRCS allele_lookup_table.cpp 
	marker.cpp 
	snp_lookup_table.cpp 
	snp_marker.cpp)

#FOREACH(SRC ${SRCS})
#	GET_FILENAME_COMPONENT(SRC_NAME ${SRC} NAME_WE)
#	#add_library(${SRC_NAME} OBJECT ${SRC})
#ENDFOREACH(SRC)

ADD_LIBRARY(marker OBJECT ${SRCS})

