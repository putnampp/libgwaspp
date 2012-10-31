SET(SRCS chromosome.cpp
	chromosome_lookup_table.cpp
	chromosome_table.cpp
	chromosomal_interval.cpp)

SET(HEADERS chromosomal_interval.h
	chromosome_exceptions.h
	chromosome.h
	chromosome_lookup_table.h
	chromosome_table.h)
#FOREACH(SRC ${SRCS})
#	GET_FILENAME_COMPONENT(SRC_NAME ${SRC} NAME_WE)
#	add_library(${SRC_NAME} OBJECT ${SRC})
#ENDFOREACH(SRC)

add_library(chromosome OBJECT ${SRCS})
