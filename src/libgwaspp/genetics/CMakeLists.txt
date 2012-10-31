SET(DIRS chromosome genotype individual marker phenotype trait)
#SET(DIRS chromosome)

include_directories(${PROJECT_SOURCE_DIR}/src)
#link_directories(${PROJECT_BINARY_DIR}/lib)

SET(SRCS genetic_data.cpp)

FOREACH(DIR ${DIRS})
	ADD_SUBDIRECTORY(${DIR})
	LIST(APPEND TAROBJS $<TARGET_OBJECTS:${DIR}>)
ENDFOREACH(DIR)

add_library(genetics SHARED ${SRCS} ${TAROBJS})
target_link_libraries(genetics util gzstream)
