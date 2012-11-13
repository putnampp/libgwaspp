include_directories(${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIR})

# Utility classes
LIST(APPEND SRCS interval/interval.cpp)
LIST(APPEND SRCS lookup_table/lookup_table.cpp)
LIST(APPEND SRCS category/category.h)
LIST(APPEND SRCS index_set/index_set.cpp)
LIST(APPEND SRCS time/timing.cpp)

ADD_LIBRARY(util SHARED ${SRCS})
#TARGET_LINK_LIBRARIES(util ${OBJS})
