set(oceanography_lib_src oceanography.c)
link_libraries(m)

add_library(oceanography SHARED ${oceanography_lib_src})

set_target_properties(oceanography
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 1)

install(TARGETS oceanography
        DESTINATION lib
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)

install(FILES oceanography.h
        DESTINATION include
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
