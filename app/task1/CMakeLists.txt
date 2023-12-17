#Compile options

add_executable(${TASK_1_NAME}
        ${TASK_1_MAIN}
        ${TASK_1_SOURCES}
)
set_target_properties(${TASK_1_NAME} PROPERTIES ENABLE_EXPORTS 1)

target_sources(${TASK_1_NAME}
        PUBLIC
        ${TASK_1_SOURCES}
)

target_include_directories(${TASK_1_NAME}
        PUBLIC
        ${TASK_1_INCLUDE_PATH}
)

target_link_libraries(${TASK_1_NAME}
        PRIVATE
        eigen
)