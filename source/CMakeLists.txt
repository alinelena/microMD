
set(microMD_SRCS 
microMD.F90
useful.F90
types.F90
parser.F90
control.F90
constants.F90
config.F90
field.F90
io.F90
potentials.F90
neighbours.F90
particles.F90
forces.F90
sampler.F90
units.F90
)

if(BUILDER)
  set_property(SOURCE setup.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "BUILDER='${BUILDER}'")
endif()
set_property(SOURCE setup.F90 APPEND PROPERTY
  COMPILE_DEFINITIONS "HOST='${HOST}'")

add_executable(${target_name} ${microMD_SRCS} )
set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${OMP_LINK_FLAGS}")
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
target_link_libraries(${target_name} ${LIBS})

