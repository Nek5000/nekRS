set(NO_INST_PROGS
  triv net_send net_recv multi_send derived_send derived_recv 
  transform_recv derived_send2 transform_recv2 transform_recv3
  dfg_master dfg_client dfg_master2 dfg_client3 dfg_master3 dfg_master4
)

foreach(PROG ${NO_INST_PROGS})
  add_executable(${PROG} ${PROG}.c)
  target_link_libraries(${PROG} evpath)
endforeach()

