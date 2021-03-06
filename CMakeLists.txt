macro(vc_add_run_target _target)
   if("${_target}" MATCHES "_mic$")
      if(MIC_NATIVELOAD)
         get_target_property(_exe "${_target}" OUTPUT_NAME)
         add_custom_target(run_${_target}
            ${MIC_NATIVELOAD} "${_exe}"
            DEPENDS ${_target}
            COMMENT "Execute ${_target} example"
            VERBATIM
            )
      endif()
   else()
      add_custom_target(run_${_target}
         ${_target}
         DEPENDS ${_target}
         COMMENT "Execute ${_target} example"
         VERBATIM
         )
   endif()
endmacro()

macro(build_example name)
   set(_SRCS)
   set(_LIBS)
   set(_disabled_impl)
   set(_state 1)
   foreach(ARG ${ARGN})
      if(ARG STREQUAL "LIBS")
         set(_state 2)
      elseif(ARG STREQUAL "DISABLE")
         set(_state 3)
      elseif(_state EQUAL 1)
         set(_SRCS ${_SRCS} ${ARG})
      elseif(_state EQUAL 2)
         set(_LIBS ${_LIBS} ${ARG})
      elseif(_state EQUAL 3)
         list(APPEND _disabled_impl ${ARG})
      endif()
   endforeach()
   add_executable("${name}_default" ${_SRCS})
   add_target_property("${name}_default" COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS}")
   target_link_libraries("${name}_default" Vc ${_LIBS})

   list(FIND _disabled_impl "Scalar" _index)
   if(_index EQUAL -1)
      add_executable("${name}_scalar" ${_SRCS})
      add_target_property("${name}_scalar" COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS} -DVC_IMPL=Scalar")
      add_target_property("${name}_scalar" LABELS "Scalar")
      add_dependencies(Scalar "${name}_scalar")
      target_link_libraries("${name}_scalar" Vc ${_LIBS})
      vc_add_run_target("${name}_scalar")
   endif()

   list(FIND _disabled_impl "SSE" _index)
   if(USE_SSE2 AND _index EQUAL -1)
      add_executable("${name}_sse" ${_SRCS})
      add_target_property("${name}_sse" COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS} -DVC_IMPL=SSE")
      add_target_property("${name}_sse" LABELS "SSE")
      add_dependencies(SSE "${name}_sse")
      target_link_libraries("${name}_sse" Vc ${_LIBS})
      vc_add_run_target("${name}_sse")
   endif()

   list(FIND _disabled_impl "AVX" _index)
   if(USE_AVX AND _index EQUAL -1)
      add_executable("${name}_avx" ${_SRCS})
      add_target_property("${name}_avx" COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS} -DVC_IMPL=AVX")
      add_target_property("${name}_avx" LABELS "AVX")
      add_dependencies(AVX "${name}_avx")
      target_link_libraries("${name}_avx" Vc ${_LIBS})
      vc_add_run_target("${name}_avx")
   endif()

   list(FIND _disabled_impl "AVX2" _index)
   if(USE_AVX2 AND _index EQUAL -1)
      add_executable("${name}_avx2" ${_SRCS})
      add_target_property("${name}_avx2" COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS} -DVC_IMPL=AVX2")
      add_target_property("${name}_avx2" LABELS "AVX2")
      add_dependencies(AVX2 "${name}_avx2")
      target_link_libraries("${name}_avx2" Vc ${_LIBS})
      vc_add_run_target("${name}_avx2")
   endif()

   list(FIND _disabled_impl "MIC" _index)
   if(MIC_NATIVE_FOUND AND "${_LIBS}" STREQUAL "" AND _index EQUAL -1)
      mic_add_executable("${name}_mic"
         #COMPILE_FLAGS -no-opt-prefetch
         SOURCES ${_SRCS}
         LINK_LIBRARIES Vc_MIC
         )
      add_target_property("${name}_mic" LABELS "MIC")
      add_dependencies(MIC "${name}_mic")
      vc_add_run_target("${name}_mic")
   endif()

   if(USE_AVX2)
      add_target_property("${name}_default" LABELS "AVX2")
      add_dependencies(AVX2 "${name}_default")
   elseif(USE_AVX)
      add_target_property("${name}_default" LABELS "AVX")
      add_dependencies(AVX "${name}_default")
   elseif(USE_SSE2)
      add_target_property("${name}_default" LABELS "SSE")
      add_dependencies(SSE "${name}_default")
   else()
      add_target_property("${name}_default" LABELS "Scalar")
      add_dependencies(Scalar "${name}_default")
   endif()
   vc_add_run_target("${name}_default")
endmacro(build_example)

build_example(kdtree main.cpp)

add_executable(simdize test.cpp)
add_target_property(simdize COMPILE_FLAGS "${Vc_ARCHITECTURE_FLAGS}")
add_target_property(simdize LABELS "other")
target_link_libraries(simdize Vc)
add_dependencies(other simdize)
