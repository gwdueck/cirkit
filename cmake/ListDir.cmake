
function(subdirlist result curdir)
  file(GLOB children ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${child})
      get_filename_component(base_name ${child} NAME_WE)
      set(dirlist ${dirlist} ${base_name})
    endif()
  endforeach()
  set(${result} ${dirlist} PARENT_SCOPE)
endfunction()

function(srcdirlist result curdir)
  if("${cirkit_PACKAGES}" STREQUAL "")
    subdirlist(tmp ${curdir})
  else()
    set(tmp ${cirkit_PACKAGES})
  endif()
  set(${result} ${tmp} PARENT_SCOPE)
endfunction()
