FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ContinuousSubmit"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ContinuousSubmit.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
