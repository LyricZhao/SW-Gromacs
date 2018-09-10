FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/run-ctest"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/run-ctest.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
