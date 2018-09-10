FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ContinuousConfigure"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ContinuousConfigure.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
