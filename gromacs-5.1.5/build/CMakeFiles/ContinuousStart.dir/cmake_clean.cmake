FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ContinuousStart"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ContinuousStart.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
