FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ContinuousMemCheck"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ContinuousMemCheck.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
