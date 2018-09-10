FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/NightlyMemoryCheck"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/NightlyMemoryCheck.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
