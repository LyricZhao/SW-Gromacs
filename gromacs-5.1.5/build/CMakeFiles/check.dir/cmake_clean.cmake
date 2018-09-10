FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/check"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/check.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
