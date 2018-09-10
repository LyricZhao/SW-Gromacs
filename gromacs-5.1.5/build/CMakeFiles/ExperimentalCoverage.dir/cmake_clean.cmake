FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ExperimentalCoverage"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ExperimentalCoverage.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
