FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/ExperimentalBuild"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/ExperimentalBuild.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
