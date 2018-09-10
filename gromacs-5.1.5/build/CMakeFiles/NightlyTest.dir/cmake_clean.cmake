FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/NightlyTest"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/NightlyTest.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
