FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/NightlySubmit"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/NightlySubmit.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
