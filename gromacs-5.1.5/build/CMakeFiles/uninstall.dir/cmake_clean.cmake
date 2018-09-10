FILE(REMOVE_RECURSE
  "install_manifest.txt"
  "CMakeFiles/uninstall"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/uninstall.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
