# include "gmxpre.h"
# include "gromacs/commandline/cmdlinemodulemanager.h"
# include "mdrun/mdrun_main.h"

# include "athread.h"

extern "C" {
  void athread_INIT();
  void athread_LEAVE();
}

int main(int argc, char *argv[]) {

    athread_INIT();
    int ret = gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);
    athread_LEAVE();

    return ret;
}
