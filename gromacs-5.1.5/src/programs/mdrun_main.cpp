# include "gmxpre.h"
# include "gromacs/commandline/cmdlinemodulemanager.h"
# include "mdrun/mdrun_main.h"

# ifdef USE_VIC_OPT
# include "athread.h"
# endif

int main(int argc, char *argv[]) {
# ifdef USE_VIC_OPT
    athread_init();
# endif

    int ret = gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);

# ifdef USE_VIC_OPT
    athread_halt();
# endif
    
    return ret;
}
