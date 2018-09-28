# include "gmxpre.h"
# include "gromacs/commandline/cmdlinemodulemanager.h"
# include "mdrun/mdrun_main.h"

# define USE_VIC_OPT

# ifdef USE_VIC_OPT
# include "athread.h"
# endif

int main(int argc, char *argv[]) {
# ifdef USE_VIC_OPT
    athread_init();
    athread_enter64();
# endif

    int ret = gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);

# ifdef USE_VIC_OPT
    athread_leave64();
    athread_halt();
# endif

    return ret;
}
