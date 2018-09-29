typedef struct tMPI_Spinlock {
    volatile unsigned int lock /*__attribute__ ((aligned(64)))*/;
} tMPI_Spinlock_t;

#define TMPI_SPINLOCK_INITIALIZER   { 0 }

#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK



static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x) {
    x->lock = 0;
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x) {
#if 1
    while (__sync_lock_test_and_set(&(x->lock), 1) == 1) {
        /* this is nicer on the system bus: */
        while (x->lock == 1) { }
    }
#else
    do
    {
    }
    while (__sync_lock_test_and_set(&(x->lock), 1) == 1);
#endif
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x) {
    return __sync_lock_test_and_set(&(x->lock), 1);
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x) {
    //__sync_lock_release(&(x->lock));
    x->lock = 0;
}

static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *x) {
    __sync_synchronize();
    return ( x->lock == 1 );
}

static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x) {
    do { }
    while (x->lock == 1);
    __sync_synchronize();
}
