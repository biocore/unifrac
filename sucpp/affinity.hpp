#define _GNU_SOURCE
#include <pthread.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#ifdef __LINUX__
#include <sched.h>
#endif


#ifdef __APPLE__
#include <mach/thread_policy.h>
#include <mach/task_info.h>
#include <sys/sysctl.h>
#include <mach/thread_policy.h>
#include <mach/thread_act.h>

// OSX code adapted from
// http://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
// these macros and methods don't exist on OSX

#define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"

typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline void
CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

static inline int
CPU_COUNT(cpu_set_t *cs) { return __builtin_popcount(cs->count); }

#define CPU_SETSIZE 32

static int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
{
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

static int pthread_setaffinity_np(pthread_t thread, size_t cpu_size,
                           cpu_set_t *cpu_set)
{
  thread_port_t mach_thread;
  int core = 0;

  for (core = 0; core < 8 * cpu_size; core++) {
    if (CPU_ISSET(core, cpu_set)) break;
  }
  thread_affinity_policy_data_t policy = { core };
  mach_thread = pthread_mach_thread_np(thread);
  thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY,
                    (thread_policy_t)&policy, 1);
  return 0;
}

#endif

#define handle_error_en(en, msg) \
       do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

static int bind_to_core(int core) {
    /* bind the calling thread to the requested core
     *
     * The use of this method is for better NUMA utilization. The 
     * default NUMA policy is local, where memory is allocated on the NUMA node
     * relative to the core if possible. The intention with this method is to
     * bind to a core first, and then allocate memory. A beneficial side effect
     * is that threads should not hop between cores either. 
     *
     * This method is cgroup safe.
     */
    // https://stackoverflow.com/a/11583550/19741
    // http://blog.saliya.org/2015/07/get-and-set-process-affinity-in-c.html
    pthread_t thread = pthread_self();
    pid_t pid = getpid();

    cpu_set_t current_set, new_set;
    int j, ret;

    CPU_ZERO(&current_set);
    CPU_ZERO(&new_set);

    ret = sched_getaffinity(pid, sizeof(current_set), &current_set);

    // find which core in our cpu_set corresponds to the callers
    // request
    int target = -1;
    for(j = 0; j < CPU_SETSIZE; j++) {
        if(CPU_ISSET(j, &current_set)) {
            target++;
        }
        if(target == core) 
            break;
    }

    if(target != core) {
        fprintf(stderr, "Unable to bind this thread to core %d. Are sufficient processors available?", thread);
        return -1;
    }

    CPU_SET(j, &new_set);
    int serr = pthread_setaffinity_np(thread, sizeof(new_set), &new_set);
    return serr;
}
