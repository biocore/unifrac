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

int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
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

int pthread_setaffinity_np(pthread_t thread, size_t cpu_size,
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

int bind_to_core(int tid) {
	// https://stackoverflow.com/a/11583550/19741
	// http://blog.saliya.org/2015/07/get-and-set-process-affinity-in-c.html
    pthread_t thread = pthread_self();
    pid_t pid = getpid();

    cpu_set_t current_set, new_set;
    int ret;

    CPU_ZERO(&current_set);
    CPU_ZERO(&new_set);

    ret = sched_getaffinity(pid, sizeof(current_set), &current_set);

    int target = -1;
	for(unsigned int j = 0; j < CPU_COUNT(&current_set); j++) {
        if(CPU_ISSET(j, &current_set)) {
            target++;
        }
        if(target == tid)
            break;
    }

    if(target != tid) {
        return -1;
    }

    CPU_SET(target, &new_set);
    return pthread_setaffinity_np(thread, sizeof(new_set), &new_set);
}
