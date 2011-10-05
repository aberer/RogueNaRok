#ifdef PARALLEL
volatile int numberOfThreads = 0; 
volatile int numberOfJobs;
volatile int jobCycle;
volatile int threadJob;
volatile char *barrierBuffer;
pthread_mutex_t mutex;
#endif

extern char *infoFileName,
  run_id[128],
  *workdir, 
  *programName,
  *programVersion, 
  *programReleaseDate;

