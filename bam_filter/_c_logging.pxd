cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_log(const char *tag, const char *msg)
    void bf_nogil_log_level(int level, const char *tag, const char *msg)
    void bf_set_verbosity(int v)
    int bf_get_verbosity()
    int bf_should_log(int level)
