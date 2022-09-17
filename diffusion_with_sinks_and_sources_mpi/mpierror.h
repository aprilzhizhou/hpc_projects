#include <stdio.h>
#include <stdlib.h>

#ifdef DO_ERROR_CHECKING

static void CHECK_ERROR( int err, const char* file, int line) {
  if (err != MPI_SUCCESS) {
    int len;
    char errstring[MPI_MAX_ERROR_STRING];
    MPI_Error_string(err, errstring, &len);
    printf("%s in %s at line %d\n", errstring, file, line);
    fflush(stdout);
    MPI_Finalize();
    exit(1);
  }
}
#define CheckError( err ) (CHECK_ERROR(err, __FILE__, __LINE__))

#else

#define CheckError( err ) err

#endif
