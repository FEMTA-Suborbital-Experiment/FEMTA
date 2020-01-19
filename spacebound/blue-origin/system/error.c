
#include "../include/program.h"

void exit_printing(int code, ...) {
  // prints a colored message to the console,
  // and then exits safely by closing pigpio.
  
  va_list args;
  va_start(args, code);
  
  char * format = va_arg(args, char *);
  
  if (code == SUCCESS) printf(GREEN);
  else                 printf(RED);
  vprintf(format, args);
  printf(RESET "\n");
  
  va_end(args);
  
  gpioTerminate();
  exit(code);
}
