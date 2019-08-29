

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <pigpio.h>

#include "../sensors/sensor.h"
#include "../system/clock.h"
#include "../system/color.h"
#include "../system/gpio.h"
#include "../system/i2c.h"
#include "../structures/list.h"
#include "../structures/selector.h"
#include "../types/types.h"
#include "../parser/y.tab.h"

FILE * yyin;
void print_config();

void parse_args(int argc, char ** argv) {

  if (argc == 1) {

    yyin = fopen("/home/noah/FEMTA/spacebound/blue-origin/experiments/default.e", "r");
    
    if (!yyin) {
      printf(RED "Experiment file does not exist\n" RESET);
      exit(1);
    }
      
    yyparse();
    
    return;

    
    // default, use all sensors
    
    for (iterate(proto_sensors -> all, ProtoSensor *, proto))
      if (proto -> hertz)
	proto -> requested = true;
    
    return;
  }
  
  if (!strncmp(argv[1], "file", 4)) {

    char * filename = argv[1] + 5;
    
    printf("%s\n", filename);
    
    yyin = fopen(filename, "r");
    
    if (!yyin) {
      printf(RED "Experiment file %s does not exist\n" RESET, filename);
      exit(1);
    }
      
    yyparse();
    
    return;
  }
  
  for (int i = 1; i < argc; i++) {
    
    char code_name[32];
    code_name[0] = '\0';    // need to protect buffer from previous iteration
    
    int  hertz = 0;
    bool print = false;
    
    sscanf(argv[i], "%[^*,],%d", code_name, &hertz);      
    
    if (!code_name[0]) {
      sscanf(argv[i], "*%[^,],%d", code_name, &hertz);
      print = true;
    }
    
    ProtoSensor * proto = hashmap_get(proto_sensors, code_name);
    
    //printf("%s at %d\n", code_name, hertz);
    
    if (!proto) {
      printf(RED "%s is not a sensor\n" RESET, code_name);
      exit(1);
    }
    
    proto -> requested = true;
    proto -> print     = print;
    
    if (hertz) proto -> hertz = hertz;
  }
}



int main(int argc, char ** argv) {
  
  // start pigpio library
  if (gpioInitialise() < 0) {
    printf(RED "pigpio unable to start\n" RESET);
    exit(2);
  }
  
  init_pins();       // set up gpio data structure
  
  schedule = calloc(1, sizeof(*schedule));
  
  init_one();        // set up the 1-wire data structures
  init_i2c();        // set up the i2c data structures
  init_sensors();    // set up sensor info and actions
  
  parse_args(argc, argv);
  
  print_config();
  
  start_sensors();
  start_one();       // start reading the 1-wire bus
  start_i2c();       // start reading the i2c bus
  
  Selector * selector = create_selector(NULL);
  
  add_selector_command(selector, 'q', "quit",  flip_bool, &reading_user_input);
  add_selector_command(selector, 'c', "char", output_str, (void *) "quit");
  
  reading_user_input = true;
  
  char input;
  while (reading_user_input) {
    input = getc(stdin);
    
    execute_selector(selector, input);
  }
  
  // tell threads to terminate
  schedule -> term_signal = true;
  
  // join with threads
  pthread_join(*schedule -> thread, NULL);
  
  terminate_i2c();
  terminate_one();
  
  free(schedule);
  
  // terminate pigpio library
  gpioTerminate();
  
  return EXIT_SUCCESS;
}