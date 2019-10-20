#ifndef HEADER_GAURD_SENSOR
#define HEADER_GAURD_SENSOR

/**
 * @Invariant_0 A sensor shall always be called "proto" until it is initialized.
 * Hence, all sensor initialization functions are "pass-through" since one may
 * never refer to an initialized sensor through a pointer called "proto".
 * Moreover, one never refers to a "ds32" until it in its "proto" form has been
 * initialized.
 * 
 * @Invariant_1 If a sensor produces data, structures must exist for triggers
 * 
 * @Invariant_2 After parsing completes, all triggers are in the same units that
 * output streams and log files use.
 */



#include <stdbool.h>
#include <stdio.h>

#include "../math/units.h"
#include "../structures/list.h"
#include "../structures/hashmap.h"
#include "../system/i2c.h"
#include "../system/one.h"
#include "../system/state.h"
#include "../types/types.h"

#define I2C_BUS 0x1
#define ONE_BUS 0x2

typedef struct Sensor Sensor;
typedef struct i2c_device i2c_device;
typedef struct one_device one_device;
typedef struct Numeric Numeric;

typedef void (* sensor_free)(Sensor * sensor);

typedef struct Trigger {
  
  char * id;                           // target (like A01)
  Numeric * threshold;                 // threshold for condition
  List * precondition;                 // states the system must be in to consider this trigger
  
  List * wires_low;                    // wires to set low
  List * wires_high;                   // wires to set high
  List * enter_set;                    // sub-states to enter
  List * leave_set;                    // sub-states to leave
  
  bool less;                           // whether comparison is '<' or '>'
  bool fired;                          // whether trigger has ever fired
  bool singular;                       // whether trigger should only ever be fired once
  bool reverses;                       // whether trigger was source for its chiral opposite
  
  // parser only
  Numeric * threshold_as_specified;    // threshold as specified in E
  
} Trigger;

typedef struct Output {
  
  float  measure;     // measurement last taken
  List * series;      // path of conversions required to get final measurement
  List * triggers;    // sensor triggers for this particular output stream
  bool   enabled;     // whether output is enabled
  
} Output;

typedef struct Sensor {
  
  char * name;              // component name
  char * code_name;         // abbreviated name
  bool   print;             // whether sensor prints
  char * print_code;        // color to use while printing
  int    print_hertz;       // frequency for printing
  
  union {
    i2c_device * i2c;       // i2c communications info
    one_device * one;       // 1-wire communications info
  };
  int bus;                  // which bus is used
  
  int data_streams;         // number of output data axes
  Output * outputs;         // everything this sensor produces
  Hashmap * targets;        // that which can be triggered (target str -> stream_index)
  
  float auto_regressive;    // smoothing constant
  
  sensor_free free;         // how to free sensor
  
  
  
  int hertz;                 // bus communication frequency in hertz
  int hertz_denominator;     // engenders fractional frequency through deferrals
  
  bool requested;            // whether sensor has actually been specified during parsing
  
  uint8 address;             // i2c address
  
} Sensor;

typedef struct Schedule {
  
  List * i2c_devices;          // list of all i2c device pointers
  List * one_devices;          // list of all 1-wire device pointers
  List * pulse_pins;           // list of all pins that may pulse
  
  long i2c_interval;           // scheduler spacing for the i2c protocol
  long one_interval;           // scheduler spacing for the 1-wire protocol
  
  bool i2c_active;             // whether experiment uses i2c
  bool one_active;             // whether experiment uses 1-wire
  
  int   interrupts;            // interrupts since schedule creation
  float interrupt_interval;    // time between each interrupt
  
  Thread * i2c_thread;         // i2c thread
  Thread * one_thread;         // 1-wire thread (deprioritized)
  bool term_signal;            // when set to true, schedule terminates
  
} Schedule;

Schedule * schedule;


List    * active_sensors;      // every active sensor on craft
Hashmap * all_sensors;         // every sensor that could be specified

int n_triggers;                // number of triggers

float time_passed();    // time since experiment start

void init_sensors();
void start_sensors();
void terminate_sensors();
Sensor * sensor_from_proto(ProtoSensor *);

void flip_print(void * nil, char * raw_text);

#endif
