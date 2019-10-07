
#define _GNU_SOURCE

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <pigpio.h>
#include <sys/prctl.h>
#include <time.h>

#include "clock.h"
#include "color.h"
#include "i2c.h"

#include "../sensors/sensor.h"
#include "../structures/list.h"
#include "../system/gpio.h"
#include "../types/types.h"
#include "../types/thread-types.h"

int8 handles[0x7F];

void * i2c_main();

i2c_device * create_i2c_device(Sensor * sensor, ProtoSensor * proto, i2c_reader reader) {
  // creates an i2c device, adding it to the device list
  
  i2c_device * i2c = calloc(1, sizeof(*i2c));
  
  i2c -> sensor   = sensor;
  i2c -> read     = reader;
  
  i2c -> address           = proto -> address;
  i2c -> hertz             = proto -> hertz;
  i2c -> hertz_denominator = proto -> hertz_denominator;
  
  if (i2c -> hertz)
    i2c -> interval = 1000 / (i2c -> hertz);
  
  if (!i2c -> hertz_denominator)
    i2c -> hertz_denominator = 1;
  
  
  // see if we need to open another i2c instance
  if (handles[i2c -> address] == -1) {
    i2c -> handle = i2cOpen(1, i2c -> address, 0);
    handles[i2c -> address] = i2c -> handle;
  }
  else {
    i2c -> handle = handles[i2c -> address];
  }
  
  return i2c;
}

uint8 i2c_read_byte(i2c_device * dev, uint8 reg) {
  // reads a single byte from an i2c device
  
  return i2cReadByteData(dev -> handle, reg);
}

bool i2c_read_bytes(i2c_device * dev, uint8 reg, uint8 * buf, char n) {
  // reads up to 32 bytes from an i2c device
  
  if (i2cReadI2CBlockData(dev -> handle, reg, buf, n) < 0) {
    printf(RED "Could not read bytes from " YELLOW "%s\n" RESET, dev -> sensor -> name);
    return false;
  }
  return true;
}

bool i2c_raw_read(i2c_device * dev, uint8 * buf, char n) {
  // reads up to 32 bytes from an i2c device, without asking for a particular register
  
  if (i2cReadDevice(dev -> handle, buf, n)) {
    printf(RED "Could not read raw bytes from " YELLOW "%s\n" RESET, dev -> sensor -> name);
    return false;
  }
  return true;
}

bool i2c_raw_write(i2c_device * dev, uint8 * buf, char n) {
  // writes up to 32 bytes from an i2c device, without specifying a particular register
  
  if (i2cWriteDevice(dev -> handle, buf, n)) {
    printf(RED "Could not write raw bytes from " YELLOW "%s\n" RESET, dev -> sensor -> name);
    return false;
  }
  return true;
}

bool i2c_write_byte(i2c_device * dev, uint8 reg, uint8 value) {
  // writes a byte to the i2c device
  
  if (i2cWriteByteData(dev -> handle, reg, value) < 0) {
    printf(RED "Could not write byte to " YELLOW "%s\n" RESET, dev -> sensor -> name);
    return false;
  }
  return true;
}

bool i2c_write_bytes(i2c_device * dev, uint8 reg, uint8 * buf, char n) {
  // writes up to 32 bytes to the i2c device
  
  if (i2cWriteI2CBlockData(dev -> handle, reg, buf, n)) {
    printf(RED "Could not write bytes to " YELLOW "%s\n" RESET, dev -> sensor -> name);
    return false;
  }
  return true;
}


void i2c_close(i2c_device * i2c) {
  // closes and frees the i2c device
  
  fclose(i2c -> log);
  i2cClose(i2c -> handle);
  free(i2c);
}

void init_i2c() {
  // initialize i2c data structures
  
  schedule -> i2c_devices = list_create();
  schedule -> i2c_thread  = malloc(sizeof(*schedule -> i2c_thread));
  
  // prepare handle array
  for (int i = 0; i < 0x7F; i++)
    handles[i] = -1;
}

void start_i2c() {
  
  if (!schedule -> i2c_active) return;
  
  printf("\nStarting i2c schedule with " MAGENTA "%d " RESET "events\n", schedule -> i2c_devices -> size);
  
  // create i2c thread
  if (pthread_create(schedule -> i2c_thread, NULL, i2c_main, NULL)) {
    printf(RED "Could not start i2c thread\n" RESET);
    return;
  }
}

void terminate_i2c() {
  // frees everything associated with the i2c system
  
  list_destroy(schedule -> i2c_devices);      // note that this kills
  
  free(schedule -> i2c_thread);
  
  schedule -> i2c_devices = NULL;
  schedule -> i2c_thread  = NULL;
}

void * i2c_main() {
  
  prctl(PR_SET_NAME, "i2c sched", NULL, NULL, NULL);
  
  FILE * i2c_log = fopen("logs/i2c.log", "a");
  fprintf(i2c_log, GRAY "Read duration [ns]\n" RESET);
  
  long bus_interval    = schedule -> i2c_interval;
  long bus_interval_ms = schedule -> i2c_interval / (long) 1E6;
  
  long last_read_duration = 0;    // tracks time taken to read i2c bus
  
  while (!schedule -> term_signal) {
    
    // get time before we perform the read
    struct timespec pre_read_time;
    clock_gettime(CLOCK_REALTIME, &pre_read_time);
    
    fprintf(i2c_log, "%ld\n", last_read_duration);
    
    
    // pulse pins
    for (iterate(schedule -> pulse_pins, Pin *, pin)) {
      
      //printf("PULES COUNTDOWN: %d " YELLOW "%d" RESET "\n", pin -> broadcom, pin -> ms_until_pulse_completes);
      
      if (!pin -> ms_until_pulse_completes) continue;
      
      pin -> ms_until_pulse_completes -= bus_interval_ms;
      
      if (pin -> ms_until_pulse_completes <= 0) {
	pin_set(pin -> broadcom, pin -> pulse_final_state);
	pin -> ms_until_pulse_completes = 0;
      }
    }
    
    
    // read sensors
    for (iterate(schedule -> i2c_devices, i2c_device *, i2c)) {
      
      i2c -> count += bus_interval_ms;
      
      if (i2c -> count == (i2c -> interval) * (i2c -> hertz_denominator) || i2c -> reading) {
	
	(i2c -> read)(i2c);
        
	i2c -> count = 0;
      }
    }    
    
    // figure out how long to sleep
    long read_duration = real_time_diff(&pre_read_time);
    
    long time_remaining = bus_interval - read_duration;
    
    if (time_remaining < 0)
      time_remaining = 0;               // taking too long to read!
    
    last_read_duration = read_duration;
    
    real_nano_sleep(time_remaining);   // interval minus time it took to read sensors
  }
  
  fclose(i2c_log);  
  return NULL;
}
