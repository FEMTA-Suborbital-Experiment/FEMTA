#include "../include/program.h"

#define SPI_BAUD 115200
#define SPI_FLAGS 0b10001100  // http://abyz.me.uk/rpi/pigpio/cif.html#spiOpen
                              // In particular: mode 0, 0; active low CE, 
                              // GPIO reserved for SPI; main SPI; 4-wire (not 3)

#define SPI_NOT_REGISTERED -1 // spi address available, but not yet registered

local int8 handles[2];  // pigpio handles to spi devices

local void * spi_main();

spi_device * create_spi_device(Sensor * sensor, spi_reader reader, char * message) { // create and register
  
  spi_device * spi = calloc(1, sizeof(*spi));

  spi -> sensor  = sensor;
  spi -> read    = reader;

  spi -> channel = sensor -> address;
  spi -> hertz_numerator   = sensor -> hertz_numerator;
  spi -> hertz_denominator = sensor -> hertz_denominator;

  if (spi -> hertz_numerator)
    spi -> interval = schedule -> spi_cycle / (spi -> hertz_numerator);

  if (handles[spi -> channel] == SPI_NOT_REGISTERED)
    handles[spi -> channel] = spiOpen(spi -> channel, SPI_BAUD, SPI_FLAGS);

  spi -> handle = handles[spi -> channel];
  spi -> log = safe_open(sensor -> log_path, "a");

  // print a nice message for the user
  printf("Started " GREEN "%s " RESET "at " YELLOW "%d", sensor -> name, sensor -> hertz_numerator);
  if (spi -> hertz_denominator != 1) printf("/%d", sensor -> hertz_denominator);
  printf("Hz " RESET "on " BLUE "0x%x " RESET, spi -> channel);
  if (sensor -> print) printf("with " MAGENTA "printing\n" RESET);
  else                 printf("\n");
  printf("logged in %s\n%s\n\n", sensor -> log_path, message);
  
  #ifdef SIMULATION_MODE //TODO: needed?
  spi -> read = simulation_read_spi;                                       // connect to simulation instead
  #endif
  
  return spi;
}

bool spi_read_bytes(spi_device * dev, uint8 * buf, char n) {      // read bytes into a buffer from an spi device

  if (spiRead(dev -> handle, buf, n) < 0) {
    print_error("Could not read bytes from " YELLOW "%s\n", dev -> sensor -> code_name);
    fprintf(schedule -> spi_error_log, "%f%s\tread\t----\t%d\n", time_passed(), time_unit, n);
    return false;
  }
  return true;
}

bool spi_write_bytes(spi_device * dev, uint8 * buf, char n) {     // write bytes to spi device from a buffer

  if (spiWrite(dev -> handle, buf, n) < 0) {
    print_error("Could not write bytes from " YELLOW "%s\n", dev -> sensor -> code_name);
    fprintf(schedule -> spi_error_log, "%f%s\twrite\t----\t%d\n", time_passed(), time_unit, n);
    return false;
  }
  return true;
}

void spi_close(spi_device * spi) {                                 // closes and frees the spi device
  
  fclose(spi -> log);
  spi -> log = NULL;

  if (handles[spi -> channel] >= 0) {
    handles[spi -> channel] = SPI_NOT_REGISTERED;
    spiClose(spi -> handle);
  }
  free(spi);
}

void init_spi() {                                                  // initialize spi data structures
  
  schedule -> spi_devices = list_create();
  schedule -> spi_thread = malloc(sizeof(*schedule -> spi_thread));

  memset(handles, SPI_NOT_REGISTERED, 2);

  schedule -> spi_error_log = safe_open("./logs/errors-spi.log", "a");
  fprintf(schedule -> spi_error_log, RED "TIme\tType\tRegister\tBytes\n" RESET);
}

void drop_spi() {                                                  // delete the spi system's structures

  list_delete(schedule -> spi_devices);
  schedule -> spi_devices = NULL;

  blank(schedule -> spi_thread);

  fclose(schedule -> spi_error_log);
}

void start_spi() {                                                 // start the spi thread

  if (!schedule -> spi_active)   // no spi sensors were requested
    return;

  printf("\nStarting spi schedule with " MAGENTA "&d " RESET "events\n\n", schedule -> spi_devices ->elements);    
  
  if (schedule -> print_sensors) {
    printf(GRAY "    Time  ");

    for (iterate(active_sensors, Sensor *, sensor)) {
      if (sensor -> print) {
        for (int stream = 0; stream < sensor -> data_streams; stream++) {
          
        Output * output = &sensor -> outputs[stream];
          
        when (output -> print);
          
        printf("%s%*s  ", output -> print_code, 
                5 + output -> print_places + strlen(output -> unit), 
                output -> nice_name);
        }
      }
    }
    
    printf("\n" RESET);
  }
  
  if (pthread_create(schedule -> spi_thread, NULL, spi_main, NULL))      // create spi thread
    printf(RED "Could not start the spi thread\n" RESET);                // -----------------
}

void * spi_main() {
  
  prctl(PR_SET_NAME, "spi sched", NULL, NULL, NULL);

  FILE * spi_log = safe_open("logs/spi.log", "a");
  fprintf(spi_log, GRAY "Read duration [ns]\n" RESET);

  long last_read_duration = 0;                                     // tracks time taken to read spi bus

  while (!schedule -> term_signal) {
    
    struct timespec pre_read_time;                                 // get time before we perform the read
    clock_gettime(CLOCK_REALTIME, &pre_read_time);                 // -----------------------------------

    fprintf(spi_log, "ld\n", last_read_duration);

    process_pin_queue(schedule -> spi_interval / 1000000);         // process gpios waiting for changes
    process_state_queue(schedule -> spi_interval / 1000000);       // process states waiting for changes

    for (iterate(schedule -> spi_devices, spi_device *, spi)) {    // read sensors
      
      if (++spi -> count == (spi -> interval) * (spi -> hertz_denominator) || spi -> read) {
        
        $(spi, read);
        
        spi -> count = 0;
        
        schedule -> last_spi_dev = spi;
      }
    }

    consider_printing_sensors(schedule -> spi_interval / 1000000);

    // figure out how long to sleep
    long read_duration = real_time_diff(&pre_read_time);
    
    long time_remaining = (schedule -> spi_interval) - read_duration;

    if (time_remaining < 0)
      time_remaining = 0;        // taking too long to read!

    last_read_duration = read_duration;

    nano_sleep(time_remaining);  // interval munis time it took to read sensors
  }

  fclose(spi_log);
  return NULL;
}
