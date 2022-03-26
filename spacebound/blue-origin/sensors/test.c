
#include "../include/program.h"

local bool read_test(i2c_device * test_i2c);

Sensor * init_test(Sensor * test) {

  test -> name = "TEST_SENSOR";
  test -> i2c = create_i2c_device(test, read_test, "A virtual test sensor");

  sensor_log_header(test, GRAY);
  return test;
}

bool read_test(i2c_device * test_i2c) {

  Sensor * test = test_i2c -> sensor;

  static double rand_ary[TEST_SMOOTHING_SAMPLES] = {0};

  //for (int i = 9; i > 0; i--) {
  //	rand_ary[i] = rand_ary[i-1];
  //}

  memmove(&rand_ary[1], &rand_ary[0], sizeof(double)*9);

  rand_ary[0] = rand();
  // memmove(&rand_ary[1], &rand_ary[0], sizeof(int)*9);

  double total = 0;
  double filled = 0;

  for (int i = 0; i < TEST_SMOOTHING_SAMPLES; i++) {
	if (rand_ary[i] > 0) {
		filled += (TEST_SMOOTHING_SAMPLES - i);
		total += ((TEST_SMOOTHING_SAMPLES - i) * rand_ary[i]);
	}
	  // total += rand_ary[i];
  }

  float identity = test -> outputs[TEST_MEASURE_IDENTITY].measure;

  bind_stream(test,                       0, TEST_MEASURE_ZERO    );
  bind_stream(test,                identity, TEST_MEASURE_IDENTITY);
  bind_stream(test, sin(time_passed() / 10), TEST_MEASURE_SINE    );
  bind_stream(test, cos(time_passed() / 10), TEST_MEASURE_COSINE  );
  bind_stream(test,             rand_ary[0], TEST_MEASURE_RANDOM  );
  bind_stream(test,                total/filled, TEST_MEASURE_RANDSMOOTH );


  //bind_stream(test, !( ((int)time_passed()) % 4), TEST_MEASURE_PULSE  );

  sensor_log_outputs(test, test_i2c -> log, NULL);
  sensor_process_triggers(test);
  return true;
}
