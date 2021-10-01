import RPi.GPIO as GPIO
import time


GPIO.setmode(GPIO.BCM)

left_sol_open = 12
left_sol_close = 16

right_sol_open = 20
right_sol_close = 21

GPIO.setup(left_sol_open, GPIO.OUT)
GPIO.setup(left_sol_close, GPIO.OUT)

GPIO.setup(right_sol_open, GPIO.OUT)
GPIO.setup(right_sol_close, GPIO.OUT)



# Open Left Sol
GPIO.output(left_sol_open, GPIO.HIGH)
GPIO.output(left_sol_close, GPIO.LOW)
time.sleep(350)
GPIO.output(left_sol_open, GPIO.LOW)


# Close Left  Sol
GPIO.output(left_sol_close, GPIO.HIGH)
GPIO.output(left_sol_open, GPIO.LOW)
time.sleep(350)
GPIO.output(left_sol_close, GPIO.LOW)


# Open Right Sol
GPIO.output(right_sol_open, GPIO.HIGH)
GPIO.output(right_sol_close, GPIO.LOW)
time.sleep(350)
GPIO.output(right_sol_open, GPIO.LOW)


# Close Right Sol
GPIO.output(right_sol_close, GPIO.HIGH)
GPIO.output(right_sol_open, GPIO.LOW)
time.sleep(350)
GPIO.output(right_sol_close, GPIO.LOW)

GPIO.cleanup()
