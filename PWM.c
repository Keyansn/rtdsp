#include "mbed.h"
#include "rtos.h"
#include "PID.h"
#include "stdint.h"
#include <limits>
#include <algorithm>
#include <RawSerial.h>
#include <cstdlib>

//Photointerrupter input pins
#define I1pin D2
#define I2pin D11
#define I3pin D12
 
//Incremental encoder input pins
#define CHA   D7
#define CHB   D8  
 
//Motor Drive output pins   //Mask in output byte
#define L1Lpin D4           //0x01
#define L1Hpin D5           //0x02
#define L2Lpin D3           //0x04
#define L2Hpin D6           //0x08
#define L3Lpin D9           //0x10
#define L3Hpin D10          //0x20
 
//Mapping from sequential drive states to motor phase outputs
/*
State   L1  L2  L3
0       H   -   L
1       -   H   L
2       L   H   -
3       L   -   H
4       -   L   H
5       H   L   -
6       -   -   -
7       -   -   -
*/

/*
typedef enum 
{

    osPriorityLow           = -2,
    osPriorityBelowNormal   = -1,
    osPriorityNormal        =  0,
    osPriorityAboveNormal   =  1,
    osPriorityHigh          =  2,
    osPriorityError        =  0x84
} osPriority;
*/





double dbl_max = std::numeric_limits<double>::max();




//Drive state to output table
const int8_t driveTable[] = {0x12,0x18,0x09,0x21,0x24,0x06,0x00,0x00};
 
//Mapping from interrupter inputs to sequential rotor states. 0x00 and 0x07 are not valid
const int8_t stateMap[] = {0x07,0x05,0x03,0x04,0x01,0x00,0x02,0x07};  
//const int8_t stateMap[] = {0x07,0x01,0x03,0x02,0x05,0x00,0x04,0x07}; //Alternative if phase order of input or drive is reversed
 
//Phase lead to make motor spin
const int8_t lead = -2;  //2 for forwards, -2 for backwards
 
//Status LED
DigitalOut led1(LED1);
 
//Photointerrupter inputs
DigitalIn I1(I1pin);
DigitalIn I2(I2pin);
DigitalIn I3(I3pin);
 
//Motor Drive outputs
PwmOut L1L(L1Lpin);
DigitalOut L1H(L1Hpin);
PwmOut L2L(L2Lpin);
DigitalOut L2H(L2Hpin);
PwmOut L3L(L3Lpin);
DigitalOut L3H(L3Hpin);

double delta = 1;
double OGdelta =0;
 
//Set a given drive state
void motorOut(int8_t driveState){
    
    //Lookup the output byte from the drive state.
    int8_t driveOut = driveTable[driveState & 0x07];
      
    //Turn off first
    if (~driveOut & 0x01) L1L = 0;
    if (~driveOut & 0x02) L1H = 1;
    if (~driveOut & 0x04) L2L = 0;
    if (~driveOut & 0x08) L2H = 1;
    if (~driveOut & 0x10) L3L = 0;
    if (~driveOut & 0x20) L3H = 1;
    
    //Then turn on
    if (driveOut & 0x01) L1L = delta;
    if (driveOut & 0x02) L1H = 0;
    if (driveOut & 0x04) L2L = delta;
    if (driveOut & 0x08) L2H = 0;
    if (driveOut & 0x10) L3L = delta;
    if (driveOut & 0x20) L3H = 0;
    }
    
    //Convert photointerrupter inputs to a rotor state
inline int8_t readRotorState(){
    return stateMap[I1 + 2*I2 + 4*I3];
    }
 
//Basic synchronisation routine    
int8_t motorHome() {
    //Put the motor in drive state 0 and wait for it to stabilise
    motorOut(0);
    wait(3.0);
    
    //Get the rotor state
    return readRotorState();
}


//Initialise the serial port
//Serial pc(SERIAL_TX, SERIAL_RX);
RawSerial pc(USBTX,USBRX);

    
/*************************************
             Interrupts
*************************************/ 
    
/* Read Photointerrupters */
void check_photo();

Ticker photo_checker_interrupt;
volatile int8_t intState;
    

/* Count Revolutions */
void increment_revolutions();

InterruptIn photointerrupter(D12);
Timer rev_count_timer;
double time_passed = 0;
double revolutions = 0;
uint32_t reference_time = 0;
double ang_velocity = 0;




/*************************************
             Threads
*************************************/ 


// General Global Variables
double wait_time; 
double prev_error = 0;
uint64_t prev_time = 0;
 

/* PID Revolutions */
void pid_rev();
 
Thread* th_pid_rev;
Timer pid_timer_rev;
double target_revolutions;
double max_ang_velocity = 100000;
double kp_rev;
double ki_rev;
double kd_rev;


/*   PID Velocity   */
void pid_vel();

Thread* th_pid_vel;
Timer pid_timer_vel;
double target_ang_velocity; 
double kp_vel;
double ki_vel;
double kd_vel;



/* Motor Control */
void move_field();

Thread* th_motor_control;
int8_t orState;




/*************************************
                Utility
*************************************/

/* Terminate, delete, & respawn threads */
void reset_threads();


/*************************************
                Main
*************************************/

int main() 
{
    orState = motorHome();
    intState = orState;

    
    char notes[32][2]; // For Melody

    pc.printf("Hello\n\rEnter command to begin.\n\r");
    
    /* PID Tuning */
    
    // pid_vel
    kp_vel = 0.5;
    ki_vel = 0.00000;
    kd_vel = 0.0;
    
    // pid_rev  
    kp_rev = 0.001;
    ki_rev = 0;
    kd_rev = 0;
    
    // Initial Wait Time
    wait_time = 100; // us

    /* Start Timers */
    rev_count_timer.start();
    pid_timer_rev.start();
    pid_timer_vel.start();
     
    /* Start Interrupts */
    photo_checker_interrupt.attach(&check_photo, 0.0001);
    photointerrupter.rise(&increment_revolutions);
    
    /* Create Threads */
    th_pid_rev = new Thread(osPriorityNormal, 2048);
    th_pid_vel = new Thread(osPriorityNormal, 2048);
    th_motor_control = new Thread(osPriorityNormal, 1024);
    
    while (true) 
    {        
        //pc.printf("%f, %f\n\r", ang_velocity, revolutions);
        //pc.printf("%f\n\r", revolutions);
        //pc.printf("%f, %f, %f\n\r", target_ang_velocity, ang_velocity, wait_time);
      
        
        // On keypress read in chars and put into string to be processed into commands
        if(pc.readable())
        {     

            unsigned index = 0; 
            unsigned size = 32;
            char cmd[size];
            char ch = ' ';
            
            // Read in entire command into `cmd`
            while(ch != '\r')
            {           
                ch = pc.getc();
                
                //Detect backspace
                if(ch == 127 && index > 0)
                {
                    index--;
                    pc.printf("\b");    
                }
                else
                {
                    cmd[index++] = ch;
                    pc.printf("%c", ch);
                }
            }


            // Set target revolutions
            if (cmd[0] == 'R')
            {
                char* next_cmd;
                
                // Read revolution target and convert to double
                target_revolutions = strtod(cmd+1, &next_cmd);
                
                // Set max velocity
                if(*next_cmd == 'V')
                {
                    // abs() to ignore sign
                    max_ang_velocity = abs(strtod(next_cmd+1, NULL));
                    
                    // Terminate and recreate Thread objects
                    reset_threads();
                    
                    // Reset count & wait time
                    revolutions = 0;
                    wait_time = 100;
                    
                    pc.printf("Target revolutions is %f revs\n\rMax velocity is %f\n\r",target_revolutions, max_ang_velocity);
                    
                    // Respawn threads
                    th_motor_control->start(&move_field);
                    th_pid_rev->start(&pid_rev);
                    
                }
                
                // If no velocity given
                else
                {
                    // Terminate and recreate Thread objects
                    reset_threads();
                    
                    // Reset count & wait time
                    revolutions = 0;
                    wait_time = 100;
                    
                    // No velocity limit given so set very high
                    max_ang_velocity = 100000;
                    
                    pc.printf("Target revolutions is %f revs\n\r",target_revolutions);
                    
                    // Respawn threads
                    th_motor_control->start(&move_field);
                    th_pid_rev->start(&pid_rev);

                    
                }
            }
            
            // Set target velocity
            else if (cmd[0] == 'V')
            {
                //pc.printf("V: %f\n\r",strtod(cmd+1, NULL));
                target_ang_velocity = strtod(cmd+1, NULL);
                
                // Terminate and recreate Thread objects
                reset_threads();

                // Reset wait time
                wait_time = 100;
                
                pc.printf("Target velocity is %f rev/s\n\r",target_ang_velocity);
                
                // Respawn threads
                th_motor_control->start(&move_field);
                th_pid_vel->start(&pid_vel);
                
            }
            
            // Melody
            else if (cmd[0] == 'T')
            {
                      
            }

        }
     Thread::wait(200);
    }
}
 
 
 void check_photo() 
{
    intState = readRotorState(); 
}

void increment_revolutions()
{
    revolutions++;
    
    // Calculate time for previous revolution
    double current_time = (double)rev_count_timer.read_us();
    time_passed = current_time - reference_time;
    
    ang_velocity = 1000000.0/time_passed;
    reference_time = current_time;
}


void pid_rev()
{

    while(true)
    {
       
        //pc.printf("PID \n\r");
        
        // Get time since last PID
        uint64_t time = pid_timer_rev.read_us();
        uint64_t time_since_last_pid = time - prev_time;
    
        // Calculate errors
        double error = target_revolutions - revolutions   ; // Proportional
        double error_sum = error * (double)time_since_last_pid; // Integral
        double error_deriv = (error - prev_error) / ((double)time_since_last_pid*1000000); // Derivative
        
        // Weight errors to calculate output
        double output = kp_rev*error + ki_rev*error_sum + kd_rev*error_deriv;
        
        if(OGdelta ==0) OGdelta = output;
        
        delta = output/OGdelta;
    
        if(delta>0.5) delta =0.5;
        if(delta<0) delta =0;
 
        pc.printf("%f,%f\n\r", delta, revolutions);

        
        // Lower limit output to small non-zero value to avoid divby0 error
        double output_angular_velocity = (output > 0) ? output : 0.00000001; 

        // Cap at max velocity
        //output_angular_velocity = (output_angular_velocity < max_ang_velocity) ? output_angular_velocity : max_ang_velocity;

        // Convert velocity to wait time
        //wait_time = 1000000.0/((output_angular_velocity)*6.0);

        
        // Store values for next iteration
        prev_error = error;
        prev_time = time;

        //pc.printf("%f, %f, %f, %f, %f, %f\n\r",target_revolutions, revolutions, error, error_deriv, max_ang_velocity, ang_velocity);

        Thread::wait(100); // ms    
    }
}

void pid_vel()
{
    while(true)
    {       
        // Get time since last PID
        uint64_t time = pid_timer_vel.read_us();
        uint64_t time_since_last_pid = time - prev_time;
    
        // Calculate errors
        double error = target_ang_velocity - ang_velocity   ; // Proportional
        double error_sum = error * (double)time_since_last_pid; // Integral
        double error_deriv = (error - prev_error) / ((double)time_since_last_pid*1000000.0); // Derivative
        
        // Weight errors to calculate output
        //double output = kp_vel*error + ki_vel*error_sum + kd_vel*error_deriv;
        delta = kp_vel*error + ki_vel*error_sum + kd_vel*error_deriv;
        pc.printf("%f\n\r", delta);
        
        // Lower limit output to small non-zero value to avoid divby0 error
       // double output_angular_velocity = (output > 0) ? output : 0.00000001; 
        
        // Convert velocity to wait time
        //wait_time = 1000000.0/((output_angular_velocity)*6.0);
        
        
        // Store values for next iteration
        prev_error = error;
        prev_time = time;
        
        //Debugging
        pc.printf("%f, %f, %f, %f, %f, %f\n\r",target_ang_velocity, ang_velocity, error, error_sum, error_deriv, wait_time);
        
        Thread::wait(100); // ms    
    }
}

void move_field()
{
    while(true)
    {
        //pc.printf("%f\n\r", wait_time);
        motorOut((intState-orState+lead+6)%6); //+6 to make sure the remainder is positive 
        //Thread::wait(wait_time/1000.0);  
        
    }
}


void reset_threads()
{
    th_pid_vel->terminate();
    th_pid_rev->terminate();
    th_motor_control->terminate();
    
    delete th_pid_vel;
    delete th_pid_rev;
    delete th_motor_control;
    
    th_pid_rev = new Thread(osPriorityNormal, 2048);
    th_pid_vel = new Thread(osPriorityNormal, 2048);
    th_motor_control = new Thread(osPriorityNormal, 1024);
}


