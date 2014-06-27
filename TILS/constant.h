#ifndef CONSTANT_H
#define CONSTANT_H

//Global constant
/* file names */
#define FILE_M      ".\\testdata\\zipf\\M1000000.fnl"
#define FILE_F      ".\\testdata\\zipf\\F50000.fnl" 
#define FILE_C      ".\\testdata\\zipf\\uC10000.fnl"
/* file size  */
#define FILE_M_SIZE 1000000 
#define FILE_F_SIZE  50000 
#define FILE_C_SIZE  10000
/* fanout for Rtrees */
#define FANOUT 102
/* variable k */
#define K  10
/* if a number small than EPSILON, we consider it equal to 0 */
#define EPSILON 0.00001



#endif CONSTANT_H