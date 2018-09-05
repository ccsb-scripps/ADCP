/**********************************************************************/
/* 16-bit multiply with carry generator x(n)=a*x(n-1)+carry mod 2^16, */
/* number and carry packed within the same 32-bit integer.            */
/**********************************************************************/

/* adopted from Glenn Rhoads, http://remus.rutgers.edu/~rhoads/ */

#define RAND16_MAX 0xFFFF

extern unsigned int rand16(void);	/* returns a random 16-bit integer */
extern void srand16(unsigned int);	/* seed the generator */
