/**********************************************************************/
/* 16-bit multiply with carry generator x(n)=a*x(n-1)+carry mod 2^16, */
/* number and carry packed within the same 32-bit integer.            */
/**********************************************************************/

/* adopted from Glenn Rhoads, http://remus.rutgers.edu/~rhoads/ */

#define RAND16_MAX 0xFFFF

unsigned int rand16(void);	/* returns a random 16-bit integer */
void srand16(unsigned int);	/* seed the generator */

static unsigned int SEED = 521288629;

unsigned int rand16()
{
/* Use any number from this list for "a"
    18000 18030 18273 18513 18879 19074 19098 19164 19215 19584       
    19599 19950 20088 20508 20544 20664 20814 20970 21153 21243       
    21423 21723 21954 22125 22188 22293 22860 22938 22965 22974       
    23109 23124 23163 23208 23508 23520 23553 23658 23865 24114       
    24219 24660 24699 24864 24948 25023 25308 25443 26004 26088       
    26154 26550 26679 26838 27183 27258 27753 27795 27810 27834       
    27960 28320 28380 28689 28710 28794 28854 28959 28980 29013       
    29379 29889 30135 30345 30459 30714 30903 30963 31059 31083
*/
	static unsigned int a = 30903;

	SEED = a * (SEED & RAND16_MAX) + (SEED >> 16);

	return SEED & RAND16_MAX;
}

void srand16(unsigned int seed)
{
	if (seed)
		SEED = seed;	/* use default seeds if parameter is 0 */
}
