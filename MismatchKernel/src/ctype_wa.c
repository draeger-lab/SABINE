#include <ctype.h>


__const unsigned short int *__ctype_b;
__const __int32_t *__ctype_tolower;
__const __int32_t *__ctype_toupper;


void ctSetup()
{
__ctype_b = *(__ctype_b_loc());
__ctype_toupper = *(__ctype_toupper_loc());
__ctype_tolower = *(__ctype_tolower_loc());
}
