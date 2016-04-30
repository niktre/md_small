#include <string.h>

/* reverse:  reverse string s in place */
void reverse(char temp_s[])
{
    int c, i, j;

    for (i = 0, j = strlen(temp_s)-1; i<j; i++, j--) {
        c = temp_s[i];
        temp_s[i] = temp_s[j];
        temp_s[j] = c;
    }
}
