void itoa(int temp_n, char temp_s[])
{
    int i;
    i = 0;
    do {       /* generate digits in reverse order */
        temp_s[i++] = temp_n % 10 + '0';   /* get next digit */
    } while ((temp_n /= 10) > 0);     /* delete it */
    temp_s[i] = '\0';
    reverse(temp_s);
}
