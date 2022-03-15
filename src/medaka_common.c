#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "medaka_common.h"


/** Allocates zero-initialised memory with a message on failure.
 *
 *  @param num number of elements to allocate.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xalloc(size_t num, size_t size, char* msg){
    void *res = calloc(num, size);
    if (res == NULL){
        fprintf(stderr, "Failed to allocate mem for %s\n", msg);
        exit(1);
    }
    return res;
}


/** Reallocates memory with a message on failure.
 *
 *  @param ptr pointer to realloc.
 *  @param size size of each element.
 *  @param msg message to describe allocation on failure.
 *  @returns pointer to allocated memory
 *
 */
void *xrealloc(void *ptr, size_t size, char* msg){
    void *res = realloc(ptr, size);
    if (res == NULL){
        fprintf(stderr, "Failed to reallocate mem for %s\n", msg);
        exit(1);
    }
    return res;
}


/** Retrieves a substring.
 *
 *  @param string input string.
 *  @param postion start position of substring.
 *  @param length length of substring required.
 *  @returns string pointer.
 *
 */
char *substring(char *string, int position, int length) {
   char *ptr;
   size_t i;

   ptr = malloc(length + 1);

   for (i = 0 ; i < length ; i++) {
      *(ptr + i) = *(string + position);
      string++;
   }

   *(ptr + i) = '\0';
   return ptr;
}


/** Format a uint32_t to a string
 *
 * @param value to format.
 * @param dst destination char.
 * @returns length of string.
 *
 */
size_t uint8_to_str(uint8_t value, char *dst) {
    static char* digits[] = {
        "0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",
        "21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40",
        "41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60",
        "61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80",
        "81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100",
        "101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120",
        "121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140",
        "141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160",
        "161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180",
        "181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200",
        "201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220",
        "221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240",
        "241","242","243","244","245","246","247","248","249","250","251","252","253","254","255"};
    static const uint8_t TEN = 10;
    static const uint8_t HUNDRED = 100;
    strcpy(dst, digits[value]);
    if (value < TEN) return 1;
    if (value < HUNDRED) return 2;
    else return 3;
}

