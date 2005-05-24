
#ifndef WIN_LINUX_AUX_H
#define WIN_LINUX_AUX_H 1

int getTIFFinfo(char *filename, 
		long unsigned int *length, long unsigned int *width);

int getTIFFimage(char *filename, long unsigned int length,
		 long unsigned int width, int *r);

#endif /* WIN_LINUX_AUX_H */


