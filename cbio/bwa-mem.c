#include "cbio.h"

/*
#include <libbwa.h>

void index_example(void)
{
    char *db = "path/to/reference.fa";
    char *prefix = "path/to/reference.fa";

    // Equivalent to `bwa index` command
    libbwa_index(db, prefix, LIBBWA_INDEX_ALGO_AUTO, 0);
}

void mem_example(void)
{
    char *db = "path/to/reference.fa";
    char *read = "path/to/read.fq";
    char *out = "path/to/out.sam";

    // Create option struct
    libbwa_mem_opt *opt = libbwa_mem_opt_init();

    // Equivalent to `bwa mem` command
    libbwa_mem(db, read, NULL, out, opt);

    // Release option
    libbwa_mem_opt_destroy(opt);
}
*/

int main(int argc, char const *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "USAGE: %s db reads out\n", argv[0]);
		return 1;
	}

	char *db = argv[1];
	char *read = argv[2];
	char *out = argv[3];

	// Create option struct
	libbwa_mem_opt *opt = libbwa_mem_opt_init();

	// Equivalent to `bwa mem` command
	libbwa_mem(db, read, NULL, out, opt);

	// Release option
	libbwa_mem_opt_destroy(opt);

	return 0;
}