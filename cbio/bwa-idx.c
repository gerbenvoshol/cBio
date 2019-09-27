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
	if (argc < 2) {
		fprintf(stderr, "USAGE: %s fasta\n", argv[0]);
		return 1;
	}

	libbwa_index(argv[1], argv[1], LIBBWA_INDEX_ALGO_AUTO, 0);

	return 0;
}