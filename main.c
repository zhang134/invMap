#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

static ko_longopt_t long_options[] = {
	{ "bucket-bits",    ko_required_argument, 300 },
	{ "mb-size",        ko_required_argument, 'K' },
	{ "seed",           ko_required_argument, 302 },
	{ "no-kalloc",      ko_no_argument,       303 },
	{ "print-qname",    ko_no_argument,       304 },
	{ "no-self",        ko_no_argument,       'D' },
	{ "print-seeds",    ko_no_argument,       306 },
	{ "max-chain-skip", ko_required_argument, 307 },
	{ "min-dp-len",     ko_required_argument, 308 },
	{ "print-aln-seq",  ko_no_argument,       309 },
	{ "splice",         ko_no_argument,       310 },
	{ "cost-non-gt-ag", ko_required_argument, 'C' },
	{ "no-long-join",   ko_no_argument,       312 },
	{ "sr",             ko_no_argument,       313 },
	{ "frag",           ko_required_argument, 314 },
	{ "secondary",      ko_required_argument, 315 },
	{ "cs",             ko_optional_argument, 316 },
	{ "end-bonus",      ko_required_argument, 317 },
	{ "no-pairing",     ko_no_argument,       318 },
	{ "splice-flank",   ko_required_argument, 319 },
	{ "idx-no-seq",     ko_no_argument,       320 },
	{ "end-seed-pen",   ko_required_argument, 321 },
	{ "for-only",       ko_no_argument,       322 },
	{ "rev-only",       ko_no_argument,       323 },
	{ "heap-sort",      ko_required_argument, 324 },
	{ "all-chain",      ko_no_argument,       'P' },
	{ "dual",           ko_required_argument, 326 },
	{ "max-clip-ratio", ko_required_argument, 327 },
	{ "min-occ-floor",  ko_required_argument, 328 },
	{ "MD",             ko_no_argument,       329 },
	{ "lj-min-ratio",   ko_required_argument, 330 },
	{ "score-N",        ko_required_argument, 331 },
	{ "eqx",            ko_no_argument,       332 },
	{ "paf-no-hit",     ko_no_argument,       333 },
	{ "split-prefix",   ko_required_argument, 334 },
	{ "no-end-flt",     ko_no_argument,       335 },
	{ "hard-mask-level",ko_no_argument,       336 },
	{ "cap-sw-mem",     ko_required_argument, 337 },
	{ "max-qlen",       ko_required_argument, 338 },
	{ "max-chain-iter", ko_required_argument, 339 },
	{ "junc-bed",       ko_required_argument, 340 },
	{ "junc-bonus",     ko_required_argument, 341 },
	{ "sam-hit-only",   ko_no_argument,       342 },
	{ "chain-gap-scale",ko_required_argument, 343 },
	{ "alt",            ko_required_argument, 344 },
	{ "alt-drop",       ko_required_argument, 345 },
	{ "mask-len",       ko_required_argument, 346 },
	{ "rmq",            ko_optional_argument, 347 },
	{ "qstrand",        ko_no_argument,       348 },
	{ "cap-kalloc",     ko_required_argument, 349 },
	{ "q-occ-frac",     ko_required_argument, 350 },
	{ "chain-skip-scale",ko_required_argument,351 },
	{ "print-chains",   ko_no_argument,       352 },
	{ "no-hash-name",   ko_no_argument,       353 },
	{ "secondary-seq",  ko_no_argument,       354 },
	{ "help",           ko_no_argument,       'h' },
	{ "max-intron-len", ko_required_argument, 'G' },
	{ "version",        ko_no_argument,       'V' },
	{ "min-count",      ko_required_argument, 'n' },
	{ "min-chain-score",ko_required_argument, 'm' },
	{ "mask-level",     ko_required_argument, 'M' },
	{ "min-dp-score",   ko_required_argument, 's' },
	{ "sam",            ko_no_argument,       'a' },
	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t mm_parse_num(const char *str)
{
	return mm_parse_num2(str, 0);
}

static inline void yes_or_no(mm_mapopt_t *opt, int64_t flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int main(int argc, char *argv[])
{
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yYPo:e:U:J:";
	ketopt_t o = KETOPT_INIT;
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = get_nprocs_conf(), n_parts, old_best_n = -1;
	char *fnw = 0, *rg = 0, *junc_bed = 0, *s, *alt_list = 0;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(o.arg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}
	}
	o = KETOPT_INIT;
	opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg);
		else if (c == 'k') ipt.k = atoi(o.arg);		
		else if (c == 't') n_threads = atoi(o.arg);
	}
	if ((opt.flag & MM_F_SPLICE) && (opt.flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR))
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;
	if (opt.best_n == 0) {
		fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
		opt.best_n = old_best_n, opt.flag |= MM_F_NO_PRINT_2ND;
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "---------------------------------------------------------------\n");
		fprintf(fp_help, "Usage: invmap [options] target.fa query.fa >output.sam\n");
		fprintf(fp_help, "\nExample: invmap genome.fa reads.fa >reads.sam\n");
		fprintf(fp_help, "\nOptions:\n");		
		fprintf(fp_help, "    -k        k-mer size (should <= 25), default: %d.\n", ipt.k);
		fprintf(fp_help, "    -w        sample window size, default: %d.\n", ipt.w);		
		fprintf(fp_help, "    -t        number of threads, default: %d (your computer has).\n", n_threads);
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - o.ind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	int n_threads_total;
	n_threads_total = get_nprocs_conf();
	time_t t;
	char buf[1024];
	time(&t);
	ctime_r(&t,buf);
	fprintf(stderr, "===================================================================\n");
	fprintf(stderr, "[INFO] program started at: %s", buf);
	fprintf(stderr, "[INFO] k-mer size: %d\n", ipt.k);
	fprintf(stderr, "[INFO] total threads: %d\n", n_threads_total);
	fprintf(stderr, "[INFO] threads you set: %d\n", n_threads);
	fprintf(stderr, "[INFO] sample k-mer window length: %d\n", ipt.w);
	fprintf(stderr, "[INFO] references details:\n");
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		int ret;
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				if (opt.split_prefix == 0)
					ret = mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
				else
					ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
			} else {
				ret = mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (opt.split_prefix == 0 && mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\033[0m\n");
			}
			if (ret != 0) {
				mm_idx_destroy(mi);
				mm_idx_reader_close(idx_rdr);
				return 1;
			}
		}
		if (mm_verbose >= 3){
			//fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
			//		__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);

		}
		//fprintf(stderr, "[INFO] build the hash table for %d target sequence(s)\n", mi->n_seq);
		if (argc != o.ind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		time(&t);
		ctime_r(&t,buf);
		fprintf(stderr, "[INFO] mapping reads at: %s\n", buf);
		if (junc_bed) mm_idx_bed_read(mi, junc_bed, 1);
		if (alt_list) mm_idx_alt_read(mi, alt_list);
		if (argc - (o.ind + 1) == 0) {
			mm_idx_destroy(mi);
			continue; // no query files
		}
		ret = 0;
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = o.ind + 1; i < argc; ++i) {
				ret = mm_map_file(mi, argv[i], &opt, n_threads);
				if (ret < 0) break;
			}
		} else {
			ret = mm_map_file_frag(mi, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
		if (ret < 0) {
			fprintf(stderr, "ERROR: failed to map the query file\n");
			exit(EXIT_FAILURE);
		}
	}
	n_parts = idx_rdr->n_parts;
	mm_idx_reader_close(idx_rdr);

	if (opt.split_prefix)
		mm_split_merge(argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_parts);

	if (fflush(stdout) == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}

	if (mm_verbose >= 3) {
		//fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		//fprintf(stderr, "[M::%s] CMD:", __func__);
		//for (i = 0; i < argc; ++i)
		//	fprintf(stderr, " %s", argv[i]);
		//fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mm_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	time(&t);
	ctime_r(&t,buf);
	fprintf(stderr, "[INFO] program ended at: %s", buf);
	return 0;
}
