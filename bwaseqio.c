#include <zlib.h>
#include <ctype.h>
#include "bwtaln.h"
#include "utils.h"
#include "bamlite.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];
static char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

struct __bwa_seqio_t {
	// for BAM input
	int is_bam, which; // 1st bit: read1, 2nd bit: read2, 3rd: SE
	bamFile fp;
	// for fastq input
	kseq_t *ks;
};

bwa_seqio_t *bwa_bam_open(const char *fn, int which)
{
	bwa_seqio_t *bs;
	bam_header_t *h;
	bs = (bwa_seqio_t*)calloc(1, sizeof(bwa_seqio_t));
	bs->is_bam = 1;
	bs->which = which;
	bs->fp = bam_open(fn, "r");
	h = bam_header_read(bs->fp);
	bam_header_destroy(h);
	return bs;
}

bwa_seqio_t *bwa_seq_open(const char *fn)
{
	gzFile fp;
	bwa_seqio_t *bs;
	bs = (bwa_seqio_t*)calloc(1, sizeof(bwa_seqio_t));
	fp = xzopen(fn, "r");
	bs->ks = kseq_init(fp);
	return bs;
}

void bwa_seq_close(bwa_seqio_t *bs)
{
	if (bs == 0) return;
	if (bs->is_bam) bam_close(bs->fp);
	else {
		gzclose(bs->ks->f->f);
		kseq_destroy(bs->ks);
	}
	free(bs);
}

void seq_reverse(int len, ubyte_t *seq, int is_comp)
{
	int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}

int bwa_trim_read(int trim_qual, bwa_seq_t *p)
{
	int s = 0, l, max = 0, max_l = p->len;
	if (trim_qual < 1 || p->qual == 0) return 0;
	for (l = p->len - 1; l >= BWA_MIN_RDLEN; --l) {
		s += trim_qual - (p->qual[l] - 33);
		if (s < 0) break;
		if (s > max) max = s, max_l = l;
	}
	p->clip_len = p->len = max_l;
	return p->full_len - p->len;
}

static bwa_seq_t *bwa_read_bam(bwa_seqio_t *bs, int n_needed, int *n, int is_comp, int trim_qual)
{
	bwa_seq_t *seqs, *p;
	int n_seqs, l, i;
	long n_trimmed = 0, n_tot = 0;
	bam1_t *b;

	b = bam_init1();
	n_seqs = 0;
	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	while (bam_read1(bs->fp, b) >= 0) {
		uint8_t *s, *q;
		int go = 0;
		if ((bs->which & 1) && (b->core.flag & BAM_FREAD1)) go = 1;
		if ((bs->which & 2) && (b->core.flag & BAM_FREAD2)) go = 1;
		if ((bs->which & 4) && !(b->core.flag& BAM_FREAD1) && !(b->core.flag& BAM_FREAD2))go = 1;
		if (go == 0) continue;
		l = b->core.l_qseq;
		p = &seqs[n_seqs++];
		p->tid = -1; // no assigned to a thread
		p->qual = 0;
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
		s = bam1_seq(b); q = bam1_qual(b);
		p->seq = (ubyte_t*)calloc(p->len + 1, 1);
		p->qual = (ubyte_t*)calloc(p->len + 1, 1);
		for (i = 0; i != p->full_len; ++i) {
			p->seq[i] = bam_nt16_nt4_table[(int)bam1_seqi(s, i)];
			p->qual[i] = q[i] + 33 < 126? q[i] + 33 : 126;
		}
		if (bam1_strand(b)) { // then reverse 
			seq_reverse(p->len, p->seq, 1);
			seq_reverse(p->len, p->qual, 0);
		}
		if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		p->rseq = (ubyte_t*)calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*)bam1_qname(b));
		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
	if (n_seqs == 0) {
		free(seqs);
		bam_destroy1(b);
		return 0;
	}
	bam_destroy1(b);
	return seqs;
}

#define BARCODE_LOW_QUAL 13

   static int ks_getuntil1(kstream_t *ks, int delimiter, kstring_t *str, int *dret)
   {
      if (dret) *dret = 0;
      str->l = 0;
      if (ks->begin >= ks->end && ks->is_eof) return -1;
      for (;;) {
         int i;
         if (ks->begin >= ks->end) {
            if (!ks->is_eof) {
               ks->begin = 0;
               ks->end = __read(ks->f, ks->buf, 4096);
               if (ks->end < 4096) ks->is_eof = 1;
               //ks->end = __read(ks->f, ks->buf, __bufsize);
               //if (ks->end < __bufsize) ks->is_eof = 1;
               if (ks->end == 0) break;
            } else break;
         }
         if (delimiter) {
            for (i = ks->begin; i < ks->end; ++i)
               if (ks->buf[i] == delimiter) break;
         } else {
            for (i = ks->begin; i < ks->end; ++i)
               if (isspace(ks->buf[i])) break;
         }
         if (str->m - str->l < i - ks->begin + 1) {
            str->m = str->l + (i - ks->begin) + 1;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
         }
         memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
         str->l = str->l + (i - ks->begin);
         ks->begin = i + 1;
         if (i < ks->end) {
            if (dret) *dret = ks->buf[i];
            break;
         }
      }
      str->s[str->l] = '\0';
      return str->l;
   }


/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
static int kseq_read1(kseq_t *seq)
{
   int c;
   kstream_t *ks = seq->f;
   if (seq->last_char == 0) { /* then jump to the next header line */
      while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
      if (c == -1) return -1; /* end of file */
      seq->last_char = c;
   } /* the first header char has been read */
   seq->comment.l = seq->seq.l = seq->qual.l = 0;
   if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1;
   if (c != '\n') ks_getuntil(ks, '\n', &seq->comment, 0);
   while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
      if (isgraph(c)) { /* printable non-space character */
         if (seq->seq.l + 1 >= seq->seq.m) { /* double the memory */
            seq->seq.m = seq->seq.l + 2;
            kroundup32(seq->seq.m); /* rounded to next closest 2^k */
            seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);
         }
         seq->seq.s[seq->seq.l++] = (char)c;
      }
   }
   if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
   seq->seq.s[seq->seq.l] = 0;   /* null terminated string */
   if (c != '+') return seq->seq.l; /* FASTA */
   if (seq->qual.m < seq->seq.m) {  /* allocate enough memory */
      seq->qual.m = seq->seq.m;
      seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);
   }
   while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
   if (c == -1) return -2; /* we should not stop here */
   while ((c = ks_getc(ks)) != -1 && seq->qual.l < seq->seq.l)
      if (c >= 33 && c <= 127) seq->qual.s[seq->qual.l++] = (unsigned char)c;
   seq->qual.s[seq->qual.l] = 0; /* null terminated string */
   seq->last_char = 0;  /* we have not come to the next header line */
   if (seq->seq.l != seq->qual.l) return -2; /* qual string is shorter than seq string */
   return seq->seq.l;
}

inline static void init_bwa_seq_t(bwa_seq_t *p)
{
   if(p->cigar) free(p->cigar);
   if(p->multi) free(p->multi);
   p->n_multi = p->n_cigar = 0;
   if(p->md) free(p->md);
   p->md = 0;
   p->sa = p->pos = 0;
   p->c1 = p->c2 = p->seQ = 0;
   p->nm = 0;
   p->strand = p->type = p->dummy = p->extra_flag = 0;
   p->n_mm = p->n_gapo = p->n_gape = p->mapQ = 0;
   p->score = 0;
   return;
}

int bwa_read_seq1(bwa_seqio_t *bs, bwa_seq_t **_seqs, int *n_avail, int n_needed, int mode, int trim_qual)
{
	bwa_seq_t *p;
   bwa_seq_t *seqs = *_seqs;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode>>24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > BWA_MAX_BCLEN) {
		fprintf(stderr, "[%s] the maximum barcode length is %d.\n", __func__, BWA_MAX_BCLEN);
		return 0;
	}
	if (bs->is_bam) {
      fprintf (stderr, "IS BAM! --- Port bwa_read_bam function\n");
      bwa_free_read_seq(*n_avail, seqs);
      // return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
      exit(0);
   }
   if (*n_avail < n_needed) {
 //fprintf(stderr, "Allocating memory\n");
      if (*n_avail) bwa_free_read_seq(*n_avail, seqs);
	   seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
      *_seqs = seqs;
      *n_avail = n_needed;
   }
	n_seqs = 0;
	while ((l = kseq_read1(seq)) >= 0) {
		if ((mode & BWA_MODE_CFY) && (seq->comment.l != 0)) {
			// skip reads that are marked to be filtered by Casava
			char *s = index(seq->comment.s, ':');
			if (s && *(++s) == 'Y') {
				continue;
			}
		}
		if (is_64 && seq->qual.l)
			for (i = 0; i < seq->qual.l; ++i) seq->qual.s[i] -= 31;
		if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
		p = &seqs[n_seqs++];
      init_bwa_seq_t(p);
		if (l_bc) { // then trim barcode
			for (i = 0; i < l_bc; ++i)
				p->bc[i] = (seq->qual.l && seq->qual.s[i]-33 < BARCODE_LOW_QUAL)? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
			p->bc[i] = 0;
			for (; i < seq->seq.l; ++i)
				seq->seq.s[i - l_bc] = seq->seq.s[i];
			seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
			if (seq->qual.l) {
				for (i = l_bc; i < seq->qual.l; ++i)
					seq->qual.s[i - l_bc] = seq->qual.s[i];
				seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
			}
			l = seq->seq.l;
		} else p->bc[0] = 0;
		p->tid = -1; // no assigned to a thread
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
      if (p->llen < p->len) {
         if(p->seq) free(p->seq);
         p->llen = p->len;
		   p->seq = (ubyte_t*)calloc(p->len, 1);
		   p->rseq = (ubyte_t*)calloc(p->full_len, 1);
      }
		for (i = 0; i != p->full_len; ++i)
			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
		if (seq->qual.l) { // copy quality
         if(p->qual) free(p->qual);
			p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		}
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
      if(p->name) free(p->name);
		p->name = strdup((const char*)seq->name.s);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) p->name[t-2] = '\0';
		}
		if (n_seqs == n_needed) break;
	}
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
	return n_seqs;
}

bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual)
{
	bwa_seq_t *seqs, *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode>>24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > BWA_MAX_BCLEN) {
		fprintf(stderr, "[%s] the maximum barcode length is %d.\n", __func__, BWA_MAX_BCLEN);
		return 0;
	}
	if (bs->is_bam) return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
	n_seqs = 0;
	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	while ((l = kseq_read1(seq)) >= 0) {
		if ((mode & BWA_MODE_CFY) && (seq->comment.l != 0)) {
			// skip reads that are marked to be filtered by Casava
			char *s = index(seq->comment.s, ':');
			if (s && *(++s) == 'Y') {
				continue;
			}
		}
		if (is_64 && seq->qual.l)
			for (i = 0; i < seq->qual.l; ++i) seq->qual.s[i] -= 31;
		if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
		p = &seqs[n_seqs++];
		if (l_bc) { // then trim barcode
			for (i = 0; i < l_bc; ++i)
				p->bc[i] = (seq->qual.l && seq->qual.s[i]-33 < BARCODE_LOW_QUAL)? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
			p->bc[i] = 0;
			for (; i < seq->seq.l; ++i)
				seq->seq.s[i - l_bc] = seq->seq.s[i];
			seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
			if (seq->qual.l) {
				for (i = l_bc; i < seq->qual.l; ++i)
					seq->qual.s[i - l_bc] = seq->qual.s[i];
				seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
			}
			l = seq->seq.l;
		} else p->bc[0] = 0;
		p->tid = -1; // no assigned to a thread
		p->qual = 0;
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
		p->seq = (ubyte_t*)calloc(p->len, 1);
		for (i = 0; i != p->full_len; ++i)
			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
		if (seq->qual.l) { // copy quality
			p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		}
		p->rseq = (ubyte_t*)calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*)seq->name.s);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t-2] == '/' && (p->name[t-1] == '1' || p->name[t-1] == '2')) p->name[t-2] = '\0';
		}
		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed/n_tot);
	if (n_seqs == 0) {
		free(seqs);
		return 0;
	}
	return seqs;
}

void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs)
{
	int i, j;
 //fprintf(stderr, "Freeing memory\n");
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		for (j = 0; j < p->n_multi; ++j)
			if (p->multi[j].cigar) free(p->multi[j].cigar);
		free(p->name);
		free(p->seq); free(p->rseq); free(p->qual); free(p->aln); free(p->md); free(p->multi);
		free(p->cigar);
	}
	free(seqs);
}
