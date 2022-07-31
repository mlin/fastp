/*
SDUST complexity-filtering algorithm

Morgulis A, Gertz EM, Schäffer AA, Agarwala R. A fast and symmetric DUST implementation to mask
low-complexity DNA sequences. J Comput Biol. 2006 Jun;13(5):1028-40.
doi: 10.1089/cmb.2006.13.1028. PMID: 16796549.

This implementation is extracted from minimap2 by Heng Li:
https://github.com/lh3/minimap2/blob/master/sdust.c
the MIT license of which follows.

--

The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "sdust.h"


/* from kalloc.h */
#define kmalloc(km, size) (malloc(size))
#define kcalloc(km, count, size) (calloc(count, size))
#define krealloc(km, ptr, size) (realloc(ptr, size))
#define kfree(km, ptr) (free(ptr))


/* from kvec.h */
#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#define kvec_t(type) struct { size_t n, m; type *a; }
#define kv_init(v) ((v).n = (v).m = 0, (v).a = 0)

#define kv_resize(type, km, v, s) do { \
		if ((v).m < (s)) { \
			(v).m = (s); \
			kv_roundup32((v).m); \
			(v).a = (type*)krealloc((km), (v).a, sizeof(type) * (v).m); \
		} \
	} while (0)

#define kv_push(type, km, v, x) do { \
		if ((v).n == (v).m) { \
			(v).m = (v).m? (v).m<<1 : 2; \
			(v).a = (type*)krealloc((km), (v).a, sizeof(type) * (v).m); \
		} \
		(v).a[(v).n++] = (x); \
	} while (0)


/* from kdq.h */
#define __KDQ_TYPE(type) \
	typedef struct { \
		uint64_t front:58, bits:6, count, mask; \
		type *a; \
		void *km; \
	} kdq_##type##_t;

#define kdq_t(type) kdq_##type##_t
#define kdq_size(q) ((q)->count)
#define kdq_at(q, i) ((q)->a[((q)->front + (i)) & (q)->mask])

#define __KDQ_IMPL(type, SCOPE) \
	SCOPE kdq_##type##_t *kdq_init_##type(void *km) \
	{ \
		kdq_##type##_t *q; \
		q = (kdq_##type##_t*)kcalloc(km, 1, sizeof(kdq_##type##_t)); \
		q->bits = 2, q->mask = (1ULL<<q->bits) - 1; \
		q->a = (type*)kmalloc(km, (1<<q->bits) * sizeof(type)); \
		q->km = km; \
		return q; \
	} \
	SCOPE void kdq_destroy_##type(kdq_##type##_t *q) \
	{ \
		if (q == 0) return; \
		kfree(q->km, q->a); kfree(q->km, q); \
	} \
	SCOPE int kdq_resize_##type(kdq_##type##_t *q, int new_bits) \
	{ \
		size_t new_size = 1ULL<<new_bits, old_size = 1ULL<<q->bits; \
		if (new_size < q->count) { /* not big enough */ \
			int i; \
			for (i = 0; i < 64; ++i) \
				if (1ULL<<i > q->count) break; \
			new_bits = i, new_size = 1ULL<<new_bits; \
		} \
		if (new_bits == q->bits) return q->bits; /* unchanged */ \
		if (new_bits > q->bits) q->a = (type*)krealloc(q->km, q->a, (1ULL<<new_bits) * sizeof(type)); \
		if (q->front + q->count <= old_size) { /* unwrapped */ \
			if (q->front + q->count > new_size) /* only happens for shrinking */ \
				memmove(q->a, q->a + new_size, (q->front + q->count - new_size) * sizeof(type)); \
		} else { /* wrapped */ \
			memmove(q->a + (new_size - (old_size - q->front)), q->a + q->front, (old_size - q->front) * sizeof(type)); \
			q->front = new_size - (old_size - q->front); \
		} \
		q->bits = new_bits, q->mask = (1ULL<<q->bits) - 1; \
		if (new_bits < q->bits) q->a = (type*)krealloc(q->km, q->a, (1ULL<<new_bits) * sizeof(type)); \
		return q->bits; \
	} \
	SCOPE void kdq_push_##type(kdq_##type##_t *q, type v) \
	{ \
		if (q->count == 1ULL<<q->bits) kdq_resize_##type(q, q->bits + 1); \
		q->a[((q->count++) + q->front) & (q)->mask] = v; \
	} \
	SCOPE type *kdq_shift_##type(kdq_##type##_t *q) \
	{ \
		type *d = 0; \
		if (q->count == 0) return 0; \
		d = &q->a[q->front++]; \
		q->front &= q->mask; \
		--q->count; \
		return d; \
	}

#define KDQ_INIT2(type, SCOPE) \
	__KDQ_TYPE(type) \
	__KDQ_IMPL(type, SCOPE)

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

#define KDQ_INIT(type) KDQ_INIT2(type, static inline klib_unused)

#define kdq_init(type, km) kdq_init_##type(km)
#define kdq_destroy(type, q) kdq_destroy_##type(q)
#define kdq_resize(type, q, new_bits) kdq_resize_##type(q, new_bits)
#define kdq_push(type, q, v) kdq_push_##type(q, v)
#define kdq_shift(type, q) kdq_shift_##type(q)


/* from sdust.c */
#define SD_WLEN 3
#define SD_WTOT (1<<(SD_WLEN<<1))
#define SD_WMSK (SD_WTOT - 1)

typedef struct {
	int start, finish;
	int r, l;
} perf_intv_t;

typedef kvec_t(perf_intv_t) perf_intv_v;
typedef kvec_t(uint64_t) uint64_v;

KDQ_INIT(int)

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

struct sdust_buf_s {
	kdq_t(int) *w;
	perf_intv_v P; // the list of perfect intervals for the current window, sorted by descending start and then by ascending finish
	uint64_v res;  // the result
	void *km;      // memory pool
};

sdust_buf_t *sdust_buf_init(void *km)
{
	sdust_buf_t *buf;
	buf = (sdust_buf_t*)kcalloc(km, 1, sizeof(sdust_buf_t));
	buf->km = km;
	buf->w = kdq_init(int, buf->km);
	kdq_resize(int, buf->w, 8);
	return buf;
}

void sdust_buf_destroy(sdust_buf_t *buf)
{
	if (buf == 0) return;
	kdq_destroy(int, buf->w);
	kfree(buf->km, buf->P.a); kfree(buf->km, buf->res.a); kfree(buf->km, buf);
}

static inline void shift_window(int t, kdq_t(int) *w, int T, int W, int *L, int *rw, int *rv, int *cw, int *cv)
{
	int s;
	if ((int)kdq_size(w) >= W - SD_WLEN + 1) { // TODO: is this right for SD_WLEN!=3?
		s = *kdq_shift(int, w);
		*rw -= --cw[s];
		if (*L > (int)kdq_size(w))
			--*L, *rv -= --cv[s];
	}
	kdq_push(int, w, t);
	++*L;
	*rw += cw[t]++;
	*rv += cv[t]++;
	if (cv[t] * 10 > T<<1) {
		do {
			s = kdq_at(w, kdq_size(w) - *L);
			*rv -= --cv[s];
			--*L;
		} while (s != t);
	}
}

static inline void save_masked_regions(void *km, uint64_v *res, perf_intv_v *P, int start)
{
	int i, saved = 0;
	perf_intv_t *p;
	if (P->n == 0 || P->a[P->n - 1].start >= start) return;
	p = &P->a[P->n - 1];
	if (res->n) {
		int s = res->a[res->n - 1]>>32, f = (uint32_t)res->a[res->n - 1];
		if (p->start <= f) // if overlapping with or adjacent to the previous interval
			saved = 1, res->a[res->n - 1] = (uint64_t)s<<32 | (f > p->finish? f : p->finish);
	}
	if (!saved) kv_push(uint64_t, km, *res, (uint64_t)p->start<<32|p->finish);
	for (i = P->n - 1; i >= 0 && P->a[i].start < start; --i); // remove perfect intervals that have falled out of the window
	P->n = i + 1;
}

static void find_perfect(void *km, perf_intv_v *P, const kdq_t(int) *w, int T, int start, int L, int rv, const int *cv)
{
	int c[SD_WTOT], r = rv, i, max_r = 0, max_l = 0;
	memcpy(c, cv, SD_WTOT * sizeof(int));
	for (i = (long)kdq_size(w) - L - 1; i >= 0; --i) {
		int j, t = kdq_at(w, i), new_r, new_l;
		r += c[t]++;
		new_r = r, new_l = kdq_size(w) - i - 1;
		if (new_r * 10 > T * new_l) {
			for (j = 0; j < (int)P->n && P->a[j].start >= i + start; ++j) { // find insertion position
				perf_intv_t *p = &P->a[j];
				if (max_r == 0 || p->r * max_l > max_r * p->l)
					max_r = p->r, max_l = p->l;
			}
			if (max_r == 0 || new_r * max_l >= max_r * new_l) { // then insert
				max_r = new_r, max_l = new_l;
				if (P->n == P->m) kv_resize(perf_intv_t, km, *P, P->n + 1);
				memmove(&P->a[j+1], &P->a[j], (P->n - j) * sizeof(perf_intv_t)); // make room
				++P->n;
				P->a[j].start = i + start, P->a[j].finish = kdq_size(w) + (SD_WLEN - 1) + start;
				P->a[j].r = new_r, P->a[j].l = new_l;
			}
		}
	}
}

const uint64_t *sdust_core(const uint8_t *seq, int l_seq, int T, int W, int *n, sdust_buf_t *buf)
{
	int rv = 0, rw = 0, L = 0, cv[SD_WTOT], cw[SD_WTOT];
	int i, start, l; // _start_: start of the current window; _l_: length of a contiguous A/C/G/T (sub)sequence
	unsigned t; // current word

	buf->P.n = buf->res.n = 0;
	buf->w->front = buf->w->count = 0;
	memset(cv, 0, SD_WTOT * sizeof(int));
	memset(cw, 0, SD_WTOT * sizeof(int));
	if (l_seq < 0) l_seq = strlen((const char*)seq);
	for (i = l = t = 0; i <= l_seq; ++i) {
		int b = i < l_seq? seq_nt4_table[seq[i]] : 4;
		if (b < 4) { // an A/C/G/T base
			++l, t = (t<<2 | b) & SD_WMSK;
			if (l >= SD_WLEN) { // we have seen a word
				start = (l - W > 0? l - W : 0) + (i + 1 - l); // set the start of the current window
				save_masked_regions(buf->km, &buf->res, &buf->P, start); // save intervals falling out of the current window?
				shift_window(t, buf->w, T, W, &L, &rw, &rv, cw, cv);
				if (rw * 10 > L * T)
					find_perfect(buf->km, &buf->P, buf->w, T, start, L, rv, cv);
			}
		} else { // N or the end of sequence; N effectively breaks input into pieces of independent sequences
			start = (l - W + 1 > 0? l - W + 1 : 0) + (i + 1 - l);
			while (buf->P.n) save_masked_regions(buf->km, &buf->res, &buf->P, start++); // clear up unsaved perfect intervals
			l = t = 0;
		}
	}
	*n = buf->res.n;
	return buf->res.a;
}

uint64_t *sdust(void *km, const uint8_t *seq, int l_seq, int T, int W, int *n)
{
	uint64_t *ret;
	sdust_buf_t *buf;
	buf = sdust_buf_init(km);
	ret = (uint64_t*)sdust_core(seq, l_seq, T, W, n, buf);
	buf->res.a = 0;
	sdust_buf_destroy(buf);
	return ret;
}
