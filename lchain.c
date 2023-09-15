#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "mmpriv.h"
#include "kalloc.h"
#include "krmq.h"

static int64_t mg_chain_bk_end(int32_t max_drop, const mm128_t *z, const int32_t *f, const int64_t *p, int32_t *t, int64_t k)
{
	int64_t i = z[k].y, end_i = -1, max_i = i;
	int32_t max_s = 0;
	if (i < 0 || t[i] != 0) return i;
	do {
		int32_t s;
		t[i] = 2;
		end_i = i = p[i];
		s = i < 0? z[k].x : (int32_t)z[k].x - f[i];
		if (s > max_s) max_s = s, max_i = i;
		else if (max_s - s > max_drop) break;
	} while (i >= 0 && t[i] == 0);
	for (i = z[k].y; i >= 0 && i != end_i; i = p[i]) // reset modified t[]
		t[i] = 0;
	return max_i;
}

uint64_t *mg_chain_backtrack(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt1, int32_t min_sc1, int32_t max_drop, int32_t *n_u_, int32_t *n_v_)
{
	mm128_t *z;
	uint64_t *u;
	int64_t i, k, n_z, n_v;
	int32_t n_u;
	int32_t min_cnt = 3, min_sc = 30; //min_sc1


	*n_u_ = *n_v_ = 0;
	for (i = 0, n_z = 0; i < n; ++i) {// precompute n_z
		if (f[i] >= min_sc) {
			++n_z;
		}
	}
	if (n_z == 0) return 0;
	z = Kmalloc(km, mm128_t, n_z);
	for (i = 0, k = 0; i < n; ++i) {// populate z[]
		if (f[i] >= min_sc) {
			z[k].x = f[i], z[k++].y = i;
		}
	}
	radix_sort_128x(z, z + n_z);

	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				++n_v, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				++n_u;
			else n_v = n_v0;
		}
	}
	u = Kmalloc(km, uint64_t, n_u);
	memset(t, 0, n * 4);
	for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
		if (t[z[k].y] == 0) {
			int64_t n_v0 = n_v, end_i;
			int32_t sc;
			end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
			for (i = z[k].y; i != end_i; i = p[i])
				v[n_v++] = i, t[i] = 1;
			sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
			if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
				u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
			else n_v = n_v0;
		}
	}
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}
uint64_t *mg_chain_backtrack_my(void *km, int64_t n, const int32_t *f, const int64_t *p, int32_t *v, int32_t *t, int32_t min_cnt1, 
	int32_t min_sc1, int32_t max_drop, int32_t *n_u_, int32_t *n_v_, int64_t index_start, int64_t index_middle, int64_t index_end)
{
	mm128_t *z;
	uint64_t *u, curr, sc;//, score_end, score_middle, score_start;
	int64_t i, k, n_z, n_v, n_v0;
	int32_t n_u;
	int32_t min_cnt = 2, min_sc = 16; //min_sc1


	*n_u_ = *n_v_ = 0;
	//-for (i = 0, n_z = 0; i < n; ++i) {// precompute n_z
	//-	if (f[i] >= min_sc) {
	//-		++n_z;
	//-	}
	//-}
	n_z = 0;
	n_u = 0;
	curr = index_end;
	if (p[curr] > -1){
		n_u++;
	}
	while(p[curr] > -1){
		//fprintf(stderr, "Curr =  %10d\n", curr);
		++n_z;
		curr = p[curr];
	}
	//score_end = f[index_end] - f[curr];
	curr = index_middle;
	if (p[curr] > -1){
		n_u++;
	}
	while(p[curr] > -1){
		++n_z;
		curr = p[curr];
	}
	//++n_z;
	curr = index_start;
	if (p[curr] > -1){
		n_u++;
	}
	while(p[curr] > -1){
		++n_z;
		curr = p[curr];
	}
	//++n_z;
	if (n_z == 0) {
		return 0;
	}
	z = Kmalloc(km, mm128_t, n_z);
	k = 0;
	i = index_end;
	while(p[i] > -1){
		z[k].x = f[i];
		z[k++].y = i;
		i = p[i];
	}
	i = index_middle;
	while(p[i] > -1){
		z[k].x = f[i];
		z[k++].y = i;
		i = p[i];
	}
	i = index_start;
	while(p[i] > -1){
		z[k].x = f[i];
		z[k++].y = i;
		i = p[i];
	}
	//-for (i = 0, k = 0; i < n; ++i) {// populate z[]
	//-	if (f[i] >= min_sc) {
	//-		z[k].x = f[i], z[k++].y = i;
	//-	}
	//-}

	radix_sort_128x(z, z + n_z);

	//memset(t, 0, n * 4);
	// for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // precompute n_u
	// 	if (t[z[k].y] == 0) {
	// 		int64_t n_v0 = n_v, end_i;
	// 		int32_t sc;
	// 		end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
	// 		for (i = z[k].y; i != end_i; i = p[i])
	// 			++n_v, t[i] = 1;
	// 		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
	// 		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
	// 			++n_u;
	// 		else n_v = n_v0;
	// 	}
	// }
	//n_u = 3;
	u = Kmalloc(km, uint64_t, n_u);
	memset(t, 0, n * 4);
	n_v = 0;
	n_v0 = 0;
	n_u = 0;
	//if(f[index_end] > f[index_middle] && f[index_end] > f[index_start]){
	//}
	curr = index_end;
	while(p[curr] > -1){
		v[n_v++] = curr;
		t[curr] = 1;
		curr = p[curr];
	}
	sc = f[index_end];// - f[curr];
	u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
	n_v0 = n_v;
	curr = index_middle;
	while(p[curr] > -1){
		v[n_v++] = curr;
		t[curr] = 1;
		curr = p[curr];
	}
	sc = f[index_middle];// - f[curr];
	u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
	n_v0 = n_v;
	curr = index_start;
	while(p[curr] > -1){
		v[n_v++] = curr;
		t[curr] = 1;
		curr = p[curr];
	}
	sc = f[index_start];// - f[curr];
	u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
	// for (k = n_z - 1, n_v = n_u = 0; k >= 0; --k) { // populate u[]
	// 	if (t[z[k].y] == 0) {
	// 		int64_t n_v0 = n_v, end_i;
	// 		int32_t sc;
	// 		end_i = mg_chain_bk_end(max_drop, z, f, p, t, k);
	// 		for (i = z[k].y; i != end_i; i = p[i])
	// 			v[n_v++] = i, t[i] = 1;
	// 		sc = i < 0? z[k].x : (int32_t)z[k].x - f[i];
	// 		if (sc >= min_sc && n_v > n_v0 && n_v - n_v0 >= min_cnt)
	// 			u[n_u++] = (uint64_t)sc << 32 | (n_v - n_v0);
	// 		else n_v = n_v0;
	// 	}
	// }
	kfree(km, z);
	assert(n_v < INT32_MAX);
	*n_u_ = n_u, *n_v_ = n_v;
	return u;
}
static mm128_t *compact_a(void *km, int32_t n_u, uint64_t *u, int32_t n_v, int32_t *v, mm128_t *a)
{
	mm128_t *b, *w;
	uint64_t *u2;
	int64_t i, j, k;

	// write the result to b[]
	b = Kmalloc(km, mm128_t, n_v);
	for (i = 0, k = 0; i < n_u; ++i) {
		int32_t k0 = k, ni = (int32_t)u[i];
		for (j = 0; j < ni; ++j)
			b[k++] = a[v[k0 + (ni - j - 1)]];
	}
	kfree(km, v);

	// sort u[] and a[] by the target position, such that adjacent chains may be joined
	w = Kmalloc(km, mm128_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		w[i].x = b[k].x, w[i].y = (uint64_t)k<<32|i;
		k += (int32_t)u[i];
	}
	radix_sort_128x(w, w + n_u);
	u2 = Kmalloc(km, uint64_t, n_u);
	for (i = k = 0; i < n_u; ++i) {
		int32_t j = (int32_t)w[i].y, n = (int32_t)u[j];
		u2[i] = u[j];
		memcpy(&a[k], &b[w[i].y>>32], n * sizeof(mm128_t));
		k += n;
	}
	memcpy(u, u2, n_u * 8);
	memcpy(b, a, k * sizeof(mm128_t)); // write _a_ to _b_ and deallocate _a_ because _a_ is oversized, sometimes a lot
	kfree(km, a); kfree(km, w); kfree(km, u2);
	return b;
}

static inline int32_t comput_sc(const mm128_t *ai, const mm128_t *aj, int32_t max_dist_x, int32_t max_dist_y, int32_t bw, float chn_pen_gap, float chn_pen_skip, int is_cdna, int n_seg)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	int32_t sidi = (ai->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	int32_t sidj = (aj->y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
	if (dq <= 0 || dq > max_dist_x) return INT32_MIN;
	dr = (int32_t)(ai->x - aj->x);
	if (sidi == sidj && (dr == 0 || dq > max_dist_y)) return INT32_MIN;
	dd = dr > dq? dr - dq : dq - dr;
	if (sidi == sidj && dd > bw) return INT32_MIN;
	if (n_seg > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) return INT32_MIN;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (dd || dg > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		if (is_cdna || sidi != sidj) {
			if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
			else if (dr > dq || sidi != sidj) sc -= (int)(lin_pen < log_pen? lin_pen : log_pen); // deletion or jump between paired ends
			else sc -= (int)(lin_pen + .5f * log_pen);
		} else sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

/* Input:
 *   a[].x: rev<<63 | tid<<32 | tpos
 *   a[].y: flags<<40 | q_span<<32 | q_pos
 * Output:
 *   n_u: #chains
 *   u[]: score<<32 | #anchors (sum of lower 32 bits of u[] is the returned length of a[])
 * input a[] is deallocated on return
 */
mm128_t *mg_lchain_dp(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw;
	int64_t *p, i, j, max_ii, st = 0, n_iter = 0;
	uint64_t *u;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw && !is_cdna) max_dist_y = bw;
	if (is_cdna) max_drop = INT32_MAX;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);

	// fill the score and backtrack arrays
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			++n_iter;
			if (sc == INT32_MIN) continue;
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		end_j = j;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) max = f[j], max_ii = j;
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) mmax_f = max_f;
	}

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}

mm128_t *mg_lchain_dp_qlen(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					  int is_cdna, int n_seg, int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km, const int qlen, int32_t *pos_inv_left, int32_t *pos_inv_right)
{ // TODO: make sure this works when n has more than 32 bits
	int32_t *f, *f1, *t, *v, n_u, n_v, mmax_f = 0, max_drop = bw, max_plus_score = INT32_MIN, max_reverse_score = INT32_MIN;
	int32_t rid, read_main_strand, anchors_num = 0;//pos_inv_left, pos_inv_right
	int32_t score_increased_max, max_score_i, left_ref, right_ref, d_read, d_ref, sc, max_score = 0, start_i, score_increased, final_start, final_end;
	int64_t *p, *p1, i, j, max_ii, st = 0, n_iter = 0, left_i, max_j, end_j, sv_index_start = 0, sv_index_end = 0;
	int64_t start1, start2, end1, end2; // the 2 pair index seperated by the inversion if sv exists
	uint64_t *u, max_index_plus = 0, max_index_reverse = 0;
	int min_sc1 = min_sc;
	

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist_x < bw) max_dist_x = bw;
	if (max_dist_y < bw && !is_cdna) max_dist_y = bw;
	if (is_cdna) max_drop = INT32_MAX;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
	// first determine the left and right position that max the increased the score
	score_increased_max = 0;
	//p1 = Kmalloc(km, int64_t, n); // store the kmer where max score comes from, -1 means does not have the previous anchor
	//f1 = Kmalloc(km, int32_t, n); // store the score of the i-th anchor

	// fill the score and backtrack arrays to determine the max increased score
	// for (i = 0; i < n; ++i){
	// 	p1[i] = -1;
	// 	f1[i] = 1;
	// 	left_i = i;
	// 	max_score_i = 0;		
	// 	max_j = -1;
	// 	////// int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
	// 	for (j = i - 1; j >= 0; --j){
	// 		if (a[i].x>>63 != a[j].x>>63){
	// 			break;
	// 		}
	// 		d_read = (int32_t)a[i].y - (int32_t)a[j].y;
	// 		d_ref  = (int32_t)(a[i].x - a[j].x);
	// 		if(d_read <= 0 || d_ref <= 0){
	// 			break;
	// 		}
	// 		//d_minn = d_read <= d_ref ? d_read : d_ref;
	// 		if (d_ref <= qlen){
	// 			left_i = j;
	// 			sc = f1[j] + 1;
	// 			if(sc > max_score_i){
	// 				max_score_i = sc;
	// 				max_j = j;
	// 				f1[i] = sc;
	// 				p1[i] = j;
	// 			}
	// 		}
	// 		else{
	// 			break;
	// 		}
	// 		// strand "+-"[a[i].x>>63]
	// 	}
	// 	if(max_score_i > max_score){
	// 		max_score = max_score_i;
	// 		start_i = i;
	// 		while((p1[start_i] != -1) && (p1[start_i] >= left_i)){
	// 			start_i = p1[start_i];
	// 		}
	// 		score_increased = f1[i] - f1[start_i];
	// 		if(score_increased > score_increased_max){
	// 			score_increased_max = score_increased;
	// 			final_start = start_i;
	// 			final_end   = i;
	// 		}
	// 	}
	// }
	// left_ref  = (int32_t)a[final_start].x - 1.1 * (int32_t)a[final_start].y <= 0 ? 0 : (int32_t)a[final_start].x - 1.1 * (int32_t)a[final_start].y;
	// right_ref = (int32_t)a[final_end].x + 1.1 * (qlen - (int32_t)a[final_end].y);
	//read_main_strand = a[final_start].x >> 63;
	// fill the score and backtrack arrays
	for (i = 0, max_ii = -1; i < n; ++i) {
		int64_t max_j = -1, end_j;
		int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) {
			++st;
		}
		if (i - st > max_iter) {
			st = i - max_iter;
		}
		for (j = i - 1; j >= st; --j) {
			int32_t sc;
			sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			++n_iter;
			//if((int32_t)a[i].x < left_ref || (int32_t)a[i].x > right_ref){
			//	sc = INT32_MIN;
			//}
			//else if((int32_t)a[j].x < left_ref || (int32_t)a[j].x > right_ref){
			//	sc = INT32_MIN;
			//}
			if (sc == INT32_MIN) {
				continue;
			}
			sc += f[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
				if (n_skip > 0) {
					--n_skip;
				}
			} else if (t[j] == (int32_t)i) {
				if (++n_skip > max_skip)
					break;
			}
			if (p[j] >= 0) {
				t[p[j]] = i;
			}
		}
		end_j = j;
		if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			int32_t max = INT32_MIN;
			max_ii = -1;
			for (j = i - 1; j >= st; --j)
				if (max < f[j]) {
					max = f[j], max_ii = j;
				}
		}
		if (max_ii >= 0 && max_ii < end_j) {
			int32_t tmp;
			tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
				max_f = tmp + f[max_ii], max_j = max_ii;
		}
		f[i] = max_f, p[i] = max_j;
		if((a[i].x>>63 == 0) && (max_f > max_plus_score)){
			max_plus_score = max_f;
			max_index_plus = i;
		}
		else if((a[i].x>>63 == 1) && (max_f > max_reverse_score)){
			max_reverse_score = max_f;
			max_index_reverse = i;
		}
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			max_ii = i;
		if (mmax_f < max_f) {
			mmax_f = max_f;
		}
		//if(f[i] > max_score){
		//	max_score = f[i];
		//	final_end = i;
		//}
	}
	if(max_plus_score < max_reverse_score){
		//min_sc1 = max_plus_score;
		//sv_index_end = max_index_plus;
		final_end = max_index_reverse;
	}
	else{
		//min_sc1 = max_reverse_score;
		//sv_index_end = max_index_reverse;
		final_end = max_index_plus;
	}
	//sv_index_start = sv_index_end;
	final_start = final_end;
	//while(p[sv_index_start] > -1) {
	//	sv_index_start = p[sv_index_start];
	//	anchors_num++;
	//}
	while(p[final_start] > - 1){
		final_start = p[final_start];
	}
	left_ref  = (int32_t)a[final_start].x;
	right_ref = (int32_t)a[final_end].x;
	// find the inversion in the loacted region
	max_score = 0;
	for(i = 0; i < n; ++i){
		if(a[i].x<<1>>33 == a[final_end].x<<1>>33 && a[i].x>>63 != a[final_end].x>>63){
			if((int32_t)a[i].x < right_ref && (int32_t)a[i].x > left_ref){
				if(f[i] > max_score){
					max_score = f[i];
					sv_index_end = i;
				}
			}
		}
	}
	sv_index_start = sv_index_end;
	while(p[sv_index_start] > -1) {
		sv_index_start = p[sv_index_start];
	}
	*pos_inv_left  = (int32_t)a[sv_index_start].x;
	*pos_inv_right = (int32_t)a[sv_index_end].x;
	uint32_t sv_existt = 1;
	//if((*pos_inv_right - *pos_inv_left < 30) || (1.0 * anchors_num / (*pos_inv_right - *pos_inv_left) < 0.1)){
	//if(*pos_inv_right - *pos_inv_left < 30 || f[sv_index_end] - f[sv_index_start] < 30){	
	if(*pos_inv_left >= right_ref){
		sv_existt = 0;
	}
	else if(*pos_inv_right <= left_ref){
		sv_existt = 0;
	}
	else if(*pos_inv_left > left_ref && *pos_inv_right <= right_ref){
		sv_existt = 1;
	}
	if(f[sv_index_end] - f[sv_index_start] < 30){
		sv_existt = 0;
	}	
	if (sv_existt == 0){
		left_ref  = (int32_t)a[final_start].x - 1.1 * (int32_t)a[final_start].y <= 0 ? 0 : (int32_t)a[final_start].x - 1.1 * (int32_t)a[final_start].y;
		right_ref = (int32_t)a[final_end].x + 1.1 * (qlen - (int32_t)a[final_end].y);
		for(i = 0; i < n; ++i){
			if(a[i].x<<1>>33 == a[final_end].x<<1>>33){
				if((int32_t)a[i].x < left_ref || (int32_t)a[i].x > right_ref){
					t[i] = 0;
					f[i] = 0;
					p[i] = -1;
				}
			}
			else{
				t[i] = 0;
				f[i] = 0;
				p[i] = -1;
			}
		}
	}
	else{
		start1 = final_start;
		end2   = final_end;	
	/////////////////////////////////////   re-fill the score and backtrack arrays if sv exists
		st = 0;
		n_iter = 0;
		mmax_f = 0;
		for(i = 0; i < n; ++i){
			if ((i >= final_start && i <= final_end) || (i >= sv_index_start && i <= sv_index_end)){
				continue;
			}
			else{
				t[i] = 0;
		 		f[i] = 0;
		 		p[i] = -1;
			}
		}
		// for (i = 0; i < final_start; ++i) {
		//  	t[i] = 0;
		//  	f[i] = 0;
		//  	p[i] = -1;
		// }
		// for (i = final_end + 1; i < n; ++i) {
		//  	t[i] = 0;
		//  	f[i] = 0;
		//  	p[i] = -1;
		// }
		//for (i = 0, max_ii = -1; i < n; ++i) {
		for (i = final_start, max_ii = -1; i <= final_end; ++i) {
			int64_t max_j = -1, end_j;
			int32_t max_f = a[i].y>>32&0xff, n_skip = 0;
			while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist_x)) {
				++st;
			}
			if (i - st > max_iter) {
				st = i - max_iter;
			}
			for (j = i - 1; j >= st; --j) {
				int32_t sc;
				sc = comput_sc(&a[i], &a[j], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
				//++n_iter;
				//if((a[j].x >> 32) != rid){
				//	sc = INT32_MIN;
				//}
				//if((int32_t)a[i].x < left_ref || (int32_t)a[i].x > right_ref){
				//	sc = INT32_MIN;
				//}
				//else if((int32_t)a[j].x < left_ref || (int32_t)a[j].x > right_ref){
				//	sc = INT32_MIN;
				//}
				if((int32_t)a[i].x < *pos_inv_right && (int32_t)a[i].x > *pos_inv_left && (int32_t)a[j].x < *pos_inv_left){
					sc = INT32_MIN;
				}
				else if((int32_t)a[i].x > *pos_inv_right && (int32_t)a[j].x > *pos_inv_left && (int32_t)a[j].x < *pos_inv_right){
					sc = INT32_MIN;
				}
				else if((int32_t)a[i].x > *pos_inv_right && (int32_t)a[j].x < *pos_inv_left){
					sc = INT32_MIN;
				}
				if (sc == INT32_MIN) {
					continue;
				}
				sc += f[j];
				if (sc > max_f) {
					max_f = sc, max_j = j;
					if (n_skip > 0) {
						--n_skip;
					}
				} else if (t[j] == (int32_t)i) {
					if (++n_skip > max_skip)
						break;
				}
				if (p[j] >= 0) {
					t[p[j]] = i;
				}
			}
			//end_j = j;
			// if (max_ii < 0 || a[i].x - a[max_ii].x > (int64_t)max_dist_x) {
			// 	int32_t max = INT32_MIN;
			// 	max_ii = -1;
			// 	for (j = i - 1; j >= st; --j)
			// 		if (max < f[j]) {
			// 			max = f[j], max_ii = j;
			// 		}
			// }
			// if (max_ii >= 0 && max_ii < end_j) {
			// 	int32_t tmp;
			// 	tmp = comput_sc(&a[i], &a[max_ii], max_dist_x, max_dist_y, bw, chn_pen_gap, chn_pen_skip, is_cdna, n_seg);
			// 	if (tmp != INT32_MIN && max_f < tmp + f[max_ii])
			// 		max_f = tmp + f[max_ii], max_j = max_ii;
			// }
			f[i] = max_f, p[i] = max_j;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
			// if (max_ii < 0 || (a[i].x - a[max_ii].x <= (int64_t)max_dist_x && f[max_ii] < f[i]))
			// 	max_ii = i;
			// if (mmax_f < max_f) {
			// 	mmax_f = max_f;
			// }
		}
	}
	// determine the left and right backtracking index based on the inversion
	// int64_t score_left = 0, score_right = 0, main_index_end = -1, main_index_start = -1;
	// for(i = 0; i < n; ++i){
	// 	if((int32_t)a[i].x < left_ref || (int32_t)a[i].x > right_ref){
	// 		continue;
	// 	}
	// 	else if ((int32_t)a[i].x < *pos_inv_left){
	// 		if(f[i] > score_left){
	// 			score_left = f[i];
	// 			main_index_start = i;
	// 		}
	// 	}
	// 	else if ((int32_t)a[i].x > *pos_inv_right){
	// 		if(f[i] > score_right){
	// 			score_right = f[i];
	// 			main_index_end = i;
	// 		}
	// 	}		
	// } 
	//u = mg_chain_backtrack_my(km, n, f, p, v, t, min_cnt, min_sc1, max_drop, &n_u, &n_v, main_index_start, sv_index_end, main_index_end);
	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc1, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); 
	//kfree(km, p1);
	kfree(km, f);
	//kfree(km, f1); 
	kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}

typedef struct lc_elem_s {
	int32_t y;
	int64_t i;
	double pri;
	KRMQ_HEAD(struct lc_elem_s) head;
} lc_elem_t;

#define lc_elem_cmp(a, b) ((a)->y < (b)->y? -1 : (a)->y > (b)->y? 1 : ((a)->i > (b)->i) - ((a)->i < (b)->i))
#define lc_elem_lt2(a, b) ((a)->pri < (b)->pri)
KRMQ_INIT(lc_elem, lc_elem_t, head, lc_elem_cmp, lc_elem_lt2)

KALLOC_POOL_INIT(rmq, lc_elem_t)

static inline int32_t comput_sc_simple(const mm128_t *ai, const mm128_t *aj, float chn_pen_gap, float chn_pen_skip, int32_t *exact, int32_t *width)
{
	int32_t dq = (int32_t)ai->y - (int32_t)aj->y, dr, dd, dg, q_span, sc;
	dr = (int32_t)(ai->x - aj->x);
	*width = dd = dr > dq? dr - dq : dq - dr;
	dg = dr < dq? dr : dq;
	q_span = aj->y>>32&0xff;
	sc = q_span < dg? q_span : dg;
	if (exact) *exact = (dd == 0 && dg <= q_span);
	if (dd || dq > q_span) {
		float lin_pen, log_pen;
		lin_pen = chn_pen_gap * (float)dd + chn_pen_skip * (float)dg;
		log_pen = dd >= 1? mg_log2(dd + 1) : 0.0f; // mg_log2() only works for dd>=2
		sc -= (int)(lin_pen + .5f * log_pen);
	}
	return sc;
}

mm128_t *mg_lchain_rmq(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km)
{
	int32_t *f,*t, *v, n_u, n_v, mmax_f = 0, max_rmq_size = 0, max_drop = bw;
	int64_t *p, i, i0, st = 0, st_inner = 0;
	uint64_t *u;
	lc_elem_t *root = 0, *root_inner = 0;
	void *mem_mp = 0;
	kmp_rmq_t *mp;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist < bw) max_dist = bw;
	if (max_dist_inner <= 0 || max_dist_inner >= max_dist) max_dist_inner = 0;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	mem_mp = km_init2(km, 0x10000);
	mp = kmp_init_rmq(mem_mp);

	// fill the score and backtrack arrays
	for (i = i0 = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span;
		lc_elem_t s, *q, *r, lo, hi;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			int64_t j;
			for (j = i0; j < i; ++j) {
				q = kmp_alloc_rmq(mp);
				q->y = (int32_t)a[j].y, q->i = j, q->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				krmq_insert(lc_elem, &root, q, 0);
				if (max_dist_inner > 0) {
					r = kmp_alloc_rmq(mp);
					*r = *q;
					krmq_insert(lc_elem, &root_inner, r, 0);
				}
			}
			i0 = i;
		}
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || krmq_size(head, root) > cap_rmq_size)) {
			s.y = (int32_t)a[st].y, s.i = st;
			if ((q = krmq_find(lc_elem, root, &s, 0)) != 0) {
				q = krmq_erase(lc_elem, &root, q, 0);
				kmp_free_rmq(mp, q);
			}
			++st;
		}
		if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || krmq_size(head, root_inner) > cap_rmq_size)) {
				s.y = (int32_t)a[st_inner].y, s.i = st_inner;
				if ((q = krmq_find(lc_elem, root_inner, &s, 0)) != 0) {
					q = krmq_erase(lc_elem, &root_inner, q, 0);
					kmp_free_rmq(mp, q);
				}
				++st_inner;
			}
		}
		// RMQ
		lo.i = INT32_MAX, lo.y = (int32_t)a[i].y - max_dist;
		hi.i = 0, hi.y = (int32_t)a[i].y;
		if ((q = krmq_rmq(lc_elem, root, &lo, &hi)) != 0) {
			int32_t sc, exact, width, n_skip = 0;
			int64_t j = q->i;
			assert(q->y >= lo.y && q->y <= hi.y);
			sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, &exact, &width);
			if (width <= bw && sc > max_f) {
				max_f = sc, max_j = j;
			}
			if (!exact && root_inner && (int32_t)a[i].y > 0) {
				lc_elem_t *lo, *hi;
				s.y = (int32_t)a[i].y - 1, s.i = n;
				krmq_interval(lc_elem, root_inner, &s, &lo, &hi);
				if (lo) {
					const lc_elem_t *q;
					int32_t width, n_rmq_iter = 0;
					krmq_itr_t(lc_elem) itr;
					krmq_itr_find(lc_elem, root_inner, lo, &itr);
					while ((q = krmq_at(&itr)) != 0) {
						if (q->y < (int32_t)a[i].y - max_dist_inner) break;
						++n_rmq_iter;
						j = q->i;
						sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
						if (width <= bw) {
							if (sc > max_f) {
								max_f = sc, max_j = j;
								if (n_skip > 0) --n_skip;
							} else if (t[j] == (int32_t)i) {
								if (++n_skip > max_chn_skip)
									break;
							}
							if (p[j] >= 0) t[p[j]] = i;
						}
						if (!krmq_itr_prev(lc_elem, &itr)) break;
					}
				}
			}
		}
		// set max
		assert(max_j < 0 || (a[max_j].x < a[i].x && (int32_t)a[max_j].y < (int32_t)a[i].y));
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (mmax_f < max_f) mmax_f = max_f;
		if (max_rmq_size < krmq_size(head, root)) max_rmq_size = krmq_size(head, root);
	}
	km_destroy(mem_mp);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}

mm128_t *mg_lchain_rmq_my(int max_dist, int max_dist_inner, int bw, int max_chn_skip, int cap_rmq_size, int min_cnt, int min_sc, float chn_pen_gap, float chn_pen_skip,
					   int64_t n, mm128_t *a, int *n_u_, uint64_t **_u, void *km, const int32_t pos_inv_left, const int32_t pos_inv_right)
{
	int32_t *f,*t, *v, n_u, n_v, mmax_f = 0, max_rmq_size = 0, max_drop = bw;
	int64_t *p, i, i0, st = 0, st_inner = 0;
	uint64_t *u;
	lc_elem_t *root = 0, *root_inner = 0;
	void *mem_mp = 0;
	kmp_rmq_t *mp;

	if (_u) *_u = 0, *n_u_ = 0;
	if (n == 0 || a == 0) {
		kfree(km, a);
		return 0;
	}
	if (max_dist < bw) max_dist = bw;
	if (max_dist_inner <= 0 || max_dist_inner >= max_dist) max_dist_inner = 0;
	p = Kmalloc(km, int64_t, n);
	f = Kmalloc(km, int32_t, n);
	t = Kcalloc(km, int32_t, n);
	v = Kmalloc(km, int32_t, n);
	mem_mp = km_init2(km, 0x10000);
	mp = kmp_init_rmq(mem_mp);

	// fill the score and backtrack arrays
	for (i = i0 = 0; i < n; ++i) {
		int64_t max_j = -1;
		int32_t q_span = a[i].y>>32&0xff, max_f = q_span;
		lc_elem_t s, *q, *r, lo, hi;
		// add in-range anchors
		if (i0 < i && a[i0].x != a[i].x) {
			int64_t j;
			for (j = i0; j < i; ++j) {
				q = kmp_alloc_rmq(mp);
				q->y = (int32_t)a[j].y, q->i = j, q->pri = -(f[j] + 0.5 * chn_pen_gap * ((int32_t)a[j].x + (int32_t)a[j].y));
				krmq_insert(lc_elem, &root, q, 0);
				if (max_dist_inner > 0) {
					r = kmp_alloc_rmq(mp);
					*r = *q;
					krmq_insert(lc_elem, &root_inner, r, 0);
				}
			}
			i0 = i;
		}
		// get rid of active chains out of range
		while (st < i && (a[i].x>>32 != a[st].x>>32 || a[i].x > a[st].x + max_dist || krmq_size(head, root) > cap_rmq_size)) {
			s.y = (int32_t)a[st].y, s.i = st;
			if ((q = krmq_find(lc_elem, root, &s, 0)) != 0) {
				q = krmq_erase(lc_elem, &root, q, 0);
				kmp_free_rmq(mp, q);
			}
			++st;
		}
		if (max_dist_inner > 0)  { // similar to the block above, but applied to the inner tree
			while (st_inner < i && (a[i].x>>32 != a[st_inner].x>>32 || a[i].x > a[st_inner].x + max_dist_inner || krmq_size(head, root_inner) > cap_rmq_size)) {
				s.y = (int32_t)a[st_inner].y, s.i = st_inner;
				if ((q = krmq_find(lc_elem, root_inner, &s, 0)) != 0) {
					q = krmq_erase(lc_elem, &root_inner, q, 0);
					kmp_free_rmq(mp, q);
				}
				++st_inner;
			}
		}
		// RMQ
		lo.i = INT32_MAX, lo.y = (int32_t)a[i].y - max_dist;
		hi.i = 0, hi.y = (int32_t)a[i].y;
		if ((q = krmq_rmq(lc_elem, root, &lo, &hi)) != 0) {
			int32_t sc, exact, width, n_skip = 0;
			int64_t j = q->i;
			assert(q->y >= lo.y && q->y <= hi.y);
			sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, &exact, &width);
			if((int32_t)a[i].x < pos_inv_right && (int32_t)a[i].x > pos_inv_left && (int32_t)a[j].x < pos_inv_left){
				sc = 0;
			}
			else if((int32_t)a[i].x > pos_inv_right && (int32_t)a[j].x > pos_inv_left && (int32_t)a[j].x < pos_inv_right){
				sc = 0;
			}
			else if((int32_t)a[i].x > pos_inv_right && (int32_t)a[j].x < pos_inv_left){
				sc = 0;
			}
			if (width <= bw && sc > max_f) {
				max_f = sc, max_j = j;
			}
			if (!exact && root_inner && (int32_t)a[i].y > 0) {
				lc_elem_t *lo, *hi;
				s.y = (int32_t)a[i].y - 1, s.i = n;
				krmq_interval(lc_elem, root_inner, &s, &lo, &hi);
				if (lo) {
					const lc_elem_t *q;
					int32_t width, n_rmq_iter = 0;
					krmq_itr_t(lc_elem) itr;
					krmq_itr_find(lc_elem, root_inner, lo, &itr);
					while ((q = krmq_at(&itr)) != 0) {
						if (q->y < (int32_t)a[i].y - max_dist_inner) break;
						++n_rmq_iter;
						j = q->i;
						sc = f[j] + comput_sc_simple(&a[i], &a[j], chn_pen_gap, chn_pen_skip, 0, &width);
						if((int32_t)a[i].x < pos_inv_right && (int32_t)a[i].x > pos_inv_left && (int32_t)a[j].x < pos_inv_left){
							sc = 0;
						}
						else if((int32_t)a[i].x > pos_inv_right && (int32_t)a[j].x > pos_inv_left && (int32_t)a[j].x < pos_inv_right){
							sc = 0;
						}
						else if((int32_t)a[i].x > pos_inv_right && (int32_t)a[j].x < pos_inv_left){
							sc = 0;
						}
						if (width <= bw) {
							if (sc > max_f) {
								max_f = sc, max_j = j;
								if (n_skip > 0) --n_skip;
							} else if (t[j] == (int32_t)i) {
								if (++n_skip > max_chn_skip)
									break;
							}
							if (p[j] >= 0) t[p[j]] = i;
						}
						if (!krmq_itr_prev(lc_elem, &itr)) break;
					}
				}
			}
		}
		// set max
		assert(max_j < 0 || (a[max_j].x < a[i].x && (int32_t)a[max_j].y < (int32_t)a[i].y));
		f[i] = max_f, p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		if (mmax_f < max_f) mmax_f = max_f;
		if (max_rmq_size < krmq_size(head, root)) max_rmq_size = krmq_size(head, root);
	}
	km_destroy(mem_mp);

	u = mg_chain_backtrack(km, n, f, p, v, t, min_cnt, min_sc, max_drop, &n_u, &n_v);
	*n_u_ = n_u, *_u = u; // NB: note that u[] may not be sorted by score here
	kfree(km, p); kfree(km, f); kfree(km, t);
	if (n_u == 0) {
		kfree(km, a); kfree(km, v);
		return 0;
	}
	return compact_a(km, n_u, u, n_v, v, a);
}
