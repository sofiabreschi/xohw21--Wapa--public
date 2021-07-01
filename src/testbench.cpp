#include "testbench_defines.h"
#define DIMREF 475
#define DIMSEQ 400


void reverseStr(char str[],const int n){
    for (int i = 0; i < n / 2; i++)
        swap(str[i], str[n - i - 1]);
}

char random_char(){
    char alphabet[4] = {'A', 'C', 'G', 'T'};
    return alphabet[rand() % 4];
}

int random_int(){
    int number[10] = {20, 34, 57, 68, 72, 12, 9, 28, 43, 62};
    return number[rand() % 10];
}

char* init_array(char*m,const int DIM){
    m=new char[DIM];
    if(m==NULL) exit(1);
    return m;
}

int** init_mat(int**m,const int DIM){
    m= new int*[DIM];
    for(int i=0;i<DIM;i++) m[i]= new int[DIM];
    return m;
}

char** init_mat_char(char**m,const int DIM){
    m= new char*[DIM];
    for(int i=0;i<DIM;i++) m[i]= new char[DIM];
    return m;
}


void print1(const char s[], const int DIM){
    cout<<endl<<"Reference: ";
    for(int i=0;i<DIM;i++) cout<<s[i];
    //cout<<endl;
}


void print2(const char s[], const int DIM){
    cout<<"\nSequence:  ";
    for(int i=0;i<DIM;i++) cout<<s[i];
    cout<<endl;
}

void free_string(string*s){
    delete[] s;
}

void free_array(char*s){
    delete [] s;
}

void free_mat(int**s, const int DIM){
    for(int i=0;i<DIM;i++) delete [] s[i];
}

void free_mat_char(char**s, const int DIM){
    for(int i=0; i<DIM; i++) delete[] s[i];
}

/*

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;


static const bntseq_t *global_bns = 0; // for debugging only

 mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = (mem_opt_t *)calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}


typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	smem_aux_t **aux;
	bseq1_t *seqs;
	mem_alnreg_v *regs;
	int64_t n_processed;
} worker_t;




typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int score;
} mem_seed_t; // unaligned memory

typedef struct {
	int n, m, first, rid;
	uint32_t w:29, kept:2, is_alt:1;
	float frac_rep;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain_t;




static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

#define MAX_BAND_TRY  2

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

	if (c->n == 0) return;
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	assert(c->rid == rid);

	srt = (uint64_t*)malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];


		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = (uint8_t*)malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = (uint8_t*)malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(srt); free(rseq);
}


*/



int main(int argc, char *argv[]){

    char ref[DIMREF];
    char seq[DIMSEQ];
    int seed_length=0;
    int start_r=0;
    int start_s=0;
    int end_r=0;
    int end_s=0;


    srand((unsigned)time(NULL));

    //genera stringa e sequenza
    for (int i = 0; i < DIMREF; i++) {
        ref[i]= random_char();
    }
    for (int i = 0; i < DIMSEQ; i++) {
        seq[i] = random_char();
    }

    print1(ref, DIMREF);
    print2(seq, DIMSEQ);

    //genera casualmente la lunghezza e gli indici del seed
    seed_length=random_int();
    cout<<"\nseed length: "<<seed_length<<endl;
    start_r=random_int();
    start_s=random_int();
    end_r=start_r+seed_length;
    end_s=start_s+seed_length;

    char ref_dx[DIMREF-end_r];

   // cout<<"Dimensione stringa"<<DIMREF-seed_length-end_r<<"  ";//avra' al massimo questa dimensione
    char seq_dx[DIMSEQ-end_s];


    //parte dx
    for (int i=0; i<DIMREF-end_r; i++) ref_dx[i]=ref[i+end_r];

    for (int i=0; i<DIMSEQ-end_s; i++) seq_dx[i]=seq[i+end_s];

    cout<<"\nstringa e sequenza DX:";
    print1(ref_dx, DIMREF-end_r);
    print2(seq_dx, DIMREF-end_s);


    char ref_sx[start_r];
    char seq_sx[start_s];

    // parte sx
    for (int i=0; i<start_r; i++) ref_sx[i]=ref[i];
    for (int i=0; i<start_s; i++) seq_sx[i]=seq[i];

    reverseStr(ref_sx, start_r);
    reverseStr(seq_sx, start_s);

    cout<<"\nstringa e sequenza SX girate:";
    print1(ref_sx, start_r);
    print2(seq_sx, start_s);

    int bandwidth = 500;
    int res_dx, res_sx;
    // parte dx
//    cout <<  (DIMREF-end_r)+1 << endl;
    cout<<"\nseed extension DX:"<<endl;
    ap_uint<512> *ref_512, *seq_512;
    ref_512 = (ap_uint<512> *)ref_dx;
    seq_512 = (ap_uint<512> *)seq_dx;
    SW(ref_512, seq_512,DIMSEQ-end_s+1 ,DIMREF-end_r+1,bandwidth ,seed_length, &res_dx);

    ref_512 = (ap_uint<512> *)ref_sx;
    seq_512 = (ap_uint<512> *)seq_sx;
    // parte sx
    cout<<"\nseed extension SX:"<<endl;
  //  SW(ref_512, seq_512,start_s+1 ,start_r+1, bandwidth,seed_length, &res_sx);

/*
    // test funzioni repo pubblica
    *mem_opt_init();
    int l_query = DIMSEQ;
    worker_t worker;
    mem_seed_t seed;
    mem_chain_t chain;
    //mem_chain2aln( *opt,*bns, *pac,l_query, *query,  *c, *av);  //gli input non sono corretti, da cercare in .h come sono inizializzati
*/
    return 0;

}

