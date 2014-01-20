/*
 * =====================================================================================
 *
 *       Filename:  multiset.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  15-01-14 20:31:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  (c) Roel Kluin 2013
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdint.h>
#include <assert.h>
#include <iostream>
#include <stack>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <iomanip>
#define SIZEOF_ARRAY(a)  (sizeof(a) / sizeof(a[0]))
using namespace std;

//#define IN_PRODUCTION
#define ITERATIONS 1000000
#define ELEMENTS 50

#if not defined(DEBUG) || defined(IN_PRODUCTION)
# define DEBUG 1
#endif

#if DEBUG > 1
struct dbgr{template<typename T> dbgr& operator ,(const T& v){std::cerr<<v<<" ";return *this;}}dbg;
#define debug(...) {dbg,__VA_ARGS__;std::cerr<<std::endl;}

#define DEC(x) std::dec << #x << ": " << (x)
#define HEX(x) std::hex << #x << ": " << (x)

#define CERR(x) do { std::cerr << x << std::endl; } while (0)
#else
#define DEC(x)
#define HEX(x)
#define CERR(x)
#define debug
#endif

/* init seed with a real random integer (from /dev/random) */
static void seed_rnd()
{
	uint8_t buf[4] = {0};
	int j = open("/dev/random", O_RDONLY);
	j = read(j, buf, 4);
	long int seedi = 0l;
	for(j = 0; j != 4; ++j) {
		seedi |= buf[j];
		seedi <<= 8;
	}
	srandom(seedi);
}

template<unsigned size>
class Ring;

/* e.g. in a window of 50 bp, below is shoen the number of deviants <n> that
 * can be calculated at each offset. Which is at most a deviant size of 6:
 *
 * [0] 3   7   11  15  19  23  <5> <4> <3> <2> <1>
 *  ---+---+---+---+---+---+--+---+---+---+---+---+---
 *  <0> <1> <2> <3> <4> <5> ^ 26  30  34  38  42  46 [50]
 *                         <6(2)>
 * 
 */
#define MAX_HDEV_SIZE(z) ((z >> 3) - !!(z & 3))
#define B2S_PER_QUAD 4


#define B2(x)   (   (x) | (   (x) >> 1) )
#define B4(x)   ( B2(x) | ( B2(x) >> 2) )
#define B8(x)   ( B4(x) | ( B4(x) >> 4) )
#define B16(x)  ( B8(x) | ( B8(x) >> 8) )
#define B32(x)  (B16(x) | (B16(x) >>16) )
#define NEXT_POWER_OF_2(x)      (B32(x-1)+1)

/* TODO: dynamic sized container of Deviants for assembly? */

/* TODO: implement adapter clipping? for paired-end: paired deviants? */


template<unsigned size>
class Deviant
{
	uint32_t val;
	char* hdev, *seq, *rc;
public:
	Deviant() : is_returned(false) {
		//static_assert(size > 4, "Error: Deviant template size <= 4!");
		assert(size > B2S_PER_QUAD);
		/* the maximum hdev length is 1/4 readlength */
		hdev = new char[MAX_HDEV_SIZE(size)]; // 4 deviants per char.
	}
	~Deviant() { delete[] hdev; }
	void update(void* s)
	{
		seq = (char*)s;/* set pointer to right quad */
		rc = (char*)--s;
		/* XXX: early quads need abbreviation marking, later need
		 * extension.
		 * XXX: deviants accross breakpoint */
	}
	void init(void* s, unsigned ct)
	{
		ct >>= 2; /* now count per quad */
		seq = (char*)s;/* set pointer to right quad */
		rc = (char*)--s;
		/*offset = ct & 3; // offset within quad, XXX: but do we care? */


		// were some deviants completed? create hdev part 1.
		for (unsigned i = 0; i != ct; ++i)
			hdev[i] = *seq ^ *rc;
	}
	uint32_t get_val() { return val; }

	Deviant* cmp(Deviant* other)
	{
		/* for deviants there will be an extension if possible in the
		 * equal case, until we're certain that one of the two is less
		 */
		if (val < other->get_val()) {
			min = this;
		} else {
			min = other;
		}
		return min;
	}
	Deviant* min;
	/* unsigned offset: 2; // offset XXX: can go? */
	unsigned is_returned: 1;
	/* yes this works: template size based bitfield size */
	unsigned index: NEXT_POWER_OF_2(MAX_HDEV_SIZE(size)); // see MAX_HDEV_SIZE
	friend class Ring<size>;
};

/**
 * see  http://stackoverflow.com/questions/4802038/implement-a-queue-in-which-push-rear-pop-front-and-get-min-are-all-consta
 * To find the minimum element of the queue, look at the smallest two elements
 * of the subarrays, and pick their minimum. using single array instead of two
 * stacks.
 */
template<unsigned size>
class Ring {
	Deviant<size> *chain, *at, *in_min, *out_min;
public:
	Ring()
	{
		//static_assert(size > 0u, "Error: Ring template size == 0!");
		assert(size > 0);
		chain = new Deviant<size>[size];
		in_min = at = &chain[0];
		at->min = at;
	}
	~Ring(){ delete[] chain; }

	void init(void* data)
	{
		at->init(data, at - &chain[0]);
		in_min = at++->cmp(in_min);
		/* caller must init from 0 to size - 1, the last is the first
		 * insertion. */
		assert(at < &chain[size]);
	}

	Deviant<size>* turn(void* s)
	{
		assert(at != NULL); assert(at < &chain[size]); assert(at >= &chain[0]);
		Deviant<size>* ret;

		at->update(s);
		at->is_returned = false;
		in_min = at->cmp(in_min);

		CERR(""<< size - (&chain[size] - at) <<"/"<< size <<"]: pushed (" <<
				HEX(at->val) <<", "<< HEX(at->min->get_val()) <<")");

		if (at != &chain[size - 1]) {
			if (at++ == out_min) {

				CERR(HEX(out_min->get_val()) << " => " << HEX(at->min->get_val()));
				out_min = at->min;
			}
			ret = in_min->cmp(out_min);
		} else {
			/* At end of array order all elements in reverse. */

			at->min = out_min = at;
			do {
				out_min = (--at)->cmp(out_min);

			} while (at != &chain[0]);

			CERR("\t*** SWAP ***\n" << HEX(out_min->get_val()));

			ret = out_min;
			in_min = at;
		}
		/* XXX: can determine offset from deviant here, using
		 * pointers. something like this? to be used in assembly.
		 */
		//ret->offset = (ret - ret->min) % size;

#if DEBUG > 1
		print();
#endif
#if defined(DEBUG)
		unsigned i, found = 0;
		Deviant<size>* p;
		for (p = &chain[0]; p != &chain[size]; ++p) {
			if (p == ret)
				found = 1;
			else if (p->val < ret->val) {
				CERR(HEX(p->get_val()) <<" is less!! then supposed minimum:"<<
					HEX(ret->get_val()));
				found = 0;
				break;
			}

		}
		if (found == 0) {
			if (p == &chain[size])
				CERR("\t*** WARNING *** No "<< HEX(ret->get_val())
					<<" in these: ");
			CERR("at ["<< size - (&chain[size] - at) <<"/"<< size <<"]:"<<
				"ASSERTION FAILED\n\nTotal min: "<< ret);
			CERR("popped ("<< at->get_val() <<", "<< at->min->get_val() <<")");
			return NULL;
		}
		CERR("popped: ("<< ret->get_val() <<","<< ret->min <<")");
#endif
		return ret;
	}
#if defined(DEBUG)
	void print() {
		Deviant<size>* p;
		cerr << "fst: ";
		for (p = &chain[0]; p != at; ++p)
			cerr << " (" << p->get_val() <<","<< p->min->get_val() <<")";

		cerr << endl << "snd: ";
		for (p = &chain[size - 1]; p >= at; --p)
			cerr << " (" << p->get_val() <<","<< p->min->get_val() <<")";

		cerr << endl;
	}
#endif
};

int main()
{
	seed_rnd();
	Ring<ELEMENTS> ring();
	/* four rotations for the 2bits in two complements */
	char seq[(MAX_HDEV_SIZE(ELEMENTS) * B2S_PER_QUAD) * 2] = {0};
	char *s = seq;

	unsigned i = 0, r = 0;
	char b2, rcb2;
	do {
		assert(s < &seq[(MAX_HDEV_SIZE(ELEMENTS) * B2S_PER_QUAD) * 2]);
		if ((i & 8) == 0) {
			r >>= 2;

			if ((i & 32) == 0) /* simulate 2bit receival */
				r = random();

			b2 = r & 0x3;
			rcb2 = (b2 << 6) ^ 0x2; // reverse complement b2
			s -= (4 * MAX_HDEV_SIZE(ELEMENTS)) - 1; /* reset and increment */
		}

		/* insert 2bit for all rotations and complements */
		*s |= b2;
		*++s |= rcb2;
		b2 <<= 2;
		rcb2 >>= 2;
		ring.init(s);
		s += MAX_HDEV_SIZE(ELEMENTS) * 2; // move to next rotation

	} while (++i != ELEMENTS - 1); /* leave last for replace */

	unsigned uniq = 0u;

	/* FIXME: the sequence is rotated once 50 bp are read;
	 * the deviants must be calculated across the breakpoints */

	unsigned frame = i;
	char m[4] = {0x3, 0xc, 0x30, 0xc0};
	do {
		if (s == &seq[(MAX_HDEV_SIZE(ELEMENTS) * B2S_PER_QUAD) * 2]) {
			s = seq;
			frame = 0;
		}
		assert(s < &seq[(MAX_HDEV_SIZE(ELEMENTS) * B2S_PER_QUAD) * 2]);

		CERR("");
		if ((i & 8) == 0) {
			r >>= 2;

			if ((i & 32) == 0) /* simulate 2bit receival */
				r = random();

			b2 = r & 0x3;
			rcb2 = (b2 << 6) ^ 0x2; // reverse complement b2
			s -= (4 * MAX_HDEV_SIZE(ELEMENTS)) - 1; /* reset and increment */
		}

		/* replace 2bit for all rotations */

		*s ^= (*s & m[frame & 3]) ^ b2;
		*++s ^= (*s & m[(3 - frame) & 3]) ^ rcb2; // reverse complement b2
		b2 <<= 2;
		rcb2 >>= 2;

		Deviant<ELEMENTS>* b = ring.turn(s);
		if (b == NULL) break; // debug test failed

		if (b->is_returned == false) { /* or was already observed as minimum */
			++uniq;
			b->is_returned = true; // after test (for assert)
			CERR("-----[new min: " << b->min->get_val() <<" ]------");
		}
		s += MAX_HDEV_SIZE(ELEMENTS) * 2;
		++frame;

	} while (++i < ITERATIONS);

	if (i == ITERATIONS) {
		cout <<setiosflags(ios::fixed) <<setprecision(2)<<
		100.0 * uniq / ITERATIONS << "% of the returns were unique.\n";
	} else {
		cerr << "Aborted." << endl;
	}

	return 0;
}


