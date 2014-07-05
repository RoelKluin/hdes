
#ifndef RK_B6_HPP
# define RK_B6_HPP
# define BITS_PER_LONG 64

# define B6_MASK			0x6

# define B6_UC				0x40
# define B6_LC				0x60

# define B6_DNA				0x10
# define B6_RNA				0x11

# define B6_A2B				0x75
# define B6_B2A				0x35

# define B6_ALT_CASE			0x20

# define ischar_ignore_cs(c, c2) (((c ^ c2) & ~B6_ALT_CASE) == 0)

/* convert ascii uc/lc DNA/RNA to 2bit format (on bits 2 and 3 -> 0x6) or
 * back to ascii the conversion is reversible, even for invalid characters
 *    twobit	+- deoxy -+	+-- oxy -+	[ribose: no_u]
 * hex	binary	uc	lc	uc	lc	[case: cs]
 * 0x00	0000	A	a	A	a
 * 0x02	0010	C	c	C	c
 * 0x04	0100	T	t	U	u
 * 0x06	0110	G	G	G	g
 *    b2a   <--->   a2b [conversion: flips 0x40 [TDE] / 0x41 [U] /0x1 (others)]
 * in two successive calls 'U' converts to b6 (0x2) and to 'T'.
 */
static inline unsigned b6(unsigned c)
{	/* don't expect first two to change very often, so may move them out */
	unsigned no_u = !ischar_ignore_cs(c, 'U');
	unsigned cs = (c | B6_UC) & B6_LC;
	unsigned TUDE_conv = (c & 0x8e) == 0x4; // conversion if qr/[TUDE]/
	return c ^ cs ^ (-TUDE_conv & B6_RNA) ^ no_u;
}

/* if uc DNA is expected use with cs = B6_UC, no_u = 1 */
static inline unsigned b6_spec(unsigned c, unsigned cs, unsigned no_u)
{
	unsigned TUDE4_conv = (c & 0x8e) == 0x4; // if qr/[TUDE4??]/
	return c ^ cs ^ (-TUDE4_conv & B6_RNA) ^ no_u;
}

/* With the conversion above only the specified characters are converted to
 * 2bits, left shifted by one, with all other bits zeroed. With the function
 * below one can test whether the ascii character prior to conversion was a
 * valid convertable character.
 */
#define isb6(b6) (((b6) | B6_MASK) == B6_MASK)

/* A similar, not entirely reversible conversion of ascii to b6. It ignores case
 * and works for both DNA and RNA, but all characters in the range [@-GP-W`-gp-w]
 * produce values interpreted as b6, including 'R', which occurs on hg19, chr3.
 */
# define __qb6(c) (c ^ ((c & (B6_ALT_CASE|B6_RNA)) | (B6_A2B ^ B6_B2A)))

static inline uint32_t revseq(uint32_t x)
{
    uint32_t m = 0x33333333;
    x = ((x & m) << 2) | ((x & ~m) >> 2);
    m = 0x0f0f0f0f;
    x = ((x & m) << 4) | ((x & ~m) >> 4);
    m = 0x00ff00ff;
    x = ((x & m) << 8) | ((x & ~m) >> 8);
    m = 0x0000ffff;
    return ((x & m) << 16) | ((x & (m << 16)) >> 16);
}

unsigned get_twisted(unsigned rev, unsigned sq)
{
    rev = rev ^ 0xaaaaaaaau;
    unsigned dev = sq ^ rev;

    unsigned m = dev & -dev; // leftmost deviant bit. NB: palindromes have none
    sq &= m;                 // test whether it's set for first sequence
    rev ^= (sq ^ -sq) & dev; // flip revcmp deviant above devbit, if set

    m -= !!m; // mask below devbit, if any - none for palindromes
    m &= rev; // revnxt below devbit, if any
    rev ^= m ^ (m << 1) ^ sq ^ !sq; // swap devbit to lowest and flip
    return ((rev & 0xffff) << 16) | (dev & 0xffff);
}

#endif // RK_B6_H


