#define ONLY_FIELDS
/* if set, assume only fields are to be generated */

#if 0
#define COUNT_MULT
#endif
/* if set, count discriminant (or cluster) multiplicities, else count # of
 * fields */

#if 0
#  define CHECK_CLUSTER /* if set check a cluster */
#  define SHIFT 7       /* ... of 2^SHIFT discriminants */
#endif

#if 0
#  define ONLY_UNRAM 
/* if set, generate unramified fields (over the quadratic subfield) only */

#  define FULL_CHECK
/* No effect if ONLY_UNRAM is not set. Otherwise, if FULL_CHECK is set, check
 * for fields rigorously. If not, use heuristics instead of full check (> 98%
 * success). More precisely, IF a form output has fundamental discriminant,
 * it corresponds to a field unramified over its quadratic subfied;
 * otherwise, it's junk. */
#endif
