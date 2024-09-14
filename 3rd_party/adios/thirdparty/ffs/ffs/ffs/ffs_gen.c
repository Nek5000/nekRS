#include "config.h"
#include "stdio.h"
#include "stdlib.h"
#include "ffs.h"

#ifdef DO_DCG
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#include "assert.h"
#include "fm_internal.h"
#include "ffs_internal.h"
#include "dill.h"
#define static_ctx c 
#define _vrr(x) x
#include "ffs_gen.h"

/* #define REG_DEBUG(x) printf x ; */
#define REG_DEBUG(x)
#define gen_fatal(str) do {fprintf(stderr, "%s\n", str); exit(0);} while (0)

iogen_oprnd
gen_operand(dill_reg src_reg, size_t offset, int size, FMdata_type data_type, int aligned, int byte_swap)
{
    iogen_oprnd ret_val;
    ret_val.address = 1;
    ret_val.size = size;
    ret_val.data_type = data_type;
    ret_val.offset = offset;
    ret_val.aligned = aligned;
    ret_val.byte_swap = byte_swap;
    ret_val.vc_reg = src_reg;
    ret_val.vc_reg2 = src_reg;
    return ret_val;
}

void
gen_load(dill_stream c, iogen_oprnd_ptr src_oprnd)
{
    iogen_oprnd tmp_val;
    tmp_val = gen_fetch(c, 
			src_oprnd->vc_reg, src_oprnd->offset, src_oprnd->size,
			src_oprnd->data_type, src_oprnd->aligned,
			src_oprnd->byte_swap);
    *src_oprnd = tmp_val;
}

iogen_oprnd
gen_bswap_fetch(dill_stream c, dill_reg src_reg, size_t offset,
		int size, FMdata_type data_type, int aligned)
{
    iogen_oprnd ret_val;
    ret_val.address = 0;
    ret_val.size = size;
    ret_val.data_type = data_type;
    ret_val.offset = 0;
    ret_val.aligned = 0;
    ret_val.byte_swap = 0;
    switch (data_type) {
    case unknown_type:
	assert(FALSE);
	break;
    case integer_type:
    case boolean_type:
    case enumeration_type:
	switch (size) {
	case 1:		/* sizeof char */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_C, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsci(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_S, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbssi(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP))
		gen_fatal("gen fetch out of registers A\n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsii(c, ret_val.vc_reg, src_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_L, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    if (!aligned || ((offset & 0x7) != 0)) {
		/* misaligned */
		if ((offset & 0x3) == 0) {
		    /* 4 byte aligned */
		    dill_reg high_reg;
		    if (!ffs_getreg(c, &high_reg, DILL_L, DILL_TEMP))
			gen_fatal("gen fetch out of registers \n");
		    if (WORDS_BIGENDIAN) {
			/* vc_reg2 holds high value */
			dill_ldbsui(c, ret_val.vc_reg, src_reg, offset);
			dill_ldbsii(c, high_reg, src_reg, offset + 4);
		    } else {
			dill_ldbsii(c, high_reg, src_reg, offset);
			dill_ldbsui(c, ret_val.vc_reg, src_reg, offset + 4);
		    }
		    dill_lshli(c, high_reg, high_reg, 32);
		    dill_orl(c, ret_val.vc_reg, high_reg, ret_val.vc_reg);
		    ffs_putreg(c, high_reg, DILL_L);
		} else {
		    assert(FALSE);
		}
	    } else {
		dill_ldbsli(c, ret_val.vc_reg, src_reg, offset);
	    }
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP) ||
		!ffs_getreg(c, &ret_val.vc_reg2, DILL_I, DILL_TEMP))
		gen_fatal("gen fetch out of registers B \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg2)));

	    if (WORDS_BIGENDIAN) {
		/* vc_reg2 holds high value */
		dill_ldbsii(c, ret_val.vc_reg, src_reg, offset);
		dill_ldbsii(c, ret_val.vc_reg2, src_reg, offset + 4);
	    } else {
		dill_ldbsii(c, ret_val.vc_reg2, src_reg, offset);
		dill_ldbsii(c, ret_val.vc_reg, src_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case unsigned_type:
    case char_type:
	switch (size) {
	case 1:		/* sizeof char */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_UC, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsuci(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_US, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsusi(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_U, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsui(c, ret_val.vc_reg, src_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_UL, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldbsuli(c, ret_val.vc_reg, src_reg, offset);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    if (!ffs_getreg(c, &ret_val.vc_reg2, DILL_U, DILL_TEMP) ||
		!ffs_getreg(c, &ret_val.vc_reg, DILL_U, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg2)));
	    /* vc_reg2 holds high value */
	    if (WORDS_BIGENDIAN) {
		dill_ldbsui(c, ret_val.vc_reg, src_reg, offset);
		dill_ldbsui(c, ret_val.vc_reg2, src_reg, offset + 4);
	    } else {
		dill_ldbsui(c, ret_val.vc_reg2, src_reg, offset);
		dill_ldbsui(c, ret_val.vc_reg, src_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case float_type:
    case string_type:
	assert(FALSE);
    }
    return ret_val;
}

iogen_oprnd
gen_set(dill_stream c, int size, char* value)
{
    iogen_oprnd ret_val;
    ret_val.address = 0;
    ret_val.size = size;
    ret_val.data_type = integer_type;
    ret_val.offset = 0;
    ret_val.aligned = 0;
    ret_val.byte_swap = 0;

    switch (size) {
    case 1:		/* sizeof char */
	if (!ffs_getreg(c, &ret_val.vc_reg, DILL_C, DILL_TEMP))
	    gen_fatal("gen fetch out of registers \n");
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	dill_setc(c, ret_val.vc_reg, value[0]);
	break;
    case 2:		/* sizeof short */
	if (!ffs_getreg(c, &ret_val.vc_reg, DILL_S, DILL_TEMP))
	    gen_fatal("gen fetch out of registers \n");
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	dill_sets(c, ret_val.vc_reg, *((short*)value));
	break;
    case 4:		/* sizeof int */
	if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP))
	    gen_fatal("gen fetch out of registers C\n");
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	dill_seti(c, ret_val.vc_reg, *((int*)value));
	break;
#if SIZEOF_SIZE_T != 4
    case SIZEOF_SIZE_T:
	if (!ffs_getreg(c, &ret_val.vc_reg, DILL_L, DILL_TEMP))
	    gen_fatal("gen fetch out of registers \n");
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	dill_setl(c, ret_val.vc_reg, *((long*)value));
	break;
#else
    case 8:
	/* simulate with double reg */
	if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP) ||
	    !ffs_getreg(c, &ret_val.vc_reg2, DILL_I, DILL_TEMP))
	    gen_fatal("gen fetch out of registers D \n");
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg2)));
	/* vc_reg2 holds high value */
	if (WORDS_BIGENDIAN) {
	    dill_seti(c, ret_val.vc_reg2, *((long*)value));
	    dill_seti(c, ret_val.vc_reg, *((long*)(value + 4)));
	} else {
	    dill_seti(c, ret_val.vc_reg, *((long*)value));
	    dill_seti(c, ret_val.vc_reg2, *((long*)(value + 4)));
	}
	break;
#endif
    }
    return ret_val;
}
   
iogen_oprnd
gen_fetch(dill_stream c, dill_reg src_reg, size_t offset, size_t size,
	  FMdata_type data_type, int aligned, int byte_swap)
{
    iogen_oprnd ret_val;

#if defined(v_ldbsi) || defined(DO_DCG)
    if (dill_has_ldbs(c)) {
	/* have byte swap load extension */
	if (byte_swap && (data_type != float_type)) {
	    return gen_bswap_fetch(c, src_reg, offset, (int)size, data_type, 
				   aligned);
	}
    }
#endif
    ret_val.address = 0;
    ret_val.size = (int)size;
    ret_val.data_type = data_type;
    ret_val.offset = 0;
    ret_val.aligned = 0;
    ret_val.byte_swap = 0;
    switch (data_type) {
    case unknown_type:
	assert(FALSE);
	break;
    case integer_type:
    case boolean_type:
    case enumeration_type:
	switch (size) {
	case 1:		/* sizeof char */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_C, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldci(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_S, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldsi(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP))
		gen_fatal("gen fetch out of registers C\n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldii(c, ret_val.vc_reg, src_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_L, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldli(c, ret_val.vc_reg, src_reg, offset);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_I, DILL_TEMP) ||
		!ffs_getreg(c, &ret_val.vc_reg2, DILL_I, DILL_TEMP))
		gen_fatal("gen fetch out of registers D \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg2)));
	    /* vc_reg2 holds high value */
	    if (WORDS_BIGENDIAN) {
		dill_ldii(c, ret_val.vc_reg2, src_reg, offset);
		dill_ldii(c, ret_val.vc_reg, src_reg, offset + 4);
	    } else {
		dill_ldii(c, ret_val.vc_reg, src_reg, offset);
		dill_ldii(c, ret_val.vc_reg2, src_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case unsigned_type:
    case char_type:
	switch (size) {
	case 1:		/* sizeof char */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_UC, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_lduci(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_US, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldusi(c, ret_val.vc_reg, src_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_U, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_ldui(c, ret_val.vc_reg, src_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    if (!ffs_getreg(c, &ret_val.vc_reg, DILL_UL, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    dill_lduli(c, ret_val.vc_reg, src_reg, offset);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    if (!ffs_getreg(c, &ret_val.vc_reg2, DILL_U, DILL_TEMP) ||
		!ffs_getreg(c, &ret_val.vc_reg, DILL_U, DILL_TEMP))
		gen_fatal("gen fetch out of registers \n");
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
	    REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg2)));
	    /* vc_reg2 holds high value */
	    if (WORDS_BIGENDIAN) {
		dill_ldui(c, ret_val.vc_reg2, src_reg, offset);
		dill_ldui(c, ret_val.vc_reg, src_reg, offset + 4);
	    } else {
		dill_ldui(c, ret_val.vc_reg, src_reg, offset);
		dill_ldui(c, ret_val.vc_reg2, src_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case float_type:
	if (byte_swap) {
	    /* best to byte_swap floating point things as we're loading */
	    ret_val.offset = offset;
	    ret_val.aligned = aligned;
	    ret_val.address = TRUE;
	    ret_val.vc_reg = src_reg;
	    ret_val.byte_swap = TRUE;
	    gen_byte_swap(c, &ret_val);
	    byte_swap = FALSE;	/* taken care of */
	} else {
	    switch (size) {
	    case SIZEOF_FLOAT:	/* sizeof char */
		ffs_getreg(c, &ret_val.vc_reg, DILL_F, DILL_TEMP);
		REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
		dill_ldfi(c, ret_val.vc_reg, src_reg, offset);
		break;
	    case SIZEOF_DOUBLE:	/* sizeof short */
		ffs_getreg(c, &ret_val.vc_reg, DILL_D, DILL_TEMP);
		REG_DEBUG(("get %d in gen_Fetch\n", _vrr(ret_val.vc_reg)));
		dill_lddi(c, ret_val.vc_reg, src_reg, offset);
		break;
	    }
	}
	break;
    case string_type:
	assert(FALSE);
    }
    if (byte_swap) {
	gen_byte_swap(c, &ret_val);
    }
    return ret_val;
}

void
gen_byte_swap(dill_stream c, iogen_oprnd_ptr src_oprnd)
{
    iogen_oprnd swap_oprnd;
    if (src_oprnd->address) {
	if (src_oprnd->data_type == float_type) {
	    swap_oprnd = *src_oprnd;
	    swap_oprnd.data_type = integer_type;
	    gen_load(c, &swap_oprnd);
	    /* gen_load did everything for us */
	    goto do_float_load_store;
	} else {
	    swap_oprnd = *src_oprnd;
	    gen_load(c, &swap_oprnd);
#if defined(v_ldbsi) || defined(DO_DCG)
	    if (dill_has_ldbs(c)) {
		/* gen_load did everything for us */
		return;
	    }
#endif
	}
    } else {
	assert(src_oprnd->data_type != float_type);
	swap_oprnd = *src_oprnd;
    }
    /* no byte swap load extension, do it manually */
    switch (swap_oprnd.size) {
    case 1:
	break;
    case 2:
	dill_bswaps(c, swap_oprnd.vc_reg, swap_oprnd.vc_reg);
	break;
    case 4:
	dill_bswapi(c, swap_oprnd.vc_reg, swap_oprnd.vc_reg);
	break;
    case 8:
	if (sizeof(size_t) == 4) {
	    /* swap top and bottom */
	    dill_reg tmp_reg = swap_oprnd.vc_reg;
	    swap_oprnd.vc_reg = swap_oprnd.vc_reg2;
	    swap_oprnd.vc_reg2 = tmp_reg;
	    /* then byte swap each */
	    dill_bswapi(c, swap_oprnd.vc_reg, swap_oprnd.vc_reg);
	    dill_bswapi(c, swap_oprnd.vc_reg2, swap_oprnd.vc_reg2);
	} else {
	    dill_bswapl(c, swap_oprnd.vc_reg, swap_oprnd.vc_reg);
	}
	break;
    default:
	assert(FALSE);
    }
 do_float_load_store:
    if (src_oprnd->address) {
	if (src_oprnd->data_type == float_type) {
	    /* lose type info by storing to memory and retrieving */
	    int tmp_base = ffs_localb(c, src_oprnd->size);
	    gen_store(c, swap_oprnd, dill_lp(c), tmp_base, src_oprnd->size,
		      integer_type, TRUE);
	    free_oprnd(c, swap_oprnd);
	    swap_oprnd = gen_fetch(c, dill_lp(c), tmp_base, src_oprnd->size,
				   src_oprnd->data_type, TRUE, FALSE);
	}
    }
    *src_oprnd = swap_oprnd;
}

void
gen_store(dill_stream c, iogen_oprnd src, dill_reg dest_reg, ssize_t offset,
	  int size, FMdata_type data_type, int aligned)
{
    assert(src.size == size);

    switch (data_type) {
    case unknown_type:
	assert(FALSE);
	break;
    case integer_type:
    case boolean_type:
    case enumeration_type:
	switch (size) {
	case 1:		/* sizeof char */
	    dill_stci(c, src.vc_reg, dest_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    dill_stsi(c, src.vc_reg, dest_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    dill_stii(c, src.vc_reg, dest_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    dill_stli(c, src.vc_reg, dest_reg, offset);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    /* vc_reg2 holds high value */
	    if (WORDS_BIGENDIAN) {
		dill_stii(c, src.vc_reg2, dest_reg, offset);
		dill_stii(c, src.vc_reg, dest_reg, offset + 4);
	    } else {
		dill_stii(c, src.vc_reg, dest_reg, offset);
		dill_stii(c, src.vc_reg2, dest_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case unsigned_type:
    case char_type:
	switch (size) {
	case 1:		/* sizeof char */
	    dill_stuci(c, src.vc_reg, dest_reg, offset);
	    break;
	case 2:		/* sizeof short */
	    dill_stusi(c, src.vc_reg, dest_reg, offset);
	    break;
	case 4:		/* sizeof int */
	    dill_stui(c, src.vc_reg, dest_reg, offset);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    dill_stuli(c, src.vc_reg, dest_reg, offset);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    /* vc_reg2 holds high value */
	    if (WORDS_BIGENDIAN) {
		dill_stui(c, src.vc_reg2, dest_reg, offset);
		dill_stui(c, src.vc_reg, dest_reg, offset + 4);
	    } else {
		dill_stui(c, src.vc_reg, dest_reg, offset);
		dill_stui(c, src.vc_reg2, dest_reg, offset + 4);
	    }
	    break;
#endif
	}
	break;
    case float_type:
	switch (size) {
	case SIZEOF_FLOAT:	/* sizeof char */
	    dill_stfi(c, src.vc_reg, dest_reg, offset);
	    break;
	case SIZEOF_DOUBLE:	/* sizeof short */
	    dill_stdi(c, src.vc_reg, dest_reg, offset);
	    break;
	}
	break;
    case string_type:

	break;
    }
}

void
gen_memcpy(dill_stream c, dill_reg src, size_t src_offset, dill_reg dest,
	   size_t dest_offset, dill_reg size, int const_size)
{
    dill_reg final_src, final_dest;
    if (src_offset != 0) {
	if (ffs_getreg(c, &final_src, DILL_P, DILL_TEMP) == 0)
	    gen_fatal("gen memcpy convert out of registers \n");
	dill_addpi(c, final_src, src, src_offset);
    } else {
	final_src = src;
    }
    if (dest_offset != 0) {
	if (ffs_getreg(c, &final_dest, DILL_P, DILL_TEMP) == 0)
	    gen_fatal("gen memcpy convert out of registers \n");
	dill_addpi(c, final_dest, dest, dest_offset);
    } else {
	final_dest = dest;
    }
    if (const_size != 0) {
	dill_scalli(c, (void*) memcpy, "memcpy", "%p%p%I", final_dest, final_src, 
		    const_size);
    } else {
	dill_scalli(c, (void*) memcpy, "memcpy", "%p%p%i", final_dest, final_src, size);
    }
    if (src_offset != 0) {
	ffs_putreg(c, final_src, DILL_P);
    }
    if (dest_offset != 0) {
	ffs_putreg(c, final_dest, DILL_P);
    }
}

void
free_oprnd(dill_stream c, iogen_oprnd oprnd)
{
    REG_DEBUG(("put %d in free\n", _vrr(oprnd.vc_reg)));
    switch (oprnd.data_type) {
    case unknown_type:
	assert(FALSE);
	break;
    case integer_type:
    case boolean_type:
    case enumeration_type:
	switch (oprnd.size) {
	case 1:		/* sizeof char */
	    ffs_putreg(c, oprnd.vc_reg, DILL_C);
	    break;
	case 2:		/* sizeof short */
	    ffs_putreg(c, oprnd.vc_reg, DILL_S);
	    break;
	case 4:		/* sizeof int */
	    ffs_putreg(c, oprnd.vc_reg, DILL_I);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    ffs_putreg(c, oprnd.vc_reg, DILL_L);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    ffs_putreg(c, oprnd.vc_reg, DILL_I);
	    ffs_putreg(c, oprnd.vc_reg2, DILL_I);
	    REG_DEBUG(("put %d in free\n", _vrr(oprnd.vc_reg2)));
	    break;
#endif
	}
	break;
    case unsigned_type:
    case char_type:
	switch (oprnd.size) {
	case 1:		/* sizeof char */
	    ffs_putreg(c, oprnd.vc_reg, DILL_UC);
	    break;
	case 2:		/* sizeof short */
	    ffs_putreg(c, oprnd.vc_reg, DILL_US);
	    break;
	case 4:		/* sizeof int */
	    ffs_putreg(c, oprnd.vc_reg, DILL_U);
	    break;
#if SIZEOF_SIZE_T != 4
	case SIZEOF_SIZE_T:
	    ffs_putreg(c, oprnd.vc_reg, DILL_UL);
	    break;
#else
	case 8:
	    /* simulate with double reg */
	    ffs_putreg(c, oprnd.vc_reg, DILL_U);
	    ffs_putreg(c, oprnd.vc_reg2, DILL_U);
	    /* vc_reg2 holds high value */
	    REG_DEBUG(("put %d in free\n", _vrr(oprnd.vc_reg2)));
	    break;
#endif
	}
	break;
    case float_type:
	switch (oprnd.size) {
	case SIZEOF_FLOAT:	/* sizeof char */
	    ffs_putreg(c, oprnd.vc_reg, DILL_F);
	    break;
	case SIZEOF_DOUBLE:	/* sizeof short */
	    ffs_putreg(c, oprnd.vc_reg, DILL_D);
	    break;
	}
	break;
    case string_type:
	break;
    }
}

iogen_oprnd
gen_type_conversion(dill_stream c, iogen_oprnd src_oprnd, FMdata_type data_type)
{
    iogen_oprnd result_oprnd = src_oprnd;
    dill_reg at;  /* temporary */
    result_oprnd.data_type = data_type;
    switch (data_type) {
    case unknown_type:
    case char_type:
	assert(FALSE);
	break;
    case integer_type:
	result_oprnd.size = sizeof(size_t);
	if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_L, DILL_TEMP))
	    gen_fatal("gen type convert out of registers \n");
	REG_DEBUG(("get %d in type_convert\n", _vrr(result_oprnd.vc_reg)));
	switch (src_oprnd.data_type) {
	case integer_type:
	case unknown_type:
	    assert(FALSE);
	    break;
	case boolean_type:
	case enumeration_type:
	case unsigned_type:
	case char_type:
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvus2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvu2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
#if SIZEOF_SIZE_T == 8
		dill_cvul2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#else
		result_oprnd.size = 8;
		if (!ffs_getreg(c, &result_oprnd.vc_reg2, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert out of registers \n");
		REG_DEBUG(("get %d in type_convert\n", _vrr(result_oprnd.vc_reg2)));
		dill_cvul2l(c, result_oprnd.vc_reg2, src_oprnd.vc_reg2);
		dill_movul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#endif
		break;
	    }
	    break;
	case float_type:
	    switch (src_oprnd.size) {
	    case SIZEOF_FLOAT:
		dill_cvf2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case SIZEOF_DOUBLE:
		dill_cvd2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    default:
		assert(FALSE);
	    }
	    break;
	default:
	    assert(FALSE);
	    break;
	}
	break;
    case boolean_type:
    case enumeration_type:
    case unsigned_type:
	result_oprnd.size = sizeof(unsigned long);
	if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_UL, DILL_TEMP))
	    gen_fatal("gen type convert out of registers \n");
	REG_DEBUG(("get %d in type_convert\n", _vrr(result_oprnd.vc_reg)));
	switch (src_oprnd.data_type) {
	case boolean_type:
	case enumeration_type:
	case unsigned_type:
	case unknown_type:
	    assert(FALSE);
	    break;
	case char_type:
	case integer_type:
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
#if SIZEOF_SIZE_T == 8
		dill_cvl2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#else
		result_oprnd.size = 8;
		if (!ffs_getreg(c, &result_oprnd.vc_reg2, DILL_UL, DILL_TEMP))
		    gen_fatal("gen type convert out of registers \n");
		REG_DEBUG(("get %d in type_convert\n", _vrr(result_oprnd.vc_reg2)));
		dill_cvl2ul(c, result_oprnd.vc_reg2, src_oprnd.vc_reg2);
		dill_movul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#endif
		dill_cvul2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case float_type:
	    switch (src_oprnd.size) {
	    case SIZEOF_FLOAT:
		dill_cvf2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case SIZEOF_DOUBLE:
		dill_cvd2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    default:
		assert(FALSE);
	    }
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case float_type:
	result_oprnd.size = sizeof(double);
	if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_D, DILL_TEMP))
	    gen_fatal("gen type convert out of registers \n");
	REG_DEBUG(("get %d in type_convert\n", _vrr(result_oprnd.vc_reg)));
	switch (src_oprnd.data_type) {
	case boolean_type:
	case enumeration_type:
	case unsigned_type:
	    switch (src_oprnd.size) {
	    case 1:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvc2l(c, at, src_oprnd.vc_reg);
		dill_cvl2d(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 2:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvs2l(c, at, src_oprnd.vc_reg);
		dill_cvl2d(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 4:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvi2l(c, at, src_oprnd.vc_reg);
		dill_cvl2d(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 8:
#if SIZEOF_SIZE_T == 8
		dill_cvl2d(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#else
		{
		    dill_reg dat; /* temporary */
		    if (!ffs_getreg(c, &dat, DILL_D, DILL_TEMP))
			gen_fatal("gen type convert2 out of registers \n");
		    dill_cvu2d(c, dat, src_oprnd.vc_reg);
		    dill_setd(c, result_oprnd.vc_reg, 4294967296.0);
		    dill_muld(c, result_oprnd.vc_reg, result_oprnd.vc_reg, dat);
		    dill_cvi2d(c, dat, src_oprnd.vc_reg);
		    dill_addd(c, result_oprnd.vc_reg, result_oprnd.vc_reg, dat);
		    ffs_putreg(c, dat, DILL_D);
		}
#endif
		break;
	    }
	    break;
	case char_type:
	case integer_type:
	    switch (src_oprnd.size) {
	    case 1:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvc2l(c, at, src_oprnd.vc_reg);
		dill_cvl2d(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 2:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvc2l(c, at, src_oprnd.vc_reg);
		dill_cvl2d(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 4:
		dill_cvi2d(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
#if SIZEOF_SIZE_T == 8
		dill_cvl2d(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
#else
		{
		    dill_reg dat; /* temporary */
		    if (!ffs_getreg(c, &dat, DILL_D, DILL_TEMP))
			gen_fatal("gen type convert2 out of registers \n");
		    dill_cvu2d(c, dat, src_oprnd.vc_reg);
		    dill_setd(c, result_oprnd.vc_reg, 4294967296.0);
		    dill_muld(c, result_oprnd.vc_reg, result_oprnd.vc_reg, dat);
		    dill_cvi2d(c, dat, src_oprnd.vc_reg);
		    dill_addd(c, result_oprnd.vc_reg, result_oprnd.vc_reg, dat);
		    ffs_putreg(c, dat, DILL_D);
		}
#endif
		break;
	    }
	    break;
	case float_type:
	    assert(FALSE);
	    break;
	case unknown_type:
	default:
	    assert(FALSE);
	    break;
	}
	break;

    case string_type:
	assert(FALSE);
	break;
    }
    return result_oprnd;
}

iogen_oprnd
gen_size_conversion(dill_stream c, iogen_oprnd src_oprnd, int size)
{
    iogen_oprnd result_oprnd = src_oprnd;
    dill_reg at;  /* temporary */
    result_oprnd.size = size;
    switch (src_oprnd.data_type) {
    case integer_type:
	switch (size) {
	case 1:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_C, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 2:
		dill_cvs2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		dill_cvl2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 2:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_S, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvc2l(c, at, src_oprnd.vc_reg);
		dill_cvl2s(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 4:
		dill_cvi2s(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2s(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 4:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_I, DILL_TEMP))
		gen_fatal("gen size convert out of registers E\n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2i(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2i(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2i(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 8:
#if SIZEOF_SIZE_T == 8
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_L, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
#else
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_L, DILL_TEMP) ||
		(!ffs_getreg(c, &result_oprnd.vc_reg2, DILL_L, DILL_TEMP)))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg2)));
	    dill_setl(c, result_oprnd.vc_reg2, 0);
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
#endif
	}
	break;
    case boolean_type:
    case enumeration_type:
    case unsigned_type:
    case char_type:
	switch (size) {
	case 1:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_UC, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 2:
		dill_cvs2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		dill_cvl2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2l(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		dill_cvl2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2c(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 2:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_US, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		if (!ffs_getreg(c, &at, DILL_L, DILL_TEMP))
		    gen_fatal("gen type convert2 out of registers \n");
		dill_cvc2l(c, at, src_oprnd.vc_reg);
		dill_cvl2s(c, result_oprnd.vc_reg, at);
		ffs_putreg(c, at, DILL_L);
		break;
	    case 4:
		dill_cvi2s(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2s(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 4:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_U, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2u(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2u(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 8:
		dill_cvl2u(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
	case 8:
#if SIZEOF_SIZE_T == 8
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_UL, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvuc2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvus2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvu2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
#else
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_UL, DILL_TEMP) ||
		!ffs_getreg(c, &result_oprnd.vc_reg2, DILL_UL, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg2)));
	    dill_setl(c, result_oprnd.vc_reg2, 0);
	    switch (src_oprnd.size) {
	    case 1:
		dill_cvc2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 2:
		dill_cvs2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    case 4:
		dill_cvi2ul(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    }
	    break;
#endif
	}
	break;
    case float_type:
	switch (size) {
	case SIZEOF_FLOAT:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_F, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case SIZEOF_FLOAT:
		assert(FALSE);
		break;
	    case SIZEOF_DOUBLE:
		dill_cvd2f(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    default:
		assert(FALSE);
	    }
	    break;
	case SIZEOF_DOUBLE:
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_D, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case SIZEOF_DOUBLE:
		assert(FALSE);
		break;
	    case SIZEOF_FLOAT:
		dill_cvf2d(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    default:
		assert(FALSE);
	    }
	    break;
	default:
	    /* really should fail */
	    if (!ffs_getreg(c, &result_oprnd.vc_reg, DILL_D, DILL_TEMP))
		gen_fatal("gen size convert out of registers \n");
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    REG_DEBUG(("get %d in size convert\n", _vrr(result_oprnd.vc_reg)));
	    switch (src_oprnd.size) {
	    case SIZEOF_DOUBLE:
/*		assert(FALSE);*/
		break;
	    case SIZEOF_FLOAT:
		dill_cvf2d(c, result_oprnd.vc_reg, src_oprnd.vc_reg);
		break;
	    default:
		assert(FALSE);
	    }
	    break;
	    
	}
	break;
    case unknown_type:
    default:
	assert(FALSE);
	break;
    }
    return result_oprnd;
}
#endif
