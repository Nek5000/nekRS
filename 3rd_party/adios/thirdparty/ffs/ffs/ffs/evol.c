#include "config.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <ffs.h>
#include "fm_internal.h"
#include "ffs_internal.h"

#define MAX_DIFF	0xFFFF
#define COMPAT_THRESH	0.8


static FMStructDescList
build_struct_list(FMFormat format)
{
    FMStructDescList ret;
    int count = 0;
    while (format->subformats[count] != NULL) count++;
    count+=2;
    ret = malloc(count * sizeof(ret[0]));
    ret[0].format_name = format->format_name;
    ret[0].field_list = format->field_list;
    ret[0].struct_size = format->record_length;
    ret[0].opt_info = NULL;
    count = 0;
    while(format->subformats[count] != NULL) {
	FMFormat subformat = format->subformats[count];
	ret[count+1].format_name = subformat->format_name;
	ret[count+1].field_list = subformat->field_list;
	ret[count+1].struct_size = subformat->record_length;
	ret[count+1].opt_info = NULL;
	count++;
    }
    ret[count+1].format_name = NULL;
    ret[count+1].field_list = NULL;
    ret[count+1].struct_size = 0;
    ret[count+1].opt_info = NULL;
    return ret;
}

extern int
FMformat_compat_cmp(FMFormat format, FMFormat *formatList,
		    int listSize, FMcompat_formats * older_format);

extern void
FFS_determine_conversion(FFSContext c, FFSTypeHandle format)
{
    int i;
    FMStructDescList struct_list;

    int nearest_format = -1, j;
    FMcompat_formats older_format = NULL;
    FMFormat *formatList;

    formatList =
	(FMFormat *) malloc(c->handle_list_size * sizeof(FMFormat));

    j = 0;
    for (i = 0; i < c->handle_list_size; i++) {
	if (c->handle_list[i] && c->handle_list[i]->is_fixed_target)
	    formatList[j++] = c->handle_list[i]->body;
    }
    nearest_format = FMformat_compat_cmp(format->body, formatList,
					 j, &older_format);

    if (nearest_format == -1) {
	free(formatList);
	format->status = none_available;
	return;
    }
    struct_list = build_struct_list(formatList[nearest_format]);
    establish_conversion(c, format, struct_list);
    format->conversion_target = FFSTypeHandle_by_index(c, formatList[nearest_format]->format_index);
    format->status = conversion_set;
    free(formatList);
    free(struct_list);
    return;
}

static int
IO_field_type_eq(const char *str1, const char *str2)
{
    FMdata_type t1, t2;
    long t1_count, t2_count;

    t1 = FMarray_str_to_data_type(str1, &t1_count);
    t2 = FMarray_str_to_data_type(str2, &t2_count);

    if ((t1_count == -1) && (t2_count == -1)) {
	/* variant array */
	char *tmp_str1 = base_data_type(str1);
	char *tmp_str2 = base_data_type(str2);
	
	char *colon1 = strchr(tmp_str1, ':');
	char *colon2 = strchr(tmp_str2, ':');
	char *lparen1 = strchr(str1, '[');
	char *lparen2 = strchr(str2, '[');
	intptr_t count1 = 0;
	intptr_t count2 = 0;

	if (colon1 != NULL) {
	    count1 = colon1 - tmp_str1;
	} else {
	    count1 = strlen(tmp_str1);
	}
	if (colon2 != NULL) {
	    count2 = colon2 - tmp_str2;
	} else {
	    count2 = strlen(tmp_str2);
	}
	/*compare base type */
	if (strncmp(tmp_str1, tmp_str2,(count1>count2)?count1:count2) != 0) {
	    /* base types differ */
	    return 0;
	}
	free(tmp_str1);
	free(tmp_str2);
	if ((lparen1 == NULL) || (lparen2 == NULL)) return -1;
	return (strcmp(lparen1, lparen2) == 0);
    }
    return ((t1 == t2) && (t1_count == t2_count));
}


/* 
 * Compares the two formats.
 * Returns FFSformat_order
 * */
static FMformat_order
FMformat_cmp_diff(FMFormat format1, FMFormat format2, int *diff1, int *diff2)
{
    FMformat_order tmp_result = Format_Equal;
    FMFieldList orig_field_list1 = format1->field_list;
    FMFieldList orig_field_list2 = format2->field_list;
    FMFieldList field_list1, field_list2;
    FMFormat *subformats1 = NULL, *subformats2 = NULL;
    int field_count1, field_count2;
    int i, j, limit;

    if (format1 == format2) return Format_Equal;

    /* count fields */
    for (field_count1 = 0; orig_field_list1[field_count1].field_name != NULL;
	 field_count1++);

    /* count fields */
    for (field_count2 = 0; orig_field_list2[field_count2].field_name != NULL;
	 field_count2++);

    field_list1 = malloc(sizeof(field_list1[0]) * (field_count1 + 1));
    field_list2 = malloc(sizeof(field_list2[0]) * (field_count2 + 1));
    memcpy(field_list1, orig_field_list1, 
	   sizeof(field_list1[0]) * (field_count1 + 1));
    memcpy(field_list2, orig_field_list2, 
	   sizeof(field_list2[0]) * (field_count2 + 1));

    qsort(field_list1, field_count1, sizeof(FMField), field_name_compar);
    qsort(field_list2, field_count2, sizeof(FMField), field_name_compar);

    limit = field_count1;
    if (field_count2 < limit)
	limit = field_count2;
    i = j = 0;
    while ((i < field_count1) && (j < field_count2)) {
	int name_cmp;
	if ((name_cmp = strcmp(field_list1[i].field_name,
			       field_list2[j].field_name)) == 0) {
	    /* fields have same name */
	    if (!IO_field_type_eq(field_list1[i].field_type,
				  field_list2[j].field_type)) {
		(*diff1)++;
		(*diff2)++;
	    }
	    i++;
	    j++;
	} else if (name_cmp < 0) {
	    /* name_cmp<0 a field in field_list1 that doesn't appear in 2 */
	    (*diff1)++;
	    i++;
	} else {
	    (*diff2)++;
	    j++;
	}
    }
    (*diff1) += (field_count1 - i);
    (*diff2) += (field_count2 - j);

    /* go through subformats */
    subformats1 = format1->subformats;
    subformats2 = format2->subformats;
/* TODO: Fix for unmatched subformats 
 * -sandip */
    while (subformats1 && (*subformats1 != NULL)) {
	char *sub1_name = name_of_FMformat(*subformats1);
	int i = 0;
	if (*subformats1 == format1) {
	    /* self appears in subformat list, skip it */
	    subformats1++;
	    continue;
	}
	while (subformats2 && (subformats2[i] != NULL)) {
	    if (strcmp(sub1_name, name_of_FMformat(subformats2[i])) == 0) {
		/* same name, compare */
		FMformat_cmp_diff(*subformats1, subformats2[i], diff1,
				  diff2);
		break;
	    }
	    i++;
	}
	subformats1++;
    }

    free(field_list1);
    free(field_list2);
    if (*diff1 == 0) {
	if (*diff2 == 0)
	    tmp_result = Format_Equal;
	else
	    tmp_result = Format_Less;
    } else {
	if (*diff2 == 0)
	    tmp_result = Format_Greater;
	else
	    tmp_result = Format_Incompatible;
    }
    return tmp_result;
}

/* 
 * Basically counts the total number of nodes in a data structure tree
 * */
static int
count_total_IOfield(FMFormat format)
{
    int count = 0;
    if (format) {
	int i;
	count = format->field_count;
	for (i = 0; i < format->field_count; i++)
	    if (format->field_subformats[i] != NULL)
		count += count_total_IOfield(format->field_subformats[i]);
    }
    return count;
}

/* 
 * Implements a simple threshold-based technique to decide whether or not
 * we should do conversion at all when the comparison result between best 
 * possible pair is Format_Incompatible
 * The idea is to keep the number of missing fields below a certain threshold.
 * There is scope of smarter algorithm here.
 * */
static int
check_compat_thresh(FMFormat_Comp_result * comp_result, FMFormat format)
{
    float curr_thresh;
    int field_count = count_total_IOfield(format);

    curr_thresh = (float) comp_result->diff2 / (float) field_count;

    return (curr_thresh < (1.0 - COMPAT_THRESH));
}

/* 
 * Try to find a format in 'formatList' which best matches the 'format'
 * Returns the index of best match found or 
 * -1 if no better comparison (as specified by comp_result) found.
 * comp_result is modified accordingly to indicate the current match purity
 * */
static int
IOformat_list_cmp(FMFormat format, FMFormat *formatList, int listSize,
		  FMFormat_Comp_result * comp_result)
{
    int i, diff1, diff2, nearest_format = -1;
    FMformat_order result = Format_Incompatible;

    for (i = 0; i < listSize; i++) {
	int order;
	if (formatList[i] == NULL)
	    continue;
	order =
	    strcmp(name_of_FMformat(format),
		   name_of_FMformat(formatList[i]));
	if (order != 0) continue;
	diff1 = diff2 = 0;
	result = FMformat_cmp_diff(format, formatList[i], &diff1, &diff2);
	if (result == Format_Equal) {
	    comp_result->diff1 = comp_result->diff2 = 0;
	    nearest_format = i;
	    break;
	}
	if ((diff2 < comp_result->diff2) ||
	    (diff2 == comp_result->diff2 && diff1 < comp_result->diff1)) {
	    comp_result->diff1 = diff1;
	    comp_result->diff2 = diff2;
	    nearest_format = i;
	}
    }
    return nearest_format;
}

static int
IOformat_list_cmp2(FMFormat format, FMFormat *formatList, int listSize,
		  FMFormat_Comp_result * comp_result)
{
    int i, diff1, diff2, nearest_format = -1;
    FMformat_order result = Format_Incompatible;

    for (i = 0; i < listSize; i++) {
	if (formatList[i] == NULL)
	    continue;
	diff1 = diff2 = 0;
	result = FMformat_cmp_diff(format, formatList[i], &diff1, &diff2);
	if (result == Format_Equal) {
	    comp_result->diff1 = comp_result->diff2 = 0;
	    nearest_format = i;
	    break;
	}
	if ((diff2 < comp_result->diff2) ||
	    (diff2 == comp_result->diff2 && diff1 < comp_result->diff1)) {
	    comp_result->diff1 = diff1;
	    comp_result->diff2 = diff2;
	    nearest_format = i;
	}
    }
    return nearest_format;
}

/* 
 * First check whether 'format' matches exactly with any of 'formatList'
 * If exact match not found, try to match with compatible formats (formats 
 * that 'format' can be translated to)
 * Returns the index of best match found or 
 * -1 if no better comparison (as specified by comp_result) found.
 * */
extern int
FMformat_compat_cmp(FMFormat format, FMFormat *formatList,
		    int listSize, FMcompat_formats * older_format)
{
    FMFormat prior_format;
    int i = 0, nearest_format = -1;
    FMcompat_formats compats;
    FMFormat_Comp_result comp_result = { MAX_DIFF, MAX_DIFF };

    *older_format = NULL;
    nearest_format =
	IOformat_list_cmp(format, formatList, listSize, &comp_result);
    if (nearest_format != -1 && !comp_result.diff1 && !comp_result.diff2)
	return nearest_format;

    compats = FMget_compat_formats(format);
    if (compats == NULL)
	return -1;
    while ((prior_format = compats[i].prior_format)) {
	int tmp = IOformat_list_cmp(prior_format, formatList, listSize,
				    &comp_result);
	if (tmp != -1) {
	    nearest_format = tmp;
	    *older_format = &compats[i];
	}
	if (comp_result.diff1 == 0 && comp_result.diff2 == 0)
	    break;
	i++;
    }
    if (nearest_format != -1 &&
	!check_compat_thresh(&comp_result, formatList[nearest_format])){
	nearest_format = -1;
	*older_format = NULL;
    }
    return nearest_format;
}

/* 
 * First check whether 'format' matches exactly with any of 'formatList'
 * If exact match not found, try to match with compatible formats (formats 
 * that 'format' can be translated to)
 * Returns the index of best match found or 
 * -1 if no better comparison (as specified by comp_result) found.
 * */
extern int
FMformat_compat_cmp2(FMFormat format, FMFormat *formatList,
		    int listSize, FMcompat_formats * older_format)
{
    FMFormat prior_format;
    int i = 0, nearest_format = -1;
    FMcompat_formats compats;
    FMFormat_Comp_result saved_comp_result = { MAX_DIFF, MAX_DIFF };

    *older_format = NULL;
    nearest_format =
	IOformat_list_cmp2(format, formatList, listSize, &saved_comp_result);
    if (nearest_format != -1 && !saved_comp_result.diff1 && 
	!saved_comp_result.diff2) {
	/* exact match */
	return nearest_format;
    }

    compats = FMget_compat_formats(format);
    if (compats == NULL) {
	if (!saved_comp_result.diff2) {
	    return nearest_format;
	}
	return -1;
    }
    while ((prior_format = compats[i].prior_format)) {
	FMFormat_Comp_result comp_result={MAX_DIFF, MAX_DIFF};

	int tmp = IOformat_list_cmp2(prior_format, formatList, listSize,
				    &comp_result);
	if (tmp != -1) {
	    if (saved_comp_result.diff1 > comp_result.diff1) {
		nearest_format = tmp;
		*older_format = &compats[i];
		saved_comp_result = comp_result;
	    }
	}
	if (comp_result.diff1 == 0 && comp_result.diff2 == 0)
	    break;
	i++;
    }
    if (nearest_format != -1 &&
	!check_compat_thresh(&saved_comp_result, formatList[nearest_format])){
	nearest_format = -1;
	*older_format = NULL;
    }
    return nearest_format;
}

extern void *
FFScreate_compat_info(FMFormat prior_format, char *xform_code, size_t *len_p)
{
    char *block;
    size_t block_len = strlen(xform_code) + prior_format->server_ID.length +1;
    block = malloc(block_len);

    memcpy(block, prior_format->server_ID.value,
	   prior_format->server_ID.length);
    memcpy(block + prior_format->server_ID.length,
	   xform_code, block_len - prior_format->server_ID.length);
    *len_p = block_len;
    return block;
}
    
#ifdef NOT_DEF
/**
 * Localize the "format" and set the conversion context to convert the wire
 * format to this localized format.
 * Return the localized format list
 */
extern FMFormatList
IOlocalize_conv(FFSContext context, FMFormat format)
{
    FMFormatList local_format_list = NULL;
    FMFormat *wire_subformats = format->subformats;
    int i = 0;

    while (wire_subformats[i] != NULL) {
	int local_struct_size;
	FMFieldList wire_field_list = wire_subformats[i]->field_list;
	FMFieldList local_field_list = copy_field_list(wire_field_list);

	/* 
	 * determine an appropriate native layout for this structure
	 * and set the conversion
	 */
	local_field_list = localize_field_list(local_field_list, context);
	local_struct_size =
	    struct_size_field_list(local_field_list, sizeof(char *));
	set_conversion_IOcontext(context, wire_subformats[i],
				 local_field_list, local_struct_size);
	local_format_list =
	    realloc(local_format_list,
		    sizeof(local_format_list[0]) * (i + 2));
	local_format_list[i].format_name =
	    strdup(name_of_FMformat(wire_subformats[i]));

	local_format_list[i].field_list = local_field_list;
	i++;
    }
    local_format_list[i].format_name = NULL;
    free(wire_subformats);
    return local_format_list;
}

/**
 * Localize the "format".
 * Register this localized format to the "context".
 * Set the conversion context to convert this localized format to native 
 * format as specified by native_field_list and native_subformat_list.
 * Return the localized format list.
 */
extern FMFormatList
IOlocalize_register_conv(FMContext context, FMFormat format,
			 FMFieldList native_field_list,
			 FMFormatList native_subformat_list,
			 FMFormat *local_prior_format,
			 int *local_struct_size_out)
{
    FMFormatList local_f1_formats = NULL;
    FMFormat *wire_subformats = format->subformats;
    int i = 0, native_struct_size;

    native_struct_size = struct_size_field_list(native_field_list,
						sizeof(char *));

#ifdef NOTDEF
    while (wire_subformats[i] != NULL) {
	int local_struct_size;
	FMFieldList wire_field_list =
	    field_list_of_IOformat(wire_subformats[i]);
	FMFieldList local_field_list = copy_field_list(wire_field_list);
	char *subformat_name = name_of_IOformat(wire_subformats[i]);
	int j = 0;

	/* 
	 * determine an appropriate native layout for this structure
	 * and set the conversion
	 */
	local_field_list = localize_field_list(local_field_list, context);
	local_struct_size =
	    struct_size_field_list(local_field_list, sizeof(char *));

	/* We don't need the subformats. The last one will be the top
	 * level FMFormat */

	*local_prior_format =
	    register_IOcontext_format(subformat_name, local_field_list,
				      context);
	while (native_subformat_list
	       && (native_subformat_list[j].format_name != NULL)) {
	    if (strcmp
		(native_subformat_list[j].format_name,
		 subformat_name) == 0) {
		FMFieldList sub_field_list;
		int sub_struct_size;
		sub_field_list = native_subformat_list[j].field_list;
		sub_struct_size = struct_size_field_list(sub_field_list,
							 sizeof(char *));
		set_conversion_IOcontext(context,
					 *local_prior_format,
					 sub_field_list, sub_struct_size);
		break;
	    }
	    j++;
	}

	local_f1_formats = realloc(local_f1_formats,
				   sizeof(local_f1_formats[0]) * (i + 2));
	local_f1_formats[i].format_name =
	    strdup(name_of_IOformat(wire_subformats[i]));
	local_f1_formats[i].field_list = local_field_list;
	if (wire_subformats[i + 1] == NULL) {
	    /* last one */
	    *local_struct_size_out = local_struct_size;
	    local_f1_formats[i + 1].format_name = NULL;
	}
	i++;
    }
#endif
    set_conversion_IOcontext(context, *local_prior_format,
			     native_field_list, native_struct_size);
    free(wire_subformats);
    return local_f1_formats;
}
#endif
