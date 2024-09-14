/*
 * t6 - was created to exercise a prior bug with arrays of polygon data
 */
#include "cod.h"
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include "malloc.h"
#endif
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>
#include <string.h>

typedef struct
{
  unsigned short x;
  unsigned short y;
  unsigned short z;
  unsigned short r;
  unsigned short g;
  unsigned short b;
} PipelinedPoint;

typedef struct
{
   int num_points;
   PipelinedPoint *polygon_points;
} PolygonPoints;
typedef struct
{
  int num_points;
  PolygonPoints* image_data;
} FrameData;

static FMField PipelinedPoint_field_list[] = {
 {"x", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, x)},
 {"y", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, y)},
 {"z", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, z)},
 {"r", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, r)},
 {"g", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, g)},
 {"b", "integer", sizeof(unsigned short), FMOffset(PipelinedPoint*, b)},
 {NULL, NULL}
};
static FMField PolygonPoints_field_list[] = {
{"num_points", "integer", sizeof(int), FMOffset(PolygonPoints*, num_points)},
{"polygon_points", "PipelinedPoint[num_points]", sizeof(PipelinedPoint), FMOffset(PolygonPoints*, polygon_points)},
 {NULL, NULL}
};
static FMField FrameData_field_list[] = {
{"num_points", "integer", sizeof(int), FMOffset(FrameData*, num_points)},
{"image_data", "PolygonPoints[num_points]", sizeof(PolygonPoints), FMOffset(FrameData*, image_data)},
{NULL, NULL, 0 , 0}
};


int
main() 
{
    static char extern_string[] = "int printf(string format, ...);";
    static cod_extern_entry externs[] =
    {
        {"printf", (void*)(intptr_t)printf},
        {(void*)0, (void*)0}
    };
    static char code[] = "{\n\
    FrameData *f;\n\
        output.num_points = 1;\
        output.image_data[0].num_points = 1;\
	return 1;\n\
     }";
    FrameData data;
   
    cod_parse_context context = new_cod_parse_context();

    cod_code gen_code;
    void (*func)(void*);

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);

    cod_add_simple_struct_type("PipelinedPoint", PipelinedPoint_field_list, context);
    cod_add_simple_struct_type("PolygonPoints", PolygonPoints_field_list, context);
    cod_add_simple_struct_type("FrameData", FrameData_field_list, context);
    cod_subroutine_declaration("int proc(FrameData *output)", context);
   
    gen_code = cod_code_gen(code, context);
    func = (void (*)(void*))(intptr_t)gen_code->func;

    data.num_points = 0;
    data.image_data = NULL;
    func(&data);
    cod_code_free(gen_code);
    cod_free_parse_context(context);
    if ((data.num_points != 1) || (data.image_data[0].num_points != 1))  {
	free(data.image_data[0].polygon_points);
	free(data.image_data);
	return 1;
    }
    free(data.image_data[0].polygon_points);
    free(data.image_data);
    return 0;
}
