#include <sort-impl.h>

static int sort_field(struct array *arr, size_t usize, gs_dom t, uint off,
                      buffer *buf, int keep) {
  uint nunits = arr->n;
  void *ptr = arr->ptr;
  switch (t) {
  case gs_double:
    gslib_sortp_double(buf, keep, (double *)((char *)ptr + off), nunits, usize);
    break;
#if 0 // FIXME
    case gs_float:
        gslib_sortp_float (buf, keep, (float *)((char*)ptr+off), 
            nunits, usize);
      break;
#endif
  case gs_long: // FIXME gs_ulong
    gslib_sortp_ull(buf, keep, (ulong *)((char *)ptr + off), nunits, usize);
    break;
  case gs_int: // FIXME gs_uint
    gslib_sortp_ui(buf, keep, (uint *)((char *)ptr + off), nunits, usize);
    break;
  default:
    break;
  }

  return 0;
}

int sort_local(struct sort *s) {
  struct array *a = s->a;
  buffer *buf = s->buf;
  size_t usize = s->unit_size;
  int i = s->nfields - 1;

  sort_field(a, usize, s->t[i], s->offset[i], buf, 0), i--;
  while (i >= 0)
    sort_field(a, usize, s->t[i], s->offset[i], buf, 1), i--;
  sarray_permute_buf_(s->align, usize, a->ptr, a->n, buf);

  return 0;
}
