#include <genmap-impl.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

int GenmapCreateVector(GenmapVector *x, GenmapInt size) {
  /* Asserts:
       - size > 0
  */
  assert(size > 0);

  GenmapMalloc(1, x);
  if(*x == NULL) {
    return 1;
  }

  (*x)->size = size;
  (*x)->data = NULL;

  GenmapMalloc((size_t)size, &(*x)->data);
  if((*x)->data == NULL) {
    return 1;
  }

  return 0;
}

int GenmapVectorsEqual(GenmapVector x, GenmapVector y,
                       GenmapScalar tol) {
  /* Asserts:
       - size of y == size of x
  */
  assert(x->size == y->size);

  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    assert(!isnan(x->data[i]) && !isnan(y->data[i]));
  }

  GenmapInt n = x->size;
  for(i = 0; i < n; i++) {
    if(fabs(x->data[i] - y->data[i]) > tol) {
      return 0;
    }
  }

  return 1;
}

int GenmapSetVector(GenmapVector x, GenmapScalar *array) {
  memcpy(x->data, array, sizeof(GenmapScalar) * (size_t)x->size);
  return 0;
}

int GenmapDestroyVector(GenmapVector x) {
  if(x->data) {
    GenmapFree(x->data);
  }

  if(x) {
    GenmapFree(x);
  }

  return 0;
}

int GenmapCopyVector(GenmapVector y, GenmapVector x) {
  /* Asserts:
       - size y = size x
  */
  assert(y->size >= x->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for(i = 0; i < n; i++) {
    y->data[i] = x->data[i];
  }

  return 0;
}

GenmapScalar GenmapNormVector(GenmapVector x, GenmapInt p) {
  assert(x->size > 0);

  GenmapInt n = x->size;
  GenmapScalar norm = 0;

  GenmapInt i;
  if(p == 1) {
    for(i = 0; i < n; i++) {
      norm += fabs(x->data[i]);
    }
  } else if(p == 2) {
    for(i = 0; i < n; i++) {
      norm += x->data[i] * x->data[i];
    }
    norm = sqrt(norm);
  } else if(p == -1) {
    norm = fabs(x->data[0]);

    for(i = 1; i < n; i++) {
      if(fabs(x->data[i]) > norm) norm = fabs(x->data[i]);
    }
  }

  return norm;
}

int GenmapScaleVector(GenmapVector y, GenmapVector x,
                      GenmapScalar alpha) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapInt n = x->size;
  GenmapInt i;
  for(i = 0; i < n; i++) {
    y->data[i] = alpha * x->data[i];
  }

  return 0;
}

int GenmapCreateOnesVector(GenmapVector *x, GenmapInt size) {
  GenmapCreateVector(x, size);

  GenmapInt i;
  for(i = 0; i < size; i++) {
    (*x)->data[i] = 1.;
  }

  return 0;
}

int GenmapCreateZerosVector(GenmapVector *x, GenmapInt size) {
  GenmapCreateVector(x, size);

  GenmapInt i;
  for(i = 0; i < size; i++) {
    (*x)->data[i] = 0.;
  }

  return 0;
}

GenmapScalar GenmapDotVector(GenmapVector y, GenmapVector x) {
  /* asserts:
       - size x = size y
  */
  assert(x->size == y->size);

  GenmapScalar result = 0.0;
  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    result += x->data[i] * y->data[i];
  }

  return result;
}

GenmapScalar GenmapAbsMaxVector(GenmapVector x) {
  GenmapScalar result = 0.0;
  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    if(fabs(x->data[i]) > result) {
      result = fabs(x->data[i]);
    }
  }

  return result;
}

GenmapScalar GenmapMaxVector(GenmapVector x) {
  GenmapScalar result = -DBL_MAX;
  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    if(x->data[i] > result) {
      result = x->data[i];
    }
  }

  return result;
}

GenmapScalar GenmapAbsMinVector(GenmapVector x) {
  GenmapScalar result = DBL_MAX;
  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    if(fabs(x->data[i]) < result) {
      result = fabs(x->data[i]);
    }
  }

  return result;
}

GenmapScalar GenmapMinVector(GenmapVector x) {
  GenmapScalar result = DBL_MAX;
  GenmapInt i;
  for(i = 0; i < x->size; i++) {
    if(x->data[i] < result) {
      result = x->data[i];
    }
  }

  return result;
}

int GenmapAxpbyVector(GenmapVector z, GenmapVector x,
                      GenmapScalar alpha,
                      GenmapVector y, GenmapScalar beta) {
  assert(z->size == x->size);
  assert(z->size == y->size);

  GenmapInt n = z->size;
  GenmapInt i;
  for(i = 0; i < n; i++) {
    z->data[i] = alpha * x->data[i] + beta * y->data[i];
  }

  return 0;
}

int GenmapPrintVector(GenmapVector x) {
  /* Asserts:
       - size x > 0
  */
  assert(x->size > 0);

  printf("(%lf", x->data[0]);
  GenmapInt i;
  for(i = 1; i < x->size - 1; i++) {
    printf(", %.10lf", x->data[i]);
  }

  if(x->size > 1) {
    printf(", %.10lf)", x->data[x->size - 1]);
  } else {
    printf(")");
  }

  return 0;
}
