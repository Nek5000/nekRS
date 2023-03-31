/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

 */

#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

// structure used to encode vertices that make
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
struct face_t
{
  dlong element;
  int face;

  dlong elementNeighbor; // neighbor element
  int faceNeighbor;    // neighbor face

  int NfaceVertices;

  hlong v[4];
};

// comparison function that orders vertices
// based on their combined vertex indices
int compareVertices(const void* a,
                    const void* b)
{
  face_t* fa = (face_t*) a;
  face_t* fb = (face_t*) b;

  for(int n = 0; n < fa->NfaceVertices; ++n) {
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;
}

/* comparison function that orders element/face
   based on their indexes */
int compareFaces(const void* a,
                 const void* b)
{
  face_t* fa = (face_t*) a;
  face_t* fb = (face_t*) b;

  if(fa->element < fb->element) return -1;
  if(fa->element > fb->element) return +1;

  if(fa->face < fb->face) return -1;
  if(fa->face > fb->face) return +1;

  return 0;
}

/* routine to find EToE (Element To Element)
   and EToF (Element To Local Face) connectivity arrays */
void meshConnect(mesh_t* mesh)
{
  /* build list of faces */
  face_t* faces =
    (face_t*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(face_t));

  dlong cnt = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for(int f = 0; f < mesh->Nfaces; ++f) {
      for(int n = 0; n < mesh->NfaceVertices; ++n) {
        dlong vid = e * mesh->Nverts + mesh->faceVertices[f * mesh->NfaceVertices + n];
        faces[cnt].v[n] = mesh->EToV[vid];
      }

      mysort(faces[cnt].v, mesh->NfaceVertices, "descending");

      faces[cnt].NfaceVertices = mesh->NfaceVertices;

      faces[cnt].element = e;
      faces[cnt].face = f;

      faces[cnt].elementNeighbor = -1;
      faces[cnt].faceNeighbor = -1;

      ++cnt;
    }
  }

  /* sort faces by their vertex number pairs */
  qsort(faces,
        mesh->Nelements * mesh->Nfaces,
        sizeof(face_t),
        compareVertices);

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for (cnt = 0; cnt < mesh->Nelements * mesh->Nfaces - 1; ++cnt) {

    if (!compareVertices(faces + cnt, faces + cnt + 1)) { // match
      faces[cnt].elementNeighbor = faces[cnt + 1].element;
      faces[cnt].faceNeighbor = faces[cnt + 1].face;

      faces[cnt + 1].elementNeighbor = faces[cnt].element;
      faces[cnt + 1].faceNeighbor = faces[cnt].face;
    }
  }

  /* resort faces back to the original element/face ordering */
  qsort(faces,
        mesh->Nelements * mesh->Nfaces,
        sizeof(face_t),
        compareFaces);

  /* extract the element to element and element to face connectivity */
  mesh->EToE = (dlong*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(dlong));
  mesh->EToF = (int*)   calloc(mesh->Nelements * mesh->Nfaces, sizeof(int  ));

  cnt = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e) {
    for(int f = 0; f < mesh->Nfaces; ++f) {
      mesh->EToE[cnt] = faces[cnt].elementNeighbor;
      mesh->EToF[cnt] = faces[cnt].faceNeighbor;
      nrsCheck(mesh->EToE[cnt] >= mesh->Nelements, MPI_COMM_SELF, EXIT_FAILURE,
               "Invalid EToE(%d,%d) = %d \n", e,f, mesh->EToE[cnt]);
      nrsCheck(mesh->EToF[cnt] >= mesh->Nfaces, MPI_COMM_SELF, EXIT_FAILURE,
               "Invalid EToF(%d,%d) = %d \n", e,f, mesh->EToF[cnt]);
      ++cnt;
    }
  }
}
