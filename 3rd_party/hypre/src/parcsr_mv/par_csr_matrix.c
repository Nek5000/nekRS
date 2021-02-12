/******************************************************************************
 * Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
 * HYPRE Project Developers. See the top-level COPYRIGHT file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 ******************************************************************************/

/******************************************************************************
 *
 * Member functions for hypre_ParCSRMatrix class.
 *
 *****************************************************************************/

#include "_hypre_parcsr_mv.h"

#include "../seq_mv/HYPRE_seq_mv.h"
#include "../seq_mv/csr_matrix.h"

/* In addition to publically accessible interface in HYPRE_mv.h, the
   implementation in this file uses accessor macros into the sequential matrix
   structure, and so includes the .h that defines that structure. Should those
   accessor functions become proper functions at some later date, this will not
   be necessary. AJC 4/99 */

#ifdef HYPRE_NO_GLOBAL_PARTITION
HYPRE_Int hypre_FillResponseParToCSRMatrix(void*, HYPRE_Int, HYPRE_Int, void*, MPI_Comm, void**, HYPRE_Int*);
#endif

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixCreate
 *--------------------------------------------------------------------------*/

/* If create is called for HYPRE_NO_GLOBAL_PARTITION and row_starts and
   col_starts are NOT null, then it is assumed that they are array of length 2
   containing the start row of the calling processor followed by the start row
   of the next processor - AHB 6/05 */

hypre_ParCSRMatrix*
hypre_ParCSRMatrixCreate( MPI_Comm comm,
                          HYPRE_BigInt global_num_rows,
                          HYPRE_BigInt global_num_cols,
                          HYPRE_BigInt *row_starts,
                          HYPRE_BigInt *col_starts,
                          HYPRE_Int num_cols_offd,
                          HYPRE_Int num_nonzeros_diag,
                          HYPRE_Int num_nonzeros_offd )
{
   hypre_ParCSRMatrix  *matrix;
   HYPRE_Int  num_procs, my_id;
   HYPRE_Int local_num_rows, local_num_cols;
   HYPRE_BigInt first_row_index, first_col_diag;

   matrix = hypre_CTAlloc(hypre_ParCSRMatrix, 1, HYPRE_MEMORY_HOST);

   hypre_MPI_Comm_rank(comm,&my_id);
   hypre_MPI_Comm_size(comm,&num_procs);

   if (!row_starts)
   {

#ifdef HYPRE_NO_GLOBAL_PARTITION
      hypre_GenerateLocalPartitioning(global_num_rows, num_procs, my_id,
                                      &row_starts);
#else
      hypre_GeneratePartitioning(global_num_rows, num_procs, &row_starts);
#endif
   }
   if (!col_starts)
   {
      if (global_num_rows == global_num_cols)
      {
         col_starts = row_starts;
      }
      else
      {
#ifdef HYPRE_NO_GLOBAL_PARTITION
         hypre_GenerateLocalPartitioning(global_num_cols, num_procs, my_id,
                                         &col_starts);
#else
         hypre_GeneratePartitioning(global_num_cols, num_procs, &col_starts);
#endif
      }
   }

#ifdef HYPRE_NO_GLOBAL_PARTITION
   /* row_starts[0] is start of local rows.  row_starts[1] is start of next
      processor's rows */
   first_row_index = row_starts[0];
   local_num_rows = row_starts[1]-first_row_index ;
   first_col_diag = col_starts[0];
   local_num_cols = col_starts[1]-first_col_diag;
#else
   first_row_index = row_starts[my_id];
   local_num_rows = row_starts[my_id+1]-first_row_index;
   first_col_diag = col_starts[my_id];
   local_num_cols = col_starts[my_id+1]-first_col_diag;
#endif

   hypre_ParCSRMatrixComm(matrix) = comm;
   hypre_ParCSRMatrixDiag(matrix) =
      hypre_CSRMatrixCreate(local_num_rows, local_num_cols,num_nonzeros_diag);
   hypre_ParCSRMatrixOffd(matrix) =
      hypre_CSRMatrixCreate(local_num_rows, num_cols_offd,num_nonzeros_offd);
   hypre_ParCSRMatrixDiagT(matrix) = NULL;
   hypre_ParCSRMatrixOffdT(matrix) = NULL; // JSP: transposed matrices are optional
   hypre_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   hypre_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;
   hypre_ParCSRMatrixFirstRowIndex(matrix) = first_row_index;
   hypre_ParCSRMatrixFirstColDiag(matrix) = first_col_diag;

   hypre_ParCSRMatrixLastRowIndex(matrix) = first_row_index + local_num_rows - 1;
   hypre_ParCSRMatrixLastColDiag(matrix) = first_col_diag + local_num_cols - 1;

   hypre_ParCSRMatrixColMapOffd(matrix) = NULL;
   hypre_ParCSRMatrixDeviceColMapOffd(matrix) = NULL;
   hypre_ParCSRMatrixProcOrdering(matrix) = NULL;

   hypre_ParCSRMatrixAssumedPartition(matrix) = NULL;
   hypre_ParCSRMatrixOwnsAssumedPartition(matrix) = 1;

   /* When NO_GLOBAL_PARTITION is set we could make these null, instead
      of leaving the range.  If that change is made, then when this create
      is called from functions like the matrix-matrix multiply, be careful
      not to generate a new partition */

   hypre_ParCSRMatrixRowStarts(matrix) = row_starts;
   hypre_ParCSRMatrixColStarts(matrix) = col_starts;

   hypre_ParCSRMatrixCommPkg(matrix) = NULL;
   hypre_ParCSRMatrixCommPkgT(matrix) = NULL;

   /* set defaults */
   hypre_ParCSRMatrixOwnsData(matrix) = 1;
   hypre_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   hypre_ParCSRMatrixOwnsColStarts(matrix) = 1;
   if (row_starts == col_starts)
   {
      hypre_ParCSRMatrixOwnsColStarts(matrix) = 0;
   }
   hypre_ParCSRMatrixRowindices(matrix) = NULL;
   hypre_ParCSRMatrixRowvalues(matrix) = NULL;
   hypre_ParCSRMatrixGetrowactive(matrix) = 0;

   matrix->bdiaginv = NULL;
   matrix->bdiaginv_comm_pkg = NULL;
   matrix->bdiag_size = -1;

   return matrix;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixDestroy
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixDestroy( hypre_ParCSRMatrix *matrix )
{
   if (matrix)
   {
      if ( hypre_ParCSRMatrixOwnsData(matrix) )
      {
         hypre_CSRMatrixDestroy(hypre_ParCSRMatrixDiag(matrix));
         hypre_CSRMatrixDestroy(hypre_ParCSRMatrixOffd(matrix));

         if ( hypre_ParCSRMatrixDiagT(matrix) )
         {
            hypre_CSRMatrixDestroy(hypre_ParCSRMatrixDiagT(matrix));
         }

         if ( hypre_ParCSRMatrixOffdT(matrix) )
         {
            hypre_CSRMatrixDestroy(hypre_ParCSRMatrixOffdT(matrix));
         }

         if (hypre_ParCSRMatrixColMapOffd(matrix))
         {
            /*ASSERT_HOST(hypre_ParCSRMatrixColMapOffd(matrix));*/
            hypre_TFree(hypre_ParCSRMatrixColMapOffd(matrix), HYPRE_MEMORY_HOST);
         }

         if (hypre_ParCSRMatrixDeviceColMapOffd(matrix))
         {
            hypre_TFree(hypre_ParCSRMatrixDeviceColMapOffd(matrix), HYPRE_MEMORY_DEVICE);
         }

         if (hypre_ParCSRMatrixCommPkg(matrix))
         {
            hypre_MatvecCommPkgDestroy(hypre_ParCSRMatrixCommPkg(matrix));
         }

         if (hypre_ParCSRMatrixCommPkgT(matrix))
         {
            hypre_MatvecCommPkgDestroy(hypre_ParCSRMatrixCommPkgT(matrix));
         }
      }

      if ( hypre_ParCSRMatrixOwnsRowStarts(matrix) )
      {
         hypre_TFree(hypre_ParCSRMatrixRowStarts(matrix), HYPRE_MEMORY_HOST);
      }

      if ( hypre_ParCSRMatrixOwnsColStarts(matrix) )
      {
         hypre_TFree(hypre_ParCSRMatrixColStarts(matrix), HYPRE_MEMORY_HOST);
      }

      hypre_TFree(hypre_ParCSRMatrixRowindices(matrix), HYPRE_MEMORY_HOST);
      hypre_TFree(hypre_ParCSRMatrixRowvalues(matrix), HYPRE_MEMORY_HOST);

      if ( hypre_ParCSRMatrixAssumedPartition(matrix) && hypre_ParCSRMatrixOwnsAssumedPartition(matrix) )
      {
         hypre_AssumedPartitionDestroy(hypre_ParCSRMatrixAssumedPartition(matrix));
      }

      if ( hypre_ParCSRMatrixProcOrdering(matrix) )
      {
         hypre_TFree(hypre_ParCSRMatrixProcOrdering(matrix), HYPRE_MEMORY_HOST);
      }

      hypre_TFree(matrix->bdiaginv, HYPRE_MEMORY_HOST);
      if (matrix->bdiaginv_comm_pkg)
      {
         hypre_MatvecCommPkgDestroy(matrix->bdiaginv_comm_pkg);
      }

      hypre_TFree(matrix, HYPRE_MEMORY_HOST);
   }

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixInitialize
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixInitialize_v2( hypre_ParCSRMatrix *matrix, HYPRE_Int memory_location )
{
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   hypre_CSRMatrixInitialize_v2(hypre_ParCSRMatrixDiag(matrix), 0, memory_location);
   hypre_CSRMatrixInitialize_v2(hypre_ParCSRMatrixOffd(matrix), 0, memory_location);

   hypre_ParCSRMatrixColMapOffd(matrix) =
      hypre_CTAlloc(HYPRE_BigInt, hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(matrix)),
                    HYPRE_MEMORY_HOST);

   return hypre_error_flag;
}

HYPRE_Int
hypre_ParCSRMatrixInitialize( hypre_ParCSRMatrix *matrix )
{
   return hypre_ParCSRMatrixInitialize_v2(matrix, HYPRE_MEMORY_SHARED);
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixClone
 * Creates and returns a new copy S of the argument A
 * The following variables are not copied because they will be constructed
 * later if needed: CommPkg, CommPkgT, rowindices, rowvalues
 *--------------------------------------------------------------------------*/

hypre_ParCSRMatrix*
hypre_ParCSRMatrixClone_v2(hypre_ParCSRMatrix *A, HYPRE_Int copy_data, HYPRE_Int memory_location)
{
   hypre_ParCSRMatrix *S;

   S = hypre_ParCSRMatrixCreate( hypre_ParCSRMatrixComm(A),
                                 hypre_ParCSRMatrixGlobalNumRows(A),
                                 hypre_ParCSRMatrixGlobalNumCols(A),
                                 hypre_ParCSRMatrixRowStarts(A),
                                 hypre_ParCSRMatrixColStarts(A),
                                 hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(A)),
                                 hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixDiag(A)),
                                 hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixOffd(A)) );

   /* !!! S does not own Row/Col-Starts */
   hypre_ParCSRMatrixSetRowStartsOwner(S, 0);
   hypre_ParCSRMatrixSetColStartsOwner(S, 0);

   hypre_ParCSRMatrixNumNonzeros(S)  = hypre_ParCSRMatrixNumNonzeros(A);
   hypre_ParCSRMatrixDNumNonzeros(S) = hypre_ParCSRMatrixNumNonzeros(A);

   hypre_ParCSRMatrixInitialize_v2(S, memory_location);

   hypre_ParCSRMatrixCopy(A, S, copy_data);

   return S;
}

hypre_ParCSRMatrix*
hypre_ParCSRMatrixClone(hypre_ParCSRMatrix *A, HYPRE_Int copy_data)
{
   return hypre_ParCSRMatrixClone_v2(A, copy_data, HYPRE_MEMORY_SHARED);
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixSetNumNonzeros
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixSetNumNonzeros( hypre_ParCSRMatrix *matrix )
{
   MPI_Comm comm;
   hypre_CSRMatrix *diag;
   HYPRE_Int *diag_i;
   hypre_CSRMatrix *offd;
   HYPRE_Int *offd_i;
   HYPRE_Int local_num_rows;
   HYPRE_BigInt total_num_nonzeros;
   HYPRE_BigInt local_num_nonzeros;
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   comm = hypre_ParCSRMatrixComm(matrix);
   diag = hypre_ParCSRMatrixDiag(matrix);
   diag_i = hypre_CSRMatrixI(diag);
   offd = hypre_ParCSRMatrixOffd(matrix);
   offd_i = hypre_CSRMatrixI(offd);
   local_num_rows = hypre_CSRMatrixNumRows(diag);

   local_num_nonzeros = (HYPRE_BigInt)(diag_i[local_num_rows] + offd_i[local_num_rows]);
   hypre_MPI_Allreduce(&local_num_nonzeros, &total_num_nonzeros, 1, HYPRE_MPI_BIG_INT,
                       hypre_MPI_SUM, comm);
   hypre_ParCSRMatrixNumNonzeros(matrix) = total_num_nonzeros;
   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixSetDNumNonzeros
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixSetDNumNonzeros( hypre_ParCSRMatrix *matrix )
{
   MPI_Comm comm;
   hypre_CSRMatrix *diag;
   HYPRE_Int *diag_i;
   hypre_CSRMatrix *offd;
   HYPRE_Int *offd_i;
   HYPRE_Int local_num_rows;
   HYPRE_Real total_num_nonzeros;
   HYPRE_Real local_num_nonzeros;
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   comm = hypre_ParCSRMatrixComm(matrix);
   diag = hypre_ParCSRMatrixDiag(matrix);
   diag_i = hypre_CSRMatrixI(diag);
   offd = hypre_ParCSRMatrixOffd(matrix);
   offd_i = hypre_CSRMatrixI(offd);

   local_num_rows = hypre_CSRMatrixNumRows(diag);
   local_num_nonzeros  = diag_i[local_num_rows];
   local_num_nonzeros += offd_i[local_num_rows];

   hypre_MPI_Allreduce(&local_num_nonzeros, &total_num_nonzeros, 1,
                       HYPRE_MPI_REAL, hypre_MPI_SUM, comm);
   hypre_ParCSRMatrixDNumNonzeros(matrix) = total_num_nonzeros;
   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixSetDataOwner
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixSetDataOwner( hypre_ParCSRMatrix *matrix,
                                HYPRE_Int              owns_data )
{
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   hypre_ParCSRMatrixOwnsData(matrix) = owns_data;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixSetRowStartsOwner
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixSetRowStartsOwner( hypre_ParCSRMatrix *matrix,
                                     HYPRE_Int owns_row_starts )
{
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   hypre_ParCSRMatrixOwnsRowStarts(matrix) = owns_row_starts;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixSetColStartsOwner
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixSetColStartsOwner( hypre_ParCSRMatrix *matrix,
                                     HYPRE_Int owns_col_starts )
{
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   hypre_ParCSRMatrixOwnsColStarts(matrix) = owns_col_starts;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixRead
 *--------------------------------------------------------------------------*/

hypre_ParCSRMatrix *
hypre_ParCSRMatrixRead( MPI_Comm    comm,
                        const char *file_name )
{
   hypre_ParCSRMatrix  *matrix;
   hypre_CSRMatrix  *diag;
   hypre_CSRMatrix  *offd;
   HYPRE_Int  my_id, i, num_procs;
   char new_file_d[80], new_file_o[80], new_file_info[80];
   HYPRE_BigInt global_num_rows, global_num_cols;
   HYPRE_Int    num_cols_offd;
   HYPRE_Int    local_num_rows;
   HYPRE_BigInt *row_starts;
   HYPRE_BigInt *col_starts;
   HYPRE_BigInt *col_map_offd;
   FILE *fp;
   HYPRE_Int   equal = 1;
#ifdef HYPRE_NO_GLOBAL_PARTITION
   HYPRE_BigInt   row_s, row_e, col_s, col_e;
#endif

   hypre_MPI_Comm_rank(comm,&my_id);
   hypre_MPI_Comm_size(comm,&num_procs);
#ifdef HYPRE_NO_GLOBAL_PARTITION
   row_starts = hypre_CTAlloc(HYPRE_BigInt, 2, HYPRE_MEMORY_HOST);
   col_starts =  hypre_CTAlloc(HYPRE_BigInt, 2, HYPRE_MEMORY_HOST);
#else
   row_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
   col_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
#endif
   hypre_sprintf(new_file_d,"%s.D.%d",file_name,my_id);
   hypre_sprintf(new_file_o,"%s.O.%d",file_name,my_id);
   hypre_sprintf(new_file_info,"%s.INFO.%d",file_name,my_id);
   fp = fopen(new_file_info, "r");
   hypre_fscanf(fp, "%b", &global_num_rows);
   hypre_fscanf(fp, "%b", &global_num_cols);
   hypre_fscanf(fp, "%d", &num_cols_offd);
#ifdef HYPRE_NO_GLOBAL_PARTITION
   /* the bgl input file should only contain the EXACT range for local processor */
   hypre_fscanf(fp, "%d %d %d %d", &row_s, &row_e, &col_s, &col_e);
   row_starts[0] = row_s;
   row_starts[1] = row_e;
   col_starts[0] = col_s;
   col_starts[1] = col_e;
#else
   for (i=0; i < num_procs; i++)
      hypre_fscanf(fp, "%b %b", &row_starts[i], &col_starts[i]);
   row_starts[num_procs] = global_num_rows;
   col_starts[num_procs] = global_num_cols;
#endif

   col_map_offd = hypre_CTAlloc(HYPRE_BigInt, num_cols_offd, HYPRE_MEMORY_HOST);

   for (i=0; i < num_cols_offd; i++)
      hypre_fscanf(fp, "%b", &col_map_offd[i]);

   fclose(fp);

#ifdef HYPRE_NO_GLOBAL_PARTITION
   for (i=1; i >= 0; i--)
   {
      if (row_starts[i] != col_starts[i])
      {
         equal = 0;
         break;
      }
   }
#else
   for (i=num_procs; i >= 0; i--)
   {
      if (row_starts[i] != col_starts[i])
      {
         equal = 0;
         break;
      }
   }
#endif
   if (equal)
   {
      hypre_TFree(col_starts, HYPRE_MEMORY_HOST);
      col_starts = row_starts;
   }

   diag = hypre_CSRMatrixRead(new_file_d);
   local_num_rows = hypre_CSRMatrixNumRows(diag);

   if (num_cols_offd)
   {
      offd = hypre_CSRMatrixRead(new_file_o);
   }
   else
   {
      offd = hypre_CSRMatrixCreate(local_num_rows,0,0);
      hypre_CSRMatrixInitialize(offd);
   }

   matrix = hypre_CTAlloc(hypre_ParCSRMatrix, 1, HYPRE_MEMORY_HOST);

   hypre_ParCSRMatrixComm(matrix) = comm;
   hypre_ParCSRMatrixGlobalNumRows(matrix) = global_num_rows;
   hypre_ParCSRMatrixGlobalNumCols(matrix) = global_num_cols;
#ifdef HYPRE_NO_GLOBAL_PARTITION
   hypre_ParCSRMatrixFirstRowIndex(matrix) = row_s;
   hypre_ParCSRMatrixFirstColDiag(matrix) = col_s;
   hypre_ParCSRMatrixLastRowIndex(matrix) = row_e - 1;
   hypre_ParCSRMatrixLastColDiag(matrix) = col_e - 1;
#else
   hypre_ParCSRMatrixFirstRowIndex(matrix) = row_starts[my_id];
   hypre_ParCSRMatrixFirstColDiag(matrix) = col_starts[my_id];
   hypre_ParCSRMatrixLastRowIndex(matrix) = row_starts[my_id+1]-1;
   hypre_ParCSRMatrixLastColDiag(matrix) = col_starts[my_id+1]-1;
#endif

   hypre_ParCSRMatrixRowStarts(matrix) = row_starts;
   hypre_ParCSRMatrixColStarts(matrix) = col_starts;
   hypre_ParCSRMatrixCommPkg(matrix) = NULL;

   /* set defaults */
   hypre_ParCSRMatrixOwnsData(matrix) = 1;
   hypre_ParCSRMatrixOwnsRowStarts(matrix) = 1;
   hypre_ParCSRMatrixOwnsColStarts(matrix) = 1;
   if (row_starts == col_starts)
      hypre_ParCSRMatrixOwnsColStarts(matrix) = 0;

   hypre_ParCSRMatrixDiag(matrix) = diag;
   hypre_ParCSRMatrixOffd(matrix) = offd;
   if (num_cols_offd)
      hypre_ParCSRMatrixColMapOffd(matrix) = col_map_offd;
   else
      hypre_ParCSRMatrixColMapOffd(matrix) = NULL;

   return matrix;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixPrint
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixPrint( hypre_ParCSRMatrix *matrix,
                         const char         *file_name )
{
   MPI_Comm comm;
   HYPRE_BigInt global_num_rows;
   HYPRE_BigInt global_num_cols;
   HYPRE_BigInt *col_map_offd;
#ifndef HYPRE_NO_GLOBAL_PARTITION
   HYPRE_BigInt *row_starts;
   HYPRE_BigInt *col_starts;
#endif
   HYPRE_Int  my_id, i, num_procs;
   char   new_file_d[80], new_file_o[80], new_file_info[80];
   FILE *fp;
   HYPRE_Int num_cols_offd = 0;
#ifdef HYPRE_NO_GLOBAL_PARTITION
   HYPRE_BigInt row_s, row_e, col_s, col_e;
#endif
   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   comm = hypre_ParCSRMatrixComm(matrix);
   global_num_rows = hypre_ParCSRMatrixGlobalNumRows(matrix);
   global_num_cols = hypre_ParCSRMatrixGlobalNumCols(matrix);
   col_map_offd = hypre_ParCSRMatrixColMapOffd(matrix);
#ifndef HYPRE_NO_GLOBAL_PARTITION
   row_starts = hypre_ParCSRMatrixRowStarts(matrix);
   col_starts = hypre_ParCSRMatrixColStarts(matrix);
#endif
   if (hypre_ParCSRMatrixOffd(matrix))
      num_cols_offd = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(matrix));

   hypre_MPI_Comm_rank(comm, &my_id);
   hypre_MPI_Comm_size(comm, &num_procs);

   hypre_sprintf(new_file_d,"%s.D.%d",file_name,my_id);
   hypre_sprintf(new_file_o,"%s.O.%d",file_name,my_id);
   hypre_sprintf(new_file_info,"%s.INFO.%d",file_name,my_id);
   hypre_CSRMatrixPrint(hypre_ParCSRMatrixDiag(matrix),new_file_d);
   if (num_cols_offd != 0)
      hypre_CSRMatrixPrint(hypre_ParCSRMatrixOffd(matrix),new_file_o);

   fp = fopen(new_file_info, "w");
   hypre_fprintf(fp, "%b\n", global_num_rows);
   hypre_fprintf(fp, "%b\n", global_num_cols);
   hypre_fprintf(fp, "%d\n", num_cols_offd);
#ifdef HYPRE_NO_GLOBAL_PARTITION
   row_s = hypre_ParCSRMatrixFirstRowIndex(matrix);
   row_e = hypre_ParCSRMatrixLastRowIndex(matrix);
   col_s =  hypre_ParCSRMatrixFirstColDiag(matrix);
   col_e =  hypre_ParCSRMatrixLastColDiag(matrix);
   /* add 1 to the ends because this is a starts partition */
   hypre_fprintf(fp, "%b %b %b %b\n", row_s, row_e + 1, col_s, col_e + 1);
#else
   for (i=0; i < num_procs; i++)
      hypre_fprintf(fp, "%b %b\n", row_starts[i], col_starts[i]);
#endif
   for (i=0; i < num_cols_offd; i++)
      hypre_fprintf(fp, "%b\n", col_map_offd[i]);
   fclose(fp);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixPrintIJ
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixPrintIJ( const hypre_ParCSRMatrix *matrix,
                           const HYPRE_Int           base_i,
                           const HYPRE_Int           base_j,
                           const char               *filename )
{
   MPI_Comm          comm;
   HYPRE_BigInt      first_row_index;
   HYPRE_BigInt      first_col_diag;
   hypre_CSRMatrix  *diag;
   hypre_CSRMatrix  *offd;
   HYPRE_BigInt     *col_map_offd;
   HYPRE_Int         num_rows;
   HYPRE_BigInt     *row_starts;
   HYPRE_BigInt     *col_starts;
   HYPRE_Complex    *diag_data;
   HYPRE_Int        *diag_i;
   HYPRE_Int        *diag_j;
   HYPRE_Complex    *offd_data;
   HYPRE_Int        *offd_i;
   HYPRE_Int        *offd_j;
   HYPRE_Int         myid, num_procs, i, j;
   HYPRE_BigInt      I, J;
   char              new_filename[255];
   FILE             *file;
   HYPRE_Int         num_nonzeros_offd;
   HYPRE_BigInt      ilower, iupper, jlower, jupper;

   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   comm = hypre_ParCSRMatrixComm(matrix);
   first_row_index = hypre_ParCSRMatrixFirstRowIndex(matrix);
   first_col_diag  = hypre_ParCSRMatrixFirstColDiag(matrix);
   diag            = hypre_ParCSRMatrixDiag(matrix);
   offd            = hypre_ParCSRMatrixOffd(matrix);
   col_map_offd    = hypre_ParCSRMatrixColMapOffd(matrix);
   num_rows        = hypre_ParCSRMatrixNumRows(matrix);
   row_starts      = hypre_ParCSRMatrixRowStarts(matrix);
   col_starts      = hypre_ParCSRMatrixColStarts(matrix);
   hypre_MPI_Comm_rank(comm, &myid);
   hypre_MPI_Comm_size(comm, &num_procs);

   hypre_sprintf(new_filename,"%s.%05d", filename, myid);

   if ((file = fopen(new_filename, "w")) == NULL)
   {
      hypre_error_w_msg(HYPRE_ERROR_GENERIC,"Error: can't open output file %s\n");
      return hypre_error_flag;
   }

   num_nonzeros_offd = hypre_CSRMatrixNumNonzeros(offd);

   diag_data = hypre_CSRMatrixData(diag);
   diag_i    = hypre_CSRMatrixI(diag);
   diag_j    = hypre_CSRMatrixJ(diag);
   offd_i    = hypre_CSRMatrixI(offd);
   if (num_nonzeros_offd)
   {
      offd_data = hypre_CSRMatrixData(offd);
      offd_j    = hypre_CSRMatrixJ(offd);
   }

#ifdef HYPRE_NO_GLOBAL_PARTITION
   ilower = row_starts[0]+(HYPRE_BigInt)base_i;
   iupper = row_starts[1]+(HYPRE_BigInt)base_i - 1;
   jlower = col_starts[0]+(HYPRE_BigInt)base_j;
   jupper = col_starts[1]+(HYPRE_BigInt)base_j - 1;
#else
   ilower = row_starts[myid]  +(HYPRE_BigInt)base_i;
   iupper = row_starts[myid+1]+(HYPRE_BigInt)base_i - 1;
   jlower = col_starts[myid]  +(HYPRE_BigInt)base_j;
   jupper = col_starts[myid+1]+(HYPRE_BigInt)base_j - 1;
#endif
   hypre_fprintf(file, "%b %b %b %b\n", ilower, iupper, jlower, jupper);

   for (i = 0; i < num_rows; i++)
   {
      I = first_row_index + (HYPRE_BigInt)(i + base_i);

      /* print diag columns */
      for (j = diag_i[i]; j < diag_i[i+1]; j++)
      {
         J = first_col_diag + (HYPRE_BigInt)(diag_j[j] + base_j);
         if ( diag_data )
         {
#ifdef HYPRE_COMPLEX
            hypre_fprintf(file, "%b %b %.14e , %.14e\n", I, J,
                          hypre_creal(diag_data[j]), hypre_cimag(diag_data[j]));
#else
            hypre_fprintf(file, "%b %b %.14e\n", I, J, diag_data[j]);
#endif
         }
         else
            hypre_fprintf(file, "%b %b\n", I, J);
      }

      /* print offd columns */
      if ( num_nonzeros_offd )
      {
         for (j = offd_i[i]; j < offd_i[i+1]; j++)
         {
            J = col_map_offd[offd_j[j]] + (HYPRE_BigInt)base_j;
            if ( offd_data )
            {
#ifdef HYPRE_COMPLEX
               hypre_fprintf(file, "%b %b %.14e , %.14e\n", I, J,
                             hypre_creal(offd_data[j]), hypre_cimag(offd_data[j]));
#else
               hypre_fprintf(file, "%b %b %.14e\n", I, J, offd_data[j]);
#endif
            }
            else
               hypre_fprintf(file, "%b %b\n", I, J );
         }
      }
   }

   fclose(file);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixReadIJ
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixReadIJ( MPI_Comm             comm,
                          const char          *filename,
                          HYPRE_Int           *base_i_ptr,
                          HYPRE_Int           *base_j_ptr,
                          hypre_ParCSRMatrix **matrix_ptr)
{
   HYPRE_BigInt        global_num_rows;
   HYPRE_BigInt        global_num_cols;
   HYPRE_BigInt        first_row_index;
   HYPRE_BigInt        first_col_diag;
   HYPRE_BigInt        last_col_diag;
   hypre_ParCSRMatrix *matrix;
   hypre_CSRMatrix    *diag;
   hypre_CSRMatrix    *offd;
   HYPRE_BigInt       *col_map_offd;
   HYPRE_BigInt       *row_starts;
   HYPRE_BigInt       *col_starts;
   HYPRE_Int           num_rows;
   HYPRE_BigInt        big_base_i, big_base_j;
   HYPRE_Int           base_i, base_j;
   HYPRE_Complex      *diag_data;
   HYPRE_Int          *diag_i;
   HYPRE_Int          *diag_j;
   HYPRE_Complex      *offd_data;
   HYPRE_Int          *offd_i;
   HYPRE_Int          *offd_j;
   HYPRE_BigInt       *tmp_j;
   HYPRE_BigInt       *aux_offd_j;
   HYPRE_BigInt        I, J;
   HYPRE_Int           myid, num_procs, i, i2, j;
   char                new_filename[255];
   FILE               *file;
   HYPRE_Int           num_cols_offd, num_nonzeros_diag, num_nonzeros_offd;
   HYPRE_Int           equal, i_col, num_cols;
   HYPRE_Int           diag_cnt, offd_cnt, row_cnt;
   HYPRE_Complex       data;

   hypre_MPI_Comm_size(comm, &num_procs);
   hypre_MPI_Comm_rank(comm, &myid);

   hypre_sprintf(new_filename,"%s.%05d", filename, myid);

   if ((file = fopen(new_filename, "r")) == NULL)
   {
      hypre_error_w_msg(HYPRE_ERROR_GENERIC,"Error: can't open output file %s\n");
      return hypre_error_flag;
   }

   hypre_fscanf(file, "%b %b", &global_num_rows, &global_num_cols);
   hypre_fscanf(file, "%d %d %d", &num_rows, &num_cols, &num_cols_offd);
   hypre_fscanf(file, "%d %d", &num_nonzeros_diag, &num_nonzeros_offd);

   row_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
   col_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);

   for (i = 0; i <= num_procs; i++)
      hypre_fscanf(file, "%b %b", &row_starts[i], &col_starts[i]);

   big_base_i = row_starts[0];
   big_base_j = col_starts[0];
   base_i = (HYPRE_Int)row_starts[0];
   base_j = (HYPRE_Int)col_starts[0];

   equal = 1;
   for (i = 0; i <= num_procs; i++)
   {
      row_starts[i] -= big_base_i;
      col_starts[i] -= big_base_j;
      if (row_starts[i] != col_starts[i]) equal = 0;
   }

   if (equal)
   {
      hypre_TFree(col_starts, HYPRE_MEMORY_HOST);
      col_starts = row_starts;
   }
   matrix = hypre_ParCSRMatrixCreate(comm, global_num_rows, global_num_cols,
                                     row_starts, col_starts, num_cols_offd,
                                     num_nonzeros_diag, num_nonzeros_offd);
   hypre_ParCSRMatrixInitialize(matrix);

   diag = hypre_ParCSRMatrixDiag(matrix);
   offd = hypre_ParCSRMatrixOffd(matrix);

   diag_data = hypre_CSRMatrixData(diag);
   diag_i    = hypre_CSRMatrixI(diag);
   diag_j    = hypre_CSRMatrixJ(diag);

   offd_i    = hypre_CSRMatrixI(offd);
   if (num_nonzeros_offd)
   {
      offd_data = hypre_CSRMatrixData(offd);
      offd_j    = hypre_CSRMatrixJ(offd);
      tmp_j     = hypre_CTAlloc(HYPRE_BigInt, num_nonzeros_offd, HYPRE_MEMORY_HOST);
   }

   first_row_index = hypre_ParCSRMatrixFirstRowIndex(matrix);
   first_col_diag = hypre_ParCSRMatrixFirstColDiag(matrix);
   last_col_diag = first_col_diag+(HYPRE_BigInt)num_cols-1;

   diag_cnt = 0;
   offd_cnt = 0;
   row_cnt = 0;
   for (i = 0; i < num_nonzeros_diag+num_nonzeros_offd; i++)
   {
      /* read values */
      hypre_fscanf(file, "%b %b %le", &I, &J, &data);
      i2 = (HYPRE_Int)(I-big_base_i-first_row_index);
      J -= big_base_j;
      if (i2 > row_cnt)
      {
         diag_i[i2] = diag_cnt;
         offd_i[i2] = offd_cnt;
         row_cnt++;
      }
      if (J < first_col_diag || J > last_col_diag)
      {
         tmp_j[offd_cnt] = J;
         offd_data[offd_cnt++] = data;
      }
      else
      {
         diag_j[diag_cnt] = (HYPRE_Int)(J - first_col_diag);
         diag_data[diag_cnt++] = data;
      }
   }
   diag_i[num_rows] = diag_cnt;
   offd_i[num_rows] = offd_cnt;

   fclose(file);

   /*  generate col_map_offd */
   if (num_nonzeros_offd)
   {
      aux_offd_j = hypre_CTAlloc(HYPRE_BigInt, num_nonzeros_offd, HYPRE_MEMORY_HOST);
      for (i=0; i < num_nonzeros_offd; i++)
         aux_offd_j[i] = (HYPRE_BigInt)offd_j[i];
      hypre_BigQsort0(aux_offd_j,0,num_nonzeros_offd-1);
      col_map_offd = hypre_ParCSRMatrixColMapOffd(matrix);
      col_map_offd[0] = aux_offd_j[0];
      offd_cnt = 0;
      for (i=1; i < num_nonzeros_offd; i++)
      {
         if (aux_offd_j[i] > col_map_offd[offd_cnt])
            col_map_offd[++offd_cnt] = aux_offd_j[i];
      }
      for (i=0; i < num_nonzeros_offd; i++)
      {
         offd_j[i] = hypre_BigBinarySearch(col_map_offd, tmp_j[i], num_cols_offd);
      }
      hypre_TFree(aux_offd_j, HYPRE_MEMORY_HOST);
      hypre_TFree(tmp_j, HYPRE_MEMORY_HOST);
   }

   /* move diagonal element in first position in each row */
   for (i=0; i < num_rows; i++)
   {
      i_col = diag_i[i];
      for (j=i_col; j < diag_i[i+1]; j++)
      {
         if (diag_j[j] == i)
         {
            diag_j[j] = diag_j[i_col];
            data = diag_data[j];
            diag_data[j] = diag_data[i_col];
            diag_data[i_col] = data;
            diag_j[i_col] = i;
            break;
         }
      }
   }

   *base_i_ptr = base_i;
   *base_j_ptr = base_j;
   *matrix_ptr = matrix;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixGetLocalRange
 * returns the row numbers of the rows stored on this processor.
 * "End" is actually the row number of the last row on this processor.
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixGetLocalRange( hypre_ParCSRMatrix *matrix,
                                 HYPRE_BigInt       *row_start,
                                 HYPRE_BigInt       *row_end,
                                 HYPRE_BigInt       *col_start,
                                 HYPRE_BigInt       *col_end )
{
   HYPRE_Int my_id;

   if (!matrix)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   hypre_MPI_Comm_rank( hypre_ParCSRMatrixComm(matrix), &my_id );

#ifdef HYPRE_NO_GLOBAL_PARTITION
   *row_start = hypre_ParCSRMatrixFirstRowIndex(matrix);
   *row_end = hypre_ParCSRMatrixLastRowIndex(matrix);
   *col_start =  hypre_ParCSRMatrixFirstColDiag(matrix);
   *col_end =  hypre_ParCSRMatrixLastColDiag(matrix);
#else
   *row_start = hypre_ParCSRMatrixRowStarts(matrix)[ my_id ];
   *row_end = hypre_ParCSRMatrixRowStarts(matrix)[ my_id + 1 ]-1;
   *col_start = hypre_ParCSRMatrixColStarts(matrix)[ my_id ];
   *col_end = hypre_ParCSRMatrixColStarts(matrix)[ my_id + 1 ]-1;
#endif

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixGetRow
 * Returns global column indices and/or values for a given row in the global
 * matrix. Global row number is used, but the row must be stored locally or
 * an error is returned. This implementation copies from the two matrices that
 * store the local data, storing them in the hypre_ParCSRMatrix structure.
 * Only a single row can be accessed via this function at any one time; the
 * corresponding RestoreRow function must be called, to avoid bleeding memory,
 * and to be able to look at another row.
 * Either one of col_ind and values can be left null, and those values will
 * not be returned.
 * All indices are returned in 0-based indexing, no matter what is used under
 * the hood. EXCEPTION: currently this only works if the local CSR matrices
 * use 0-based indexing.
 * This code, semantics, implementation, etc., are all based on PETSc's hypre_MPI_AIJ
 * matrix code, adjusted for our data and software structures.
 * AJC 4/99.
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixGetRow( hypre_ParCSRMatrix  *mat,
                          HYPRE_BigInt         row,
                          HYPRE_Int           *size,
                          HYPRE_BigInt       **col_ind,
                          HYPRE_Complex      **values )
{
   HYPRE_Int my_id;
   HYPRE_BigInt row_start, row_end;
   hypre_CSRMatrix *Aa;
   hypre_CSRMatrix *Ba;

   if (!mat)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   Aa = (hypre_CSRMatrix *) hypre_ParCSRMatrixDiag(mat);
   Ba = (hypre_CSRMatrix *) hypre_ParCSRMatrixOffd(mat);

   if (hypre_ParCSRMatrixGetrowactive(mat)) return(-1);

   hypre_MPI_Comm_rank( hypre_ParCSRMatrixComm(mat), &my_id );

   hypre_ParCSRMatrixGetrowactive(mat) = 1;
#ifdef HYPRE_NO_GLOBAL_PARTITION
   row_start = hypre_ParCSRMatrixFirstRowIndex(mat);
   row_end = hypre_ParCSRMatrixLastRowIndex(mat) + 1;
#else
   row_end = hypre_ParCSRMatrixRowStarts(mat)[ my_id + 1 ];
   row_start = hypre_ParCSRMatrixRowStarts(mat)[ my_id ];
#endif
   if (row < row_start || row >= row_end) return(-1);

   /* if buffer is not allocated and some information is requested,
      allocate buffer */
   if (!hypre_ParCSRMatrixRowvalues(mat) && ( col_ind || values ))
   {
      /*
        allocate enough space to hold information from the longest row.
      */
      HYPRE_Int     max = 1,tmp;
      HYPRE_Int i;
      HYPRE_Int     m = row_end-row_start;

      for ( i=0; i<m; i++ ) {
         tmp = hypre_CSRMatrixI(Aa)[i+1] - hypre_CSRMatrixI(Aa)[i] +
            hypre_CSRMatrixI(Ba)[i+1] - hypre_CSRMatrixI(Ba)[i];
         if (max < tmp) { max = tmp; }
      }

      hypre_ParCSRMatrixRowvalues(mat) = (HYPRE_Complex *) hypre_CTAlloc( HYPRE_Complex, max , HYPRE_MEMORY_HOST);
      hypre_ParCSRMatrixRowindices(mat) = (HYPRE_BigInt *) hypre_CTAlloc( HYPRE_BigInt, max , HYPRE_MEMORY_HOST);
   }

   /* Copy from dual sequential matrices into buffer */
   {
      HYPRE_Complex     *vworkA, *vworkB, *v_p;
      HYPRE_Int        i, *cworkA, *cworkB;
      HYPRE_BigInt     cstart = hypre_ParCSRMatrixFirstColDiag(mat);
      HYPRE_Int        nztot, nzA, nzB, lrow=(HYPRE_Int)(row-row_start);
      HYPRE_BigInt     *cmap, *idx_p;

      nzA = hypre_CSRMatrixI(Aa)[lrow+1]-hypre_CSRMatrixI(Aa)[lrow];
      cworkA = &( hypre_CSRMatrixJ(Aa)[ hypre_CSRMatrixI(Aa)[lrow] ] );
      vworkA = &( hypre_CSRMatrixData(Aa)[ hypre_CSRMatrixI(Aa)[lrow] ] );

      nzB = hypre_CSRMatrixI(Ba)[lrow+1]-hypre_CSRMatrixI(Ba)[lrow];
      cworkB = &( hypre_CSRMatrixJ(Ba)[ hypre_CSRMatrixI(Ba)[lrow] ] );
      vworkB = &( hypre_CSRMatrixData(Ba)[ hypre_CSRMatrixI(Ba)[lrow] ] );

      nztot = nzA + nzB;

      cmap  = hypre_ParCSRMatrixColMapOffd(mat);

      if (values  || col_ind) {
         if (nztot) {
            /* Sort by increasing column numbers, assuming A and B already sorted */
            HYPRE_Int imark = -1;
            if (values) {
               *values = v_p = hypre_ParCSRMatrixRowvalues(mat);
               for ( i=0; i<nzB; i++ ) {
                  if (cmap[cworkB[i]] < cstart)   v_p[i] = vworkB[i];
                  else break;
               }
               imark = i;
               for ( i=0; i<nzA; i++ )     v_p[imark+i] = vworkA[i];
               for ( i=imark; i<nzB; i++ ) v_p[nzA+i]   = vworkB[i];
            }
            if (col_ind) {
               *col_ind = idx_p = hypre_ParCSRMatrixRowindices(mat);
               if (imark > -1) {
                  for ( i=0; i<imark; i++ ) {
                     idx_p[i] = cmap[cworkB[i]];
                  }
               } else {
                  for ( i=0; i<nzB; i++ ) {
                     if (cmap[cworkB[i]] < cstart)   idx_p[i] = cmap[cworkB[i]];
                     else break;
                  }
                  imark = i;
               }
               for ( i=0; i<nzA; i++ )     idx_p[imark+i] = cstart + cworkA[i];
               for ( i=imark; i<nzB; i++ ) idx_p[nzA+i]   = cmap[cworkB[i]];
            }
         }
         else {
            if (col_ind) *col_ind = 0;
            if (values)   *values   = 0;
         }
      }
      *size = nztot;

   } /* End of copy */

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixRestoreRow
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixRestoreRow( hypre_ParCSRMatrix *matrix,
                              HYPRE_BigInt        row,
                              HYPRE_Int          *size,
                              HYPRE_BigInt      **col_ind,
                              HYPRE_Complex     **values )
{
   if (!hypre_ParCSRMatrixGetrowactive(matrix))
   {
      hypre_error(HYPRE_ERROR_GENERIC);
      return hypre_error_flag;
   }

   hypre_ParCSRMatrixGetrowactive(matrix)=0;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_CSRMatrixToParCSRMatrix:
 * generates a ParCSRMatrix distributed across the processors in comm
 * from a CSRMatrix on proc 0 .
 *
 * This shouldn't be used with the HYPRE_NO_GLOBAL_PARTITON option
 *
 *--------------------------------------------------------------------------*/

hypre_ParCSRMatrix *
hypre_CSRMatrixToParCSRMatrix( MPI_Comm         comm,
                               hypre_CSRMatrix *A,
                               HYPRE_BigInt    *row_starts,
                               HYPRE_BigInt    *col_starts )
{
   HYPRE_BigInt       *global_data;
   HYPRE_BigInt        global_size;
   HYPRE_BigInt        global_num_rows;
   HYPRE_BigInt        global_num_cols;
   HYPRE_Int          *local_num_rows;

   HYPRE_Int           num_procs, my_id;
   HYPRE_Int          *local_num_nonzeros=NULL;
   HYPRE_Int           num_nonzeros;

   HYPRE_Complex      *a_data;
   HYPRE_Int          *a_i;
   HYPRE_Int          *a_j;

   hypre_CSRMatrix    *local_A;

   hypre_MPI_Request  *requests;
   hypre_MPI_Status   *status, status0;
   hypre_MPI_Datatype *csr_matrix_datatypes;

   hypre_ParCSRMatrix *par_matrix;

   HYPRE_BigInt        first_col_diag;
   HYPRE_BigInt        last_col_diag;
   HYPRE_Int           i, j, ind;

   hypre_MPI_Comm_rank(comm, &my_id);
   hypre_MPI_Comm_size(comm, &num_procs);

   global_data = hypre_CTAlloc(HYPRE_BigInt, 2*num_procs+6, HYPRE_MEMORY_HOST);
   if (my_id == 0)
   {
      global_size = 3;
      if (row_starts)
      {
         if (col_starts)
         {
            if (col_starts != row_starts)
            {
               /* contains code for what to expect,
                  if 0:  row_starts = col_starts, only row_starts given
                  if 1: only row_starts given, col_starts = NULL
                  if 2: both row_starts and col_starts given
                  if 3: only col_starts given, row_starts = NULL */
               global_data[3] = 2;
               global_size = (HYPRE_BigInt)(2*num_procs+6);
               for (i=0; i < num_procs+1; i++)
                  global_data[i+4] = row_starts[i];
               for (i=0; i < num_procs+1; i++)
                  global_data[i+num_procs+5] = col_starts[i];
            }
            else
            {
               global_data[3] = 0;
               global_size = (HYPRE_BigInt)num_procs+5;
               for (i=0; i < num_procs+1; i++)
                  global_data[i+4] = row_starts[i];
            }
         }
         else
         {
            global_data[3] = 1;
            global_size = (HYPRE_BigInt)num_procs+5;
            for (i=0; i < num_procs+1; i++)
               global_data[i+4] = row_starts[i];
         }
      }
      else
      {
         if (col_starts)
         {
            global_data[3] = 3;
            global_size = (HYPRE_BigInt)num_procs+5;
            for (i=0; i < num_procs+1; i++)
               global_data[i+4] = col_starts[i];
         }
      }
      global_data[0] = (HYPRE_BigInt)hypre_CSRMatrixNumRows(A);
      global_data[1] = (HYPRE_BigInt)hypre_CSRMatrixNumCols(A);
      global_data[2] = global_size;
      a_data = hypre_CSRMatrixData(A);
      a_i = hypre_CSRMatrixI(A);
      a_j = hypre_CSRMatrixJ(A);
   }
   hypre_MPI_Bcast(global_data,3,HYPRE_MPI_BIG_INT,0,comm);
   global_num_rows = global_data[0];
   global_num_cols = global_data[1];

   global_size = global_data[2];
   if (global_size > 3)
   {
      hypre_MPI_Bcast(&global_data[3],global_size-3,HYPRE_MPI_BIG_INT,0,comm);
      if (my_id > 0)
      {
         if (global_data[3] < 3)
         {
            row_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
            for (i=0; i< num_procs+1; i++)
            {
               row_starts[i] = global_data[i+4];
            }
            if (global_data[3] == 0)
               col_starts = row_starts;
            if (global_data[3] == 2)
            {
               col_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
               for (i=0; i < num_procs+1; i++)
               {
                  col_starts[i] = global_data[i+num_procs+5];
               }
            }
         }
         else
         {
            col_starts = hypre_CTAlloc(HYPRE_BigInt, num_procs+1, HYPRE_MEMORY_HOST);
            for (i=0; i< num_procs+1; i++)
            {
               col_starts[i] = global_data[i+4];
            }
         }
      }
   }
   hypre_TFree(global_data, HYPRE_MEMORY_HOST);

   local_num_rows = hypre_CTAlloc(HYPRE_Int,  num_procs, HYPRE_MEMORY_HOST);
   csr_matrix_datatypes = hypre_CTAlloc(hypre_MPI_Datatype,  num_procs, HYPRE_MEMORY_HOST);

   par_matrix = hypre_ParCSRMatrixCreate(
      comm, global_num_rows, global_num_cols,row_starts,col_starts,0,0,0);

   row_starts = hypre_ParCSRMatrixRowStarts(par_matrix);
   col_starts = hypre_ParCSRMatrixColStarts(par_matrix);

   for (i=0; i < num_procs; i++)
      local_num_rows[i] = (HYPRE_Int)(row_starts[i+1] - row_starts[i]);

   if (my_id == 0)
   {
      local_num_nonzeros = hypre_CTAlloc(HYPRE_Int,  num_procs, HYPRE_MEMORY_HOST);
      for (i=0; i < num_procs-1; i++)
         local_num_nonzeros[i] = a_i[(HYPRE_Int)row_starts[i+1]]
            - a_i[(HYPRE_Int)row_starts[i]];
      local_num_nonzeros[num_procs-1] = a_i[(HYPRE_Int)global_num_rows]
         - a_i[(HYPRE_Int)row_starts[num_procs-1]];
   }
   hypre_MPI_Scatter(local_num_nonzeros,1,HYPRE_MPI_INT,&num_nonzeros,1,
                     HYPRE_MPI_INT,0,comm);

   if (my_id == 0) num_nonzeros = local_num_nonzeros[0];

   local_A = hypre_CSRMatrixCreate(local_num_rows[my_id], (HYPRE_Int)global_num_cols,
                                   num_nonzeros);
   if (my_id == 0)
   {
      requests = hypre_CTAlloc(hypre_MPI_Request,  num_procs-1, HYPRE_MEMORY_HOST);
      status = hypre_CTAlloc(hypre_MPI_Status,  num_procs-1, HYPRE_MEMORY_HOST);
      j=0;
      for (i=1; i < num_procs; i++)
      {
         ind = a_i[(HYPRE_Int)row_starts[i]];
         hypre_BuildCSRMatrixMPIDataType(local_num_nonzeros[i],
                                         local_num_rows[i],
                                         &a_data[ind],
                                         &a_i[(HYPRE_Int)row_starts[i]],
                                         &a_j[ind],
                                         &csr_matrix_datatypes[i]);
         hypre_MPI_Isend(hypre_MPI_BOTTOM, 1, csr_matrix_datatypes[i], i, 0, comm,
                         &requests[j++]);
         hypre_MPI_Type_free(&csr_matrix_datatypes[i]);
      }
      hypre_CSRMatrixData(local_A) = a_data;
      hypre_CSRMatrixI(local_A) = a_i;
      hypre_CSRMatrixJ(local_A) = a_j;
      hypre_CSRMatrixOwnsData(local_A) = 0;
      hypre_MPI_Waitall(num_procs-1,requests,status);
      hypre_TFree(requests, HYPRE_MEMORY_HOST);
      hypre_TFree(status, HYPRE_MEMORY_HOST);
      hypre_TFree(local_num_nonzeros, HYPRE_MEMORY_HOST);
   }
   else
   {
      hypre_CSRMatrixInitialize(local_A);
      hypre_BuildCSRMatrixMPIDataType(num_nonzeros,
                                      local_num_rows[my_id],
                                      hypre_CSRMatrixData(local_A),
                                      hypre_CSRMatrixI(local_A),
                                      hypre_CSRMatrixJ(local_A),
                                      csr_matrix_datatypes);
      hypre_MPI_Recv(hypre_MPI_BOTTOM,1,csr_matrix_datatypes[0],0,0,comm,&status0);
      hypre_MPI_Type_free(csr_matrix_datatypes);
   }

   first_col_diag = col_starts[my_id];
   last_col_diag = col_starts[my_id+1]-1;

   GenerateDiagAndOffd(local_A, par_matrix, first_col_diag, last_col_diag);

   /* set pointers back to NULL before destroying */
   if (my_id == 0)
   {
      hypre_CSRMatrixData(local_A) = NULL;
      hypre_CSRMatrixI(local_A) = NULL;
      hypre_CSRMatrixJ(local_A) = NULL;
   }
   hypre_CSRMatrixDestroy(local_A);
   hypre_TFree(local_num_rows, HYPRE_MEMORY_HOST);
   hypre_TFree(csr_matrix_datatypes, HYPRE_MEMORY_HOST);

   return par_matrix;
}

HYPRE_Int
GenerateDiagAndOffd(hypre_CSRMatrix *A,
                    hypre_ParCSRMatrix *matrix,
                    HYPRE_BigInt first_col_diag,
                    HYPRE_BigInt last_col_diag)
{
   HYPRE_Int  i, j;
   HYPRE_Int  jo, jd;
   HYPRE_Int  num_rows = hypre_CSRMatrixNumRows(A);
   HYPRE_Int  num_cols = hypre_CSRMatrixNumCols(A);
   HYPRE_Complex *a_data = hypre_CSRMatrixData(A);
   HYPRE_Int *a_i = hypre_CSRMatrixI(A);
   HYPRE_Int *a_j = hypre_CSRMatrixJ(A);

   hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(matrix);
   hypre_CSRMatrix *offd = hypre_ParCSRMatrixOffd(matrix);

   HYPRE_BigInt  *col_map_offd;

   HYPRE_Complex *diag_data, *offd_data;
   HYPRE_Int  *diag_i, *offd_i;
   HYPRE_Int  *diag_j, *offd_j;
   HYPRE_Int  *marker;
   HYPRE_Int num_cols_diag, num_cols_offd;
   HYPRE_Int first_elmt = a_i[0];
   HYPRE_Int num_nonzeros = a_i[num_rows]-first_elmt;
   HYPRE_Int counter;

   num_cols_diag = (HYPRE_Int)(last_col_diag - first_col_diag +1);
   num_cols_offd = 0;

   if (num_cols - num_cols_diag)
   {
      hypre_CSRMatrixInitialize(diag);
      diag_i = hypre_CSRMatrixI(diag);

      hypre_CSRMatrixInitialize(offd);
      offd_i = hypre_CSRMatrixI(offd);
      marker = hypre_CTAlloc(HYPRE_Int, num_cols, HYPRE_MEMORY_HOST);

      for (i=0; i < num_cols; i++)
         marker[i] = 0;

      jo = 0;
      jd = 0;
      for (i=0; i < num_rows; i++)
      {
         offd_i[i] = jo;
         diag_i[i] = jd;

         for (j=a_i[i]-first_elmt; j < a_i[i+1]-first_elmt; j++)
            if (a_j[j] < first_col_diag || a_j[j] > last_col_diag)
            {
               if (!marker[a_j[j]])
               {
                  marker[a_j[j]] = 1;
                  num_cols_offd++;
               }
               jo++;
            }
            else
            {
               jd++;
            }
      }
      offd_i[num_rows] = jo;
      diag_i[num_rows] = jd;

      hypre_ParCSRMatrixColMapOffd(matrix) = hypre_CTAlloc(HYPRE_BigInt, num_cols_offd, HYPRE_MEMORY_HOST);
      col_map_offd = hypre_ParCSRMatrixColMapOffd(matrix);

      counter = 0;
      for (i=0; i < num_cols; i++)
         if (marker[i])
         {
            col_map_offd[counter] = (HYPRE_BigInt)i;
            marker[i] = counter;
            counter++;
         }

      hypre_CSRMatrixNumNonzeros(diag) = jd;
      hypre_CSRMatrixInitialize(diag);
      diag_data = hypre_CSRMatrixData(diag);
      diag_j = hypre_CSRMatrixJ(diag);

      hypre_CSRMatrixNumNonzeros(offd) = jo;
      hypre_CSRMatrixNumCols(offd) = num_cols_offd;
      hypre_CSRMatrixInitialize(offd);
      offd_data = hypre_CSRMatrixData(offd);
      offd_j = hypre_CSRMatrixJ(offd);

      jo = 0;
      jd = 0;
      for (i=0; i < num_rows; i++)
      {
         for (j=a_i[i]-first_elmt; j < a_i[i+1]-first_elmt; j++)
            if (a_j[j] < (HYPRE_Int)first_col_diag || a_j[j] > (HYPRE_Int)last_col_diag)
            {
               offd_data[jo] = a_data[j];
               offd_j[jo++] = marker[a_j[j]];
            }
            else
            {
               diag_data[jd] = a_data[j];
               diag_j[jd++] = (HYPRE_Int)(a_j[j]-first_col_diag);
            }
      }
      hypre_TFree(marker, HYPRE_MEMORY_HOST);
   }
   else
   {
      hypre_CSRMatrixNumNonzeros(diag) = num_nonzeros;
      hypre_CSRMatrixInitialize(diag);
      diag_data = hypre_CSRMatrixData(diag);
      diag_i = hypre_CSRMatrixI(diag);
      diag_j = hypre_CSRMatrixJ(diag);

      for (i=0; i < num_nonzeros; i++)
      {
         diag_data[i] = a_data[i];
         diag_j[i] = a_j[i];
      }
      offd_i = hypre_CTAlloc(HYPRE_Int,  num_rows+1, HYPRE_MEMORY_HOST);

      for (i=0; i < num_rows+1; i++)
      {
         diag_i[i] = a_i[i];
         offd_i[i] = 0;
      }

      hypre_CSRMatrixNumCols(offd) = 0;
      hypre_CSRMatrixI(offd) = offd_i;
   }

   return hypre_error_flag;
}

hypre_CSRMatrix *
hypre_MergeDiagAndOffd(hypre_ParCSRMatrix *par_matrix)
{
   hypre_CSRMatrix  *diag = hypre_ParCSRMatrixDiag(par_matrix);
   hypre_CSRMatrix  *offd = hypre_ParCSRMatrixOffd(par_matrix);
   hypre_CSRMatrix  *matrix;

   HYPRE_BigInt       num_cols = hypre_ParCSRMatrixGlobalNumCols(par_matrix);
   HYPRE_BigInt       first_col_diag = hypre_ParCSRMatrixFirstColDiag(par_matrix);
   HYPRE_BigInt      *col_map_offd = hypre_ParCSRMatrixColMapOffd(par_matrix);
   HYPRE_Int          num_rows = hypre_CSRMatrixNumRows(diag);

   HYPRE_Int          *diag_i = hypre_CSRMatrixI(diag);
   HYPRE_Int          *diag_j = hypre_CSRMatrixJ(diag);
   HYPRE_Complex      *diag_data = hypre_CSRMatrixData(diag);
   HYPRE_Int          *offd_i = hypre_CSRMatrixI(offd);
   HYPRE_Int          *offd_j = hypre_CSRMatrixJ(offd);
   HYPRE_Complex      *offd_data = hypre_CSRMatrixData(offd);

   HYPRE_Int          *matrix_i;
   HYPRE_BigInt       *matrix_j;
   HYPRE_Complex      *matrix_data;

   HYPRE_Int          num_nonzeros, i, j;
   HYPRE_Int          count;
   HYPRE_Int          size, rest, num_threads, ii;

   num_nonzeros = diag_i[num_rows] + offd_i[num_rows];

   matrix = hypre_CSRMatrixCreate(num_rows,num_cols,num_nonzeros);
   hypre_CSRMatrixBigInitialize(matrix);

   matrix_i = hypre_CSRMatrixI(matrix);
   matrix_j = hypre_CSRMatrixBigJ(matrix);
   matrix_data = hypre_CSRMatrixData(matrix);
   num_threads = hypre_NumThreads();
   size = num_rows/num_threads;
   rest = num_rows - size*num_threads;

#ifdef HYPRE_USING_OPENMP
#pragma omp parallel for private(ii, i, j, count) HYPRE_SMP_SCHEDULE
#endif
   for (ii=0; ii < num_threads; ii++)
   {
      HYPRE_Int ns, ne;
      if (ii < rest)
      {
         ns = ii*size+ii;
         ne = (ii+1)*size+ii+1;
      }
      else
      {
         ns = ii*size+rest;
         ne = (ii+1)*size+rest;
      }
      count = diag_i[ns]+offd_i[ns];;
      for (i=ns; i < ne; i++)
      {
         matrix_i[i] = count;
         for (j=diag_i[i]; j < diag_i[i+1]; j++)
         {
            matrix_data[count] = diag_data[j];
            matrix_j[count++] = (HYPRE_BigInt)diag_j[j]+first_col_diag;
         }
         for (j=offd_i[i]; j < offd_i[i+1]; j++)
         {
            matrix_data[count] = offd_data[j];
            matrix_j[count++] = col_map_offd[offd_j[j]];
         }
      }
   } /* end parallel region */
   matrix_i[num_rows] = num_nonzeros;

   return matrix;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixToCSRMatrixAll:
 * generates a CSRMatrix from a ParCSRMatrix on all processors that have
 * parts of the ParCSRMatrix
 * Warning: this only works for a ParCSRMatrix that is smaller than 2^31-1
 *--------------------------------------------------------------------------*/

hypre_CSRMatrix *
hypre_ParCSRMatrixToCSRMatrixAll(hypre_ParCSRMatrix *par_matrix)
{
   MPI_Comm comm = hypre_ParCSRMatrixComm(par_matrix);
   hypre_CSRMatrix *matrix;
   hypre_CSRMatrix *local_matrix;
   HYPRE_Int num_rows = (HYPRE_Int)hypre_ParCSRMatrixGlobalNumRows(par_matrix);
   HYPRE_Int num_cols = (HYPRE_Int)hypre_ParCSRMatrixGlobalNumCols(par_matrix);
#ifndef HYPRE_NO_GLOBAL_PARTITION
   HYPRE_BigInt *row_starts = hypre_ParCSRMatrixRowStarts(par_matrix);
#endif
   HYPRE_Int *matrix_i;
   HYPRE_Int *matrix_j;
   HYPRE_Complex *matrix_data;

   HYPRE_Int *local_matrix_i;
   HYPRE_Int *local_matrix_j;
   HYPRE_Complex *local_matrix_data;

   HYPRE_Int i, j;
   HYPRE_Int local_num_rows;
   HYPRE_Int local_num_nonzeros;
   HYPRE_Int num_nonzeros;
   HYPRE_Int num_data;
   HYPRE_Int num_requests;
   HYPRE_Int vec_len, offset;
   HYPRE_Int start_index;
   HYPRE_Int proc_id;
   HYPRE_Int num_procs, my_id;
   HYPRE_Int num_types;
   HYPRE_Int *used_procs;

   hypre_MPI_Request *requests;
   hypre_MPI_Status *status;

#ifdef HYPRE_NO_GLOBAL_PARTITION

   HYPRE_Int *new_vec_starts;

   HYPRE_Int num_contacts;
   HYPRE_Int contact_proc_list[1];
   HYPRE_Int contact_send_buf[1];
   HYPRE_Int contact_send_buf_starts[2];
   HYPRE_Int max_response_size;
   HYPRE_Int *response_recv_buf=NULL;
   HYPRE_Int *response_recv_buf_starts = NULL;
   hypre_DataExchangeResponse response_obj;
   hypre_ProcListElements send_proc_obj;

   HYPRE_Int *send_info = NULL;
   hypre_MPI_Status  status1;
   HYPRE_Int count, tag1 = 11112, tag2 = 22223, tag3 = 33334;
   HYPRE_Int start;

#endif

   hypre_MPI_Comm_size(comm, &num_procs);
   hypre_MPI_Comm_rank(comm, &my_id);

#ifdef HYPRE_NO_GLOBAL_PARTITION

   local_num_rows = (HYPRE_Int)(hypre_ParCSRMatrixLastRowIndex(par_matrix)  -
      hypre_ParCSRMatrixFirstRowIndex(par_matrix) + 1);


   local_matrix = hypre_MergeDiagAndOffd(par_matrix); /* creates matrix */
   hypre_CSRMatrixBigJtoJ(local_matrix); /* copies big_j to j */
   local_matrix_i = hypre_CSRMatrixI(local_matrix);
   local_matrix_j = hypre_CSRMatrixJ(local_matrix);
   local_matrix_data = hypre_CSRMatrixData(local_matrix);


   /* determine procs that have vector data and store their ids in used_procs */
   /* we need to do an exchange data for this.  If I own row then I will contact
      processor 0 with the endpoint of my local range */

   if (local_num_rows > 0)
   {
      num_contacts = 1;
      contact_proc_list[0] = 0;
      contact_send_buf[0] =  (HYPRE_Int)hypre_ParCSRMatrixLastRowIndex(par_matrix);
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 1;

   }
   else
   {
      num_contacts = 0;
      contact_send_buf_starts[0] = 0;
      contact_send_buf_starts[1] = 0;
   }
   /*build the response object*/
   /*send_proc_obj will  be for saving info from contacts */
   send_proc_obj.length = 0;
   send_proc_obj.storage_length = 10;
   send_proc_obj.id = hypre_CTAlloc(HYPRE_Int,  send_proc_obj.storage_length, HYPRE_MEMORY_HOST);
   send_proc_obj.vec_starts =
      hypre_CTAlloc(HYPRE_Int,  send_proc_obj.storage_length + 1, HYPRE_MEMORY_HOST);
   send_proc_obj.vec_starts[0] = 0;
   send_proc_obj.element_storage_length = 10;
   send_proc_obj.elements =
      hypre_CTAlloc(HYPRE_BigInt,  send_proc_obj.element_storage_length, HYPRE_MEMORY_HOST);

   max_response_size = 0; /* each response is null */
   response_obj.fill_response = hypre_FillResponseParToCSRMatrix;
   response_obj.data1 = NULL;
   response_obj.data2 = &send_proc_obj; /*this is where we keep info from contacts*/


   hypre_DataExchangeList(num_contacts,
                          contact_proc_list, contact_send_buf,
                          contact_send_buf_starts, sizeof(HYPRE_Int),
                          sizeof(HYPRE_Int), &response_obj,
                          max_response_size, 1,
                          comm, (void**) &response_recv_buf,
                          &response_recv_buf_starts);

   /* now processor 0 should have a list of ranges for processors that have rows -
      these are in send_proc_obj - it needs to create the new list of processors
      and also an array of vec starts - and send to those who own row*/
   if (my_id)
   {
      if (local_num_rows)
      {
         /* look for a message from processor 0 */
         hypre_MPI_Probe(0, tag1, comm, &status1);
         hypre_MPI_Get_count(&status1, HYPRE_MPI_INT, &count);

         send_info = hypre_CTAlloc(HYPRE_Int,  count, HYPRE_MEMORY_HOST);
         hypre_MPI_Recv(send_info, count, HYPRE_MPI_INT, 0, tag1, comm, &status1);

         /* now unpack */
         num_types = send_info[0];
         used_procs =  hypre_CTAlloc(HYPRE_Int,  num_types, HYPRE_MEMORY_HOST);
         new_vec_starts = hypre_CTAlloc(HYPRE_Int,  num_types+1, HYPRE_MEMORY_HOST);

         for (i=1; i<= num_types; i++)
         {
            used_procs[i-1] = send_info[i];
         }
         for (i=num_types+1; i< count; i++)
         {
            new_vec_starts[i-num_types-1] = send_info[i] ;
         }
      }
      else /* clean up and exit */
      {
         hypre_TFree(send_proc_obj.vec_starts, HYPRE_MEMORY_HOST);
         hypre_TFree(send_proc_obj.id, HYPRE_MEMORY_HOST);
         hypre_TFree(send_proc_obj.elements, HYPRE_MEMORY_HOST);
         if(response_recv_buf)        hypre_TFree(response_recv_buf, HYPRE_MEMORY_HOST);
         if(response_recv_buf_starts) hypre_TFree(response_recv_buf_starts, HYPRE_MEMORY_HOST);


         if (hypre_CSRMatrixOwnsData(local_matrix))
            hypre_CSRMatrixDestroy(local_matrix);
         else
            hypre_TFree(local_matrix, HYPRE_MEMORY_HOST);


         return NULL;
      }
   }
   else /* my_id ==0 */
   {
      num_types = send_proc_obj.length;
      used_procs =  hypre_CTAlloc(HYPRE_Int,  num_types, HYPRE_MEMORY_HOST);
      new_vec_starts = hypre_CTAlloc(HYPRE_Int,  num_types+1, HYPRE_MEMORY_HOST);

      new_vec_starts[0] = 0;
      for (i=0; i< num_types; i++)
      {
         used_procs[i] = send_proc_obj.id[i];
         new_vec_starts[i+1] = send_proc_obj.elements[i]+1;
      }
      hypre_qsort0(used_procs, 0, num_types-1);
      hypre_qsort0(new_vec_starts, 0, num_types);
      /*now we need to put into an array to send */
      count =  2*num_types+2;
      send_info = hypre_CTAlloc(HYPRE_Int,  count, HYPRE_MEMORY_HOST);
      send_info[0] = num_types;
      for (i=1; i<= num_types; i++)
      {
         send_info[i] = (HYPRE_BigInt)used_procs[i-1];
      }
      for (i=num_types+1; i< count; i++)
      {
         send_info[i] = new_vec_starts[i-num_types-1];
      }
      requests = hypre_CTAlloc(hypre_MPI_Request,  num_types, HYPRE_MEMORY_HOST);
      status =  hypre_CTAlloc(hypre_MPI_Status,  num_types, HYPRE_MEMORY_HOST);

      /* don't send to myself  - these are sorted so my id would be first*/
      start = 0;
      if (num_types && used_procs[0] == 0)
      {
         start = 1;
      }

      for (i=start; i < num_types; i++)
      {
         hypre_MPI_Isend(send_info, count, HYPRE_MPI_INT, used_procs[i], tag1,
                         comm, &requests[i-start]);
      }
      hypre_MPI_Waitall(num_types-start, requests, status);

      hypre_TFree(status, HYPRE_MEMORY_HOST);
      hypre_TFree(requests, HYPRE_MEMORY_HOST);
   }
   /* clean up */
   hypre_TFree(send_proc_obj.vec_starts, HYPRE_MEMORY_HOST);
   hypre_TFree(send_proc_obj.id, HYPRE_MEMORY_HOST);
   hypre_TFree(send_proc_obj.elements, HYPRE_MEMORY_HOST);
   hypre_TFree(send_info, HYPRE_MEMORY_HOST);
   if(response_recv_buf)        hypre_TFree(response_recv_buf, HYPRE_MEMORY_HOST);
   if(response_recv_buf_starts) hypre_TFree(response_recv_buf_starts, HYPRE_MEMORY_HOST);

   /* now proc 0 can exit if it has no rows */
   if (!local_num_rows)
   {
      if (hypre_CSRMatrixOwnsData(local_matrix))
         hypre_CSRMatrixDestroy(local_matrix);
      else
         hypre_TFree(local_matrix, HYPRE_MEMORY_HOST);

      hypre_TFree(new_vec_starts, HYPRE_MEMORY_HOST);
      hypre_TFree(used_procs, HYPRE_MEMORY_HOST);

      return NULL;
   }

   /* everyone left has rows and knows: new_vec_starts, num_types, and used_procs */

   /* this matrix should be rather small */
   matrix_i = hypre_CTAlloc(HYPRE_Int,  num_rows+1, HYPRE_MEMORY_HOST);

   num_requests = 4*num_types;
   requests = hypre_CTAlloc(hypre_MPI_Request,  num_requests, HYPRE_MEMORY_HOST);
   status = hypre_CTAlloc(hypre_MPI_Status,  num_requests, HYPRE_MEMORY_HOST);

   /* exchange contents of local_matrix_i - here we are sending to ourself also*/

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = (HYPRE_Int)(new_vec_starts[i+1] - new_vec_starts[i]);
      hypre_MPI_Irecv(&matrix_i[new_vec_starts[i]+1], vec_len, HYPRE_MPI_INT,
                      proc_id, tag2, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      hypre_MPI_Isend(&local_matrix_i[1], local_num_rows, HYPRE_MPI_INT,
                      proc_id, tag2, comm, &requests[j++]);
   }

   hypre_MPI_Waitall(j, requests, status);

   /* generate matrix_i from received data */
   /* global numbering?*/
   offset = matrix_i[new_vec_starts[1]];
   for (i=1; i < num_types; i++)
   {
      for (j = new_vec_starts[i]; j < new_vec_starts[i+1]; j++)
         matrix_i[j+1] += offset;
      offset = matrix_i[new_vec_starts[i+1]];
   }

   num_nonzeros = matrix_i[num_rows];

   matrix = hypre_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);

   hypre_CSRMatrixMemoryLocation(matrix) = HYPRE_MEMORY_HOST;

   hypre_CSRMatrixI(matrix) = matrix_i;
   hypre_CSRMatrixInitialize(matrix);
   matrix_j = hypre_CSRMatrixJ(matrix);
   matrix_data = hypre_CSRMatrixData(matrix);

   /* generate datatypes for further data exchange and exchange remaining
      data, i.e. column info and actual data */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      start_index = matrix_i[(HYPRE_Int)new_vec_starts[i]];
      num_data = matrix_i[(HYPRE_Int)new_vec_starts[i+1]] - start_index;
      hypre_MPI_Irecv(&matrix_data[start_index], num_data, HYPRE_MPI_COMPLEX,
                      used_procs[i], tag1, comm, &requests[j++]);
      hypre_MPI_Irecv(&matrix_j[start_index], num_data, HYPRE_MPI_INT,
                      used_procs[i], tag3, comm, &requests[j++]);
   }
   local_num_nonzeros = local_matrix_i[local_num_rows];
   for (i=0; i < num_types; i++)
   {
      hypre_MPI_Isend(local_matrix_data, local_num_nonzeros, HYPRE_MPI_COMPLEX,
                      used_procs[i], tag1, comm, &requests[j++]);
      hypre_MPI_Isend(local_matrix_j, local_num_nonzeros, HYPRE_MPI_INT,
                      used_procs[i], tag3, comm, &requests[j++]);
   }


   hypre_MPI_Waitall(num_requests, requests, status);

   hypre_TFree(new_vec_starts, HYPRE_MEMORY_HOST);

#else

   local_num_rows = (HYPRE_Int)(row_starts[my_id+1] - row_starts[my_id]);

   /* if my_id contains no data, return NULL */

   if (!local_num_rows)
      return NULL;

   local_matrix = hypre_MergeDiagAndOffd(par_matrix);
   hypre_CSRMatrixBigJtoJ(local_matrix); /* copies big_j to j */
   local_matrix_i = hypre_CSRMatrixI(local_matrix);
   local_matrix_j = hypre_CSRMatrixJ(local_matrix);
   local_matrix_data = hypre_CSRMatrixData(local_matrix);

   matrix_i = hypre_CTAlloc(HYPRE_Int,  num_rows+1, HYPRE_MEMORY_HOST);

   /* determine procs that have vector data and store their ids in used_procs */

   num_types = 0;
   for (i=0; i < num_procs; i++)
      if (row_starts[i+1]-row_starts[i] && i-my_id)
         num_types++;
   num_requests = 4*num_types;

   used_procs = hypre_CTAlloc(HYPRE_Int,  num_types, HYPRE_MEMORY_HOST);
   j = 0;
   for (i=0; i < num_procs; i++)
      if (row_starts[i+1]-row_starts[i] && i-my_id)
         used_procs[j++] = i;

   requests = hypre_CTAlloc(hypre_MPI_Request,  num_requests, HYPRE_MEMORY_HOST);
   status = hypre_CTAlloc(hypre_MPI_Status,  num_requests, HYPRE_MEMORY_HOST);
   /* data_type = hypre_CTAlloc(hypre_MPI_Datatype, num_types+1); */

   /* exchange contents of local_matrix_i */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      vec_len = (HYPRE_Int)(row_starts[proc_id+1] - row_starts[proc_id]);
      hypre_MPI_Irecv(&matrix_i[(HYPRE_Int)row_starts[proc_id]+1], vec_len, HYPRE_MPI_INT,
                      proc_id, 0, comm, &requests[j++]);
   }
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      hypre_MPI_Isend(&local_matrix_i[1], local_num_rows, HYPRE_MPI_INT,
                      proc_id, 0, comm, &requests[j++]);
   }

   vec_len = (HYPRE_Int)(row_starts[my_id+1] - row_starts[my_id]);
   for (i=1; i <= vec_len; i++)
      matrix_i[(HYPRE_Int)row_starts[my_id]+i] = local_matrix_i[i];

   hypre_MPI_Waitall(j, requests, status);

   /* generate matrix_i from received data */

   offset = matrix_i[(HYPRE_Int)row_starts[1]];
   for (i=1; i < num_procs; i++)
   {
      for (j = (HYPRE_Int)row_starts[i]; j < (HYPRE_Int)row_starts[i+1]; j++)
         matrix_i[j+1] += offset;
      offset = matrix_i[(HYPRE_Int)row_starts[i+1]];
   }

   num_nonzeros = matrix_i[num_rows];

   matrix = hypre_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);

   hypre_CSRMatrixMemoryLocation(matrix) = HYPRE_MEMORY_HOST;

   hypre_CSRMatrixI(matrix) = matrix_i;
   hypre_CSRMatrixInitialize(matrix);
   matrix_j = hypre_CSRMatrixJ(matrix);
   matrix_data = hypre_CSRMatrixData(matrix);

   /* generate datatypes for further data exchange and exchange remaining
      data, i.e. column info and actual data */

   j = 0;
   for (i = 0; i < num_types; i++)
   {
      proc_id = used_procs[i];
      start_index = matrix_i[(HYPRE_Int)row_starts[proc_id]];
      num_data = matrix_i[(HYPRE_Int)row_starts[proc_id+1]] - start_index;
      hypre_MPI_Irecv(&matrix_data[start_index], num_data, HYPRE_MPI_COMPLEX,
                      used_procs[i], 0, comm, &requests[j++]);
      hypre_MPI_Irecv(&matrix_j[start_index], num_data, HYPRE_MPI_INT,
                      used_procs[i], 0, comm, &requests[j++]);
   }
   local_num_nonzeros = local_matrix_i[local_num_rows];
   for (i=0; i < num_types; i++)
   {
      hypre_MPI_Isend(local_matrix_data, local_num_nonzeros, HYPRE_MPI_COMPLEX,
                      used_procs[i], 0, comm, &requests[j++]);
      hypre_MPI_Isend(local_matrix_j, local_num_nonzeros, HYPRE_MPI_INT,
                      used_procs[i], 0, comm, &requests[j++]);
   }

   start_index = matrix_i[(HYPRE_Int)row_starts[my_id]];
   for (i=0; i < local_num_nonzeros; i++)
   {
      matrix_j[start_index+i] = local_matrix_j[i];
      matrix_data[start_index+i] = local_matrix_data[i];
   }
   hypre_MPI_Waitall(num_requests, requests, status);
   start_index = matrix_i[(HYPRE_Int)row_starts[my_id]];
   for (i=0; i < local_num_nonzeros; i++)
   {
      matrix_j[start_index+i] = local_matrix_j[i];
      matrix_data[start_index+i] = local_matrix_data[i];
   }
   hypre_MPI_Waitall(num_requests, requests, status);

#endif

   if (hypre_CSRMatrixOwnsData(local_matrix))
      hypre_CSRMatrixDestroy(local_matrix);
   else
      hypre_TFree(local_matrix, HYPRE_MEMORY_HOST);

   if (num_requests)
   {
      hypre_TFree(requests, HYPRE_MEMORY_HOST);
      hypre_TFree(status, HYPRE_MEMORY_HOST);
      hypre_TFree(used_procs, HYPRE_MEMORY_HOST);
   }

   return matrix;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixCopy,
 * copies B to A,
 * if copy_data = 0, only the structure of A is copied to B
 * the routine does not check whether the dimensions of A and B are compatible
 *--------------------------------------------------------------------------*/

HYPRE_Int
hypre_ParCSRMatrixCopy( hypre_ParCSRMatrix *A,
                        hypre_ParCSRMatrix *B,
                        HYPRE_Int copy_data )
{
   hypre_CSRMatrix *A_diag;
   hypre_CSRMatrix *A_offd;
   HYPRE_BigInt *col_map_offd_A;
   hypre_CSRMatrix *B_diag;
   hypre_CSRMatrix *B_offd;
   HYPRE_BigInt *col_map_offd_B;
   HYPRE_Int num_cols_offd_A;
   HYPRE_Int num_cols_offd_B;

   if (!A)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }
   if (!B)
   {
      hypre_error_in_arg(1);
      return hypre_error_flag;
   }

   A_diag = hypre_ParCSRMatrixDiag(A);
   A_offd = hypre_ParCSRMatrixOffd(A);
   B_diag = hypre_ParCSRMatrixDiag(B);
   B_offd = hypre_ParCSRMatrixOffd(B);

   num_cols_offd_A = hypre_CSRMatrixNumCols(A_offd);
   num_cols_offd_B = hypre_CSRMatrixNumCols(B_offd);

   hypre_assert(num_cols_offd_A == num_cols_offd_B);

   col_map_offd_A = hypre_ParCSRMatrixColMapOffd(A);
   col_map_offd_B = hypre_ParCSRMatrixColMapOffd(B);

   hypre_CSRMatrixCopy(A_diag, B_diag, copy_data);
   hypre_CSRMatrixCopy(A_offd, B_offd, copy_data);

   /* should not happen if B has been initialized */
   if (num_cols_offd_B && col_map_offd_B == NULL)
   {
      col_map_offd_B = hypre_TAlloc(HYPRE_BigInt, num_cols_offd_B, HYPRE_MEMORY_HOST);
      hypre_ParCSRMatrixColMapOffd(B) = col_map_offd_B;
   }

   hypre_TMemcpy(col_map_offd_B, col_map_offd_A, HYPRE_BigInt, num_cols_offd_B,
                 HYPRE_MEMORY_HOST, HYPRE_MEMORY_HOST);

   return hypre_error_flag;
}

/*--------------------------------------------------------------------
 * hypre_FillResponseParToCSRMatrix
 * Fill response function for determining the send processors
 * data exchange
 *--------------------------------------------------------------------*/

HYPRE_Int
hypre_FillResponseParToCSRMatrix( void       *p_recv_contact_buf,
                                  HYPRE_Int   contact_size,
                                  HYPRE_Int   contact_proc,
                                  void       *ro,
                                  MPI_Comm    comm,
                                  void      **p_send_response_buf,
                                  HYPRE_Int *response_message_size )
{
   HYPRE_Int    myid;
   HYPRE_Int    i, index, count, elength;

   HYPRE_BigInt    *recv_contact_buf = (HYPRE_BigInt * ) p_recv_contact_buf;

   hypre_DataExchangeResponse  *response_obj = (hypre_DataExchangeResponse*)ro;

   hypre_ProcListElements      *send_proc_obj = (hypre_ProcListElements*)response_obj->data2;

   hypre_MPI_Comm_rank(comm, &myid );

   /*check to see if we need to allocate more space in send_proc_obj for ids*/
   if (send_proc_obj->length == send_proc_obj->storage_length)
   {
      send_proc_obj->storage_length +=10; /*add space for 10 more processors*/
      send_proc_obj->id = hypre_TReAlloc(send_proc_obj->id, HYPRE_Int,
                                         send_proc_obj->storage_length, HYPRE_MEMORY_HOST);
      send_proc_obj->vec_starts =
         hypre_TReAlloc(send_proc_obj->vec_starts, HYPRE_Int,
                        send_proc_obj->storage_length + 1, HYPRE_MEMORY_HOST);
   }

   /*initialize*/
   count = send_proc_obj->length;
   index = send_proc_obj->vec_starts[count]; /*this is the number of elements*/

   /*send proc*/
   send_proc_obj->id[count] = contact_proc;

   /*do we need more storage for the elements?*/
   if (send_proc_obj->element_storage_length < index + contact_size)
   {
      elength = hypre_max(contact_size, 10);
      elength += index;
      send_proc_obj->elements = hypre_TReAlloc(send_proc_obj->elements,
                                               HYPRE_BigInt,  elength, HYPRE_MEMORY_HOST);
      send_proc_obj->element_storage_length = elength;
   }
   /*populate send_proc_obj*/
   for (i=0; i< contact_size; i++)
   {
      send_proc_obj->elements[index++] = recv_contact_buf[i];
   }
   send_proc_obj->vec_starts[count+1] = index;
   send_proc_obj->length++;

   /*output - no message to return (confirmation) */
   *response_message_size = 0;

   return hypre_error_flag;
}

/*--------------------------------------------------------------------------
 * hypre_ParCSRMatrixUnion
 * Creates and returns a new matrix whose elements are the union of A and B.
 * Data is not copied, only structural information is created.
 * A and B must have the same communicator, numbers and distributions of rows
 * and columns (they can differ in which row-column pairs are nonzero, thus
 * in which columns are in a offd block)
 *--------------------------------------------------------------------------*/

hypre_ParCSRMatrix * hypre_ParCSRMatrixUnion( hypre_ParCSRMatrix * A,
                                              hypre_ParCSRMatrix * B )
{
   hypre_ParCSRMatrix * C;
   HYPRE_BigInt * col_map_offd_C = NULL;
   HYPRE_Int  num_procs, my_id, p;
   MPI_Comm comm = hypre_ParCSRMatrixComm( A );

   hypre_MPI_Comm_rank(comm,&my_id);
   hypre_MPI_Comm_size(comm,&num_procs);

   C = hypre_CTAlloc( hypre_ParCSRMatrix,  1 , HYPRE_MEMORY_HOST);
   hypre_ParCSRMatrixComm( C ) = hypre_ParCSRMatrixComm( A );
   hypre_ParCSRMatrixGlobalNumRows( C ) = hypre_ParCSRMatrixGlobalNumRows( A );
   hypre_ParCSRMatrixGlobalNumCols( C ) = hypre_ParCSRMatrixGlobalNumCols( A );
   hypre_ParCSRMatrixFirstRowIndex( C ) = hypre_ParCSRMatrixFirstRowIndex( A );
   hypre_assert( hypre_ParCSRMatrixFirstRowIndex( B )
                 == hypre_ParCSRMatrixFirstRowIndex( A ) );
   hypre_ParCSRMatrixRowStarts( C ) = hypre_ParCSRMatrixRowStarts( A );
   hypre_ParCSRMatrixOwnsRowStarts( C ) = 0;
   hypre_ParCSRMatrixColStarts( C ) = hypre_ParCSRMatrixColStarts( A );
   hypre_ParCSRMatrixOwnsColStarts( C ) = 0;
   for ( p=0; p<=num_procs; ++p )
      hypre_assert( hypre_ParCSRMatrixColStarts(A)
                    == hypre_ParCSRMatrixColStarts(B) );
   hypre_ParCSRMatrixFirstColDiag( C ) = hypre_ParCSRMatrixFirstColDiag( A );
   hypre_ParCSRMatrixLastRowIndex( C ) = hypre_ParCSRMatrixLastRowIndex( A );
   hypre_ParCSRMatrixLastColDiag( C ) = hypre_ParCSRMatrixLastColDiag( A );

   hypre_ParCSRMatrixDiag( C ) =
      hypre_CSRMatrixUnion( hypre_ParCSRMatrixDiag(A), hypre_ParCSRMatrixDiag(B),
                            0, 0, 0 );
   hypre_ParCSRMatrixOffd( C ) =
      hypre_CSRMatrixUnion( hypre_ParCSRMatrixOffd(A), hypre_ParCSRMatrixOffd(B),
                            hypre_ParCSRMatrixColMapOffd(A),
                            hypre_ParCSRMatrixColMapOffd(B), &col_map_offd_C );
   hypre_ParCSRMatrixColMapOffd( C ) = col_map_offd_C;
   hypre_ParCSRMatrixCommPkg( C ) = NULL;
   hypre_ParCSRMatrixCommPkgT( C ) = NULL;
   hypre_ParCSRMatrixOwnsData( C ) = 1;
   /*  SetNumNonzeros, SetDNumNonzeros are global, need hypre_MPI_Allreduce.
       I suspect, but don't know, that other parts of hypre do not assume that
       the correct values have been set.
       hypre_ParCSRMatrixSetNumNonzeros( C );
       hypre_ParCSRMatrixSetDNumNonzeros( C );*/
   hypre_ParCSRMatrixNumNonzeros( C ) = 0;
   hypre_ParCSRMatrixDNumNonzeros( C ) = 0.0;
   hypre_ParCSRMatrixRowindices( C ) = NULL;
   hypre_ParCSRMatrixRowvalues( C ) = NULL;
   hypre_ParCSRMatrixGetrowactive( C ) = 0;

   return C;
}


/* drop the entries that are not on the diagonal and smaller than
 * its row norm: type 1: 1-norm, 2: 2-norm, -1: infinity norm */
HYPRE_Int
hypre_ParCSRMatrixDropSmallEntries( hypre_ParCSRMatrix *A,
                                    HYPRE_Real tol,
                                    HYPRE_Int type)
{
   HYPRE_Int i, j, k, nnz_diag, nnz_offd, A_diag_i_i, A_offd_i_i;

   MPI_Comm         comm     = hypre_ParCSRMatrixComm(A);
   /* diag part of A */
   hypre_CSRMatrix *A_diag   = hypre_ParCSRMatrixDiag(A);
   HYPRE_Real      *A_diag_a = hypre_CSRMatrixData(A_diag);
   HYPRE_Int       *A_diag_i = hypre_CSRMatrixI(A_diag);
   HYPRE_Int       *A_diag_j = hypre_CSRMatrixJ(A_diag);
   /* off-diag part of A */
   hypre_CSRMatrix *A_offd   = hypre_ParCSRMatrixOffd(A);
   HYPRE_Real      *A_offd_a = hypre_CSRMatrixData(A_offd);
   HYPRE_Int       *A_offd_i = hypre_CSRMatrixI(A_offd);
   HYPRE_Int       *A_offd_j = hypre_CSRMatrixJ(A_offd);

   HYPRE_Int  num_cols_A_offd = hypre_CSRMatrixNumCols(A_offd);
   HYPRE_BigInt *col_map_offd_A  = hypre_ParCSRMatrixColMapOffd(A);
   HYPRE_Int *marker_offd = NULL;

   HYPRE_BigInt first_row  = hypre_ParCSRMatrixFirstRowIndex(A);
   HYPRE_Int nrow_local = hypre_CSRMatrixNumRows(A_diag);
   HYPRE_Int my_id, num_procs;
   /* MPI size and rank*/
   hypre_MPI_Comm_size(comm, &num_procs);
   hypre_MPI_Comm_rank(comm, &my_id);

   if (tol <= 0.0)
   {
      return hypre_error_flag;
   }

   marker_offd = hypre_CTAlloc(HYPRE_Int, num_cols_A_offd, HYPRE_MEMORY_HOST);

   nnz_diag = nnz_offd = A_diag_i_i = A_offd_i_i = 0;
   for (i = 0; i < nrow_local; i++)
   {
      /* compute row norm */
      HYPRE_Real row_nrm = 0.0;
      for (j = A_diag_i_i; j < A_diag_i[i+1]; j++)
      {
         HYPRE_Complex v = A_diag_a[j];
         if (type == 1)
         {
            row_nrm += fabs(v);
         }
         else if (type == 2)
         {
            row_nrm += v*v;
         }
         else
         {
            row_nrm = hypre_max(row_nrm, fabs(v));
         }
      }
      if (num_procs > 1)
      {
         for (j = A_offd_i_i; j < A_offd_i[i+1]; j++)
         {
            HYPRE_Complex v = A_offd_a[j];
            if (type == 1)
            {
               row_nrm += fabs(v);
            }
            else if (type == 2)
            {
               row_nrm += v*v;
            }
            else
            {
               row_nrm = hypre_max(row_nrm, fabs(v));
            }
         }
      }

      if (type == 2)
      {
         row_nrm = sqrt(row_nrm);
      }

      /* drop small entries based on tol and row norm */
      for (j = A_diag_i_i; j < A_diag_i[i+1]; j++)
      {
         HYPRE_Int     col = A_diag_j[j];
         HYPRE_Complex val = A_diag_a[j];
         if (i == col || fabs(val) >= tol * row_nrm)
         {
            A_diag_j[nnz_diag] = col;
            A_diag_a[nnz_diag] = val;
            nnz_diag ++;
         }
      }
      if (num_procs > 1)
      {
         for (j = A_offd_i_i; j < A_offd_i[i+1]; j++)
         {
            HYPRE_Int     col = A_offd_j[j];
            HYPRE_Complex val = A_offd_a[j];
            /* in normal cases: diagonal entry should not
             * appear in A_offd (but this can still be possible) */
            if (i + first_row == col_map_offd_A[col] || fabs(val) >= tol * row_nrm)
            {
               if (0 == marker_offd[col])
               {
                  marker_offd[col] = 1;
               }
               A_offd_j[nnz_offd] = col;
               A_offd_a[nnz_offd] = val;
               nnz_offd ++;
            }
         }
      }
      A_diag_i_i = A_diag_i[i+1];
      A_offd_i_i = A_offd_i[i+1];
      A_diag_i[i+1] = nnz_diag;
      A_offd_i[i+1] = nnz_offd;
   }

   hypre_CSRMatrixNumNonzeros(A_diag) = nnz_diag;
   hypre_CSRMatrixNumNonzeros(A_offd) = nnz_offd;
   hypre_ParCSRMatrixSetNumNonzeros(A);
   hypre_ParCSRMatrixDNumNonzeros(A) = (HYPRE_Real) hypre_ParCSRMatrixNumNonzeros(A);

   for (i = 0, k = 0; i < num_cols_A_offd; i++)
   {
      if (marker_offd[i])
      {
         col_map_offd_A[k] = col_map_offd_A[i];
         marker_offd[i] = k++;
      }
   }
   /* num_cols_A_offd = k; */
   hypre_CSRMatrixNumCols(A_offd) = k;
   for (i = 0; i < nnz_offd; i++)
   {
      A_offd_j[i] = marker_offd[A_offd_j[i]];
   }

   if ( hypre_ParCSRMatrixCommPkg(A) )
   {
      hypre_MatvecCommPkgDestroy( hypre_ParCSRMatrixCommPkg(A) );
   }
   hypre_MatvecCommPkgCreate(A);

   hypre_TFree(marker_offd, HYPRE_MEMORY_HOST);

   return hypre_error_flag;
}

/*
#ifdef HYPRE_USING_UNIFIED_MEMORY
hypre_int hypre_ParCSRMatrixIsManaged(hypre_ParCSRMatrix *a){
  if (hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(a)))
    return ((hypre_CSRMatrixIsManaged(hypre_ParCSRMatrixDiag(a))) && (hypre_CSRMatrixIsManaged(hypre_ParCSRMatrixOffd(a))));
  else
    return hypre_CSRMatrixIsManaged(hypre_ParCSRMatrixDiag(a));
}
#endif
*/
