#ifndef __COD__H__
#define __COD__H__
/*! \file */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <stdint.h>

/*!
 * cod_parse_context is the basic handle controlling generation of
 * subroutines and functions.  It is used in nearly all calls to setup
 * parameter profiles and other information and to generate the code.  Using
 * a separate handle rather than static variables to store this information
 * over the multiple calls needed for code generation ensures thread
 * safety. 
 */
typedef struct parse_struct *cod_parse_context;

/*!
 * cod_exec_context is used to provide an execution context for 
 * COD-generated procedures.  The full set of uses for this have not been 
 * determined, but it will support execution emulation on platforms where 
 * we do not have dynamic code generation, extraction of static state from
 * functions, and the ability to pass information "through" a generated
 * function from the calling context and make it available to a subroutine
 * that might be called by the generated code.
 */
typedef struct cod_exec_struct *cod_exec_context;

/*!
 * sm_ref is a pointer to an internal COD data structure representing a
 * semantic element of a program fragment (E.G. a parameter definition).
 */
typedef struct sm_struct *sm_ref;

#ifndef FMOffset
#define FMOffset(p_type,field) \
	((int) (((char *) (&(((p_type)0)->field))) - ((char *) 0)))
#if defined(__STDC__) || defined(__ANSI_CPP__) || defined(_MSC_VER)
#define FMstr(s) #s
#else
#define FMstr(s) "s"
#endif
#define FMArrayDecl(type, size) FMstr(type[size])
#define FMArrayDecl2(type, size, size2) FMstr(type[size][size2])

#define FMDefaultDecl(name, val) FMstr(name(val))

typedef struct _FMField {
    const char *field_name;	/* Field name */
    const char *field_type;	/* Representation type desired */
    int field_size;		/* Size in bytes of representation */
    int field_offset;		/* Offset from base to put field value */
} FMField, *FMFieldList;

typedef struct _FMOptInfo {
    int info_type;
    int info_len;
    char *info_block;
} FMOptInfo;

/*!
 * A structure to hold Format Name / Field List associations.
 *
 *
 *  This is used to associate names with field lists.  Together these define 
 *  a structure that can be composed into larger structures.
 */
typedef struct _FMformat_list {
    /*! the name to be associated with this structure */
    char *format_name;
    /*! the FFS-style list of fields within this structure */
    FMFieldList field_list;
    int struct_size;
    FMOptInfo *opt_info;
}FMStructDescRec, *FMStructDescList;
#endif

/*!
  A structure to hold external entries.

  This is used to associate textual names with addresses.
*/
typedef struct extern_entry {
    /*! the textual name of the external entry */
    char *extern_name;
    /*! the address of the external entry */
    void *extern_value;
} cod_extern_entry;

/*!
 * A list of cod_extern_entry structures.  This is used to specify the
 addresses of external  functions or variables that a function may access.
*/
typedef cod_extern_entry *cod_extern_list;

/*! Create a handle to a new cod_parse_context.
 *
 * This is generally the first call in dynamically generating a
 * subroutine/function.  The cod_parse_context created here is used in all
 * future calls, setting up the parameter profile and generating the code.
 * \return Will return a new initialized cod_parse_context unless there is
 * no available memory.
 */
extern cod_parse_context new_cod_parse_context(void);

/*!
 * Free a handle to an cod_parse_context.
 *
 * \param context the cod_parse_context to be free'd.
 * Calling this routine frees all memory associated with the parse context,
 * but not that of code that has been generated from this context.
 */
extern void cod_free_parse_context(cod_parse_context context);

/*!
 * Associate a set of "name, external address" pairs with a parse context 
 *
 *  These names and addresses are used to resolve (external references in
 *  generated code).  In addition to associating the addresses, the full
 *  prototypes of the functions or variables must be added to the parse
 *  context, usually with cod_parse_for_context().
 *
 * \param context the cod_parse_context to which the externs are associated.
 * \param externs the list of "name, external address" pairs to be
 *  associated.  This list should be terminated with a {NULL, 0} pair.
 */
extern void cod_assoc_externs(cod_parse_context context,
				    cod_extern_list externs);

/*!
 * \brief This is used to establish the parameter profile and return type 
 * of the subroutine to be generated.
 *
 * \param decl  the textual declaration for the subroutine
 * \param context  the context for subroutine generation
 * 
 * 
 *  The 'decl' parameter is a C-style declaration, without a 
 *  terminating semicolon.  I.E. something like:
 * 	"int proc(double d, int *i)"
 *  
 */
extern void
cod_subroutine_declaration(const char *decl, cod_parse_context context);

/*!
 * \brief Make a structured type, represented by a FFS-style FMFieldList,
 *  available in a particular context.
 *
 * \param name the name to be associated with the new structured type.
 * \param field_list the FFS-style FMFieldList that describes the layout of
 *  the structure.
 * \param context the context in which the type is to be made available.
 */
extern void cod_add_simple_struct_type(const char *name, FMFieldList field_list, 
				       cod_parse_context context);

/*!
 * \brief Make a set of structured types, represented by a FFS-style FMStructDescList,
 *  available in a particular context.
 *
 * \param format_list the FFS-style FMStructDescList that describes the layout of
 *  the structures.
 * \param context the context in which the type is to be made available.
 */
extern void cod_add_struct_type(FMStructDescList format_list, 
				cod_parse_context context);

/*!
 * \brief As an alternative to cod_subroutine_declaration(), add declaration for the Nth parameter.
 *
 * \param name the name to be associated with the new parameter.
 * \param typ the data type of the new parameter.
 * \param param_num the numeral of the new parameter (0 is first)
 * \param context the context in which the subroutine is being declared.
 */
extern void
cod_add_param(const char *id, const char *typ, int param_num, 
	      cod_parse_context context);

/*!
 * \brief add a declaration for the Nth parameter as an FFS-encoded structure
 *
 * \param name the name to be associated with the new parameter.
 * \param data an encoded data block 
 * \param param_num the numeral of the new parameter (0 is first)
 * \param c an FMContext value in which to format information for the data is available
 * \param context the context in which the subroutine is being declared.
 */
#ifdef __FM__
extern void
cod_add_encoded_param(const char *id, char *data, int param_num, 
		      FMContext c, cod_parse_context context);
#endif

/*!
 * \brief As an alternative to cod_subroutine_declaration(), set the return type of a subroutine
 *
 * \param typ the data type of the return value.
 * \param context the context in which the subroutine is being declared.
 */
extern void
cod_set_return_type(char *typ, cod_parse_context context);

/*!
 * cod_code is a handle to the generated code.  In addition to the generated
 * code block, it holds external data pointers and supplemental information.
 */
typedef struct _cod_code_struct {
    /*! the entry point of the generated function */
    void (*func)(void);
    /*! set to true if the code requires an cod_exec_ctx to run */
    int has_exec_ctx;
    /*! the memory block containing the code of the generated function */
    void *code_memory_block;
    /*! the static data block referenced by the generated function */
    void *data;
    /*! the size of the static data block referenced by the generated function */
    int static_size_required;
    /*! the virtual register that holds the address of the static data block  */
    int static_block_address_register;
    /*! the IOFormats that describe the static data block */
    FMStructDescList static_formats;
    /*! a pointer to the DRISC context used to generate the function \internal 
    */
    void *drisc_context;
    /*! a pointer to the DRISC execution handle used to free the function \internal 
    */
    void *execution_handle;
} *cod_code;

/*!
 * generate code for a function body in a given context
 *
 * \param code the string representing the function body.
 * \param context the context in which the function body is to be generated.
 */
cod_code cod_code_gen(char *code, cod_parse_context context);

/*!
 * perform syntactical and semantic checking of a function body without
 * actually generating code.
 *
 * \param code the string representing the function body.
 * \param context the context in which the function body is to be checked.
 */
int cod_code_verify(char *code, cod_parse_context context);

/*!
 * Free all resources associated with the generated code associated with the
 * cod_code handle.
 *
 * \param code the handle to the resources that will be free'd.
 */
extern void cod_code_free(cod_code code);

/*!
 * create an execution context associated with a code block
 *
 * \param code the handle to the code bloc
 * \return the created execution context
 */
extern cod_exec_context cod_create_exec_context(cod_code code);

/*!
 * Free all resources associated with the generated code associated with the
 * cod_exec_context handle.
 *
 * \param code the handle to the resources that will be free'd.
 */
extern void cod_exec_context_free(cod_exec_context ec);

/*!
 * Associate application-level data with an execution context.  This is
 * typically done so that it can be retrieved by a called subprocess.  The
 * key value is any 32-bit value, presumably chosen to minimize the
 * possibility of conflict.  The value that can be associated is a long.
 *
 * \param ec the cod_exec_context value to which to associate the data
 * \param key the value that will serve as a key to retrieve the data
 * \param value the 'long' data that will be associated with the key
 */
extern void cod_assoc_client_data(cod_exec_context ec, int key, intptr_t value);

/*!
 * Retrieve application-level data with an execution context.  This is
 * typically done by a called subprocess to retrieve data that was available
 * to the broader application.  The key value is any 32-bit value,
 * presumably chosen to minimize the possibility of conflict.  The value
 * returned is one which had been associated via cod_assoc_client_data().
 * If no value had been associated, -1 is returned. * 
 *
 * \param ec the cod_exec_context value to which to associate the data
 * \param key the value that will serve as a key to retrieve the data
 * \return the 'long' data that was associated with the key
 */
extern intptr_t cod_get_client_data(cod_exec_context ec, int key);

/*!
 * Extract static state from an execution context.  
 *
 * \param ec the cod_exec_context value from which to extract static state
 * \param length_pointer a pointer to an integer which will be set to 
 *  the length of the returned state block
 * \return a pointer to the extracted state
 */
extern void *cod_extract_state(cod_exec_context ec, int *length_p);

/*!
 * Install static state into an execution context.  
 *
 * \param ec the cod_exec_context value into which to install static state
 * \param state the state block to install
 * \param state_size the size of the state block
 * \return 1 on success, 0 on failure
 */
extern int cod_install_state(cod_exec_context ec, void *state, int length);

/*!
 * \brief This parses a string to determine what external
 *  functions and variables are visible during code generation.
 *
 * The input string is of standard C syntax for external variable and
 * function declarations.  All previously registered types are available.
 * \param code The input string containing declarations.
 * \param context The parse context in which the declarations should be
 * visible.
*/
int cod_parse_for_context(char *code, cod_parse_context context);

/*!
 * \brief This parses a string to setup global variables that are 
 *  visible during code generation.
 *
 * The input string is of standard C syntax for external variable and
 * function declarations.  All previously registered types are available.
 * \param code The input string containing declarations.
 * \param context The parse context in which the declarations should be
 * visible.
*/
int cod_parse_for_globals(char *code, cod_parse_context context);

/*!
 * Duplicate a handle to an cod_parse_context.
 *
 * \param context the cod_parse_context to be duplicated.
 */
extern cod_parse_context cod_copy_context(cod_parse_context context);

/*!
 * Duplicate a handle to an cod_parse_context, specifically adapting the results to 
 * reuse the global declarations of the prior context.  This is a specialized 
 * function that helps run external regression tests and may not be broadly useful.
 *
 * \param context the cod_parse_context to be duplicated.
 */
extern cod_parse_context cod_copy_globals(cod_parse_context context);

/*!
 *  err_out_func_t is a function pointer type.   Functions matching this
 *  profile can be used as call-out handlers for COD errors.
 *  \param client_data an uninspected value passed in from 
 *    cod_set_error_func()
 *  \param string the textual representation of the error.
*/
typedef void (*err_out_func_t)(void *client_data, char *string);

/*!
 * cod_set_error_func establishes a new error output routine for COD.
 *
 * By default, errors are sent to stdout.  If the error output is set then
 * errors can be captured for other disposition.
 *  \param context the context in which errors are to be captured
 *  \param err_func the function to be called when errors occur
 */
void cod_set_error_func(cod_parse_context context, 
			      err_out_func_t err_func);

/*!
 * cod_set_dont_coerce_return restricts COD to more strict type matching for expressions and return values
 *
 *  \param context the context to restrict
 *  \param value True if coercion is not to be applied to return values, false by default
 */
void cod_set_dont_coerce_return(cod_parse_context context, int value);

/*!
 * This will dump (to stdout) a disassembled version of the 
 * machine code that has been generated
 *  \param code the cod_code handle containing the code to be dumped.
 */
void cod_dump(cod_code code);

/*!
 * This will generate rollback code for message morphing 
 *  \param format1 the new format. This is the format of the input message.
 *  \param format2 the old format. This is the format of the output message.
 *  \param xform_code  The COD code string that transforms data from format1 to format2.
 */
extern cod_code
gen_rollback_code(FMStructDescList format1, FMStructDescList format2, char *xform_code);

/*!
 * \brief This is used to add an integer constant available in a particular context.
 * \param id The name to be associated with the type.
 * \param value The value of the constant
 * \param context the context in which this is to be created
*/
extern void cod_add_int_constant_to_parse_context(const char *id, int value,
    cod_parse_context context); 

/*!
 * \brief Sets what value cod_closure_context will send to a subroutine.
 * \param name The name of the (external) subroutine.
 * \param value The value to send
 * \param context The context in which the subroutine has been declared.
*/
extern void cod_set_closure(char *name, void* value, cod_parse_context context);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif
