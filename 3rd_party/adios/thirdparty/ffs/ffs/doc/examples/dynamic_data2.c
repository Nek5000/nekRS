#include "ffs.h"
#include "cod.h"
#include "fcntl.h"

int main()     /* receiving program */
{
    FFSTypeHandle handle;
    FMFormat fmf;
    FMStructDescList local_struct;
    FFSContext ffsc = create_FFSContext();
    int fd, encode_size, decode_size;
    char encoded_buffer[2048];  /* hopefully big enough */
    void *unknown_record;
    cod_parse_context context = new_cod_parse_context();
    cod_code gen_code;
    int (*func)(void *);

    /* "receive" encoded record over a file */
    fd = open("enc_file", O_RDONLY, 0777);
    encode_size = read(fd, encoded_buffer, sizeof(encoded_buffer));

    handle = FFSTypeHandle_from_encode(ffsc, encoded_buffer);
    fmf = FMFormat_of_original(handle);
    local_struct = get_localized_formats(fmf);
   
    /* decode the data into an acceptable local format */
    FFSset_fixed_target(ffsc, local_struct);
    FFS_target_from_encode(ffsc, encoded_buffer);
    decode_size = FFS_est_decode_length(ffsc, encoded_buffer, encode_size);
    unknown_record = malloc(decode_size);
    FFSdecode(ffsc, encoded_buffer, unknown_record);

    /* generate code to operate on it */
    cod_add_struct_type(local_struct, context);
    cod_add_param("input", FFSTypeHandle_name(handle), 0, context);
    gen_code = cod_code_gen("{return input.j", context);
    if (!gen_code) {
        printf("The input did not have an acceptable field 'j'\n");
    } else {
        func = (long(*)()) gen_code->func;
        printf("Field 'j' in the input has value %d\n", func(unknown_record));
    }
}
