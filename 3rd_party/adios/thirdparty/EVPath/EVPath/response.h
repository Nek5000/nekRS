extern int response_data_free(CManager cm, void *resp_void);

extern void *
install_response_handler(CManager cm, int stone_id, char *response_spec, 
			 void *local_data, FMFormat **ref_ptr);

extern void
dump_mrd(void *mrd);

extern int
response_determination(CManager cm, stone_type stone, action_class stage, event_item *event);

