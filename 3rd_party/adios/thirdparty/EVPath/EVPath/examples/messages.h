typedef struct _custom_sensor {
        char IPMI_Hostname[25];
        char IPMI_Sensor_ID[25];
        int IPMI_Sensor_Number;
        char IPMI_Sensor_Status[25];
        int IPMI_Entity_ID;
        int IPMI_Instance;
        char IPMI_Sensor_Reading[25];
        char IPMI_Timestamp[100];
} custom_sensor_ctrl, *custom_sensor_ctrl_ptr;

static FMField custom_sensor_field_list[] =
{
    {"IPMI_Hostname", "char[25]", sizeof(char), FMOffset(custom_sensor_ctrl_ptr, IPMI_Hostname)},
    {"IPMI_Sensor_ID", "char[25]", sizeof(char), FMOffset(custom_sensor_ctrl_ptr, IPMI_Sensor_ID)},
    {"IPMI_Sensor_Number", "integer", sizeof(int), FMOffset(custom_sensor_ctrl_ptr, IPMI_Sensor_Number)},
    {"IPMI_Sensor_Status", "char[25]", sizeof(char), FMOffset(custom_sensor_ctrl_ptr, IPMI_Sensor_Status)},
    {"IPMI_Entity_ID", "integer", sizeof(int), FMOffset(custom_sensor_ctrl_ptr, IPMI_Entity_ID)},
    {"IPMI_Instance", "integer", sizeof(int), FMOffset(custom_sensor_ctrl_ptr, IPMI_Instance)},
    {"IPMI_Sensor_Reading", "char[25]", sizeof(char), FMOffset(custom_sensor_ctrl_ptr, IPMI_Sensor_Reading)},
    {"IPMI_Timestamp", "char[100]", sizeof(char), FMOffset(custom_sensor_ctrl_ptr, IPMI_Timestamp)},
    {NULL, NULL, 0, 0}
};
static FMStructDescRec custom_sensor_format_list[] =
{
    {"custom_sensor", custom_sensor_field_list, sizeof(custom_sensor_ctrl), NULL},
    {NULL, NULL}
};



