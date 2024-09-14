#include "evpath.h"
#include "messages.h"
#include <stdio.h>


static int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    custom_sensor_ctrl_ptr event = vevent;
    printf("%s  | ",event->IPMI_Hostname);
    printf(" %s | ",event->IPMI_Sensor_ID);
    printf(" %X | ",event->IPMI_Sensor_Number);
    printf(" %s | ",event->IPMI_Sensor_Status);
    printf(" %s | ", event->IPMI_Sensor_Reading);
    printf(" %d | ",event->IPMI_Entity_ID);
    printf(" %d | ",event->IPMI_Instance);
    printf(" %s\n\n",event->IPMI_Timestamp);

}

int main(int argc, char **argv)
{
    CManager cm;
    EVstone stone;
    char *string_list;

    cm = CManager_create();
    CMlisten(cm);

    stone = EValloc_stone(cm);
    EVassoc_terminal_action(cm, stone, custom_sensor_format_list, simple_handler, NULL);
    string_list = attr_list_to_string(CMget_contact_list(cm));
    printf("Contact list %d:%s\n", stone, string_list);
    CMrun_network(cm);
}
