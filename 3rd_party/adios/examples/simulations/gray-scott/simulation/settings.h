/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>

struct Settings
{
    size_t L;
    int steps;
    int plotgap;
    double F;
    double k;
    double dt;
    double Du;
    double Dv;
    double noise;
    std::string output;
    bool checkpoint;
    int checkpoint_freq;
    std::string checkpoint_output;
    bool restart;
    std::string restart_input;
    std::string adios_config;
    bool adios_span;
    bool adios_memory_selection;
    std::string mesh_type;

    Settings();
    static Settings from_json(const std::string &fname);
};

#endif
