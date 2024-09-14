/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include "settings.h"

#include <fstream>

#include "json.hpp"

void to_json(nlohmann::json &j, const Settings &s)
{
    j = nlohmann::json{{"L", s.L},
                       {"steps", s.steps},
                       {"plotgap", s.plotgap},
                       {"F", s.F},
                       {"k", s.k},
                       {"dt", s.dt},
                       {"Du", s.Du},
                       {"Dv", s.Dv},
                       {"noise", s.noise},
                       {"output", s.output},
                       {"checkpoint", s.checkpoint},
                       {"checkpoint_freq", s.checkpoint_freq},
                       {"checkpoint_output", s.checkpoint_output},
                       {"restart", s.restart},
                       {"restart_input", s.restart_input},
                       {"adios_config", s.adios_config},
                       {"adios_span", s.adios_span},
                       {"adios_memory_selection", s.adios_memory_selection},
                       {"mesh_type", s.mesh_type}};
}

void from_json(const nlohmann::json &j, Settings &s)
{
    j.at("L").get_to(s.L);
    j.at("steps").get_to(s.steps);
    j.at("plotgap").get_to(s.plotgap);
    j.at("F").get_to(s.F);
    j.at("k").get_to(s.k);
    j.at("dt").get_to(s.dt);
    j.at("Du").get_to(s.Du);
    j.at("Dv").get_to(s.Dv);
    j.at("noise").get_to(s.noise);
    j.at("output").get_to(s.output);
    j.at("checkpoint").get_to(s.checkpoint);
    j.at("checkpoint_freq").get_to(s.checkpoint_freq);
    j.at("checkpoint_output").get_to(s.checkpoint_output);
    j.at("restart").get_to(s.restart);
    j.at("restart_input").get_to(s.restart_input);
    j.at("adios_config").get_to(s.adios_config);
    j.at("adios_span").get_to(s.adios_span);
    j.at("adios_memory_selection").get_to(s.adios_memory_selection);
    j.at("mesh_type").get_to(s.mesh_type);
}

Settings::Settings()
{
    L = 128;
    steps = 20000;
    plotgap = 200;
    F = 0.04;
    k = 0.06075;
    dt = 0.2;
    Du = 0.05;
    Dv = 0.1;
    noise = 0.0;
    output = "foo.bp";
    checkpoint = false;
    checkpoint_freq = 2000;
    checkpoint_output = "ckpt.bp";
    restart = false;
    restart_input = "ckpt.bp";
    adios_config = "adios2.xml";
    adios_span = false;
    adios_memory_selection = false;
    mesh_type = "image";
}

Settings Settings::from_json(const std::string &fname)
{
    std::ifstream ifs(fname);
    nlohmann::json j;

    ifs >> j;

    return j.get<Settings>();
}
