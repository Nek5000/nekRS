/*
 * ioGroup.h
 *
 *  Created on: Nov 2018
 *      Author: Norbert Podhorszki
 */

#ifndef IOGROUP_H
#define IOGROUP_H

#include "adios2.h"
#include "settings.h"

#include <map>
#include <string>

class ioGroup
{
public:
    const std::string name;
    adios2::IO adiosio;
    ioGroup(const std::string &name) : name(name){};
    virtual ~ioGroup() = 0;
};

class adiosIOGroup : public ioGroup
{
public:
    adiosIOGroup(const std::string &name, adios2::ADIOS &adiosobj) : ioGroup(name)
    {
        adiosio = adiosobj.DeclareIO(name);
    };
    ~adiosIOGroup(){};
};

class hdf5IOGroup : public ioGroup
{
public:
    hdf5IOGroup(const std::string &name) : ioGroup(name){};
    ~hdf5IOGroup(){};
};

std::shared_ptr<ioGroup> createGroup(const std::string &name, IOLib iolib, adios2::ADIOS &adiosobj);

#endif /* IOGROUP_H */
